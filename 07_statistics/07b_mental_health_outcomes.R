#!/usr/bin/env Rscript
# ==============================================================================
# 07b_mental_health_outcomes.R
#
# Mental health outcomes: depressive symptoms and subjective wellbeing in the
# matched sample.
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Cognitive and mental health outcomes -> Mental health outcomes;
#          Supplementary material, "Mental health outcome measures"
#
# Depressive symptoms are the mean of four self-report items (low mood,
# anhedonia, restlessness, tiredness; 0-3 scale) and subjective wellbeing the
# mean of six items (happiness, health, work, family, friendships, finances;
# 1-6 scale, reverse-scored so that higher values indicate better wellbeing).
# Missing items are imputed by within-person means where possible, or by
# group-specific medians for total non-response. Group differences use the
# same regression framework and covariate set as the cognitive analyses:
# linear regression with Type III sums of squares, adjusting for age, sex,
# education, diastolic blood pressure and log-transformed total WMH; adjusted
# Cohen's d = (AF - controls) / residual SD; Benjamini-Hochberg FDR across the
# two outcomes.
#
# Input:  DERIV_DIR/cohort_matched.csv          (01_cohort_assembly.R)
# Output: OUTPUT_DIR/mental_health/mental_health_results.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(car); library(emmeans)
})

deriv_dir <- Sys.getenv("DERIV_DIR",  "derivatives")
out_root  <- Sys.getenv("OUTPUT_DIR", "output")
in_csv    <- file.path(deriv_dir, "cohort_matched.csv")
out_dir   <- file.path(out_root, "mental_health")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(in_csv))

core_covariates <- c("age_i2", "sex", "education_i2", "dia_mean", "log_norm_wmh_ml")

dep_items <- c("mh_dep_mood_2w_i2", "mh_anhedonia_2w_i2",
               "mh_tense_restless_2w_i2", "mh_tired_lethargic_2w_i2")
wb_items  <- c("wb_happiness_i2", "wb_health_satisfaction_i2", "wb_work_satisfaction_i2",
               "wb_family_satisfaction_i2", "wb_friendship_satisfaction_i2",
               "wb_financial_satisfaction_i2")

df <- read_csv(in_csv, show_col_types = FALSE, col_types = cols(eid = col_character())) %>%
  mutate(group = factor(af, levels = c(0, 1), labels = c("control", "AF")))
message(sprintf("Matched sample: n = %d (AF %d, controls %d)", nrow(df), sum(df$af == 1), sum(df$af == 0)))

# ---- item cleaning ---------------------------------------------------------
# Non-response codes (-1 do not know, -3 prefer not to answer; 7 = work
# satisfaction "not employed") are set to missing. Depressive items are coded
# 1-4 in UK Biobank and rescaled to 0-3; wellbeing items (1 = extremely happy /
# satisfied ... 6 = extremely unhappy / unsatisfied) are reverse-scored.
df <- df %>%
  mutate(across(all_of(c(dep_items, wb_items)), ~ ifelse(.x %in% c(-1, -3, 7), NA_real_, .x)),
         across(all_of(dep_items), ~ .x - 1),
         across(all_of(wb_items),  ~ 7 - .x))

impute_items <- function(data, items) {
  m <- as.matrix(data[, items])
  row_means <- rowMeans(m, na.rm = TRUE); row_means[is.nan(row_means)] <- NA_real_
  for (j in seq_along(items)) {                       # within-person mean
    miss <- is.na(m[, j]) & !is.na(row_means)
    m[miss, j] <- row_means[miss]
  }
  all_missing <- is.na(row_means)
  for (j in seq_along(items)) {                       # group median fallback
    for (g in c(0, 1)) {
      idx <- all_missing & data$af == g
      if (any(idx)) m[idx, j] <- median(m[data$af == g, j], na.rm = TRUE)
    }
  }
  data[, items] <- m
  data
}
df <- df %>% impute_items(dep_items) %>% impute_items(wb_items) %>%
  mutate(depressive_symptoms = rowMeans(across(all_of(dep_items))),
         subjective_wellbeing = rowMeans(across(all_of(wb_items))))

# ---- models ----------------------------------------------------------------
fit_group_model <- function(data, outcome, covariates) {
  d <- data %>% select(score = all_of(outcome), group, all_of(covariates)) %>% na.omit()
  fit <- lm(score ~ group + ., data = d)
  an  <- car::Anova(fit, type = 3)
  Fv  <- an["group", "F value"]; p <- an["group", "Pr(>F)"]; df2 <- an["Residuals", "Df"]
  ctr <- summary(emmeans::contrast(emmeans::emmeans(fit, "group"), method = "revpairwise"))
  beta <- ctr$estimate[1]; dval <- beta / sigma(fit)
  n1 <- sum(d$group == "AF"); n2 <- sum(d$group == "control")
  se_d <- sqrt((n1 + n2) / (n1 * n2) + dval^2 / (2 * (n1 + n2)))
  tibble(n_AF = n1, n_control = n2,
         AF_mean = mean(d$score[d$group == "AF"]), AF_sd = sd(d$score[d$group == "AF"]),
         Control_mean = mean(d$score[d$group == "control"]), Control_sd = sd(d$score[d$group == "control"]),
         beta_adj = beta, cohens_d = dval,
         d_ci_low = dval - 1.96 * se_d, d_ci_high = dval + 1.96 * se_d,
         t = sign(beta) * sqrt(Fv), df_resid = df2, P = p)
}

outcomes <- c("Depressive symptoms (0-3 scale)" = "depressive_symptoms",
              "Subjective wellbeing (1-6 scale)" = "subjective_wellbeing")
res <- map_dfr(names(outcomes), function(lab)
  bind_cols(tibble(Outcome = lab), fit_group_model(df, outcomes[[lab]], core_covariates))) %>%
  mutate(P_FDR = p.adjust(P, method = "BH"))

write_csv(res, file.path(out_dir, "mental_health_results.csv"))
print(res %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
        select(Outcome, n_AF, n_control, cohens_d, d_ci_low, d_ci_high, t, df_resid, P, P_FDR),
      n = Inf, width = Inf)
message("Done. Outputs in: ", out_dir)
