#!/usr/bin/env Rscript
# ==============================================================================
# 07a_cognitive_outcomes.R
#
# Primary cognitive outcome analysis: group differences across the four
# cognitive domains (Table 3) and the ten individual tests (Supplementary
# Table 3) in the matched sample.
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Cognitive and mental health outcomes -> Cognitive outcomes;
#          Statistical analysis
#
# Ten UK Biobank tasks spanning four domains. Time-based measures (Trail
# Making A and B, Reaction Time) are sign-reversed so that higher values
# indicate better performance; Symbol Digit Substitution and Matrix Pattern
# Completion are expressed as accuracy ratios. All variables are z-scored
# within the analytic (matched) sample and domain scores are the unweighted
# mean of the constituent tests.
#
# Group differences are evaluated by linear regression with Type III sums of
# squares, with AF status as predictor and the core covariates: age, sex,
# education, diastolic blood pressure and log-transformed total WMH. The
# additional matching variables are not entered as regression covariates, and
# mean framewise displacement is not included because no imaging predictor
# enters these models. Effect sizes are covariate-adjusted Cohen's d from
# estimated marginal means, computed as (AF - controls) / residual SD, with
# large-sample confidence intervals based on group sizes and Benjamini-Hochberg
# FDR correction within the family (four domains; ten tests).
#
# Input:   DERIV_DIR/cohort_matched.csv                (01_cohort_assembly.R)
# Outputs: OUTPUT_DIR/cognition/domain_results.csv     (Table 3)
#          OUTPUT_DIR/cognition/test_results.csv       (Supplementary Table 3)
#          OUTPUT_DIR/cognition/domains_barplot.svg
#          DERIV_DIR/cohort_matched_cognition.csv      z-scored tests and
#                                                      domain scores, used by
#                                                      07c and 07d
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(purrr)
  library(car); library(emmeans); library(ggplot2)
})

# ---- configuration ---------------------------------------------------------
deriv_dir <- Sys.getenv("DERIV_DIR",  "derivatives")
out_root  <- Sys.getenv("OUTPUT_DIR", "output")
in_csv    <- file.path(deriv_dir, "cohort_matched.csv")
out_dir   <- file.path(out_root, "cognition")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(in_csv))

core_covariates <- c("age_i2", "sex", "education_i2", "dia_mean", "log_norm_wmh_ml")

test_groups <- list(
  "Attention / executive function" = c("tower_rearr_i2", "trail_mak_2_i2"),
  "Processing speed"               = c("reaction_time_i2", "sym_dig_accuracy", "trail_mak_1_i2"),
  "Memory"                         = c("num_mem_nmbr_i2", "paired_assc_learn_i2", "pros_mem_cog_summ_i2"),
  "Reasoning"                      = c("fluid_intl_i2", "mat_pat_compl_accuracy")
)
domain_keys <- c("Attention / executive function" = "domain_executive_std",
                 "Processing speed"               = "domain_processing_std",
                 "Memory"                         = "domain_memory_std",
                 "Reasoning"                      = "domain_reasoning_std")
test_labels <- c(
  tower_rearr_i2         = "Tower Rearranging",
  trail_mak_2_i2         = "Trail Making Test B",
  reaction_time_i2       = "Reaction Time",
  sym_dig_accuracy       = "Symbol Digit Substitution",
  trail_mak_1_i2         = "Trail Making Test A",
  num_mem_nmbr_i2        = "Numeric Memory",
  paired_assc_learn_i2   = "Paired Associate Learning",
  pros_mem_cog_summ_i2   = "Prospective Memory",
  fluid_intl_i2          = "Fluid Intelligence",
  mat_pat_compl_accuracy = "Matrix Pattern Completion"
)
all_tests <- unique(unlist(test_groups))

# ---- data ------------------------------------------------------------------
df <- read_csv(in_csv, show_col_types = FALSE, col_types = cols(eid = col_character()))

if (!"sym_dig_accuracy" %in% names(df))
  df$sym_dig_accuracy <- ifelse(df$sym_dig_attempt_i2 > 0, df$sym_dig_corr_i2 / df$sym_dig_attempt_i2, NA_real_)
if (!"mat_pat_compl_accuracy" %in% names(df))
  df$mat_pat_compl_accuracy <- ifelse(df$mat_pat_compl_viewed_i2 > 0, df$mat_pat_compl_corr_i2 / df$mat_pat_compl_viewed_i2, NA_real_)

df <- df %>%
  drop_na(all_of(c(all_tests, core_covariates))) %>%
  mutate(across(c(trail_mak_1_i2, trail_mak_2_i2, reaction_time_i2), ~ -.x),   # higher = better
         group = factor(af, levels = c(0, 1), labels = c("control", "AF")))

zscore <- function(x) (x - mean(x)) / sd(x)
for (v in all_tests) df[[v]] <- zscore(df[[v]])          # within the analytic sample
for (dn in names(test_groups)) df[[domain_keys[[dn]]]] <- rowMeans(df[, test_groups[[dn]]])

message(sprintf("Matched analytic sample: n = %d (AF %d, controls %d)",
                nrow(df), sum(df$af == 1), sum(df$af == 0)))

write_csv(df %>% select(eid, af, all_of(all_tests), all_of(unname(domain_keys))),
          file.path(deriv_dir, "cohort_matched_cognition.csv"))

# ---- model -----------------------------------------------------------------
# Linear regression, Type III sums of squares; adjusted Cohen's d from the
# estimated marginal means (AF - control) divided by the residual SD;
# t = sign(beta) * sqrt(F) on F(1, df_resid).
fit_group_model <- function(data, outcome, covariates) {
  d <- data %>% select(score = all_of(outcome), group, all_of(covariates))
  fit <- lm(score ~ group + ., data = d)
  an  <- car::Anova(fit, type = 3)
  Fv  <- an["group", "F value"]; p <- an["group", "Pr(>F)"]; df2 <- an["Residuals", "Df"]
  em  <- emmeans::emmeans(fit, specs = "group")
  ctr <- summary(emmeans::contrast(em, method = "revpairwise"))     # AF - control
  beta <- ctr$estimate[1]
  dval <- beta / sigma(fit)
  n1 <- sum(d$group == "AF"); n2 <- sum(d$group == "control")
  se_d <- sqrt((n1 + n2) / (n1 * n2) + dval^2 / (2 * (n1 + n2)))
  tibble(
    AF_mean = mean(d$score[d$group == "AF"]),      AF_sd = sd(d$score[d$group == "AF"]),
    Control_mean = mean(d$score[d$group == "control"]), Control_sd = sd(d$score[d$group == "control"]),
    n_AF = n1, n_control = n2,
    beta_adj = beta, cohens_d = dval, d_ci_low = dval - 1.96 * se_d, d_ci_high = dval + 1.96 * se_d,
    t = sign(beta) * sqrt(Fv), df_resid = df2, F = Fv, P = p
  )
}

domain_tbl <- map_dfr(names(test_groups), function(dn)
  bind_cols(tibble(Domain = dn), fit_group_model(df, domain_keys[[dn]], core_covariates))) %>%
  mutate(P_FDR = p.adjust(P, method = "BH"))                      # across the four domains

test_tbl <- map_dfr(names(test_groups), function(dn)
  map_dfr(test_groups[[dn]], function(tv)
    bind_cols(tibble(Domain = dn, Test = test_labels[[tv]]),
              fit_group_model(df, tv, core_covariates)))) %>%
  mutate(P_FDR = p.adjust(P, method = "BH"))                      # across the ten tests

write_csv(domain_tbl, file.path(out_dir, "domain_results.csv"))
write_csv(test_tbl,   file.path(out_dir, "test_results.csv"))

cat("\n--- Table 3: cognitive domains (matched sample) ---\n")
print(domain_tbl %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
        select(Domain, n_AF, n_control, cohens_d, d_ci_low, d_ci_high, t, df_resid, P, P_FDR),
      n = Inf, width = Inf)
cat("\n--- Supplementary Table 3: individual tests ---\n")
print(test_tbl %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
        select(Domain, Test, cohens_d, P, P_FDR), n = Inf, width = Inf)

# ---- figure ----------------------------------------------------------------
domain_long <- df %>%
  select(group, all_of(unname(domain_keys))) %>%
  pivot_longer(-group, names_to = "key", values_to = "score") %>%
  mutate(Domain = factor(names(domain_keys)[match(key, domain_keys)], levels = names(test_groups)))

domain_sum <- domain_long %>%
  group_by(Domain, group) %>%
  summarise(mean = mean(score), se = sd(score) / sqrt(n()), .groups = "drop")

anno <- domain_tbl %>%
  transmute(Domain = factor(Domain, levels = names(test_groups)),
            label = sprintf("d = %+.2f, P[FDR] = %s", cohens_d,
                            ifelse(P_FDR < 0.001, "< 0.001", sprintf("%.3f", P_FDR))))

p <- ggplot(domain_sum, aes(Domain, mean, colour = group)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15,
                position = position_dodge(0.4)) +
  geom_point(size = 2.5, position = position_dodge(0.4)) +
  geom_text(data = anno, aes(x = Domain, y = 0.22, label = label),
            inherit.aes = FALSE, size = 2.8) +
  scale_colour_manual(values = c(control = "#6E7B8B", AF = "#9FD3D3")) +
  labs(x = NULL, y = "Domain score (z, mean ± SE)", colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", panel.grid.major.x = element_blank())
ggsave(file.path(out_dir, "domains_barplot.svg"), p, width = 7, height = 3.6)

message("Done. Outputs in: ", out_dir)
