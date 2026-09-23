#!/usr/bin/env Rscript
# ==============================================================================
# 07f_sensitivity_cognition_wmh.R
#
# White matter hyperintensity location sensitivity for the cognitive domain
# results (Supplementary Table 2).
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Brain-behaviour analyses -> White matter hyperintensity location
#          sensitivity
#
# The cognitive domain models of 07a_cognitive_outcomes.R are repeated in the
# same matched sample with the log-transformed total WMH covariate replaced by
#   (A) log(1 + x)-transformed deep WMH,
#   (B) log(1 + x)-transformed periventricular WMH,
#   (C) both compartments entered simultaneously.
# All other covariates (age, sex, education, diastolic blood pressure) are held
# constant. FDR correction is applied across the four domains within each
# specification.
#
# Inputs:  DERIV_DIR/cohort_matched.csv, DERIV_DIR/cohort_matched_cognition.csv
# Output:  OUTPUT_DIR/cognition/domain_results_wmh_sensitivity.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(car); library(emmeans)
})

deriv_dir <- Sys.getenv("DERIV_DIR",  "derivatives")
out_root  <- Sys.getenv("OUTPUT_DIR", "output")
out_dir   <- file.path(out_root, "cognition")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cohort_csv <- file.path(deriv_dir, "cohort_matched.csv")
cog_csv    <- file.path(deriv_dir, "cohort_matched_cognition.csv")   # from 07a
stopifnot(file.exists(cohort_csv), file.exists(cog_csv))

base_covariates <- c("age_i2", "sex", "education_i2", "dia_mean")
wmh_specs <- list(
  "Model A: deep WMH"                   = "log_norm_deep_wmh_ml",
  "Model B: periventricular WMH"        = "log_norm_peri_wmh_ml",
  "Model C: periventricular + deep WMH" = c("log_norm_peri_wmh_ml", "log_norm_deep_wmh_ml")
)
domain_keys <- c("Processing speed"               = "domain_processing_std",
                 "Attention / executive function" = "domain_executive_std",
                 "Memory"                         = "domain_memory_std",
                 "Reasoning"                      = "domain_reasoning_std")

cohort <- read_csv(cohort_csv, show_col_types = FALSE, col_types = cols(eid = col_character()))
cog    <- read_csv(cog_csv,    show_col_types = FALSE, col_types = cols(eid = col_character()))
df <- cohort %>%
  select(eid, af, all_of(base_covariates), all_of(unique(unlist(wmh_specs)))) %>%
  inner_join(cog %>% select(eid, all_of(unname(domain_keys))), by = "eid") %>%
  mutate(group = factor(af, levels = c(0, 1), labels = c("control", "AF")))
message(sprintf("Matched sample: n = %d (AF %d, controls %d)", nrow(df), sum(df$af == 1), sum(df$af == 0)))

fit_group_model <- function(data, outcome, covariates) {
  d <- data %>% select(score = all_of(outcome), group, all_of(covariates)) %>% na.omit()
  fit <- lm(score ~ group + ., data = d)
  an  <- car::Anova(fit, type = 3)
  Fv  <- an["group", "F value"]; p <- an["group", "Pr(>F)"]; df2 <- an["Residuals", "Df"]
  ctr <- summary(emmeans::contrast(emmeans::emmeans(fit, "group"), method = "revpairwise"))
  beta <- ctr$estimate[1]; dval <- beta / sigma(fit)
  n1 <- sum(d$group == "AF"); n2 <- sum(d$group == "control")
  se_d <- sqrt((n1 + n2) / (n1 * n2) + dval^2 / (2 * (n1 + n2)))
  tibble(n = nrow(d), cohens_d = dval, d_ci_low = dval - 1.96 * se_d, d_ci_high = dval + 1.96 * se_d,
         t = sign(beta) * sqrt(Fv), df_resid = df2, P = p)
}

res <- map_dfr(names(wmh_specs), function(spec) {
  covs <- c(base_covariates, wmh_specs[[spec]])
  map_dfr(names(domain_keys), function(dn)
    bind_cols(tibble(WMH_specification = spec, Domain = dn),
              fit_group_model(df, domain_keys[[dn]], covs))) %>%
    mutate(P_FDR = p.adjust(P, method = "BH"))          # across the four domains, within spec
})

write_csv(res, file.path(out_dir, "domain_results_wmh_sensitivity.csv"))
print(res %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = Inf, width = Inf)
message("Done. Output in: ", out_dir)
