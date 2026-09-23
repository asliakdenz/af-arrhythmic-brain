#!/usr/bin/env Rscript
# ==============================================================================
# 07e_sensitivity_unmatched.R
#
# Sensitivity analysis: AF effect on the cognitive domains in the unmatched
# pool of all eligible participants with complete cognitive testing
# (Supplementary Fig. 3).
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Cognitive and mental health outcomes -> Cognitive outcomes
#          ("To assess generalizability beyond the matched cohort, analyses
#          were repeated without matching in all eligible participants with
#          complete cognitive testing. This pool was defined by cognitive-data
#          completeness and was not restricted by connectome quality control.")
#
# Same five-covariate adjustment as the primary analysis (age, sex, education,
# diastolic blood pressure, log-transformed total WMH). Cohen's d is the AF
# coefficient divided by the residual SD; 95% CIs on the d scale from the
# coefficient standard error. FDR across the four domains.
#
# Input:  DERIV_DIR/cohort_unmatched_cogcomplete.csv   (01_cohort_assembly.R)
# Output: OUTPUT_DIR/cognition/sensitivity_unmatched.csv,
#         OUTPUT_DIR/cognition/supp_sensitivity_forest.{svg,pdf,png}
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(purrr); library(ggplot2); library(forcats)
})

deriv_dir <- Sys.getenv("DERIV_DIR",  "derivatives")
out_root  <- Sys.getenv("OUTPUT_DIR", "output")
in_csv    <- file.path(deriv_dir, "cohort_unmatched_cogcomplete.csv")
out_dir   <- file.path(out_root, "cognition")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(in_csv))

covariates <- c("age_i2", "sex", "education_i2", "dia_mean", "log_norm_wmh_ml")
domains <- list(
  "Attention / executive function" = c("tower_rearr_i2", "trail_mak_2_i2"),
  "Processing speed"               = c("reaction_time_i2", "sym_dig_accuracy", "trail_mak_1_i2"),
  "Memory"                         = c("num_mem_nmbr_i2", "paired_assc_learn_i2", "pros_mem_cog_summ_i2"),
  "Reasoning"                      = c("fluid_intl_i2", "mat_pat_compl_accuracy")
)
all_tests <- unique(unlist(domains))

df <- read_csv(in_csv, show_col_types = FALSE, col_types = cols(eid = col_character())) %>%
  drop_na(all_of(c(all_tests, covariates, "af"))) %>%
  mutate(across(c(trail_mak_1_i2, trail_mak_2_i2, reaction_time_i2), ~ -.x))
zscore <- function(x) (x - mean(x)) / sd(x)
for (v in all_tests) df[[v]] <- zscore(df[[v]])                     # within the pool
for (dn in names(domains)) df[[dn]] <- rowMeans(df[, domains[[dn]]])

n_af <- sum(df$af == 1); n_ct <- sum(df$af == 0)
message(sprintf("Unmatched pool: n = %d (AF %d, controls %d)", nrow(df), n_af, n_ct))

res <- map_dfr(names(domains), function(dn) {
  d   <- df %>% select(score = all_of(dn), af, all_of(covariates))
  fit <- lm(score ~ af + ., data = d)
  co  <- summary(fit)$coefficients["af", ]
  s   <- sigma(fit)
  tibble(Domain = dn, n_AF = n_af, n_control = n_ct,
         beta = co[["Estimate"]], SE = co[["Std. Error"]],
         cohens_d = co[["Estimate"]] / s,
         d_ci_low  = (co[["Estimate"]] - 1.96 * co[["Std. Error"]]) / s,
         d_ci_high = (co[["Estimate"]] + 1.96 * co[["Std. Error"]]) / s,
         t = co[["t value"]], df_resid = fit$df.residual, P = co[["Pr(>|t|)"]])
}) %>% mutate(P_FDR = p.adjust(P, method = "BH"))

write_csv(res, file.path(out_dir, "sensitivity_unmatched.csv"))
print(res %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = Inf, width = Inf)

# ---- forest plot (Supplementary Fig. 3) ------------------------------------
plot_df <- res %>%
  mutate(Domain = fct_rev(factor(Domain, levels = names(domains))),
         label  = ifelse(P < 0.001, "P < 0.001", sprintf("P = %.3f", P)))
p <- ggplot(plot_df, aes(cohens_d, Domain)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = d_ci_low, xmax = d_ci_high), height = 0.18, linewidth = 0.6, colour = "#3C5488") +
  geom_point(size = 3, colour = "#3C5488") +
  geom_text(aes(label = label), vjust = -1.1, size = 3.2, colour = "grey25") +
  labs(x = "AF effect on cognitive domain (Cohen's d, 95% CI)", y = NULL,
       subtitle = sprintf("Unmatched pool: n = %d (%d AF, %d controls); same five-covariate adjustment",
                          nrow(df), n_af, n_ct)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
for (ext in c("svg", "pdf", "png"))
  ggsave(file.path(out_dir, paste0("supp_sensitivity_forest.", ext)), p, width = 7.5, height = 4, dpi = 300)
message("Done. Outputs in: ", out_dir)
