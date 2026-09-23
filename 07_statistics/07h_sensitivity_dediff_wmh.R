#!/usr/bin/env Rscript
# ==============================================================================
# 07h_sensitivity_dediff_wmh.R
#
# White matter hyperintensity location sensitivity for network segregation
# and integration (Supplementary Table 10).
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Network segregation and integration; Brain-behaviour analyses ->
#          White matter hyperintensity location sensitivity
#
# For each of the four measures (functional and structural system
# segregation and participation coefficient), the covariate-adjusted group
# model of 08a_figure3_dedifferentiation.R is refitted with the log-transformed
# total WMH covariate replaced by (a) deep WMH, (b) periventricular WMH, or
# (c) both entered simultaneously:
#
#   z(metric) ~ group + age + sex + education + diastolic BP
#               + mean framewise displacement (functional only) + <WMH>
#
# OLS with HC3 robust standard errors; the metric is z-scored so the
# coefficient is the AF effect in SD units; Benjamini-Hochberg FDR across the
# two measures within each modality and WMH specification.
#
# Inputs (DERIV_DIR): cohort_matched.csv, segregation_functional.csv (05a),
#   participation_functional.csv (05b), segregation_participation_structural.csv (05c)
# Output: OUTPUT_DIR/dedifferentiation/segregation_participation_by_wmh.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(tibble); library(lmtest); library(sandwich)
})

deriv_dir <- Sys.getenv("DERIV_DIR",  "derivatives")
out_dir   <- file.path(Sys.getenv("OUTPUT_DIR", "output"), "dedifferentiation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cohort_csv <- file.path(deriv_dir, "cohort_matched.csv")
seg_csv    <- file.path(deriv_dir, "segregation_functional.csv")
pc_csv     <- file.path(deriv_dir, "participation_functional.csv")
str_csv    <- file.path(deriv_dir, "segregation_participation_structural.csv")
stopifnot(file.exists(cohort_csv), file.exists(seg_csv), file.exists(pc_csv))

base_covariates <- c("age_std", "sex", "education_std", "dia_bp_std", "mean_fd_std")
fd_covariate    <- "mean_fd_std"
wmh_specs <- list("Deep WMH" = "log_deep_wmh_std", "Periventricular WMH" = "log_peri_wmh_std",
                  "Peri + Deep" = c("log_peri_wmh_std", "log_deep_wmh_std"))

measures <- list(
  list(file = seg_csv, col = "sys_cort",     modality = "Functional", measure = "System segregation"),
  list(file = pc_csv,  col = "pc_mean_cort", modality = "Functional", measure = "Participation coefficient"))
if (file.exists(str_csv)) {
  measures <- c(measures, list(
    list(file = str_csv, col = "sys_cort",     modality = "Structural", measure = "System segregation"),
    list(file = str_csv, col = "pc_mean_cort", modality = "Structural", measure = "Participation coefficient")))
} else {
  message("segregation_participation_structural.csv not found (run 05c); structural rows skipped.")
}

cov_frame <- read_csv(cohort_csv, show_col_types = FALSE, col_types = cols(eid = col_character())) %>%
  select(eid, af, all_of(c(base_covariates, unique(unlist(wmh_specs)))))

run_one <- function(m, covariates) {
  if (identical(m$modality, "Structural")) covariates <- setdiff(covariates, fd_covariate)
  d <- read_csv(m$file, show_col_types = FALSE, col_types = cols(eid = col_character())) %>%
    select(eid, y = all_of(m$col)) %>%
    inner_join(cov_frame, by = "eid") %>%
    select(y, af, all_of(covariates)) %>% na.omit()
  d$y <- as.numeric(scale(d$y))
  fit <- lm(reformulate(c("af", covariates), response = "y"), data = d)
  ct  <- coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
  tibble(Modality = m$modality, Measure = m$measure, N = nrow(d),
         N_AF = sum(d$af == 1), N_Control = sum(d$af == 0),
         Beta_SD = ct["af", "Estimate"], SE = ct["af", "Std. Error"], t = ct["af", "t value"],
         df_resid = fit$df.residual, Cohens_d = ct["af", "Estimate"] / sigma(fit), P = ct["af", "Pr(>|t|)"])
}

res <- map_dfr(names(wmh_specs), function(spec) {
  map_dfr(measures, ~ run_one(.x, c(base_covariates, wmh_specs[[spec]]))) %>%
    group_by(Modality) %>% mutate(P_FDR = p.adjust(P, method = "BH")) %>% ungroup() %>%
    mutate(WMH_specification = spec, .before = 1)
})
write_csv(res, file.path(out_dir, "segregation_participation_by_wmh.csv"))
print(res %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = Inf, width = Inf)
message("Done. Output in: ", out_dir)
