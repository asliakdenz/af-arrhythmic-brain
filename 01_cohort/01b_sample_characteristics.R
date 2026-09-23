#!/usr/bin/env Rscript
# ==============================================================================
# 01b_sample_characteristics.R
#
# Sample characteristics of the matched sample: Table 1 (matching covariates
# with effect sizes and FDR-corrected P-values), Table 2 (cardiac and other
# comorbidities by group) and Supplementary Fig. 2 (propensity score matching
# balance).
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Participants and matching; Results, "Sample characteristics and
#          matching"
#
# Table 1: continuous variables as mean +/- SD with Cohen's d
#          (AF - controls) / pooled SD and Welch t-test; categorical variables
#          as n (%) with the standardized difference in proportions and a
#          chi-squared test; Benjamini-Hochberg FDR across all matching
#          covariates.
# Table 2: ICD-10 comorbidity counts (hospital-episode diagnoses) with
#          Fisher's exact test, or the chi-squared test where all expected
#          counts are at least 5.
# Supplementary Fig. 2: standardized mean differences after matching for the
#          matching covariates (circles) and the WMH measures (triangles), with
#          the |SMD| < 0.10 balance threshold.
#
# Inputs:  DERIV_DIR/cohort_matched.csv, DERIV_DIR/hesin_diag_cohort.csv
# Outputs: OUTPUT_DIR/sample/table1_matching_covariates.csv
#          OUTPUT_DIR/sample/table2_comorbidities.csv
#          OUTPUT_DIR/sample/supp_fig2_balance.{svg,pdf}
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(tibble); library(ggplot2)
})

deriv_dir  <- Sys.getenv("DERIV_DIR",  "derivatives")
out_root   <- Sys.getenv("OUTPUT_DIR", "output")
cohort_csv <- file.path(deriv_dir, "cohort_matched.csv")
hesin_csv  <- file.path(deriv_dir, "hesin_diag_cohort.csv")
out_dir    <- file.path(out_root, "sample")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(cohort_csv))

df <- read_csv(cohort_csv, show_col_types = FALSE, col_types = cols(eid = col_character()))
af <- df$af == 1; ct <- df$af == 0
message(sprintf("Matched sample: n = %d (AF %d, controls %d)", nrow(df), sum(af), sum(ct)))

pooled_sd <- function(x, g1, g2) sqrt((var(x[g1]) * (sum(g1) - 1) + var(x[g2]) * (sum(g2) - 1)) /
                                        (sum(g1) + sum(g2) - 2))
smd_continuous <- function(x) (mean(x[af]) - mean(x[ct])) / pooled_sd(x, af, ct)
smd_binary     <- function(x) { p1 <- mean(x[af]); p0 <- mean(x[ct])
                                (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2) }

# ---- Table 1 -----------------------------------------------------------------
continuous <- c("Age, years" = "age_i2", "Education (ISCED)" = "education_i2",
                "BMI, kg/m2" = "bmi_i2", "Systolic BP, mmHg" = "sys_mean",
                "Diastolic BP, mmHg" = "dia_mean", "Total cholesterol, mmol/L" = "cholesterol",
                "Blood glucose, mmol/L" = "glucose", "WMH load, mL/L" = "norm_wmh_ml")
binary     <- c("Female, n (%)" = "sex", "Current smokers, n (%)" = "smo")
# UK Biobank field 31: 0 = female, 1 = male
df$female <- as.integer(df$sex == 0)
binary[["Female, n (%)"]] <- "female"

fmt_cont <- function(x) sprintf("%.2f ± %.2f", mean(x), sd(x))
fmt_bin  <- function(x) sprintf("%d (%.1f%%)", sum(x), 100 * mean(x))

table1 <- bind_rows(
  map_dfr(names(continuous), function(lab) {
    x <- df[[continuous[[lab]]]]
    tibble(Measure = lab, AF = fmt_cont(x[af]), Controls = fmt_cont(x[ct]),
           effect_size = smd_continuous(x), statistic = "Cohen's d",
           P = t.test(x[af], x[ct])$p.value)
  }),
  map_dfr(names(binary), function(lab) {
    x <- df[[binary[[lab]]]]
    tibble(Measure = lab, AF = fmt_bin(x[af]), Controls = fmt_bin(x[ct]),
           effect_size = smd_binary(x), statistic = "SMD",
           P = suppressWarnings(chisq.test(table(df$af, x))$p.value))
  })
) %>% mutate(P_FDR = p.adjust(P, method = "BH"))
write_csv(table1, file.path(out_dir, "table1_matching_covariates.csv"))
cat("\n--- Table 1 ---\n"); print(table1, n = Inf, width = Inf)

# ---- Table 2 -----------------------------------------------------------------
if (file.exists(hesin_csv)) {
  hes <- read_csv(hesin_csv, show_col_types = FALSE, col_types = cols(.default = col_character())) %>%
    mutate(eid = as.character(eid), code = toupper(gsub("\\.", "", trimws(diag_icd10)))) %>%
    filter(eid %in% df$eid)
  comorbidities <- list(
    "Atherosclerotic heart disease (I25.1)" = "I251",
    "Congestive heart failure (I50.0)"      = "I500",
    "Alcohol abuse (F10.1/F10.2)"           = c("F101", "F102"),
    "Sleep apnoea (G47.3)"                  = "G473",
    "Hyperthyroidism (E05.0/.1/.2)"         = c("E050", "E051", "E052")
  )
  table2 <- map_dfr(names(comorbidities), function(lab) {
    pref <- comorbidities[[lab]]
    has  <- df$eid %in% unique(hes$eid[Reduce(`|`, lapply(pref, function(p) startsWith(hes$code, p)))])
    tab  <- table(factor(df$af, levels = c(1, 0)), factor(has, levels = c(TRUE, FALSE)))
    expected_ok <- all(suppressWarnings(chisq.test(tab)$expected) >= 5)
    p <- if (expected_ok) suppressWarnings(chisq.test(tab)$p.value) else fisher.test(tab)$p.value
    tibble(Diagnosis = lab, AF = fmt_bin(has[af]), Controls = fmt_bin(has[ct]),
           test = if (expected_ok) "chi-squared" else "Fisher's exact", P = p)
  })
  write_csv(table2, file.path(out_dir, "table2_comorbidities.csv"))
  cat("\n--- Table 2 ---\n"); print(table2, n = Inf, width = Inf)
} else {
  message("hesin_diag_cohort.csv not found; Table 2 skipped.")
}

# ---- Supplementary Fig. 2: balance after matching -----------------------------
bal <- bind_rows(
  map_dfr(names(continuous)[continuous != "norm_wmh_ml"], function(lab)
    tibble(Covariate = lab, SMD = smd_continuous(df[[continuous[[lab]]]]), type = "Matching covariate")),
  map_dfr(names(binary), function(lab)
    tibble(Covariate = lab, SMD = smd_binary(df[[binary[[lab]]]]), type = "Matching covariate")),
  tibble(Covariate = "Total WMH (TIV-normalized; matching variable)", SMD = smd_continuous(df$norm_wmh_ml),   type = "WMH measure"),
  tibble(Covariate = "Deep WMH",                                        SMD = smd_continuous(df$norm_deep_wmh_ml), type = "WMH measure"),
  tibble(Covariate = "Periventricular WMH",                             SMD = smd_continuous(df$norm_peri_wmh_ml), type = "WMH measure"),
  tibble(Covariate = "Total WMH (unnormalized)",                        SMD = smd_continuous(df$total_wmh_i2),     type = "WMH measure")
) %>% mutate(Covariate = factor(Covariate, levels = rev(Covariate)))
write_csv(bal, file.path(out_dir, "supp_fig2_balance.csv"))

p <- ggplot(bal, aes(SMD, Covariate, shape = type)) +
  geom_vline(xintercept = c(-0.10, 0.10), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(size = 2.8, colour = "#3C5488") +
  scale_shape_manual(values = c("Matching covariate" = 16, "WMH measure" = 17)) +
  labs(x = "Standardized mean difference (AF - controls) after matching", y = NULL, shape = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "top")
for (ext in c("svg", "pdf"))
  ggsave(file.path(out_dir, paste0("supp_fig2_balance.", ext)), p, width = 6.5, height = 4)
message("Done. Outputs in: ", out_dir)
