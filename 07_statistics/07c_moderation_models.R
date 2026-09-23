#!/usr/bin/env Rscript
# ==============================================================================
# 07c_moderation_models.R
#
# Brain-behaviour analyses in the matched sample: moderation models, simple
# slopes, hierarchical regression, specificity of the moderation, the
# exploratory structure-function interaction within AF, the white matter
# hyperintensity location sensitivity of all of these, and Figure 4.
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Brain-behaviour analyses -> Moderation models; Hierarchical
#          regression and nested cross-validation (in-sample part);
#          Specificity of the moderation; White matter hyperintensity location
#          sensitivity. Results: Figure 4; Supplementary Tables 11, 12 and 14.
#
# Models (OLS with HC3 heteroscedasticity-consistent standard errors):
#
#   Moderation        processing speed ~ predictor x AF + core covariates
#                     + mean framewise displacement, for the functional
#                     deviation score, whole-brain functional system
#                     segregation and the whole-brain participation
#                     coefficient (the last two as alternative moderators,
#                     testing specificity). Simple slopes per group from the
#                     same model.
#
#   Hierarchical      baseline (core covariates + mean FD) versus baseline +
#                     functional deviation score; F-test; delta R2 as the
#                     change in adjusted R2 (in-sample). Full sample under each
#                     WMH specification, and within AF patients and controls
#                     with total WMH. The leakage-controlled out-of-fold
#                     estimates are produced by 07d_nested_cross_validation.py.
#
#   Structure-function   within AF only: processing speed ~ functional
#                     deviation x structural system segregation (and x
#                     structural participation) + covariates; predictors
#                     z-scored.
#
#   WMH specifications   log-transformed total WMH (primary), periventricular,
#                     deep, and periventricular + deep entered simultaneously.
#
# Inputs (DERIV_DIR):
#   cohort_matched.csv                          01_cohort_assembly.R
#   cohort_matched_cognition.csv                07a_cognitive_outcomes.R
#   fc_deviation.csv                            06a_fc_deviation_score.py
#   segregation_functional.csv                  05a_segregation_functional.py
#   participation_functional.csv                05b_participation_functional.py
#   segregation_participation_structural.csv    05c_segregation_participation_structural.py
#
# Outputs (OUTPUT_DIR/brain_behaviour): moderation_results.csv,
#   simple_slopes.csv, hierarchical_results.csv, structure_function_af.csv,
#   Figure4.{svg,pdf,png}
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(purrr); library(tibble)
  library(sandwich); library(lmtest); library(ggplot2)
})

deriv_dir <- Sys.getenv("DERIV_DIR",  "derivatives")
out_root  <- Sys.getenv("OUTPUT_DIR", "output")
out_dir   <- file.path(out_root, "brain_behaviour")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rd <- function(name, required = TRUE) {
  f <- file.path(deriv_dir, name)
  if (!file.exists(f)) { if (required) stop("missing input: ", f); return(NULL) }
  read_csv(f, show_col_types = FALSE, col_types = cols(eid = col_character()))
}
cohort <- rd("cohort_matched.csv")
cog    <- rd("cohort_matched_cognition.csv")
fcdev  <- rd("fc_deviation.csv")
segf   <- rd("segregation_functional.csv")
pcf    <- rd("participation_functional.csv")
str_m  <- rd("segregation_participation_structural.csv", required = FALSE)

outcome   <- "domain_processing_std"
base_covs <- c("age_std", "sex", "education_std", "dia_bp_std", "mean_fd_std")
wmh_specs <- list(
  "Total WMH (primary)" = "log_total_wmh_std",
  "Periventricular WMH" = "log_peri_wmh_std",
  "Deep WMH"            = "log_deep_wmh_std",
  "Peri + Deep"         = c("log_peri_wmh_std", "log_deep_wmh_std")
)

z <- function(x) as.numeric(scale(x))
df <- cohort %>%
  select(eid, af, all_of(base_covs), all_of(unique(unlist(wmh_specs)))) %>%
  inner_join(cog   %>% select(eid, all_of(outcome)),                   by = "eid") %>%
  inner_join(fcdev %>% select(eid, fc_deviation),                      by = "eid") %>%
  inner_join(segf  %>% select(eid, fc_segregation = sys_cort),         by = "eid") %>%
  inner_join(pcf   %>% select(eid, fc_participation = pc_mean_cort),   by = "eid")
if (!is.null(str_m))
  df <- df %>% left_join(str_m %>% select(eid, sc_segregation = sys_cort,
                                          sc_participation = pc_mean_cort), by = "eid")
df <- df %>%
  mutate(fc_deviation_z     = z(fc_deviation),
         fc_segregation_z   = z(fc_segregation),
         fc_participation_z = z(fc_participation))
message(sprintf("Analysis sample: n = %d (AF %d, controls %d)", nrow(df), sum(df$af == 1), sum(df$af == 0)))

hc3_row <- function(fit, term) {
  V  <- sandwich::vcovHC(fit, type = "HC3")
  ct <- lmtest::coeftest(fit, vcov = V)
  b <- ct[term, "Estimate"]; se <- ct[term, "Std. Error"]
  tibble(beta = b, SE = se, t = ct[term, "t value"], df_resid = fit$df.residual,
         ci_low = b - 1.96 * se, ci_high = b + 1.96 * se, P = ct[term, "Pr(>|t|)"])
}

# ==============================================================================
# 1. MODERATION MODELS AND SIMPLE SLOPES
# ==============================================================================
predictors <- c("FC deviation"                 = "fc_deviation",
                "FC deviation (per SD)"        = "fc_deviation_z",
                "Functional system segregation" = "fc_segregation_z",
                "Functional participation"     = "fc_participation_z")

moderation <- list(); slopes <- list()
for (spec in names(wmh_specs)) {
  covs <- c(base_covs, wmh_specs[[spec]])
  for (pn in names(predictors)) {
    pv  <- predictors[[pn]]
    d   <- df %>% select(all_of(c(outcome, pv, "af", covs))) %>% na.omit()
    fit <- lm(reformulate(c(sprintf("%s * af", pv), covs), response = outcome), data = d)
    int_term <- sprintf("%s:af", pv)
    moderation[[length(moderation) + 1]] <- bind_cols(
      tibble(WMH_specification = spec, Model = sprintf("%s x AF", pn), n = nrow(d)),
      hc3_row(fit, int_term))
    # simple slopes: controls = beta_pred; AF = beta_pred + beta_int (HC3 covariance)
    V <- sandwich::vcovHC(fit, type = "HC3"); b <- coef(fit)
    s_ctrl <- b[[pv]];             se_ctrl <- sqrt(V[pv, pv])
    s_af   <- b[[pv]] + b[[int_term]]
    se_af  <- sqrt(V[pv, pv] + V[int_term, int_term] + 2 * V[pv, int_term])
    slopes[[length(slopes) + 1]] <- tibble(
      WMH_specification = spec, Predictor = pn,
      Group = c("AF", "Controls"), beta = c(s_af, s_ctrl), SE = c(se_af, se_ctrl),
      ci_low = beta - 1.96 * SE, ci_high = beta + 1.96 * SE,
      t = beta / SE, P = 2 * pt(abs(t), df = fit$df.residual, lower.tail = FALSE))
  }
}
moderation <- bind_rows(moderation); slopes <- bind_rows(slopes)
write_csv(moderation, file.path(out_dir, "moderation_results.csv"))
write_csv(slopes,     file.path(out_dir, "simple_slopes.csv"))

cat("\n--- Moderation (interaction terms, HC3; Supplementary Table 11) ---\n")
print(moderation %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = Inf, width = Inf)
cat("\n--- Simple slopes, primary specification ---\n")
print(slopes %>% filter(WMH_specification == "Total WMH (primary)", Predictor == "FC deviation") %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))), width = Inf)

# ==============================================================================
# 2. HIERARCHICAL REGRESSION (in-sample; Supplementary Table 12)
# ==============================================================================
hierarchical <- function(d, covs, label) {
  d <- d %>% select(all_of(c(outcome, "fc_deviation", covs))) %>% na.omit()
  base <- lm(reformulate(covs, response = outcome), data = d)
  full <- lm(reformulate(c(covs, "fc_deviation"), response = outcome), data = d)
  an   <- anova(base, full)
  tibble(Model = label, n = nrow(d),
         baseline_adj_R2 = summary(base)$adj.r.squared,
         full_adj_R2     = summary(full)$adj.r.squared,
         delta_adj_R2    = summary(full)$adj.r.squared - summary(base)$adj.r.squared,
         F = an$F[2], df1 = an$Df[2], df2 = an$Res.Df[2], P = an$`Pr(>F)`[2])
}
hier <- bind_rows(
  map_dfr(names(wmh_specs), function(spec)
    hierarchical(df, c(base_covs, wmh_specs[[spec]]), sprintf("Full matched sample: %s", spec))),
  hierarchical(df %>% filter(af == 1), c(base_covs, wmh_specs[[1]]), "AF patients (total WMH)"),
  hierarchical(df %>% filter(af == 0), c(base_covs, wmh_specs[[1]]), "Controls (total WMH)")
)
write_csv(hier, file.path(out_dir, "hierarchical_results.csv"))
cat("\n--- Hierarchical regression (Supplementary Table 12) ---\n")
print(hier %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = Inf, width = Inf)

# ==============================================================================
# 3. STRUCTURE-FUNCTION MODERATION WITHIN AF (Supplementary Table 14)
# ==============================================================================
if (!is.null(str_m)) {
  covs <- c(base_covs, wmh_specs[[1]])
  d_af <- df %>% filter(af == 1) %>%
    mutate(sc_segregation_z = z(sc_segregation), sc_participation_z = z(sc_participation))
  sf <- map_dfr(c("sc_segregation_z", "sc_participation_z"), function(sv) {
    d   <- d_af %>% select(all_of(c(outcome, "fc_deviation_z", sv, covs))) %>% na.omit()
    fit <- lm(reformulate(c(sprintf("fc_deviation_z * %s", sv), covs), response = outcome), data = d)
    map_dfr(c("fc_deviation_z", sv, sprintf("fc_deviation_z:%s", sv)), function(term)
      bind_cols(tibble(Model = sprintf("FC deviation x %s", sv), Term = term, n = nrow(d)),
                hc3_row(fit, term)))
  })
  write_csv(sf, file.path(out_dir, "structure_function_af.csv"))
  cat("\n--- Structure-function moderation within AF (Supplementary Table 14) ---\n")
  print(sf %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = Inf, width = Inf)
} else {
  message("segregation_participation_structural.csv not found; Supplementary Table 14 skipped.")
}

# ==============================================================================
# 4. FIGURE 4: partial regression plot of the FC deviation x AF interaction
# ==============================================================================
covs <- c(base_covs, wmh_specs[[1]])
d4 <- df %>% select(all_of(c(outcome, "fc_deviation", "af", covs))) %>% na.omit()
d4$y_adj <- resid(lm(reformulate(covs, response = outcome), data = d4))
d4$x_adj <- resid(lm(reformulate(covs, response = "fc_deviation"), data = d4))
d4$Group <- factor(d4$af, levels = c(0, 1), labels = c("Matched controls", "AF patients"))
int_row  <- moderation %>% filter(WMH_specification == "Total WMH (primary)", Model == "FC deviation x AF")

p4 <- ggplot(d4, aes(x_adj, y_adj, colour = Group, fill = Group)) +
  geom_point(alpha = 0.45, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.9, alpha = 0.2) +
  scale_colour_manual(values = c("Matched controls" = "#6E7B8B", "AF patients" = "#2A9D8F")) +
  scale_fill_manual(values   = c("Matched controls" = "#6E7B8B", "AF patients" = "#2A9D8F")) +
  labs(x = "Functional deviation score (covariate-adjusted)",
       y = "Processing speed, z (covariate-adjusted)",
       subtitle = sprintf("FC deviation x AF: beta = %.3f, P = %.3f", int_row$beta, int_row$P),
       colour = NULL, fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
for (ext in c("svg", "pdf", "png"))
  ggsave(file.path(out_dir, paste0("Figure4.", ext)), p4, width = 89 / 25.4, height = 76 / 25.4,
         dpi = 600, bg = "white")
message("Done. Outputs in: ", out_dir)
