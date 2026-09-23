#!/usr/bin/env Rscript
# =============================================================================
# 08a_figure3_dedifferentiation.R
#
# Figure 3 and the covariate-adjusted group models for network segregation
# and integration.
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Network segregation and integration; Statistical analysis
# Figure:  Main-text Figure 3
#
# This script is the single point of origin for the adjusted effects of
# Figure 3 and the Results (functional and structural system segregation and
# participation coefficient; signed within- and between-network
# connectivity). The extraction scripts 05a-05c compute the metrics; the
# models below produce the reported betas.
#
# Every number shown comes from the models described in the Statistical
# analysis section:
#   y ~ group + age + sex + education + diastolic BP + log total WMH (+ mean FD)
#   OLS with HC3 heteroscedasticity-consistent standard errors.
#
# Mean framewise displacement is dropped for the structural models, matching
# the Methods (FD indexes functional acquisition motion).
#
# Panel A  points = covariate-adjusted group means in Fisher-z units
#                   NOTE: the two facets have INDEPENDENT y-axes; this must
#                   be stated in the figure legend in the manuscript file
# Panel B  points = covariate-adjusted group means, z-scored to a common scale
#
# Brain formatting decisions:
#   - Rendered at 185 mm (double column) so text keeps its true point size
#   - Panels stacked A over B; side-by-side gives each facet ~45 mm
#   - No embedded captions, titles or subtitles: legends live at the end of
#     the main text document, and interpretation belongs in the text
#   - Panel B uses points rather than bars, matching Panel A and avoiding
#     the bar-chart-for-continuous-data issue
#   - Asterisks use <= thresholds per Brain style
#
# Outputs: Figure3.png / .svg / .pdf + figure3_adjusted_stats.csv
# =============================================================================

suppressMessages({
  library(ggplot2); library(dplyr); library(patchwork)
  library(sandwich); library(lmtest)
})

# ---- configuration ---------------------------------------------------------
# Edit these to point at your own input and output locations. Defaults are
# relative to the repository root.
DERIV   <- Sys.getenv("DERIV_DIR",  "derivatives")   # extraction outputs
OUTDIR  <- Sys.getenv("OUTPUT_DIR", "output")

SEG_CSV <- file.path(DERIV, "segregation_functional.csv")   # from 05a
PC_CSV  <- file.path(DERIV, "participation_functional.csv") # from 05b
STR_CSV <- file.path(DERIV, "segregation_participation_structural.csv")  # from 05c
# If STR_CSV is absent (05c not yet run), panel B is drawn with the functional
# measures only and the structural rows are skipped.
MASTER  <- file.path(DERIV, "cohort_matched.csv")            # one row per subject (01)

OUTBASE <- file.path(OUTDIR, "Figure3")
OUT_CSV <- file.path(OUTDIR, "figure3_adjusted_stats.csv")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- columns ---------------------------------------------------------------
# Functional and structural extractions share a column naming scheme, so the
# same names serve both modalities.
ID_COL  <- "eid"
SEG_COL <- "sys_cort";        PC_COL  <- "pc_mean_cort"
ZW      <- "zw_cort_signed";  ZB      <- "zb_cort_signed"
STR_SEG <- "sys_cort";        STR_PC  <- "pc_mean_cort"
COVARS  <- c("age_std", "sex", "education_std",
             "log_total_wmh_std", "dia_bp_std", "mean_fd_std")

# ---- figure geometry and type sizes ----------------------------------------
FIG_W  <- 185 / 25.4          # 7.283 in = 185 mm, Brain double column
FIG_H  <- FIG_W * 0.78
BASE   <- 8                   # ggplot base_size, points
SZ_BIG <- 2.5                 # geom_text size ~ 7 pt
SZ_SML <- 2.1                 # geom_text size ~ 6 pt

col_ctrl <- "#6E7B8B"; col_af <- "#9FD3D3"
GRP_COL <- "af"      # binary group column, 1 = AF, 0 = control

# Brain style: *P <= 0.05, **P <= 0.01, ***P <= 0.001
star <- function(p) ifelse(p <= 1e-3, "***",
                           ifelse(p <= 0.01, "**",
                                  ifelse(p <= 0.05, "*", "n.s.")))

ma <- read.csv(MASTER); ma[[ID_COL]] <- as.character(ma[[ID_COL]])
ma_cov <- ma[, c(ID_COL, intersect(COVARS, names(ma)))]

# =============================================================================
# Core: fit one adjusted model, return effect + adjusted group means
#   standardise = TRUE  -> y z-scored, means on a common SD scale (Panel B)
#   standardise = FALSE -> y kept in native units, means in Fisher-z (Panel A)
# Adjusted means are model predictions with all covariates held at their
# sample means; their SEs are taken from the HC3 covariance matrix.
# =============================================================================
fit_adj <- function(file, col, drop_fd = FALSE, standardise = TRUE, label = col) {
  df <- read.csv(file); df[[ID_COL]] <- as.character(df[[ID_COL]])
  m  <- merge(df[, c(ID_COL, GRP_COL, col)], ma_cov, by = ID_COL)   # ma_cov holds covariates only
  names(m)[names(m) == GRP_COL] <- "grp"
  names(m)[names(m) == col] <- "y"
  m <- m[stats::complete.cases(m[, c("y", "grp", intersect(COVARS, names(m)))]), ]
  
  cv <- intersect(COVARS, names(m))
  if (drop_fd) cv <- setdiff(cv, "mean_fd_std")   # structural models: no head-motion covariate
  
  sd_y <- sd(m$y, na.rm = TRUE)
  if (standardise) m$y <- as.vector(scale(m$y))
  
  fit <- lm(reformulate(c("grp", cv), "y"), data = m)
  V   <- vcovHC(fit, type = "HC3")
  ct  <- coeftest(fit, vcov = V)["grp", ]
  
  # AF effect in SD units of the outcome, regardless of the plotting scale
  beta_sd <- if (standardise) unname(ct[1]) else unname(ct[1]) / sd_y
  se_sd   <- if (standardise) unname(ct[2]) else unname(ct[2]) / sd_y
  
  # covariate-adjusted group means at mean covariate values
  nd <- data.frame(grp = c(0, 1))
  for (v in cv) nd[[v]] <- mean(m[[v]], na.rm = TRUE)
  X   <- model.matrix(delete.response(terms(fit)), nd)
  est <- as.vector(X %*% coef(fit))
  sem <- sqrt(diag(X %*% V %*% t(X)))
  
  stats <- data.frame(
    label = label, beta_sd = beta_sd, se_sd = se_sd,
    t = unname(ct[3]), df = df.residual(fit), p = unname(ct[4]),
    n = nobs(fit), n_case = sum(m$grp == 1), n_ctrl = sum(m$grp == 0),
    fd_included = !drop_fd, stringsAsFactors = FALSE
  )
  means <- data.frame(
    label = label,
    group = factor(c("Controls", "AF"), levels = c("Controls", "AF")),
    m = est, se = sem, stringsAsFactors = FALSE
  )
  list(stats = stats, means = means)
}

# =============================================================================
# PANEL A -- within vs between-network signed connectivity (Fisher-z retained)
# =============================================================================
aW <- fit_adj(SEG_CSV, ZW, standardise = FALSE, label = "Within-network")
aB <- fit_adj(SEG_CSV, ZB, standardise = FALSE, label = "Between-network")

dA <- bind_rows(aW$means, aB$means)
dA$component <- factor(dA$label, levels = c("Within-network", "Between-network"))

sA <- bind_rows(aW$stats, aB$stats)
sA$component <- factor(sA$label, levels = c("Within-network", "Between-network"))
sA$star  <- star(sA$p)
sA$blab  <- sprintf("\u03b2 = %+.2f SD", sA$beta_sd)

rngA <- dA %>% group_by(component) %>%
  summarise(top = max(m + 1.96 * se), bot = min(m - 1.96 * se), .groups = "drop")
brkA <- left_join(sA[, c("component", "star", "blab")], rngA, by = "component")
brkA$y <- brkA$top + 0.18 * (brkA$top - brkA$bot + 1e-6)

zline <- data.frame(
  component = factor("Between-network",
                     levels = c("Within-network", "Between-network")), y = 0)

pA <- ggplot(dA, aes(group, m, color = group)) +
  geom_hline(data = zline, aes(yintercept = y), linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.08, linewidth = 0.5) +
  geom_point(size = 2.6) +
  geom_segment(data = brkA, aes(x = 1, xend = 2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.4, color = "grey30") +
  geom_text(data = brkA, aes(x = 1.5, y = y, label = star),
            inherit.aes = FALSE, vjust = -0.3, size = SZ_BIG, color = "grey15") +
  # beta label sits directly under its bracket, not on the axis floor
  geom_text(data = brkA, aes(x = 1.5, y = y, label = blab),
            inherit.aes = FALSE, vjust = 1.6, size = SZ_SML, color = "grey35") +
  facet_wrap(~component, scales = "free_y") +
  scale_color_manual(values = c(Controls = col_ctrl, AF = col_af)) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.18))) +
  labs(x = NULL, y = "Mean connectivity (Fisher-z)", color = NULL) +
  theme_minimal(base_size = BASE) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 7),
        axis.text.x = element_blank(),      # group names are in the legend
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7.5),
        strip.text = element_text(face = "bold", size = 8))

# =============================================================================
# PANEL B -- functional vs structural dissociation (common SD scale)
# =============================================================================
specs <- list(
  list(SEG_CSV, SEG_COL, FALSE, "Functional\nsegregation"),
  list(PC_CSV,  PC_COL,  FALSE, "Functional\nparticipation")
)
if (file.exists(STR_CSV)) {
  specs <- c(specs, list(
    list(STR_CSV, STR_SEG, TRUE, "Structural\nsegregation"),
    list(STR_CSV, STR_PC,  TRUE, "Structural\nparticipation")))
} else {
  message("STR_CSV not found (run 05c); drawing panel B with functional measures only.")
}
resB <- lapply(specs, function(sp)
  fit_adj(sp[[1]], sp[[2]], drop_fd = sp[[3]], standardise = TRUE, label = sp[[4]]))

dB <- bind_rows(lapply(resB, `[[`, "means"))
sB <- bind_rows(lapply(resB, `[[`, "stats"))

lv <- vapply(specs, `[[`, character(1), 4)
dB$facet <- factor(dB$label, levels = lv)
sB$facet <- factor(sB$label, levels = lv)
sB$star  <- star(sB$p)

# All significance brackets at one common height so they read as a row
gtopB   <- max(dB$m + dB$se)
sgB     <- sB[, c("facet", "star")]
sgB$xi  <- as.integer(sgB$facet)
sgB$y   <- gtopB + 0.10

# Points with SE bars rather than bars: Brain asks that continuously
# distributed data be shown as data points or box-and-whisker plots.
pB <- ggplot(dB, aes(facet, m, color = group)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se),
                position = position_dodge(0.5), width = 0.12, linewidth = 0.5) +
  geom_point(position = position_dodge(0.5), size = 2.6) +
  geom_segment(data = sgB, aes(x = xi - 0.17, xend = xi + 0.17, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.4, color = "grey30") +
  geom_text(data = sgB, aes(x = xi, y = y + 0.015, label = star),
            inherit.aes = FALSE, vjust = 0, size = SZ_BIG, color = "grey15") +
  scale_color_manual(values = c(Controls = col_ctrl, AF = col_af)) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.15))) +
  labs(x = NULL, y = "Metric value (SD units)", color = NULL) +
  theme_minimal(base_size = BASE) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7.5))

# =============================================================================
# ASSEMBLE -- stacked, no embedded caption; the legend lives in the manuscript
# Panel tags are 9 pt bold capitals per Brain figure requirements.
# =============================================================================
fig <- (pA / pB) +
  plot_layout(guides = "collect", heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top",
        plot.tag = element_text(size = 9, face = "bold"))

print(fig)

ggsave(paste0(OUTBASE, ".png"), fig, width = FIG_W, height = FIG_H,
       dpi = 600, bg = "white")
ggsave(paste0(OUTBASE, ".svg"), fig, width = FIG_W, height = FIG_H,
       bg = "white")
ggsave(paste0(OUTBASE, ".pdf"), fig, width = FIG_W, height = FIG_H,
       device = cairo_pdf, bg = "white")

# =============================================================================
# STATS OUT -- for the Results text and the Supplementary Table 10 primary rows
# NOTE: sA carries a `component` column and sB a `facet` column, so bind them
# on an explicit shared column set rather than names(sB).
# =============================================================================
keep <- c("label","beta_sd","se_sd","t","df","p","n","n_case","n_ctrl","fd_included")
allstats <- bind_rows(sA[, keep], sB[, keep])
allstats$label <- gsub("\n", " ", allstats$label)
allstats$partial_r2 <- allstats$t^2 / (allstats$t^2 + allstats$df)
write.csv(allstats, OUT_CSV, row.names = FALSE)

cat("\n--- Figure 3, covariate-adjusted (HC3) ---\n")
for (i in seq_len(nrow(allstats))) with(allstats[i, ],
                                        cat(sprintf("%-26s beta=%+.3f SD  SE=%.3f  t(%d)=%+.2f  P=%.3g  partial R2=%.3f  n=%d\n",
                                                    label, beta_sd, se_sd, df, t, p, partial_r2, n)))

cat("\nsaved:", paste0(OUTBASE, ".png"), "\n")
cat("saved:", paste0(OUTBASE, ".svg"), "\n")
cat("saved:", paste0(OUTBASE, ".pdf"), "\n")
cat("saved:", OUT_CSV, "\n")

cat("\nNext: convert PNG to TIFF for submission --\n")
cat("python -c \"from PIL import Image; im=Image.open('", paste0(OUTBASE, ".png"),
    "').convert('RGB'); im.save('", paste0(OUTBASE, ".tif"),
    "', compression='tiff_lzw', dpi=(600,600))\"\n", sep = "")
cat("\nREMINDER for the figure legend in the manuscript file:\n")
cat("  - the two facets in A have independent y-axes\n")
cat("  - error bars are SE of covariate-adjusted means\n")
cat("  - define *P <= 0.05, **P <= 0.01, ***P <= 0.001\n")
cat("  - list the covariates\n")