# ==============================================================================
# 03a_build_design.R
#
# Design matrices and contrast files for the network-based statistics, from
# the primary matched sample.
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Network-based statistics -> Functional NBS and hub analysis;
#          Structural NBS; White matter hyperintensity location sensitivity;
#          Supplementary material, "Motion-related Quality Control" (v)
#
# Functional design (8 columns): AF, control, age, education, sex, diastolic
#   BP, log total WMH, mean rs-fMRI head motion (mean framewise displacement).
# Structural design (7 columns): as above without head motion, which indexes
#   functional acquisition motion and is not meaningful for diffusion data.
# WMH sensitivity: the network-based statistics were repeated with the deep
#   and periventricular specifications only, so three designs are written per
#   modality: design_1_total (primary), design_2_peri, design_3_deep.
# Contrasts: AF > Control [+1 -1 0 ...] and Control > AF [-1 +1 0 ...].
#
# High-motion sensitivity cohort: set FD_EXCLUDE_MM to 0.30 to exclude
#   participants with mean rs-fMRI head motion above that value; the designs
#   are then written to OUT_BASE_sens and subjects.txt lists the retained
#   participants for building the matching matrix directory (03d).
#
# Row order: participants sorted by eid as character strings. NBS-Connectome
#   reads the matrix files of a directory in alphabetical order, and the
#   matrix files are named conn_<eid>.txt (functional) or subjectNNN.txt in
#   eid order (structural), so design row k corresponds to the kth file.
#   If MATS_DIR is set, the script verifies that exactly one matrix file
#   exists per participant and no other files are present.
#
# Input:   cohort_matched.csv                (01_cohort_assembly.R)
# Outputs: <OUT_BASE>/{fc,sc}/design_<k>_<wmh>/design.txt, design_reference.csv,
#          subjects.txt, contrasts/contrast_01_af_gt_control.txt,
#          contrasts/contrast_02_control_gt_af.txt
# ==============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr) })

# ------------------------------------------------------------------------------
# Configuration  --- edit paths for your environment
# ------------------------------------------------------------------------------
COHORT_CSV    <- "/path/to/derivatives/cohort_matched.csv"
OUT_BASE      <- "/path/to/nbs/01_designs"
MATS_DIR      <- ""            # optional: directory of conn_<eid>.txt files to verify alignment
FD_EXCLUDE_MM <- NA            # NA = main cohort; 0.30 = high-motion sensitivity cohort

COL_EID <- "eid"; COL_AF <- "af"; COL_AGE <- "age_i2"; COL_EDU <- "education_i2"
COL_SEX <- "sex"; COL_DBP <- "dia_mean"; COL_WMH_TOT <- "log_norm_wmh_ml"
COL_WMH_PERI <- "log_norm_peri_wmh_ml"; COL_WMH_DEEP <- "log_norm_deep_wmh_ml"
COL_MEANFD <- "MeanFD_i2"

# ==============================================================================
# 1. LOAD AND PREPARE
# ==============================================================================
cohort <- read_csv(COHORT_CSV, show_col_types = FALSE, col_types = cols(.default = col_guess(), eid = col_character())) %>%
  arrange(.data[[COL_EID]])

needed <- c(COL_EID, COL_AF, COL_AGE, COL_EDU, COL_SEX, COL_DBP, COL_WMH_TOT, COL_WMH_PERI, COL_WMH_DEEP, COL_MEANFD)
missing <- setdiff(needed, names(cohort))
if (length(missing)) stop("Missing column(s) in cohort_matched.csv: ", paste(missing, collapse = ", "))

if (!is.na(FD_EXCLUDE_MM)) {
  n0 <- nrow(cohort)
  cohort <- cohort %>% filter(.data[[COL_MEANFD]] <= FD_EXCLUDE_MM)
  OUT_BASE <- paste0(OUT_BASE, "_sens")
  cat(sprintf("High-motion sensitivity cohort: excluded %d participants with motion > %.2f mm\n",
              n0 - nrow(cohort), FD_EXCLUDE_MM))
}

z <- function(x) as.numeric(scale(x))
cohort <- cohort %>%
  mutate(age_z = z(.data[[COL_AGE]]), edu_z = z(.data[[COL_EDU]]), dbp_z = z(.data[[COL_DBP]]),
         log_wmh_z = z(.data[[COL_WMH_TOT]]), log_peri_z = z(.data[[COL_WMH_PERI]]),
         log_deep_z = z(.data[[COL_WMH_DEEP]]), meanfd_z = z(.data[[COL_MEANFD]]))
cat(sprintf("Cohort: N = %d (AF = %d, Control = %d)\n", nrow(cohort),
            sum(cohort[[COL_AF]] == 1), sum(cohort[[COL_AF]] == 0)))

if (nzchar(MATS_DIR)) {
  files <- sort(list.files(MATS_DIR, pattern = "^conn_.*\\.txt$"))
  expected <- paste0("conn_", cohort[[COL_EID]], ".txt")
  if (!identical(files, expected))
    stop("Matrix directory does not match the cohort: ", length(files), " files vs ",
         length(expected), " participants (or different eids/order). Every participant needs ",
         "exactly one conn_<eid>.txt and the directory must contain nothing else.")
  cat("Matrix directory verified against the cohort.\n")
}

# ==============================================================================
# 2. WRITE ONE DESIGN VARIANT
# ==============================================================================
write_design <- function(out_dir, design_mat, col_labels, eids) {
  dir.create(file.path(out_dir, "contrasts"), showWarnings = FALSE, recursive = TRUE)
  write.table(design_mat, file.path(out_dir, "design.txt"), sep = " ",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  ref <- as.data.frame(design_mat); names(ref) <- col_labels
  write_csv(cbind(row_in_design = seq_len(nrow(design_mat)), eid = eids, ref),
            file.path(out_dir, "design_reference.csv"))
  writeLines(eids, file.path(out_dir, "subjects.txt"))
  k <- ncol(design_mat)
  write.table(matrix(c( 1, -1, rep(0, k - 2)), nrow = 1),
              file.path(out_dir, "contrasts", "contrast_01_af_gt_control.txt"),
              sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(matrix(c(-1,  1, rep(0, k - 2)), nrow = 1),
              file.path(out_dir, "contrasts", "contrast_02_control_gt_af.txt"),
              sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# ==============================================================================
# 3. DESIGNS
# ==============================================================================
af <- cohort[[COL_AF]]; ctrl <- 1 - af; sex <- cohort[[COL_SEX]]; eids <- cohort[[COL_EID]]
wmh_variants <- list(design_1_total = "log_wmh_z", design_2_peri = "log_peri_z", design_3_deep = "log_deep_z")

for (dn in names(wmh_variants)) {
  wcol <- wmh_variants[[dn]]
  fc <- cbind(af, ctrl, cohort$age_z, cohort$edu_z, sex, cohort$dbp_z, cohort[[wcol]], cohort$meanfd_z)
  sc <- cbind(af, ctrl, cohort$age_z, cohort$edu_z, sex, cohort$dbp_z, cohort[[wcol]])
  write_design(file.path(OUT_BASE, "fc", dn), fc,
               c("af", "control", "age_z", "edu_z", "sex", "dbp_z", wcol, "meanfd_z"), eids)
  write_design(file.path(OUT_BASE, "sc", dn), sc,
               c("af", "control", "age_z", "edu_z", "sex", "dbp_z", wcol), eids)
  cat(sprintf("  wrote fc/%s (8 cols) and sc/%s (7 cols)\n", dn, dn))
}

cat("\nSanity check, primary FC design (first row):\n")
print(read_csv(file.path(OUT_BASE, "fc", "design_1_total", "design_reference.csv"), show_col_types = FALSE)[1, ])
cat(sprintf("\nGroup totals: AF = %d, Control = %d\nDone.\n", sum(af), sum(ctrl)))
