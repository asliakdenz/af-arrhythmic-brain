# ==============================================================================
# 01_cohort_assembly.R
#
# Cohort assembly: ECG-confirmed atrial fibrillation versus propensity-matched
# controls from the UK Biobank imaging cohort.
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Participants and matching; Supplementary Fig. 1; Supplementary
#          Table 1
#
# Inputs (UK Biobank extract; never committed to this repository):
#   - cohort_final.csv        one row per participant: imaging-visit ECG report
#                             text (field 12653, instance 2 and 3), demographics,
#                             vascular risk factors, cognition, mental health
#                             items, WMH volumes (25781, 24485, 24486), TIV
#                             (26521), mean rs-fMRI head motion (25741),
#                             self-reported illness (20002) and main ICD-10
#                             diagnoses (41202)
#   - hesin_diag_cohort.csv   long table of hospital-episode ICD-10 diagnoses
#                             (eid, diag_icd10); the record-level content of
#                             field 41270
#   - subject_list.csv        identifiers (column subject_id) of participants
#                             with complete, quality-controlled functional and
#                             structural connectomes in the UK Biobank
#                             Connectome Resource
#
# Outputs (written to OUTPUT_DIR):
#   - cohort_matched.csv                 primary analysis sample: participants
#                                        with complete cognitive testing,
#                                        1:2 propensity-matched
#   - cohort_heldout_matched.csv         internal held-out sample: participants
#                                        with incomplete cognitive testing,
#                                        matched separately with the same
#                                        procedure
#   - cohort_unmatched_cogcomplete.csv   all eligible participants with complete
#                                        cognitive testing, not restricted by
#                                        connectome quality control and not
#                                        matched (unmatched-pool sensitivity)
#   - exclusion_summary.csv, flow_counts.csv
#
# Procedure (Methods, "Participants and matching"):
#   1. Exclusion of pre-existing conditions listed in Supplementary Table 1,
#      using ICD-10 hospital diagnoses (fields 41270 / 41202) and self-reported
#      illness codes (field 20002). Exclusions are applied sequentially in the
#      order of the table; each participant is counted once, at the first
#      criterion met.
#   2. AF cases: atrial fibrillation on the imaging-visit 12-lead ECG report.
#      Controls: normal sinus rhythm on the imaging-visit ECG, no AF on the ECG
#      at either imaging visit, and no AF diagnosis in ICD-based records.
#   3. Imaging-derived measures (TIV, total, periventricular and deep WMH, mean
#      rs-fMRI head motion) are required to be complete and are not imputed.
#   4. Missing clinical matching covariates are imputed by chained equations
#      (mice), using only the clinical covariates, age, sex and AF status as
#      predictors. The imputation is run once, on the eligible pool, so that
#      every branch below uses the same imputed values.
#   5. Branches: the unmatched pool keeps every eligible participant with
#      complete cognitive testing; the connectome-restricted pool is split into
#      the primary sample (complete cognitive testing) and the held-out sample
#      (incomplete cognitive testing). Each of the last two is matched 1:2
#      (AF:control, nearest neighbour on the propensity score, ATT estimand)
#      on age, sex, education, systolic and diastolic blood pressure, BMI, total
#      cholesterol, blood glucose, smoking status and total WMH load.
#
# The cohort used in the reported analyses was generated once and frozen;
# downstream scripts read the written files rather than re-running this step.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(mice)
  library(MatchIt)
  library(cobalt)
})

# ------------------------------------------------------------------------------
# Configuration  --- edit paths for your environment
# ------------------------------------------------------------------------------
INPUT_DIR           <- "/path/to/derivatives"
ECG_COHORT_CSV      <- file.path(INPUT_DIR, "cohort_final.csv")
HESIN_DIAG_CSV      <- file.path(INPUT_DIR, "hesin_diag_cohort.csv")
COMMON_SUBJECTS_CSV <- "/path/to/ukb_extract/subject_list.csv"

OUTPUT_DIR <- INPUT_DIR

# Supplementary Table 1 lists "Active malignancy (within 5 years,
# brain-affecting)" as "various C-codes" without enumerating them. Enter the
# ICD-10 code prefixes used for that criterion here. Left empty, the criterion
# is not applied.
MALIGNANCY_ICD10 <- character(0)

# ==============================================================================
# 1. LOAD AND PARSE THE EXTRACT
# ==============================================================================
frac_numeric <- function(v) {
  x <- trimws(as.character(v))
  x <- x[nzchar(x) & !is.na(x)]
  if (!length(x)) return(0)
  suppressWarnings(mean(!is.na(as.numeric(x))))
}

# columns that must never be coerced to numeric (list-like or free text)
never_numeric_patterns <- c("^eid$", "^p12653", "^p20002", "^p41202", "^p6138")

data <- read_csv(ECG_COHORT_CSV,
                 col_types = cols(.default = col_character()),
                 show_col_types = FALSE) %>%
  select(-starts_with("...")) %>%
  mutate(eid = as.character(eid))

is_num        <- sapply(data, frac_numeric)
never_numeric <- names(data)[Reduce(`|`, lapply(never_numeric_patterns,
                                                 function(rx) grepl(rx, names(data))))]
to_numeric    <- setdiff(names(data)[is_num >= 0.98], never_numeric)
data <- data %>% mutate(across(all_of(to_numeric), ~ suppressWarnings(as.numeric(.))))

hesin_diag <- read_csv(HESIN_DIAG_CSV,
                       col_types = cols(.default = col_character()),
                       show_col_types = FALSE) %>%
  select(-starts_with("...")) %>%
  mutate(eid        = as.character(eid),
         diag_nodot = toupper(gsub("\\.", "", trimws(diag_icd10)))) %>%
  filter(!is.na(diag_nodot), nzchar(diag_nodot))

flow <- tibble(step = character(), n = integer(), n_af = integer(), n_control = integer())
add_flow <- function(step, df) {
  n_af <- if ("af" %in% names(df)) sum(df$af == 1) else NA_integer_
  n_ct <- if ("af" %in% names(df)) sum(df$af == 0) else NA_integer_
  flow <<- bind_rows(flow, tibble(step = step, n = nrow(df), n_af = n_af, n_control = n_ct))
  message(sprintf("%-55s n = %d%s", step, nrow(df),
                  if (!is.na(n_af)) sprintf("  (AF %d | Control %d)", n_af, n_ct) else ""))
}
add_flow("Extract loaded", data)

# helpers for list-like fields ("['I63.9', 'G45']" or "1081,1263")
parse_list_field <- function(x) {
  if (is.null(x) || is.na(x)) return(character(0))
  s <- gsub("\\[|\\]|'|\"", "", as.character(x))
  toks <- trimws(unlist(strsplit(s, ",")))
  toks[nzchar(toks)]
}

# ==============================================================================
# 2. EXCLUSIONS  (Supplementary Table 1)
# ==============================================================================
# ICD-10 entries are code prefixes (dots removed at matching time), so "I60"
# covers I60.x and "F06.7" covers F06.7x. Self-reported codes are field 20002
# (coding 6) values. A participant is excluded if any ICD-10 code in the main
# diagnoses (41202) or in the hospital-episode diagnoses (41270), or any
# self-reported code, matches the criterion.
exclusion_criteria <- list(
  list(name = "Stroke (any subtype)",
       icd = c("I60", "I61", "I63", "I64"),
       sr  = c(1081, 1082, 1083, 1086, 1491, 1583)),
  list(name = "Transient ischaemic attack",
       icd = "G45", sr = 1082),
  list(name = "Other cerebrovascular disease",
       icd = c("I65", "I66", "I67", "I68", "I69"), sr = c(1086, 1583)),
  list(name = "Subarachnoid / intracerebral haemorrhage",
       icd = c("I60", "I61", "I62"), sr = c(1086, 1491)),
  list(name = "Dementia (all causes)",
       icd = c("F00", "F01", "F02", "F03", "G30"), sr = 1263),
  list(name = "Mild cognitive impairment",
       icd = "F06.7", sr = NULL),
  list(name = "Parkinson's disease",
       icd = "G20", sr = 1262),
  list(name = "Other extrapyramidal disease",
       icd = c("G21", "G22", "G23", "G24", "G25"), sr = NULL),
  list(name = "Multiple sclerosis",
       icd = "G35", sr = 1261),
  list(name = "Other demyelinating disease",
       icd = c("G36", "G37"), sr = NULL),
  list(name = "Motor neuron disease",
       icd = "G12", sr = 1259),
  list(name = "Epilepsy",
       icd = c("G40", "G41"), sr = 1264),
  list(name = "Migraine",
       icd = "G43", sr = 1265),
  list(name = "Brain tumour (any)",
       icd = c("C71", "D33", "D43"), sr = c(1031, 1032, 1437, 1438)),
  list(name = "Encephalitis",
       icd = c("G04", "G05"), sr = 1247),
  list(name = "Meningitis",
       icd = c("G00", "G01", "G02", "G03"), sr = 1246),
  list(name = "Traumatic brain injury",
       icd = c("S06", "S07"), sr = 1266),
  list(name = "Hydrocephalus",
       icd = "G91", sr = NULL),
  list(name = "Schizophrenia",
       icd = paste0("F2", 0:9), sr = 1289),
  list(name = "Bipolar disorder",
       icd = c("F30", "F31"), sr = 1291),
  list(name = "Major depressive disorder (recurrent severe)",
       icd = c("F33.2", "F33.3"), sr = NULL),
  list(name = "Substance use disorders (other than alcohol)",
       icd = paste0("F1", 1:9), sr = c(1408, 1409, 1410, 1411)),
  list(name = "Severe psychiatric disorders, other",
       icd = c(paste0("F0", 0:9), paste0("F2", 2:9)), sr = NULL),
  list(name = "Cardiac arrest / sudden cardiac death",
       icd = "I46", sr = NULL),
  list(name = "Severe valvular heart disease",
       icd = c("I05", "I06", "I07", "I08", paste0("I3", 4:9)), sr = c(1490, 1586, 1587)),
  list(name = "Cardiomyopathy (severe)",
       icd = c("I42", "I43"), sr = 1079),
  list(name = "Implanted cardiac device (pacemaker / ICD)",
       icd = "Z95.0", sr = NULL),
  list(name = "Active malignancy (brain-affecting)",
       icd = MALIGNANCY_ICD10, sr = NULL),
  list(name = "End-stage renal / hepatic disease",
       icd = c("N18.5", "N18.6", "K72", "K74"), sr = NULL)
)

# per-participant code sets
icd_main_list <- lapply(data$p41202, function(x) toupper(gsub("\\.", "", parse_list_field(x))))
hes_by_eid    <- split(hesin_diag$diag_nodot, hesin_diag$eid)
icd_hes_list  <- lapply(data$eid, function(e) { v <- hes_by_eid[[e]]; if (is.null(v)) character(0) else v })
icd_all_list  <- Map(c, icd_main_list, icd_hes_list)

sr_cols <- grep("^p20002", names(data), value = TRUE)
sr_mat  <- as.matrix(data[sr_cols])
sr_list <- lapply(seq_len(nrow(data)), function(i) {
  if (!length(sr_cols)) return(numeric(0))
  v <- suppressWarnings(as.numeric(unlist(lapply(sr_mat[i, ], parse_list_field))))
  v[!is.na(v)]
})

matches_criterion <- function(icd_codes, sr_codes, crit) {
  icd_prefixes <- toupper(gsub("\\.", "", crit$icd))
  hit_icd <- length(icd_prefixes) > 0 && length(icd_codes) > 0 &&
    any(vapply(icd_prefixes, function(p) any(startsWith(icd_codes, p)), logical(1)))
  hit_sr  <- length(crit$sr) > 0 && length(sr_codes) > 0 && any(sr_codes %in% crit$sr)
  hit_icd || hit_sr
}

remaining <- rep(TRUE, nrow(data))
exclusion_summary <- tibble(criterion = character(), excluded = integer(), remaining = integer())
for (crit in exclusion_criteria) {
  hit <- remaining & vapply(seq_len(nrow(data)), function(i)
    matches_criterion(icd_all_list[[i]], sr_list[[i]], crit), logical(1))
  remaining <- remaining & !hit
  exclusion_summary <- bind_rows(exclusion_summary,
                                 tibble(criterion = crit$name, excluded = sum(hit),
                                        remaining = sum(remaining)))
}
data <- data[remaining, ]
add_flow("After exclusions (Supplementary Table 1)", data)
print(exclusion_summary, n = Inf)

# ==============================================================================
# 3. AF AND CONTROL DEFINITIONS
# ==============================================================================
# 3a. ECG report text at the imaging visit (instance 2) and the repeat imaging
#     visit (instance 3)
i2_cols <- grep("^p12653_i2_a(0|[1-9]|1[0-4])$", names(data), value = TRUE)
i3_cols <- grep("^p12653_i3_a(0|[1-9]|1[0-4])$", names(data), value = TRUE)
to_low  <- function(x) stringr::str_to_lower(as.character(x))

af_rx  <- paste0("(",
                 "atrial\\s*fibrill\\w*",
                 "|",
                 "(?:a\\s*[-\\.]?\\s*fib\\w*|\\baf\\b)(?!\\s*:?\\s*1\\b|:\\s*1\\b)",
                 ")")
nsr_rx <- "(normal\\s+sinus\\s+rhythm|\\bnsr\\b|\\bsinus\\s+rhythm\\b)"

detect_any <- function(cols, rx) {
  if (!length(cols)) return(rep(FALSE, nrow(data)))
  m <- data %>% transmute(across(all_of(cols), ~ str_detect(to_low(.x), regex(rx, TRUE))))
  rowSums(m, na.rm = TRUE) > 0L
}
af_in_i2  <- detect_any(i2_cols, af_rx)
nsr_in_i2 <- detect_any(i2_cols, nsr_rx)
af_in_i3  <- detect_any(i3_cols, af_rx)

# 3b. AF in ICD-based records: any I48 code other than atrial flutter
#     (I48.3, I48.4), in main diagnoses or hospital-episode diagnoses
is_af_code <- function(codes) any(startsWith(codes, "I48") & !(codes %in% c("I483", "I484")))
af_icd <- vapply(icd_all_list[remaining], is_af_code, logical(1))

data <- data %>%
  mutate(.af_i2 = af_in_i2, .nsr_i2 = nsr_in_i2, .af_i3 = af_in_i3, .af_icd = af_icd)

af_df      <- data %>% filter(.af_i2) %>% mutate(af = 1L)
control_df <- data %>% filter(.nsr_i2, !.af_i2, !.af_i3, !.af_icd) %>% mutate(af = 0L)
stopifnot(length(intersect(af_df$eid, control_df$eid)) == 0)

data <- bind_rows(af_df, control_df) %>% select(-.af_i2, -.nsr_i2, -.af_i3, -.af_icd)
add_flow("Group assignment (ECG-defined AF; sinus-rhythm controls)", data)

# ==============================================================================
# 4. VARIABLES
# ==============================================================================
# 4a. Education: highest qualification (field 6138) mapped to ISCED level
isced_mapping <- c("1" = 5L, "2" = 3L, "3" = 2L, "4" = 2L, "5" = 5L, "6" = 4L,
                   "-7" = 0L, "-3" = NA_integer_)
convert_education <- function(x) {
  v <- suppressWarnings(as.integer(parse_list_field(x)))
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_integer_)
  mapped <- isced_mapping[as.character(v)]
  if (all(is.na(mapped))) return(NA_integer_)
  as.integer(max(mapped, na.rm = TRUE))
}

pick_num <- function(df, col) if (col %in% names(df)) suppressWarnings(as.numeric(df[[col]])) else NA_real_
rowmeans_na <- function(...) {
  m <- cbind(...)
  v <- rowMeans(m, na.rm = TRUE); v[is.nan(v)] <- NA_real_; v
}

data <- data %>%
  mutate(
    age_i2       = pick_num(., "p21003_i2"),
    sex          = pick_num(., "p31"),
    education_i2 = vapply(p6138_i2, convert_education, integer(1)),
    bmi_i2       = coalesce(pick_num(., "p23104_i2"), pick_num(., "p21001_i2")),
    sys_mean     = coalesce(rowmeans_na(pick_num(., "p4080_i2_a0"), pick_num(., "p4080_i2_a1")),
                            rowmeans_na(pick_num(., "p93_i2_a0"),   pick_num(., "p93_i2_a1"))),
    dia_mean     = coalesce(rowmeans_na(pick_num(., "p4079_i2_a0"), pick_num(., "p4079_i2_a1")),
                            rowmeans_na(pick_num(., "p94_i2_a0"),   pick_num(., "p94_i2_a1"))),
    cholesterol  = rowmeans_na(pick_num(., "p30690_i0"), pick_num(., "p30690_i1")),
    glucose      = rowmeans_na(pick_num(., "p30740_i0"), pick_num(., "p30740_i1")),
    smo          = case_when(pick_num(., "p1239_i2") %in% c(1, 2) ~ 1L,
                             pick_num(., "p1239_i2") == 0         ~ 0L,
                             TRUE                                 ~ NA_integer_),
    # imaging-derived measures (instance 2)
    tiv_i2       = pick_num(., "p26521_i2"),
    total_wmh_i2 = pick_num(., "p25781_i2"),
    peri_wmh_i2  = pick_num(., "p24485_i2"),
    deep_wmh_i2  = pick_num(., "p24486_i2"),
    MeanFD_i2    = pick_num(., "p25741_i2"),
    # cognition (instance 2)
    tower_rearr_i2          = pick_num(., "p21004_i2"),
    trail_mak_1_i2          = pick_num(., "p6348_i2"),
    trail_mak_2_i2          = pick_num(., "p6350_i2"),
    reaction_time_i2        = pick_num(., "p20023_i2"),
    sym_dig_attempt_i2      = pick_num(., "p23323_i2"),
    sym_dig_corr_i2         = pick_num(., "p23324_i2"),
    num_mem_nmbr_i2         = pick_num(., "p4282_i2"),
    paired_assc_learn_i2    = pick_num(., "p20197_i2"),
    pros_mem_cog_summ_i2    = pick_num(., "p20018_i2"),
    fluid_intl_i2           = pick_num(., "p20016_i2"),
    mat_pat_compl_corr_i2   = pick_num(., "p6373_i2"),
    mat_pat_compl_viewed_i2 = pick_num(., "p6374_i2"),
    sym_dig_accuracy        = ifelse(!is.na(sym_dig_corr_i2) & !is.na(sym_dig_attempt_i2) &
                                       sym_dig_attempt_i2 > 0,
                                     sym_dig_corr_i2 / sym_dig_attempt_i2, NA_real_),
    mat_pat_compl_accuracy  = ifelse(!is.na(mat_pat_compl_corr_i2) & !is.na(mat_pat_compl_viewed_i2) &
                                       mat_pat_compl_viewed_i2 > 0,
                                     mat_pat_compl_corr_i2 / mat_pat_compl_viewed_i2, NA_real_),
    # mental health items (instance 2)
    mh_dep_mood_2w_i2        = pick_num(., "p2050_i2"),
    mh_anhedonia_2w_i2       = pick_num(., "p2060_i2"),
    mh_tense_restless_2w_i2  = pick_num(., "p2070_i2"),
    mh_tired_lethargic_2w_i2 = pick_num(., "p2080_i2"),
    wb_happiness_i2               = pick_num(., "p4526_i2"),
    wb_health_satisfaction_i2     = pick_num(., "p4548_i2"),
    wb_work_satisfaction_i2       = pick_num(., "p4537_i2"),
    wb_family_satisfaction_i2     = pick_num(., "p4559_i2"),
    wb_friendship_satisfaction_i2 = pick_num(., "p4570_i2"),
    wb_financial_satisfaction_i2  = pick_num(., "p4581_i2")
  )

# 4b. WMH normalised by total intracranial volume (mL per litre of TIV) and
#     log(1 + x) transformed. Missing values are kept missing (no imputation).
normalize_wmh <- function(wmh, tiv) {
  v <- (wmh / tiv) * 1000
  v[!is.finite(v)] <- NA_real_
  v
}
data <- data %>%
  mutate(norm_wmh_ml          = normalize_wmh(total_wmh_i2, tiv_i2),
         norm_peri_wmh_ml     = normalize_wmh(peri_wmh_i2,  tiv_i2),
         norm_deep_wmh_ml     = normalize_wmh(deep_wmh_i2,  tiv_i2),
         log_norm_wmh_ml      = log1p(norm_wmh_ml),
         log_norm_peri_wmh_ml = log1p(norm_peri_wmh_ml),
         log_norm_deep_wmh_ml = log1p(norm_deep_wmh_ml))

# 4c. Imaging-derived measures must be complete (not imputed)
data <- data %>% drop_na(tiv_i2, norm_wmh_ml, norm_peri_wmh_ml, norm_deep_wmh_ml, MeanFD_i2)
add_flow("Complete imaging-derived measures (TIV, WMH, head motion)", data)

cognitive_tests <- c("tower_rearr_i2", "trail_mak_2_i2",
                     "reaction_time_i2", "sym_dig_accuracy", "trail_mak_1_i2",
                     "num_mem_nmbr_i2", "paired_assc_learn_i2", "pros_mem_cog_summ_i2",
                     "fluid_intl_i2", "mat_pat_compl_accuracy")
data <- data %>% mutate(cognition_complete = complete.cases(across(all_of(cognitive_tests))))

columns_to_keep <- c(
  "eid", "af", "age_i2", "sex", "education_i2", "bmi_i2", "sys_mean", "dia_mean",
  "cholesterol", "glucose", "smo",
  "tiv_i2", "total_wmh_i2", "peri_wmh_i2", "deep_wmh_i2",
  "norm_wmh_ml", "norm_peri_wmh_ml", "norm_deep_wmh_ml",
  "log_norm_wmh_ml", "log_norm_peri_wmh_ml", "log_norm_deep_wmh_ml", "MeanFD_i2",
  cognitive_tests, "sym_dig_attempt_i2", "sym_dig_corr_i2",
  "mat_pat_compl_corr_i2", "mat_pat_compl_viewed_i2", "cognition_complete",
  "mh_dep_mood_2w_i2", "mh_anhedonia_2w_i2", "mh_tense_restless_2w_i2", "mh_tired_lethargic_2w_i2",
  "wb_happiness_i2", "wb_health_satisfaction_i2", "wb_work_satisfaction_i2",
  "wb_family_satisfaction_i2", "wb_friendship_satisfaction_i2", "wb_financial_satisfaction_i2"
)
data <- data[, columns_to_keep]

# ==============================================================================
# 5. IMPUTATION OF CLINICAL MATCHING COVARIATES (chained equations)
# ==============================================================================
# Only clinical covariates, age, sex and AF status enter the imputation model;
# no imaging or cognitive variable is used as a predictor.
imputation_targets    <- c("bmi_i2", "education_i2", "sys_mean", "dia_mean",
                           "cholesterol", "glucose", "smo")
imputation_predictors <- c("af", "age_i2", "sex")

miss_rate <- colMeans(is.na(data[imputation_targets]))
message("Missingness of clinical matching covariates (fraction):")
print(round(miss_rate, 3))

if (any(miss_rate > 0)) {
  age_mean <- mean(data$age_i2, na.rm = TRUE); age_sd <- sd(data$age_i2, na.rm = TRUE)
  edu_levels <- sort(unique(na.omit(data$education_i2)))
  impute_prep <- data %>%
    select(eid, all_of(c(imputation_targets, imputation_predictors))) %>%
    mutate(age_i2       = (age_i2 - age_mean) / age_sd,
           education_i2 = factor(education_i2, levels = edu_levels, ordered = TRUE),
           smo          = factor(smo, levels = c(0, 1)))

  meth <- make.method(impute_prep)
  meth[c("bmi_i2", "cholesterol", "glucose", "sys_mean", "dia_mean")] <- "norm"
  meth["education_i2"] <- if (length(edu_levels) > 2) "polr" else "logreg"
  meth["smo"]          <- "logreg"
  meth["eid"]          <- ""

  init_mice <- mice(impute_prep, maxit = 0, method = meth, printFlag = FALSE)
  pred_mat  <- init_mice$predictorMatrix
  pred_mat[, ] <- 0
  for (tgt in imputation_targets) {
    pred_mat[tgt, c(imputation_targets, imputation_predictors)] <- 1
    pred_mat[tgt, tgt] <- 0
  }

  imp_results <- mice(impute_prep, m = 20, maxit = 10, method = meth,
                      predictorMatrix = pred_mat, ridge = 1e-4, printFlag = FALSE)

  completed <- complete(imp_results, 1) %>%
    mutate(age_i2       = age_i2 * age_sd + age_mean,
           smo          = as.integer(as.character(smo)),
           education_i2 = as.integer(as.character(education_i2))) %>%
    select(eid, all_of(imputation_targets))

  data <- data %>% select(-all_of(imputation_targets)) %>% left_join(completed, by = "eid")
}
stopifnot(!anyNA(data[imputation_targets]))

# ==============================================================================
# 6. BRANCHES
# ==============================================================================
standardise_covariates <- function(df) {
  z <- function(x) as.numeric(scale(x))
  df %>% mutate(age_std           = z(age_i2),
                education_std     = z(education_i2),
                dia_bp_std        = z(dia_mean),
                log_total_wmh_std = z(log_norm_wmh_ml),
                log_peri_wmh_std  = z(log_norm_peri_wmh_ml),
                log_deep_wmh_std  = z(log_norm_deep_wmh_ml),
                mean_fd_std       = z(MeanFD_i2))
}

# 6a. Unmatched pool: every eligible participant with complete cognitive
#     testing, not restricted by connectome quality control
unmatched_pool <- data %>% filter(cognition_complete) %>% standardise_covariates()
add_flow("Unmatched pool: complete cognitive testing (no connectome restriction)", unmatched_pool)
write_csv(unmatched_pool, file.path(OUTPUT_DIR, "cohort_unmatched_cogcomplete.csv"))

# 6b. Connectome-restricted pool
common_subjects <- read_csv(COMMON_SUBJECTS_CSV, show_col_types = FALSE,
                            col_types = cols(.default = col_character()))
connectome_pool <- data %>% filter(eid %in% common_subjects$subject_id)
add_flow("Complete, quality-controlled functional and structural connectomes", connectome_pool)

primary_pool <- connectome_pool %>% filter(cognition_complete)
heldout_pool <- connectome_pool %>% filter(!cognition_complete)
add_flow("Primary branch before matching: complete cognitive testing", primary_pool)
add_flow("Held-out branch before matching: incomplete cognitive testing", heldout_pool)

# ==============================================================================
# 7. PROPENSITY SCORE MATCHING  (1:2 nearest neighbour, ATT)
# ==============================================================================
matching_formula <- af ~ age_i2 + sex + education_i2 + bmi_i2 + sys_mean + dia_mean +
                         smo + glucose + cholesterol + norm_wmh_ml

match_branch <- function(df, label) {
  if (sum(df$af == 1) == 0 || sum(df$af == 0) < 2) {
    message(sprintf("[%s] too few participants to match; skipped", label))
    return(NULL)
  }
  mod <- matchit(matching_formula, data = df, method = "nearest", ratio = 2,
                 replace = FALSE, estimand = "ATT", m.order = "random")
  matched <- match.data(mod, data = df) %>%
    group_by(subclass) %>%
    filter(sum(af == 1) == 1, sum(af == 0) == 2) %>%   # complete 1:2 sets only
    ungroup() %>%
    standardise_covariates()
  cat(sprintf("\n--- Covariate balance after matching: %s ---\n", label))
  print(bal.tab(mod, binary = "std", thresholds = c(m = 0.10)))
  matched
}

cohort_matched <- match_branch(primary_pool, "primary sample")
if (!is.null(cohort_matched)) {
  add_flow("Primary sample after 1:2 matching", cohort_matched)
  write_csv(cohort_matched, file.path(OUTPUT_DIR, "cohort_matched.csv"))
}

cohort_heldout <- match_branch(heldout_pool, "held-out sample")
if (!is.null(cohort_heldout)) {
  add_flow("Held-out sample after 1:2 matching", cohort_heldout)
  write_csv(cohort_heldout, file.path(OUTPUT_DIR, "cohort_heldout_matched.csv"))
}

# ==============================================================================
# 8. BOOK-KEEPING
# ==============================================================================
write_csv(exclusion_summary, file.path(OUTPUT_DIR, "exclusion_summary.csv"))
write_csv(flow, file.path(OUTPUT_DIR, "flow_counts.csv"))
message("\nParticipant flow (Supplementary Fig. 1):")
print(flow, n = Inf)
message("\nOutputs written to: ", OUTPUT_DIR)
