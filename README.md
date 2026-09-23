# The arrhythmic brain

Analysis code accompanying:

> Akdeniz A, Kiakou D, Mueller K, Valk S, Villringer A.
> **The arrhythmic brain: Interoceptive overload and cognitive slowing in atrial fibrillation.**
> *Under review.*
>
> Preprint: <https://doi.org/10.21203/rs.3.rs-10662928/v1>

## Overview

This repository contains the analysis pipeline for a cross-sectional study of
resting-state functional and diffusion MRI connectomes in atrial fibrillation,
using the UK Biobank imaging cohort under application no. 37721. It covers
cohort assembly and propensity score matching, connectome construction and
motion quality control, network-based statistics and hub analysis, network
segregation and integration, the control-referenced functional deviation
score, brain–behaviour modelling with nested cross-validation, the sensitivity
analyses and figure generation.

**No data are included.** UK Biobank participant data are obtained through a
formal application. In line with UK Biobank's data-sharing policy,
individual-level derived measures generated for this study are not
redistributed here; they are returned to UK Biobank and made available to
approved researchers through the standard access procedure.

## Repository structure

```
.
├── 01_cohort/                    Cohort assembly, matching, sample characteristics
├── 02_construction/              Connectome construction and motion QC
├── 03_nbs/                       Network-based statistics pipeline
├── 04_subnetwork/                Subnetwork characterisation, labels, validation
├── 05_segregation_integration/   System segregation and participation
├── 06_deviation_scores/          Functional deviation score
├── 07_statistics/                Outcome models, brain–behaviour models, sensitivity
├── 08_figures/                   Figures 1 and 3
└── data/atlas/                   216-node atlas reference files (Schaefer 200 + Tian S1)
```

## Methods sections → scripts

| Methods subsection | Script |
|---|---|
| Participants and matching | `01_cohort/01_cohort_assembly.R` |
| — Sample characteristics (Tables 1–2, Supplementary Fig. 2) | `01_cohort/01b_sample_characteristics.R` |
| Cognitive outcomes (Table 3, Supplementary Table 3) | `07_statistics/07a_cognitive_outcomes.R` |
| — Unmatched-pool sensitivity (Supplementary Fig. 3) | `07_statistics/07e_sensitivity_unmatched.R` |
| Mental health outcomes | `07_statistics/07b_mental_health_outcomes.R` |
| MRI acquisition and connectome construction | UK Biobank pipeline and Connectome Resource; no custom code |
| — Functional connectomes | `02_construction/02a_fc_construction.py` |
| — Structural connectomes | `02_construction/02b_sc_prep.py` |
| — Residual motion benchmarking (Supplementary Figs 4–6, Supplementary Table 7) | `02_construction/02c_qc_fd_distribution.py`, `02d_qc_fc_distance.py` |
| Network-based statistics | `03_nbs/03a_build_design.R` → `03b_run_nbs.m`, `03c_run_nbs.sh`, `03d_submit_nbs_fc.sh`, `03e_submit_nbs_sc.sh` |
| — Functional NBS and hub analysis (Figures 1–2, Supplementary Tables 8–9) | `04_subnetwork/04a_extract_edges.m`, `04b_eigenvector_centrality.m`, `04e_harvard_oxford_labels.py` |
| — Threshold stability and sensitivity (Supplementary Tables 5a, 5b, 6) | `04_subnetwork/04c_threshold_summary.m`; high-motion cohort via `03a` (`FD_EXCLUDE_MM`) and `03d` (`SUBJECTS_TXT`) |
| — Within versus outside the subnetwork (Supplementary Table 4) | `04_subnetwork/04d_heldout_validation.py --mode primary` |
| — Structural NBS | `03_nbs/03a`, `03b`, `03c`, `03e` |
| Network segregation and integration | `05_segregation_integration/05a_segregation_functional.py`, `05b_participation_functional.py`, `05c_segregation_participation_structural.py` |
| — Group models for those metrics (Figure 3) | `08_figures/08a_figure3_dedifferentiation.R` |
| Functional deviation score | `06_deviation_scores/06a_fc_deviation_score.py` |
| Moderation models; specificity of the moderation; structure–function interaction (Figure 4, Supplementary Tables 11, 14) | `07_statistics/07c_moderation_models.R` |
| Hierarchical regression (Supplementary Table 12) | `07_statistics/07c_moderation_models.R` |
| Nested cross-validation (Supplementary Table 13) | `07_statistics/07d_nested_cross_validation.py` |
| White matter hyperintensity location sensitivity (Supplementary Tables 2, 5a, 10, 11, 12) | `07_statistics/07f_sensitivity_cognition_wmh.R`, `07h_sensitivity_dediff_wmh.R`, `07c_moderation_models.R`, `03a`/`03d`/`03e` (deep and periventricular designs) |
| Internal held-out validation of the functional subnetwork (Supplementary Fig. 7) | `04_subnetwork/04d_heldout_validation.py --mode heldout` |

### Figures

| Figure | Script |
|---|---|
| Figure 1, subnetwork chord diagram | `08_figures/08b_circular_connectivity_figure.R` |
| Figure 2, hub rendering | BrainNet Viewer, from `04_subnetwork/04b_eigenvector_centrality.m` and `04e_harvard_oxford_labels.py` |
| Figure 3, dedifferentiation | `08_figures/08a_figure3_dedifferentiation.R` |
| Figure 4, connectivity–processing speed interaction | `07_statistics/07c_moderation_models.R` |
| Supplementary Fig. 1, participant flow | diagram drawn from the counts printed by `01_cohort_assembly.R` (`flow_counts.csv`, `exclusion_summary.csv`) |

## Execution order

```
01_cohort/01_cohort_assembly.R                       cohort, matching, held-out and unmatched pools
01_cohort/01b_sample_characteristics.R               Tables 1-2, Supplementary Fig. 2
02_construction/02a_fc_construction.py               functional connectomes (primary + held-out)
02_construction/02b_sc_prep.py                       structural connectomes
02_construction/02c_qc_fd_distribution.py            head motion (Supplementary Fig. 4)
02_construction/02d_qc_fc_distance.py                QC-FC, network blocks, distance dependence
03_nbs/03a_build_design.R                            design matrices and contrasts
03_nbs/03d_submit_nbs_fc.sh                          functional NBS (SLURM)
03_nbs/03e_submit_nbs_sc.sh                          structural NBS (SLURM)
04_subnetwork/04a_extract_edges.m                    NBS output to edge lists
04_subnetwork/04b_eigenvector_centrality.m           hub analysis
04_subnetwork/04e_harvard_oxford_labels.py           anatomical labels of the hubs
04_subnetwork/04c_threshold_summary.m                threshold stability
04_subnetwork/04d_heldout_validation.py              within/outside subnetwork; held-out validation
05_segregation_integration/05a_segregation_functional.py
05_segregation_integration/05b_participation_functional.py
05_segregation_integration/05c_segregation_participation_structural.py
06_deviation_scores/06a_fc_deviation_score.py
07_statistics/07a_cognitive_outcomes.R
07_statistics/07b_mental_health_outcomes.R
07_statistics/07c_moderation_models.R
07_statistics/07d_nested_cross_validation.py
07_statistics/07e_sensitivity_unmatched.R
07_statistics/07f_sensitivity_cognition_wmh.R
07_statistics/07h_sensitivity_dediff_wmh.R
08_figures/08a_figure3_dedifferentiation.R
08_figures/08b_circular_connectivity_figure.R
```

Steps `03d` and `03e` are SLURM-based and produce a nested tree of result
files; the `04` scripts walk that tree. The sensitivity analyses with
alternative WMH covariates and with the high-motion exclusion re-run the same
scripts with the alternative designs written by `03a`.

## Data flow

`01_cohort_assembly.R` writes three subject tables into `DERIV_DIR`:

| File | Content |
|---|---|
| `cohort_matched.csv` | primary matched sample (complete cognitive testing); raw and standardised covariates (`age_std`, `sex`, `education_std`, `dia_bp_std`, `log_total_wmh_std`, `log_peri_wmh_std`, `log_deep_wmh_std`, `mean_fd_std`) |
| `cohort_heldout_matched.csv` | held-out matched sample (incomplete cognitive testing) |
| `cohort_unmatched_cogcomplete.csv` | all eligible participants with complete cognitive testing, not connectome-restricted, unmatched |

Later steps add one file each, keyed on `eid`: `cohort_matched_cognition.csv`
(07a; z-scored tests and domain scores), `fc_deviation.csv` (06a),
`segregation_functional.csv` (05a), `participation_functional.csv` (05b) and
`segregation_participation_structural.csv` (05c). The brain–behaviour scripts
(07c, 07d, 07h, 08a) merge these tables.

## Configuration

Input and output locations are not hard-coded. The Python analysis scripts
take them as arguments:

```bash
python 05_segregation_integration/05a_segregation_functional.py \
    --subjects  derivatives/cohort_matched.csv \
    --matrix-dir <with_gsr/primary> \
    --labels     data/atlas/labels_216.txt \
    --out        derivatives/segregation_functional.csv
```

R scripts read `DERIV_DIR` and `OUTPUT_DIR` from the environment, defaulting to
`derivatives/` and `output/` relative to the repository root:

```bash
DERIV_DIR=/path/to/derivatives OUTPUT_DIR=/path/to/output \
    Rscript 07_statistics/07a_cognitive_outcomes.R
```

The data-preparation scripts (`01_cohort_assembly.R`, `02a`–`02d`,
`03a`–`03e`, `04a`–`04c`, `08b`) use a configuration block of placeholder
paths near the top. Edit that block once per script.

## Scope

The repository covers the analyses reported in the manuscript. Steps that
are performed with third-party graphical tools are documented but not
scripted: the NBS permutation test itself runs in NBS-Connectome v1.2 (called
from `03b`), and Figure 2 is rendered in BrainNet Viewer from the hub tables.
The participant flow diagram (Supplementary Fig. 1) is drawn from the counts
written by `01_cohort_assembly.R`. Supplementary Table 1 lists one criterion
("active malignancy, brain-affecting") as "various C-codes"; the code
prefixes for that criterion are entered in `MALIGNANCY_ICD10` in
`01_cohort_assembly.R`.

## Software dependencies

### Python

See `environment.yml`. Tested with Python 3.10–3.11.

- `numpy`, `pandas`, `scipy`, `statsmodels`
- `nilearn` (time-series filtering; atlas fetching for the anatomical labels)
- `matplotlib`, `seaborn`

### R

See `requirements_r.md`. Tested with R 4.3 or later.

- `dplyr`, `readr`, `tidyr`, `purrr`, `stringr`, `tibble`, `forcats`
- `MatchIt`, `mice`, `cobalt`
- `car`, `emmeans`, `sandwich`, `lmtest`
- `ggplot2`, `patchwork`, `svglite`, `circlize`

### MATLAB

- MATLAB R2016b or later
- [NBS-Connectome v1.2](https://www.nitrc.org/projects/nbs/)
- [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/)
- BrainNet Viewer (visualisation only)

### SLURM

Steps `03d` and `03e` submit batch jobs. Edit the partition, time, CPU and
memory directives at the top of each submission script for your cluster.

## Data access

Researchers wishing to reproduce these analyses should:

1. Apply for UK Biobank access at <https://www.ukbiobank.ac.uk/>.
2. Request demographics, cognition and lifestyle measures (instance 2); the
   imaging-visit ECG report text (field 12653); the parcellated functional
   time series and structural connectomes of the UK Biobank Connectome
   Resource (fields 31017–31019 and 31024); the head-motion summary (field
   25741); white matter hyperintensity volumes (fields 25781, 24485, 24486)
   and total intracranial volume (field 26521); self-reported illness (field
   20002); and ICD-coded hospital episode statistics (fields 41270, 41202).
3. Point each script at the resulting extract using the configuration
   described above.

## License

MIT — see `LICENSE`.

## Contact

Aslı Akdeniz, Max Planck Institute for Human Cognitive and Brain Sciences,
Leipzig. Questions about the code can be raised through the repository's issue
tracker.
