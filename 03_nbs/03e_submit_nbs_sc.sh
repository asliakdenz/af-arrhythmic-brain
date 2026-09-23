#!/bin/bash
# ==============================================================================
# 03e_submit_nbs_sc.sh
#
# SLURM submission for the structural network-based statistics. One job per
# (measure x design x contrast x threshold) combination.
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Network-based statistics -> Structural NBS; White matter
#          hyperintensity location sensitivity
#
# Measures (02b):  sift2_fbc_lengthnorm_log1p  primary: log(1+x) SIFT2 fibre
#                                              bundle capacity / mean length
#                  streamline_count_log1p      sensitivity
#                  mean_FA_raw                 sensitivity
# Designs (03a):   design_1_total (primary), design_2_peri, design_3_deep
# Contrasts:       AF > Control, Control > AF
# Thresholds:      t = 1.5, 2.0, 2.5 (lower edge-level variance of
#                  log-transformed SIFT2 capacity than Fisher-z correlations)
# Component size:  Extent
#
# Total jobs: 3 measures x 3 designs x 2 contrasts x 3 thresholds = 54
# ==============================================================================
set -euo pipefail

# ------------------------------------------------------------------------------
# Configuration  --- edit paths for your environment
# ------------------------------------------------------------------------------
BASE_NBS="/path/to/nbs"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNNER="${SCRIPT_DIR}/03c_run_nbs.sh"
MATS_ROOT="${BASE_NBS}/00_input/structural_inputs"
DESIGN_ROOT="${BASE_NBS}/01_designs/sc"
RESULTS_BASE="${BASE_NBS}/02_results/sc"

MEASURES=(sift2_fbc_lengthnorm_log1p streamline_count_log1p mean_FA_raw)
DESIGNS=(design_1_total design_2_peri design_3_deep)
THRESH_LIST=("1.5" "2" "2.5")
SIZES=("Extent")
METHOD="Run NBS"

# SLURM settings
PART=standard
TIME="1-00:00:00"
CPUS=2
MEM="48G"

chmod +x "$RUNNER"
echo "[INFO] Structural NBS submission"

TOTAL=0
for MEASURE in "${MEASURES[@]}"; do
  MATS_DIR="${MATS_ROOT}/${MEASURE}"
  [[ -d "${MATS_DIR}" ]] || { echo "[WARN] Missing measure folder: ${MATS_DIR}; skipping"; continue; }
  read -r MATS_FILE < <(find "${MATS_DIR}" -maxdepth 1 -type f -name "subject*.txt" | sort | head -n 1)
  [[ -n "${MATS_FILE:-}" ]] || { echo "[WARN] No subject*.txt in ${MATS_DIR}; skipping"; continue; }
  N_MATS=$(find "${MATS_DIR}" -maxdepth 1 -type f -name "subject*.txt" | wc -l | tr -d ' ')
  echo "[INFO] MEASURE: ${MEASURE}  (${N_MATS} matrices)"

  for DESIGN_NAME in "${DESIGNS[@]}"; do
    DESIGN_DIR="${DESIGN_ROOT}/${DESIGN_NAME}"
    DESIGN_TXT="${DESIGN_DIR}/design.txt"
    CONTRAST_DIR="${DESIGN_DIR}/contrasts"
    [[ -f "${DESIGN_TXT}" ]] || { echo "[FATAL] Missing design file: ${DESIGN_TXT}" >&2; exit 1; }
    N_ROWS=$(wc -l < "${DESIGN_TXT}" | tr -d ' ')
    [[ "${N_ROWS}" == "${N_MATS}" ]] || { echo "[FATAL] ${DESIGN_NAME}: ${N_ROWS} design rows but ${N_MATS} matrix files" >&2; exit 1; }
    mapfile -t CONTRAST_FILES < <(find "${CONTRAST_DIR}" -maxdepth 1 -type f -name "*.txt" | sort)
    echo "  DESIGN: ${DESIGN_NAME}  (${#CONTRAST_FILES[@]} contrasts, ${N_ROWS} rows)"

    for CONTRAST_FILE in "${CONTRAST_FILES[@]}"; do
      CONTRAST_KEY=$(basename "${CONTRAST_FILE}" .txt)
      for THRESH in "${THRESH_LIST[@]}"; do
        for SIZE in "${SIZES[@]}"; do
          OUT_DIR="${RESULTS_BASE}/${MEASURE}/${DESIGN_NAME}/${CONTRAST_KEY}/T${THRESH}/${SIZE}"
          mkdir -p "${OUT_DIR}"
          JOB_NAME="NBSsc_${MEASURE:0:8}_${DESIGN_NAME#design_}_${CONTRAST_KEY}_T${THRESH}_${SIZE:0:3}"
          echo "    Submitting: ${JOB_NAME}"
          sbatch <<EOS
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=${OUT_DIR}/slurm_%j.out
#SBATCH --error=${OUT_DIR}/slurm_%j.err
#SBATCH --time=${TIME}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --partition=${PART}

"${RUNNER}" "${OUT_DIR}" "${MATS_FILE}" "${DESIGN_TXT}" "${CONTRAST_FILE}" "${THRESH}" "${SIZE}" "${METHOD}"
EOS
          TOTAL=$((TOTAL + 1))
        done
      done
    done
  done
done
echo "[DONE] Submitted ${TOTAL} structural NBS jobs. Results in ${RESULTS_BASE}"
