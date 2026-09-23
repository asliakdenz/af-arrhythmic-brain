#!/bin/bash
# ==============================================================================
# 03d_submit_nbs_fc.sh
#
# SLURM submission for the functional network-based statistics. One job per
# (design x contrast x threshold) combination.
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Network-based statistics -> Functional NBS and hub analysis;
#          White matter hyperintensity location sensitivity; Supplementary
#          material, "Motion-related Quality Control" (v)
#
# Designs   (03a): design_1_total (primary), design_2_peri, design_3_deep
# Contrasts:       AF > Control, Control > AF
# Thresholds:      t = 4.5 (primary), 5.0, 5.5, 6.0
# Component size:  Extent
# Permutations:    5,000 (03b)
#
# High-motion sensitivity cohort: point DESIGN_ROOT at the *_sens designs
# written by 03a with FD_EXCLUDE_MM = 0.30 and set SUBJECTS_TXT to the
# subjects.txt of a design; a matrix directory containing only those
# participants is then assembled from symbolic links, because NBS reads every
# file of a directory.
#
# Total jobs: 3 designs x 2 contrasts x 4 thresholds = 24
# ==============================================================================
set -euo pipefail

# ------------------------------------------------------------------------------
# Configuration  --- edit paths for your environment
# ------------------------------------------------------------------------------
BASE_NBS="/path/to/nbs"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNNER="${SCRIPT_DIR}/03c_run_nbs.sh"
MATS_DIR="${BASE_NBS}/00_input/matrices/with_gsr/primary"
DESIGN_ROOT="${BASE_NBS}/01_designs/fc"
RESULTS_BASE="${BASE_NBS}/02_results/fc"
SUBJECTS_TXT=""                      # optional: e.g. ${BASE_NBS}/01_designs_sens/fc/design_1_total/subjects.txt

DESIGNS=(design_1_total design_2_peri design_3_deep)
THRESH_LIST=("4.5" "5" "5.5" "6")
SIZES=("Extent")
METHOD="Run NBS"

# SLURM settings
PART=standard
TIME="1-00:00:00"
CPUS=2
MEM="48G"

# ------------------------------------------------------------------------------
# Matrix directory (optionally restricted to a subject list)
# ------------------------------------------------------------------------------
chmod +x "$RUNNER"
if [[ -n "${SUBJECTS_TXT}" ]]; then
  SUBSET_DIR="${MATS_DIR}_subset_$(basename "$(dirname "$(dirname "${SUBJECTS_TXT}")")")"
  rm -rf "${SUBSET_DIR}"; mkdir -p "${SUBSET_DIR}"
  while read -r eid; do
    [[ -z "${eid}" ]] && continue
    [[ -f "${MATS_DIR}/conn_${eid}.txt" ]] || { echo "[FATAL] missing ${MATS_DIR}/conn_${eid}.txt" >&2; exit 1; }
    ln -s "${MATS_DIR}/conn_${eid}.txt" "${SUBSET_DIR}/conn_${eid}.txt"
  done < "${SUBJECTS_TXT}"
  MATS_DIR="${SUBSET_DIR}"
  echo "[INFO] Subset matrix directory: ${MATS_DIR} ($(ls "${MATS_DIR}" | wc -l) files)"
fi

read -r MATS_FILE < <(find "${MATS_DIR}" -maxdepth 1 \( -type f -o -type l \) -name "conn_*.txt" | sort | head -n 1)
[[ -n "${MATS_FILE:-}" ]] || { echo "[FATAL] No conn_*.txt files in ${MATS_DIR}" >&2; exit 1; }
N_MATS=$(find "${MATS_DIR}" -maxdepth 1 \( -type f -o -type l \) -name "conn_*.txt" | wc -l | tr -d ' ')
echo "[INFO] Functional NBS submission: ${N_MATS} matrices in ${MATS_DIR}"

# ------------------------------------------------------------------------------
# Submission loop
# ------------------------------------------------------------------------------
TOTAL=0
for DESIGN_NAME in "${DESIGNS[@]}"; do
  DESIGN_DIR="${DESIGN_ROOT}/${DESIGN_NAME}"
  DESIGN_TXT="${DESIGN_DIR}/design.txt"
  CONTRAST_DIR="${DESIGN_DIR}/contrasts"
  [[ -f "${DESIGN_TXT}" ]] || { echo "[FATAL] Missing design file: ${DESIGN_TXT}" >&2; exit 1; }
  N_ROWS=$(wc -l < "${DESIGN_TXT}" | tr -d ' ')
  [[ "${N_ROWS}" == "${N_MATS}" ]] || { echo "[FATAL] ${DESIGN_NAME}: ${N_ROWS} design rows but ${N_MATS} matrix files" >&2; exit 1; }
  mapfile -t CONTRAST_FILES < <(find "${CONTRAST_DIR}" -maxdepth 1 -type f -name "*.txt" | sort)
  echo "[INFO] DESIGN: ${DESIGN_NAME}  (${#CONTRAST_FILES[@]} contrasts, ${N_ROWS} rows)"

  for CONTRAST_FILE in "${CONTRAST_FILES[@]}"; do
    CONTRAST_KEY=$(basename "${CONTRAST_FILE}" .txt)
    for THRESH in "${THRESH_LIST[@]}"; do
      for SIZE in "${SIZES[@]}"; do
        OUT_DIR="${RESULTS_BASE}/${DESIGN_NAME}/${CONTRAST_KEY}/T${THRESH}/${SIZE}"
        mkdir -p "${OUT_DIR}"
        JOB_NAME="NBSfc_${DESIGN_NAME#design_}_${CONTRAST_KEY}_T${THRESH}_${SIZE:0:3}"
        echo "  Submitting: ${JOB_NAME}"
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
echo "[DONE] Submitted ${TOTAL} functional NBS jobs. Results in ${RESULTS_BASE}"
