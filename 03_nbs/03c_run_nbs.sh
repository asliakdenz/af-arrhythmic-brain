#!/bin/bash
# ==============================================================================
# 03c_run_nbs.sh
#
# Bash wrapper around the MATLAB NBS runner (03b_run_nbs.m).
# Called from the SLURM submission scripts (03d_submit_nbs_fc.sh,
# 03e_submit_nbs_sc.sh) — one invocation per (design × contrast ×
# threshold × size) combination.
#
# Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
#           in atrial fibrillation" (Akdeniz et al.)
# Methods: Network-based statistics
#
# Usage:
#   bash 03c_run_nbs.sh OUT_DIR MATS_FILE DESIGN_TXT CONTRAST_PATH \
#                       THRESH_STR SIZE_STR METHOD_STR
# ==============================================================================
set -euo pipefail

if [ "$#" -ne 7 ]; then
  echo "Usage: $0 OUT_DIR MATS_FILE DESIGN_TXT CONTRAST_PATH THRESH_STR SIZE_STR METHOD_STR" >&2
  exit 1
fi

OUT_DIR="$1"
MATS_FILE="$2"
DESIGN_TXT="$3"
CONTRAST_PATH="$4"
THRESH_STR="$5"
SIZE_STR="$6"
METHOD_STR="$7"

echo "--- Starting MATLAB NBS Job ---"
echo "Output Dir: ${OUT_DIR}"
echo "Matrix    : ${MATS_FILE}"
echo "Design    : ${DESIGN_TXT}"
echo "Contrast  : ${CONTRAST_PATH}"
echo "Thresh    : ${THRESH_STR}"
echo "Size      : ${SIZE_STR}"
echo "Method    : ${METHOD_STR}"
mkdir -p "${OUT_DIR}"

# Isolate MATLAB preferences per-job so concurrent SLURM jobs don't clash
export MATLAB_PREFDIR="/path/to/matlab_tmp/${SLURM_JOB_ID:-interactive_$$}"
mkdir -p "$MATLAB_PREFDIR"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

matlab -noFigureWindows -nodesktop -nosplash -r " \
try, \
  addpath('${SCRIPT_DIR}'); \
  run_nbs('${OUT_DIR}','${MATS_FILE}','${DESIGN_TXT}','${CONTRAST_PATH}','${THRESH_STR}','${SIZE_STR}','${METHOD_STR}'); \
  exit(0); \
catch ME, \
  fprintf('MATLAB script failed.\n'); \
  disp(getReport(ME,'extended')); \
  exit(1); \
end;"
