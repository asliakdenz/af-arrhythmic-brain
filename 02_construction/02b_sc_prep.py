#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02b_sc_prep.py
==============
Structural connectome preparation from the UK Biobank Connectome Resource
(field 31024, Mansour et al. 2023): edge weighting, numerical clean-up, and
export in the per-subject text format read by NBS-Connectome v1.2 and by the
graph-metric script 05c.

Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
          in atrial fibrillation" (Akdeniz et al.)
Methods: MRI acquisition and connectome construction -> Structural
         connectomes; Network-based statistics -> Structural NBS

Edge weights are SIFT2 fibre bundle capacity normalised by mean bundle
length (edge-wise division; edges without a length estimate are set to
zero). Matrices are processed for numerical stability: non-finite and
negative values zeroed, diagonal zeroed, symmetry enforced. A log(1 + x)
transformation is applied to the NBS inputs; raw (unlogged) weights are
written for graph-metric extraction. Sensitivity metrics for the structural
NBS are log(1 + x)-transformed streamline count and raw mean fractional
anisotropy.

Inputs:
    - cohort_matched.csv                          (01_cohort_assembly.R)
    - UKB field 31024 per-subject CSVs at
        {STRUCT_DIR}/{eid}/Baseline/{eid}_31024_2_0/
            connectome_sift2_fbc_10M.csv
            connectome_mean_length_10M.csv
            connectome_streamline_count_10M.csv
            connectome_mean_FA_10M.csv

Outputs ({OUTPUT_DIR}):
    sift2_fbc_lengthnorm_log1p/subjectNNN.txt   primary NBS input
    streamline_count_log1p/subjectNNN.txt       sensitivity NBS input
    mean_FA_raw/subjectNNN.txt                  sensitivity NBS input
    sift2_fbc_lengthnorm_raw/sc_{eid}.txt       graph-metric input (05c)
    subject_list.csv                            design row -> filename -> eid
    qc_matrix_availability.csv

NBS-Connectome reads every file in a directory in alphabetical order, so the
kth subjectNNN.txt must correspond to the kth row of the design matrix.
Both this script and 03a_build_design.R order participants by eid
(character sort). The script stops if any participant of the cohort lacks a
matrix, because a silent omission would misalign files and design rows.
"""

from pathlib import Path

import numpy as np
import pandas as pd

# =============================================================================
# Configuration  --- edit paths for your environment
# =============================================================================
COHORT_CSV = Path("/path/to/derivatives/cohort_matched.csv")
STRUCT_DIR = Path("/path/to/ukb_extract/31024")
OUTPUT_DIR = Path("/path/to/nbs/00_input/structural_inputs")

N_NODES = 216
FILES = {
    "sift2_fbc":        "connectome_sift2_fbc_10M.csv",
    "mean_length":      "connectome_mean_length_10M.csv",
    "streamline_count": "connectome_streamline_count_10M.csv",
    "mean_FA":          "connectome_mean_FA_10M.csv",
}
NBS_FMT = "%.6f"


# =============================================================================
# Helpers
# =============================================================================
def load_matrix(eid: str, filename: str) -> np.ndarray:
    direct = STRUCT_DIR / eid / "Baseline" / f"{eid}_31024_2_0" / filename
    path = direct if direct.exists() else next(iter(STRUCT_DIR.glob(f"{eid}/**/{filename}")), None)
    if path is None:
        raise FileNotFoundError(f"{eid}: {filename} not found under {STRUCT_DIR}")
    mat = pd.read_csv(path, header=None).values.astype(float)
    if mat.shape != (N_NODES, N_NODES):
        raise ValueError(f"{eid} {filename}: shape {mat.shape}, expected ({N_NODES}, {N_NODES})")
    return mat


def stabilise(mat: np.ndarray) -> np.ndarray:
    """Non-finite and negative values zeroed, diagonal zeroed, symmetry enforced."""
    mat = mat.copy()
    mat[~np.isfinite(mat)] = 0.0
    mat[mat < 0] = 0.0
    np.fill_diagonal(mat, 0.0)
    mat = (mat + mat.T) / 2.0
    return mat


def length_normalise(fbc: np.ndarray, length: np.ndarray) -> np.ndarray:
    out = np.zeros_like(fbc)
    ok = np.isfinite(length) & (length > 0)
    out[ok] = fbc[ok] / length[ok]
    return out


# =============================================================================
# 1. Cohort, ordered by eid (character sort, as in 03a_build_design.R)
# =============================================================================
cohort = pd.read_csv(COHORT_CSV, dtype={"eid": str})
cohort = cohort.sort_values("eid").reset_index(drop=True)
eids = list(cohort["eid"])
pad = len(str(len(eids)))
print(f"Cohort: N = {len(eids)} (AF = {int(cohort['af'].sum())}, "
      f"Control = {int((cohort['af'] == 0).sum())})")

folders = {
    "sift2_fbc_lengthnorm_log1p": OUTPUT_DIR / "sift2_fbc_lengthnorm_log1p",
    "streamline_count_log1p":     OUTPUT_DIR / "streamline_count_log1p",
    "mean_FA_raw":                OUTPUT_DIR / "mean_FA_raw",
    "sift2_fbc_lengthnorm_raw":   OUTPUT_DIR / "sift2_fbc_lengthnorm_raw",
}
for d in folders.values():
    d.mkdir(parents=True, exist_ok=True)

# =============================================================================
# 2. Per-subject processing and export
# =============================================================================
subject_log, qc_rows = [], []
for k, eid in enumerate(eids, start=1):
    try:
        fbc    = load_matrix(eid, FILES["sift2_fbc"])
        length = load_matrix(eid, FILES["mean_length"])
        count  = load_matrix(eid, FILES["streamline_count"])
        fa     = load_matrix(eid, FILES["mean_FA"])
    except (FileNotFoundError, ValueError) as e:
        raise SystemExit(f"Missing or malformed input; the design matrix would be misaligned. {e}")

    w_primary = stabilise(length_normalise(fbc, length))     # SIFT2 FBC / mean length
    w_count   = stabilise(count)
    w_fa      = stabilise(fa)

    fname = f"subject{k:0{pad}d}.txt"
    np.savetxt(folders["sift2_fbc_lengthnorm_log1p"] / fname, np.log1p(w_primary), fmt=NBS_FMT)
    np.savetxt(folders["streamline_count_log1p"]     / fname, np.log1p(w_count),   fmt=NBS_FMT)
    np.savetxt(folders["mean_FA_raw"]                / fname, w_fa,                fmt=NBS_FMT)
    np.savetxt(folders["sift2_fbc_lengthnorm_raw"]   / f"sc_{eid}.txt", w_primary,  fmt=NBS_FMT)

    subject_log.append({"row_in_design_matrix": k, "filename": fname, "eid": eid,
                        "af": int(cohort.loc[k - 1, "af"])})
    qc_rows.append({"eid": eid, "n_nonzero_edges": int((np.triu(w_primary, 1) > 0).sum()),
                    "max_weight": float(w_primary.max()),
                    "n_edges_without_length": int(((np.triu(fbc, 1) > 0) & ~(np.triu(length, 1) > 0)).sum())})

pd.DataFrame(subject_log).to_csv(OUTPUT_DIR / "subject_list.csv", index=False)
pd.DataFrame(qc_rows).to_csv(OUTPUT_DIR / "qc_matrix_availability.csv", index=False)

# =============================================================================
# 3. Verification: read back the first file
# =============================================================================
first = np.loadtxt(folders["sift2_fbc_lengthnorm_log1p"] / f"subject{1:0{pad}d}.txt")
assert first.shape == (N_NODES, N_NODES)
assert np.allclose(first, first.T, atol=1e-5) and np.all(np.diag(first) == 0)
assert np.isfinite(first).all() and (first >= 0).all()
print(f"Wrote {len(eids)} subjects x 4 matrix sets to {OUTPUT_DIR}")
print("subject_list.csv maps design rows to filenames and eids.")
