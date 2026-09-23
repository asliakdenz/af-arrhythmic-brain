#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02a_fc_construction.py
======================
Functional connectome construction from UK Biobank parcellated rsfMRI
time-series. Produces 216-node (Schaefer 200 cortical + Tian S1 16
subcortical) Fisher-z connectomes with global signal regression applied.

Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
          in atrial fibrillation" (Akdeniz et al.)
Methods: MRI acquisition and connectome construction -> Functional connectomes

Inputs:
    - cohort_matched.csv          Primary matched sample (01_cohort_assembly.R)
    - cohort_heldout_matched.csv  Held-out matched sample (01_cohort_assembly.R)
    - Per-subject parcellated time series:
        Schaefer 200 cortical:    {ROOT_SCHAEFER}/{eid}/Baseline/.../<file>
        Tian S1 16 subcortical:   {ROOT_TIAN}/{eid}/Baseline/.../<file>
        Global signal:            {ROOT_GS}/{eid}/Baseline/.../<file>

Outputs:
    - {OUTPUT_ROOT}/with_gsr/primary/conn_{eid}.txt   primary sample
    - {OUTPUT_ROOT}/with_gsr/heldout/conn_{eid}.txt   held-out sample
    - {OUTPUT_ROOT}/qc_summary_216.csv

The two samples are written to separate directories because NBS-Connectome
reads every file in a directory; the held-out matrices must not enter the
network-based statistics that define the subnetwork.

Pipeline per subject (Methods, "Functional connectomes"):
    1. Load Schaefer 200 + Tian S1 16 parcellated time series (UK Biobank
       fields 31018, 31019) and concatenate (490 timepoints x 216 regions).
    2. Drop all-NaN parcels; column-mean impute residual NaNs.
    3. Low-pass filter at 0.1 Hz (fifth-order zero-phase Butterworth, as
       implemented in nilearn.signal.clean; TR = 0.735 s), complementing the
       highpass temporal filtering of the standard UK Biobank pipeline; no
       further detrending.
    4. Global signal regression: regress the filtered, mean-centred global
       signal (field 31017) out of every regional time series by OLS.
    5. Pearson correlation matrix -> Fisher r-to-z transform.
    6. Write matrix and QC statistics (min/max/mean Fisher-z, fraction of
       negative edges).
"""

import os
import numpy as np
import pandas as pd
from nilearn.connectome import ConnectivityMeasure
from nilearn import signal

# =============================================================================
# Configuration  --- edit paths for your environment
# =============================================================================
# Cohort files (output of 01_cohort_assembly.R); one output directory each
COHORTS = {
    "primary": "/path/to/derivatives/cohort_matched.csv",
    "heldout": "/path/to/derivatives/cohort_heldout_matched.csv",
}

# UK Biobank Connectome Resource project IDs (Mansour et al., 2023)
PID_SCHAEFER = "31018"
PID_TIAN     = "31019"
PID_GS       = "31017"

ROOT_SCHAEFER = f"/path/to/ukb_extract/{PID_SCHAEFER}/merged_output"
ROOT_TIAN     = f"/path/to/ukb_extract/{PID_TIAN}/merged_output"
ROOT_GS       = f"/path/to/ukb_extract/{PID_GS}/merged_output"

NAME_SCHAEFER = "fMRI.Schaefer7n200p.csv.gz"
NAME_TIAN     = "fMRI.Tian_Subcortex_S1_3T.csv.gz"
NAME_GS       = "fMRI.global_signal.csv.gz"

OUTPUT_ROOT     = "/path/to/nbs/00_input/matrices"
OUT_WITH_GSR    = os.path.join(OUTPUT_ROOT, "with_gsr")     # + /<cohort>/conn_{eid}.txt
QC_SUMMARY_FILE = os.path.join(OUTPUT_ROOT, "qc_summary_216.csv")

# Skip subjects whose output matrices already exist (resume-friendly)
SKIP_IF_Z_EXISTS = True

# Acquisition parameters
TR        = 0.735     # UK Biobank rsfMRI repetition time (s)
LOW_PASS  = 0.10      # Hz
HIGH_PASS = None      # already applied by UK Biobank preprocessing

# Expected parcel counts
N_SCHAEFER = 200
N_TIAN     = 16
N_TOTAL    = N_SCHAEFER + N_TIAN


# =============================================================================
# Helpers
# =============================================================================
def fisher_z(R: np.ndarray) -> np.ndarray:
    """Fisher z-transform with diagonal zeroed and edge values clipped."""
    Rz = np.arctanh(np.clip(R, -0.999999, 0.999999))
    np.fill_diagonal(Rz, 0.0)
    return Rz


def tri_stats(M: np.ndarray):
    """Min/max/mean/negative-fraction of upper triangle (excl. diagonal)."""
    iu = np.triu_indices_from(M, 1)
    v  = M[iu]
    return (float(np.nanmin(v)),
            float(np.nanmax(v)),
            float(np.nanmean(v)),
            float((v < 0).mean()))


def load_and_clean_ukb_ts(path: str) -> pd.DataFrame:
    """Load UK Biobank parcellated time series CSV (header row + ROI column);
    return (n_time, n_parcels) numeric DataFrame."""
    df = pd.read_csv(path, header=None)
    df_clean = (df.drop(index=0)
                  .drop(columns=0)
                  .apply(pd.to_numeric, errors="coerce"))
    return df_clean.T


# =============================================================================
# Main batch loop
# =============================================================================
measure = ConnectivityMeasure(kind="correlation", standardize="zscore_sample")
qc_records = []

jobs = []
for cohort_name, sample_csv in COHORTS.items():
    if not os.path.exists(sample_csv):
        print(f"[skip] {cohort_name}: {sample_csv} not found")
        continue
    out_dir = os.path.join(OUT_WITH_GSR, cohort_name)
    os.makedirs(out_dir, exist_ok=True)
    eids = pd.read_csv(sample_csv, dtype=str)["eid"].tolist()
    print(f"{cohort_name}: {len(eids)} subjects -> {out_dir}")
    jobs += [(cohort_name, eid, out_dir) for eid in eids]
print(f"Starting {N_TOTAL}-node FC construction for {len(jobs)} subject entries...")

for cohort_name, eid, out_dir in jobs:
    out_with_gsr = os.path.join(out_dir, f"conn_{eid}.txt")

    if SKIP_IF_Z_EXISTS and os.path.exists(out_with_gsr):
        print(f"  skip   {eid}  (exists)")
        qc_records.append({"cohort": cohort_name, "eid": eid, "status": "skipped_exists"})
        continue

    print(f"  run    {eid}")

    try:
        path_s = os.path.join(ROOT_SCHAEFER, eid, "Baseline",
                              f"{eid}_{PID_SCHAEFER}_2_0", NAME_SCHAEFER)
        path_t = os.path.join(ROOT_TIAN,     eid, "Baseline",
                              f"{eid}_{PID_TIAN}_2_0",     NAME_TIAN)
        path_g = os.path.join(ROOT_GS,       eid, "Baseline",
                              f"{eid}_{PID_GS}_2_0",       NAME_GS)

        if not all(os.path.exists(p) for p in (path_s, path_t, path_g)):
            print(f"         missing input file(s)")
            qc_records.append({"cohort": cohort_name, "eid": eid, "error": "missing_files"})
            continue

        # ---- 1. Load and merge cortical + subcortical time series ----
        ts_s = load_and_clean_ukb_ts(path_s)
        ts_t = load_and_clean_ukb_ts(path_t)

        if ts_s.shape[0] != ts_t.shape[0]:
            raise ValueError(
                f"Time mismatch: Schaefer={ts_s.shape[0]}, Tian={ts_t.shape[0]}")
        if ts_s.shape[1] != N_SCHAEFER:
            print(f"         warning: expected {N_SCHAEFER} Schaefer parcels, "
                  f"got {ts_s.shape[1]}")
        if ts_t.shape[1] != N_TIAN:
            print(f"         warning: expected {N_TIAN} Tian parcels, "
                  f"got {ts_t.shape[1]}")

        ts_combined = pd.concat([ts_s, ts_t], axis=1)
        n_time, n_parcels = ts_combined.shape
        if n_time < 10:
            raise ValueError("Too few timepoints")
        if n_parcels != N_TOTAL:
            print(f"         warning: expected {N_TOTAL} parcels, got {n_parcels}")

        # ---- 2. NaN handling ----
        valid_cols = ~np.all(np.isnan(ts_combined.values), axis=0)
        if not np.all(valid_cols):
            dropped = int((~valid_cols).sum())
            print(f"         dropped {dropped} all-NaN parcel(s)")
            ts_combined = ts_combined.loc[:, valid_cols]
            n_time, n_parcels = ts_combined.shape

        ts_vals = ts_combined.values
        if np.isnan(ts_vals).any():
            col_means = np.nanmean(ts_vals, axis=0)
            nan_idx   = np.where(np.isnan(ts_vals))
            ts_vals[nan_idx] = np.take(col_means, nan_idx[1])

        # ---- 3. Low-pass filter brain time series ----
        ts_filtered = signal.clean(
            ts_vals,
            t_r=TR,
            low_pass=LOW_PASS,
            high_pass=HIGH_PASS,
            detrend=False,
            standardize=False,
        )

        # ---- 4. Global signal regression ----
        gs_df = pd.read_csv(path_g)
        gs_raw = (gs_df.drop(columns=["label_name"]).values.flatten()
                  if "label_name" in gs_df.columns
                  else gs_df.values.flatten())

        if len(gs_raw) != n_time:
            raise ValueError(
                f"GS length ({len(gs_raw)}) != n_time ({n_time})")

        gs_filtered = signal.clean(
            gs_raw.reshape(-1, 1),
            t_r=TR,
            low_pass=LOW_PASS,
            high_pass=HIGH_PASS,
            detrend=False,
            standardize=False,
        ).flatten()
        gs_dm = gs_filtered - gs_filtered.mean()

        if np.allclose(gs_dm, 0):
            ts_gsr = ts_filtered.copy()
        else:
            X = np.column_stack((np.ones_like(gs_dm), gs_dm))
            beta, *_ = np.linalg.lstsq(X, ts_filtered, rcond=None)
            ts_gsr = ts_filtered - X @ beta

        # ---- 5. Connectome with GSR ----
        R_with = measure.fit_transform([ts_gsr])[0]

        # ---- 6. Fisher-z and save ----
        Z_with = fisher_z(R_with)
        np.savetxt(out_with_gsr, Z_with, fmt="%.6f")

        mn, mx, mu, neg = tri_stats(Z_with)
        qc_records.append({
            "cohort":    cohort_name,
            "eid":       eid,
            "status":    "written",
            "n_time":    n_time,
            "n_parcels": n_parcels,
            "z_min":     mn,
            "z_max":     mx,
            "z_mean":    mu,
            "z_negfrac": neg,
        })

    except Exception as e:
        print(f"         FAILED: {e}")
        qc_records.append({"cohort": cohort_name, "eid": eid, "error": str(e)})


# =============================================================================
# Save QC summary
# =============================================================================
pd.DataFrame(qc_records).to_csv(QC_SUMMARY_FILE, index=False)
print(f"\nDone. Outputs in {OUTPUT_ROOT}")
print(f"QC summary: {QC_SUMMARY_FILE}")
