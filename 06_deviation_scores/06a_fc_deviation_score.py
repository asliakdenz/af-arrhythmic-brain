#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06a_fc_deviation_score.py
=========================
Functional deviation score: a per-participant, control-referenced index of
connectivity magnitude within the NBS-identified functional subnetwork.

Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
          in atrial fibrillation" (Akdeniz et al.)
Methods: Brain-behaviour analyses -> Functional deviation score

For each subnetwork edge, an OLS model is fitted in matched controls only,
regressing the absolute Fisher-z edge weight on the core covariates plus mean
framewise displacement. Absolute values are used so that the score reflects
connectivity magnitude irrespective of direction. Covariate-adjusted
residuals are computed for all participants from the control-derived
coefficients and standardised by the control residual SD, yielding per-edge
z-scores. The functional deviation score is the mean z-score across edges: a
control-referenced deviation score in the tradition of normative modelling,
implemented as a linear W-score.

Inputs:
    --subjects     cohort_matched.csv (01_cohort_assembly.R): identifier,
                   group and standardised covariates
    --matrix-dir   per-subject functional connectivity matrices with GSR
                   (02a_fc_construction.py, primary cohort directory)
    --edges        significant-edge CSV of the primary subnetwork
                   (04a_extract_edges.m, t = 4.5, Extent)
    --labels       216-node atlas label file
Output:
    --out          CSV with eid, af, fc_deviation
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

N_NODES = 216
COVARIATES = ["age_std", "sex", "education_std", "log_total_wmh_std", "dia_bp_std", "mean_fd_std"]


def get_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--subjects", required=True, type=Path)
    p.add_argument("--matrix-dir", required=True, type=Path)
    p.add_argument("--matrix-pattern", default="conn_{sid}.txt")
    p.add_argument("--edges", required=True, type=Path)
    p.add_argument("--labels", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--id-col", default="eid")
    p.add_argument("--group-col", default="af")
    return p.parse_args()


def main():
    args = get_args()
    labels = [l.strip() for l in args.labels.read_text().splitlines() if l.strip()]
    if len(labels) != N_NODES:
        raise ValueError(f"expected {N_NODES} labels, found {len(labels)}")

    edges = pd.read_csv(args.edges)
    if "Region1" not in edges.columns:
        edges = pd.read_csv(args.edges, sep="\t")
    edge_idx = [(labels.index(a), labels.index(b)) for a, b in zip(edges["Region1"], edges["Region2"])]
    print(f"[mask] {len(edge_idx)} subnetwork edges")

    df = pd.read_csv(args.subjects, dtype={args.id_col: str})
    df = df.dropna(subset=[args.group_col] + COVARIATES).reset_index(drop=True)

    Y = np.full((len(df), len(edge_idx)), np.nan)
    for k, sid in enumerate(df[args.id_col].values):
        path = args.matrix_dir / args.matrix_pattern.format(sid=sid)
        if not path.exists():
            continue
        W = np.abs(np.loadtxt(path))                   # absolute Fisher-z edge weights
        Y[k] = [W[i, j] for i, j in edge_idx]
    have = np.isfinite(Y).all(axis=1)
    if (~have).any():
        print(f"[warn] {int((~have).sum())} participants without a matrix were dropped")
    df, Y = df.loc[have].reset_index(drop=True), Y[have]

    X = np.column_stack([np.ones(len(df))] + [df[c].values.astype(float) for c in COVARIATES])
    ctrl = (df[args.group_col] == 0).values

    B = np.linalg.pinv(X[ctrl].T @ X[ctrl]) @ X[ctrl].T @ Y[ctrl]   # control-only OLS per edge
    resid = Y - X @ B                                                # residuals for everyone
    sd_ctrl = resid[ctrl].std(axis=0, ddof=0)                        # control residual SD per edge
    sd_ctrl[sd_ctrl <= 0] = np.inf
    z = resid / sd_ctrl
    df["fc_deviation"] = z.mean(axis=1)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df[[args.id_col, args.group_col, "fc_deviation"]].to_csv(args.out, index=False)
    print(f"[saved] {args.out}  (n = {len(df)}, controls = {int(ctrl.sum())})")


if __name__ == "__main__":
    main()
