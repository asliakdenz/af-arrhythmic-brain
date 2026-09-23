#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05c_segregation_participation_structural.py
===========================================
Whole-brain structural system segregation and nodal participation coefficient.

Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
          in atrial fibrillation" (Akdeniz et al.)
Methods: Network segregation and integration ("To test whether any group
         difference was specific to functional connectivity, identical
         analyses were applied to the structural connectomes, with system
         segregation and participation computed from each participant's
         SIFT2-weighted matrix. Structural edge weights are non-negative, so
         no negative-edge handling was required.")

The same two measures as 05a and 05b, computed from the raw (unlogged)
length-normalised SIFT2 fibre bundle capacity matrices written by
02b_sc_prep.py:

    system segregation   SyS = (Zw - Zb) / Zw
    participation        P_i = 1 - sum_m (k_i(m) / k_i)^2, averaged over nodes

Networks are the seven canonical Schaefer systems (cortical, primary) and,
additionally, all regions with the subcortex as an eighth module.

Inputs:
    - Subject table with participant identifiers and group membership
    - Per-subject structural matrices (216 x 216, raw length-normalised
      SIFT2 weights), one plain-text file per subject
    - 216-node atlas label file

Output:
    - One row per subject: sys_cort, zw_cort, zb_cort, pc_mean_cort,
      sys_all, zw_all, zb_all, pc_mean_wb
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

N_NODES = 216


def get_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--subjects", required=True, type=Path,
                   help="CSV with one row per subject (e.g. cohort_matched.csv).")
    p.add_argument("--matrix-dir", required=True, type=Path,
                   help="Directory of per-subject structural matrices "
                        "(sift2_fbc_lengthnorm_raw from 02b_sc_prep.py).")
    p.add_argument("--matrix-pattern", default="sc_{sid}.txt",
                   help="Filename pattern; '{sid}' is replaced by the identifier.")
    p.add_argument("--labels", required=True, type=Path,
                   help="Text file of 216 atlas labels, one per line.")
    p.add_argument("--out", required=True, type=Path, help="Destination CSV.")
    p.add_argument("--id-col", default="eid", help="Identifier column name.")
    p.add_argument("--group-col", default="af",
                   help="Binary group column name (1 = case, 0 = control).")
    return p.parse_args()


def network_of(label: str) -> str:
    parts = label.split("_")
    return parts[2] if label.startswith("7Networks") and len(parts) > 2 \
        else "Subcortex"


def load_atlas(path: Path):
    labels = [line.strip() for line in path.read_text().splitlines()
              if line.strip()]
    if len(labels) != N_NODES:
        raise ValueError(f"expected {N_NODES} labels, found {len(labels)}")
    net = np.array([network_of(l) for l in labels])
    cortical = np.array([l.startswith("7Networks") for l in labels])
    return net, cortical


def segregation(W, net_vec, node_mask):
    idx = np.where(node_mask)[0]
    Wm = W[np.ix_(idx, idx)].astype(float).copy()
    Wm[Wm < 0] = 0.0                     # no-op for non-negative structural weights
    iu = np.triu_indices(len(idx), 1)
    weights = Wm[iu]
    within = net_vec[idx][iu[0]] == net_vec[idx][iu[1]]
    zw = weights[within].mean()
    zb = weights[~within].mean()
    return ((zw - zb) / zw if zw != 0 else np.nan), zw, zb


def participation(W, net_vec, node_mask):
    idx = np.where(node_mask)[0]
    Wm = W[np.ix_(idx, idx)].astype(float).copy()
    Wm[Wm < 0] = 0.0
    np.fill_diagonal(Wm, 0.0)
    modules = net_vec[idx]
    present = np.unique(modules)
    strength = Wm.sum(axis=1)
    pc = np.full(len(idx), np.nan)
    for node in range(len(idx)):
        k = strength[node]
        if k == 0:
            continue
        share = sum((Wm[node, modules == m].sum() / k) ** 2 for m in present)
        pc[node] = 1.0 - share
    return pc


def main():
    args = get_args()
    net, cortical = load_atlas(args.labels)
    all_nodes = np.ones(N_NODES, bool)

    subjects = pd.read_csv(args.subjects, dtype={args.id_col: str})
    rows, missing = [], 0
    for record in subjects.itertuples(index=False):
        sid = getattr(record, args.id_col)
        path = args.matrix_dir / args.matrix_pattern.format(sid=sid)
        if not path.exists():
            missing += 1
            continue
        W = np.loadtxt(path)
        if W.shape != (N_NODES, N_NODES):
            raise ValueError(f"{path.name}: expected {N_NODES} x {N_NODES}, found {W.shape}")
        if (W < 0).any() or not np.isfinite(W).all():
            raise ValueError(f"{path.name}: structural weights must be finite and non-negative")

        row = {args.id_col: sid, args.group_col: int(getattr(record, args.group_col))}
        row["sys_cort"], row["zw_cort"], row["zb_cort"] = segregation(W, net, cortical)
        row["sys_all"],  row["zw_all"],  row["zb_all"]  = segregation(W, net, all_nodes)
        row["pc_mean_cort"] = float(np.nanmean(participation(W, net, cortical)))
        row["pc_mean_wb"]   = float(np.nanmean(participation(W, net, all_nodes)))
        rows.append(row)

    if missing:
        print(f"[warn] {missing} subjects had no matrix file and were skipped")
    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[saved] {args.out.name}  ({len(out)} subjects)")


if __name__ == "__main__":
    main()
