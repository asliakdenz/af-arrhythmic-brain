#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05b_participation_functional.py
===============================
Whole-brain functional nodal participation coefficient.

Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
          in atrial fibrillation" (Akdeniz et al.)
Methods: Network segregation and integration

Weighted participation coefficient (Guimera & Amaral 2005; weighted form as
implemented in the Brain Connectivity Toolbox):

    P_i = 1 - sum_m ( k_i(m) / k_i )^2

    k_i(m) = summed connectivity of node i to network m
    k_i    = total connectivity of node i

Computed per node directly from the raw weighted Fisher z-transformed matrix
and averaged across the parcellation, so no thresholding and no graph
construction is involved. Higher participation indicates broader connectivity
across networks, i.e. lower segregation.

Modules are the seven canonical Schaefer networks. Negative edges are set to
zero, matching the handling used for the segregation index.

Two node sets are produced:
    cortical   200 Schaefer nodes, 7 modules               (primary analysis)
    all-node   216 nodes, subcortex treated as an 8th module

Inputs:
    - Per-subject functional connectivity matrices (216 x 216 Fisher z,
      with global signal regression; output of the FC construction step)
    - Subject table with participant identifiers and group membership
    - 216-node atlas label file

Outputs:
    - One row per subject with the parcellation-averaged participation
      coefficient for each node set.
    - The full nodal participation matrix (subjects x 216), saved alongside.

Group-level inference for these metrics is not performed here. The
covariate-adjusted HC3 models reported in the manuscript are fitted in
08a_figure3_dedifferentiation.R, so that every reported number has a single
point of origin.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

N_NODES = 216


# ==============================================================================
# CONFIGURATION
# ==============================================================================
def get_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--subjects", required=True, type=Path,
                   help="CSV with one row per subject; must contain the "
                        "identifier column and a binary group column.")
    p.add_argument("--matrix-dir", required=True, type=Path,
                   help="Directory of per-subject functional connectivity "
                        "matrices, one plain-text 216 x 216 file per subject.")
    p.add_argument("--matrix-pattern", default="conn_{sid}.txt",
                   help="Filename pattern for a subject matrix; '{sid}' is "
                        "replaced by the identifier. Default: 'conn_{sid}.txt'")
    p.add_argument("--labels", required=True, type=Path,
                   help="Text file of 216 atlas labels, one per line, in "
                        "matrix row order.")
    p.add_argument("--out", required=True, type=Path,
                   help="Destination CSV. The nodal matrix is written "
                        "alongside it with a '_nodal.npy' suffix.")
    p.add_argument("--id-col", default="eid", help="Identifier column name.")
    p.add_argument("--group-col", default="af",
                   help="Binary group column name (1 = case, 0 = control).")
    return p.parse_args()


# ==============================================================================
# ATLAS -> MODULE ASSIGNMENT
# ==============================================================================
def network_of(label: str) -> str:
    """Schaefer '7Networks_LH_Vis_1' -> 'Vis'; anything else -> 'Subcortex'."""
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


# ==============================================================================
# PARTICIPATION COEFFICIENT
# ==============================================================================
def participation(W, net_vec, node_mask, zero_negative=True):
    """
    Weighted participation coefficient per node over the subgraph selected by
    node_mask. Nodes outside the mask are returned as NaN.
    """
    idx = np.where(node_mask)[0]
    Wm = W[np.ix_(idx, idx)].astype(float).copy()
    if zero_negative:
        Wm[Wm < 0] = 0.0
    np.fill_diagonal(Wm, 0.0)

    modules = net_vec[idx]
    present = np.unique(modules)
    strength = Wm.sum(axis=1)

    pc = np.zeros(len(idx))
    for node in range(len(idx)):
        k = strength[node]
        if k == 0:
            continue
        share = sum((Wm[node, modules == m].sum() / k) ** 2 for m in present)
        pc[node] = 1.0 - share

    out = np.full(len(net_vec), np.nan)
    out[idx] = pc
    return out


# ==============================================================================
# MAIN
# ==============================================================================
def main():
    args = get_args()
    net, cortical = load_atlas(args.labels)
    all_nodes = np.ones(N_NODES, bool)

    print(f"[atlas] {cortical.sum()} cortical / {(~cortical).sum()} subcortical")
    print(f"[atlas] modules: {sorted(set(net))}")

    subjects = pd.read_csv(args.subjects, dtype={args.id_col: str})
    rows, nodal, missing = [], [], 0

    for record in subjects.itertuples(index=False):
        sid = getattr(record, args.id_col)
        path = args.matrix_dir / args.matrix_pattern.format(sid=sid)
        if not path.exists():
            missing += 1
            continue

        W = np.loadtxt(path)
        if W.shape != (N_NODES, N_NODES):
            raise ValueError(f"{path.name}: expected {N_NODES} x {N_NODES}, "
                             f"found {W.shape}")

        pc_all = participation(W, net, all_nodes)   # subcortex = 8th module
        pc_cort = participation(W, net, cortical)   # 200 nodes, 7 modules

        nodal.append(pc_all)
        rows.append({
            args.id_col: sid,
            args.group_col: int(getattr(record, args.group_col)),
            "pc_mean_cort": float(np.nanmean(pc_cort)),
            "pc_mean_wb": float(np.nanmean(pc_all)),
        })

    if missing:
        print(f"[warn] {missing} subjects had no matrix file and were skipped")

    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)

    nodal_path = args.out.with_name(args.out.stem + "_nodal.npy")
    np.save(nodal_path, np.vstack(nodal))

    print(f"[saved] {args.out.name}  ({len(out)} subjects)")
    print(f"[saved] {nodal_path.name}  (nodal participation matrix)")


if __name__ == "__main__":
    main()
