#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05a_segregation_functional.py
=============================
Whole-brain functional system segregation.

Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
          in atrial fibrillation" (Akdeniz et al.)
Methods: Network segregation and integration

System segregation index (Chan et al. 2014):

    SyS = (Zw - Zb) / Zw

    Zw = mean connectivity between regions within the same network
    Zb = mean connectivity between regions in different networks

Computed directly from the raw weighted Fisher z-transformed matrices, so no
thresholding and no graph construction is involved.

Networks are the seven canonical systems of the Schaefer 200-region cortical
parcellation. Following convention, negative edges are set to zero before the
ratio is formed. Zw and Zb are additionally computed with edge sign retained,
so the signed connectivity changes underlying any segregation difference can
be characterised.

Two node sets are produced:
    cortical   200 Schaefer nodes, 7 networks              (primary analysis)
    all-node   216 nodes, subcortex treated as an 8th module

Inputs:
    - Per-subject functional connectivity matrices (216 x 216 Fisher z,
      with global signal regression; output of the FC construction step)
    - Subject table with participant identifiers and group membership
    - 216-node atlas label file

Output:
    - One row per subject, with the segregation index and its signed and
      zeroed components for each node set.

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
                   help="Destination CSV.")
    p.add_argument("--id-col", default="eid", help="Identifier column name.")
    p.add_argument("--group-col", default="af",
                   help="Binary group column name (1 = case, 0 = control).")
    return p.parse_args()


# ==============================================================================
# ATLAS -> NETWORK ASSIGNMENT
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
# SEGREGATION
# ==============================================================================
def segregation(W, net_vec, node_mask, zero_negative=True):
    """
    Return (SyS, Zw, Zb) over the nodes selected by node_mask.

    With zero_negative=True the ratio follows the published convention.
    With zero_negative=False the returned Zw and Zb are signed means, used to
    describe the direction of any change rather than to form the ratio.
    """
    idx = np.where(node_mask)[0]
    Wm = W[np.ix_(idx, idx)].astype(float).copy()
    if zero_negative:
        Wm[Wm < 0] = 0.0

    iu = np.triu_indices(len(idx), 1)
    weights = Wm[iu]
    within = net_vec[idx][iu[0]] == net_vec[idx][iu[1]]

    zw = weights[within].mean()
    zb = weights[~within].mean()
    sys_index = (zw - zb) / zw if zw != 0 else np.nan
    return sys_index, zw, zb


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
    rows, missing = [], 0

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

        row = {args.id_col: sid,
               args.group_col: int(getattr(record, args.group_col))}

        # Primary: cortical parcellation, seven networks.
        row["sys_cort"], row["zw_cort"], row["zb_cort"] = \
            segregation(W, net, cortical, zero_negative=True)
        _, row["zw_cort_signed"], row["zb_cort_signed"] = \
            segregation(W, net, cortical, zero_negative=False)

        # Secondary: all regions, subcortex as an eighth module.
        row["sys_all"], row["zw_all"], row["zb_all"] = \
            segregation(W, net, all_nodes, zero_negative=True)
        _, row["zw_all_signed"], row["zb_all_signed"] = \
            segregation(W, net, all_nodes, zero_negative=False)

        rows.append(row)

    if missing:
        print(f"[warn] {missing} subjects had no matrix file and were skipped")

    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[saved] {args.out.name}  ({len(out)} subjects)")


if __name__ == "__main__":
    main()
