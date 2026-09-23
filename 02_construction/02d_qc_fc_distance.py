#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02d_qc_fc_distance.py
=====================
Edge-wise QC-FC benchmarking, network-level motion contamination and
distance dependence (Supplementary Figs 5 and 6; Supplementary Table 7).

Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
          in atrial fibrillation" (Akdeniz et al.)
Methods: MRI acquisition and connectome construction ("Residual
         motion-related variance was benchmarked with edge-wise QC-FC,
         network-level and distance-dependence analyses"); Supplementary
         material, "Motion-related Quality Control" (ii)-(iv)

For each of the 23,220 upper-triangular edges of the 216-node connectome,
the Pearson correlation between the edge weight and participants' mean
rs-fMRI head motion (UK Biobank field 25741) is computed (QC-FC). Reported:
median |QC-FC|, 95th percentile of |QC-FC|, the proportion of edges with
P < 0.05, the median signed QC-FC, the median |QC-FC| per network block
(seven Yeo/Schaefer networks plus the subcortical block) and the Spearman
correlation between Euclidean inter-node distance and QC-FC.

Set FD_EXCLUDE_MM to 0.30 to repeat the analysis in the sensitivity cohort
that excludes high-motion participants; outputs are then suffixed "_sens".

Inputs:
    - cohort_matched.csv                  (01_cohort_assembly.R; MeanFD_i2)
    - data/atlas/labels_216.txt, coords_216.txt
    - per-subject FC matrices with GSR (02a_fc_construction.py, primary)
Outputs (PLOT_DIR / TABLE_DIR):
    - heatmap_qcfc{suffix}.{png,svg,csv}      network-level median |QC-FC|
    - qcfc_distance{suffix}.{png,svg}         distance dependence
    - qcfc_summary{suffix}.csv                Supplementary Table 7 row
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr, t as t_dist

# =====================================================
# CONFIG  --- edit paths for your environment
# =====================================================
REPO_ROOT     = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
COHORT_CSV    = "/path/to/derivatives/cohort_matched.csv"
LABELS_FILE   = os.path.join(REPO_ROOT, "data", "atlas", "labels_216.txt")
COORDS_FILE   = os.path.join(REPO_ROOT, "data", "atlas", "coords_216.txt")
MATRIX_DIR    = "/path/to/nbs/00_input/matrices/with_gsr/primary"
MATRIX_PATTERN = "conn_{eid}.txt"
PLOT_DIR      = "/path/to/nbs/04_qc/qc_216_main/plots"
TABLE_DIR     = "/path/to/nbs/04_qc/qc_216_main/tables"
MOTION_COL    = "MeanFD_i2"
FD_EXCLUDE_MM = None            # None = main cohort; 0.30 = sensitivity cohort

os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)
suffix = "" if FD_EXCLUDE_MM is None else "_sens"
CORTICAL_ORDER = ["Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default"]


def save_both(basename):
    plt.savefig(os.path.join(PLOT_DIR, basename + ".png"), dpi=200, bbox_inches="tight")
    plt.savefig(os.path.join(PLOT_DIR, basename + ".svg"), format="svg", bbox_inches="tight")


def get_network_name(label):
    parts = str(label).split("_")
    return parts[2] if str(label).startswith("7Networks") and len(parts) > 2 else "Subcortical"


# ---- atlas ------------------------------------------------------------------
with open(LABELS_FILE) as f:
    labels = [ln.strip() for ln in f if ln.strip()]
xyz = np.loadtxt(COORDS_FILE)
if len(labels) != xyz.shape[0]:
    raise ValueError(f"Row mismatch: {len(labels)} labels vs {xyz.shape[0]} coord rows")
networks = np.array([get_network_name(l) for l in labels])
P = len(labels)
triu_idx = np.triu_indices(P, k=1)
edge_dist = squareform(pdist(xyz, metric="euclidean"))[triu_idx]
print(f"Atlas: {P} nodes, {len(edge_dist)} edges")

# ---- cohort + motion ---------------------------------------------------------
df = pd.read_csv(COHORT_CSV, dtype={"eid": str})
df[MOTION_COL] = pd.to_numeric(df[MOTION_COL], errors="coerce")
df = df.dropna(subset=[MOTION_COL]).reset_index(drop=True)
if FD_EXCLUDE_MM is not None:
    n0 = len(df)
    df = df[df[MOTION_COL] <= FD_EXCLUDE_MM].reset_index(drop=True)
    print(f"Sensitivity cohort: excluded {n0 - len(df)} participants with motion > {FD_EXCLUDE_MM} mm")
print(f"Cohort N = {len(df)}")

conn_list, fd_list = [], []
for _, row in df.iterrows():
    f = os.path.join(MATRIX_DIR, MATRIX_PATTERN.format(eid=row["eid"]))
    if not os.path.exists(f):
        continue
    conn_list.append(np.loadtxt(f))
    fd_list.append(row[MOTION_COL])
conn = np.stack(conn_list, axis=0)
fd = np.array(fd_list)
N = conn.shape[0]
print(f"Matrices loaded: {N}")

# ---- edge-wise QC-FC ----------------------------------------------------------
edges = conn[:, triu_idx[0], triu_idx[1]]
x = edges - edges.mean(axis=0, keepdims=True)
y = fd - fd.mean()
r_qcfc = (y @ x) / (np.sqrt((x ** 2).sum(axis=0)) * np.sqrt((y ** 2).sum()) + 1e-12)
t_stat = r_qcfc * np.sqrt((N - 2) / np.clip(1 - r_qcfc ** 2, 1e-12, None))
p_qcfc = 2 * t_dist.sf(np.abs(t_stat), df=N - 2)
rho, p_rho = spearmanr(edge_dist, r_qcfc)

summary = {
    "cohort": "sensitivity" if FD_EXCLUDE_MM is not None else "main",
    "N": N,
    "median_abs_qcfc": float(np.median(np.abs(r_qcfc))),
    "pct95_abs_qcfc": float(np.percentile(np.abs(r_qcfc), 95)),
    "prop_p_lt_0.05": float((p_qcfc < 0.05).mean()),
    "median_signed_qcfc": float(np.median(r_qcfc)),
    "distance_spearman_rho": float(rho), "distance_spearman_P": float(p_rho),
}
pd.DataFrame([summary]).to_csv(os.path.join(TABLE_DIR, f"qcfc_summary{suffix}.csv"), index=False)
for k, v in summary.items():
    print(f"  {k:<24s} {v}")

# ---- network-level heatmap (Supplementary Fig. 5) ------------------------------
qcfc_abs = np.zeros((P, P))
qcfc_abs[triu_idx] = np.abs(r_qcfc)
qcfc_abs = qcfc_abs + qcfc_abs.T
order = [n for n in CORTICAL_ORDER if n in set(networks)] + (["Subcortical"] if "Subcortical" in set(networks) else [])
idx = {n: np.where(networks == n)[0] for n in order}
heat = np.zeros((len(order), len(order)))
for i, a in enumerate(order):
    for j, b in enumerate(order):
        block = qcfc_abs[np.ix_(idx[a], idx[b])]
        vals = block[np.triu_indices_from(block, k=1)] if i == j else block.ravel()
        heat[i, j] = float(np.median(vals)) if vals.size else 0.0
heat_df = pd.DataFrame(heat, index=order, columns=order)
heat_df.to_csv(os.path.join(TABLE_DIR, f"heatmap_qcfc{suffix}.csv"))
plt.figure(figsize=(8.5, 7))
sns.heatmap(heat_df, annot=True, fmt=".3f", cmap="Reds", vmin=0, vmax=0.20, square=True,
            cbar_kws={"label": "median |QC-FC|"})
plt.title(f"Median |QC-FC| by network block (N = {N})")
plt.tight_layout()
save_both(f"heatmap_qcfc{suffix}")
plt.close()

# ---- distance dependence (Supplementary Fig. 6) ------------------------------
plt.figure(figsize=(7.5, 5.5))
hb = plt.hexbin(edge_dist, r_qcfc, gridsize=60, cmap="inferno", mincnt=1, bins="log")
plt.colorbar(hb).set_label("log10(edge count)")
sns.regplot(x=edge_dist, y=r_qcfc, scatter=False, color="cyan", line_kws={"lw": 2})
plt.axhline(0, color="white", ls="--", lw=1)
plt.xlabel("Euclidean distance between nodes (mm)")
plt.ylabel("QC-FC correlation (Pearson r)")
plt.title(f"Distance dependence (N = {N}), Spearman rho = {rho:+.3f}")
plt.ylim(-0.6, 0.6)
plt.tight_layout()
save_both(f"qcfc_distance{suffix}")
plt.close()
print(f"\nOutputs in {PLOT_DIR} and {TABLE_DIR}")
