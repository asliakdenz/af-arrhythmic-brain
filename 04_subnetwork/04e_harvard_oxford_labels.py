#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04e_harvard_oxford_labels.py
============================
Anatomical labels for the 216 atlas nodes by spatial overlap with the
Harvard-Oxford atlas, and the hub table of the NBS subnetwork
(Supplementary Table 9; node labels of Figures 1 and 2).

Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
          in atrial fibrillation" (Akdeniz et al.)
Methods: Network-based statistics -> Functional NBS and hub analysis
         ("anatomical labels were assigned by spatial overlap with the
         Harvard-Oxford atlas")

For each Schaefer 200 cortical parcel, the Harvard-Oxford cortical
maximum-probability label (thr25, 2 mm) with the largest voxel overlap is
assigned, together with the overlap fraction. Subcortical nodes of the
Melbourne Subcortical Atlas (Tian S1) are labelled by structure name.

Atlases are fetched with nilearn (Schaefer 2018, 200 parcels, 7 networks,
2 mm; Harvard-Oxford cort-maxprob-thr25-2mm) and cached in the nilearn data
directory. No participant data are involved.

Inputs:
    --labels   216-node atlas label file (data/atlas/labels_216.txt)
    --hubs     optional weighted hub CSV from 04b_eigenvector_centrality.m
               (RegionName, WeightedEigenvectorCentrality)
Output:
    --out      CSV with RegionName, Hemisphere, Network,
               Top_Overlap_HO_Label, Overlap_fraction and, if --hubs was
               given, WeightedEigenvectorCentrality and Rank
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

TIAN_NAMES = {"HIP": "Hippocampus", "AMY": "Amygdala", "pTHA": "Posterior thalamus",
              "aTHA": "Anterior thalamus", "NAc": "Nucleus accumbens", "GP": "Globus pallidus",
              "PUT": "Putamen", "CAU": "Caudate"}
NETWORK_NAMES = {"Vis": "Visual", "SomMot": "Somatomotor", "DorsAttn": "Dorsal Attention",
                 "SalVentAttn": "Salience/Ventral Attn", "Limbic": "Limbic",
                 "Cont": "Frontoparietal Control", "Default": "Default Mode"}


def get_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--labels", required=True, type=Path)
    p.add_argument("--hubs", type=Path, default=None)
    p.add_argument("--out", required=True, type=Path)
    return p.parse_args()


def hemisphere(label):
    if "_LH_" in label or label.endswith("-lh"):
        return "L"
    if "_RH_" in label or label.endswith("-rh"):
        return "R"
    return ""


def network(label):
    parts = label.split("_")
    return NETWORK_NAMES.get(parts[2], parts[2]) if label.startswith("7Networks") and len(parts) > 2 \
        else "Subcortical"


def schaefer_overlap_labels():
    from nilearn import datasets, image

    schaefer = datasets.fetch_atlas_schaefer_2018(n_rois=200, yeo_networks=7, resolution_mm=2)
    ho = datasets.fetch_atlas_harvard_oxford("cort-maxprob-thr25-2mm")
    s_img = image.load_img(schaefer.maps)
    h_img = image.resample_to_img(ho.maps, s_img, interpolation="nearest")
    s = np.asarray(s_img.dataobj).astype(int)
    h = np.asarray(h_img.dataobj).astype(int)
    s_labels = [l.decode() if isinstance(l, bytes) else str(l) for l in schaefer.labels]
    ho_labels = list(ho.labels)          # index 0 = Background

    out = {}
    for k, name in enumerate(s_labels, start=1):
        vox = h[s == k]
        vox = vox[vox > 0]
        if vox.size == 0:
            out[name] = ("(no cortical overlap)", 0.0)
            continue
        counts = np.bincount(vox, minlength=len(ho_labels))
        top = int(np.argmax(counts))
        out[name] = (ho_labels[top], float(counts[top] / (s == k).sum()))
    return out


def main():
    args = get_args()
    labels = [l.strip() for l in args.labels.read_text().splitlines() if l.strip()]
    overlap = schaefer_overlap_labels()

    rows = []
    for name in labels:
        if name.startswith("7Networks"):
            if name not in overlap:
                raise SystemExit(f"{name} not found in the nilearn Schaefer 200 / 7-network labels")
            ho_label, frac = overlap[name]
        else:
            ho_label, frac = TIAN_NAMES.get(name.split("-")[0], name), np.nan
        rows.append({"RegionName": name, "Hemisphere": hemisphere(name), "Network": network(name),
                     "Top_Overlap_HO_Label": ho_label, "Overlap_fraction": frac})
    table = pd.DataFrame(rows)

    if args.hubs is not None:
        hubs = pd.read_csv(args.hubs)
        table = hubs.merge(table, on="RegionName", how="left")
        ec_col = [c for c in hubs.columns if "Eigenvector" in c][0]
        table = table.sort_values(ec_col, ascending=False).reset_index(drop=True)
        table.insert(0, "Rank", np.arange(1, len(table) + 1))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, index=False)
    print(f"[saved] {args.out}  ({len(table)} rows)")


if __name__ == "__main__":
    main()
