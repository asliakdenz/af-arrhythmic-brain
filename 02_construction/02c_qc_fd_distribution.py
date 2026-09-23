#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02c_qc_fd_distribution.py
=========================
Head-motion distribution diagnostics in the primary matched sample
(Supplementary Fig. 4): (A) distribution of mean rs-fMRI head motion with
the 0.30 mm high-motion threshold, and (B) motion by group with a
Mann-Whitney test.

Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
          in atrial fibrillation" (Akdeniz et al.)
Methods: MRI acquisition and connectome construction; Supplementary
         material, "Motion-related Quality Control" (i)

Motion measure: UK Biobank field 25741 ("Mean rs-fMRI head motion, averaged
across space and time points"), column MeanFD_i2 of cohort_matched.csv.
Participants with motion > 0.30 mm are retained in the main analysis and
excluded in the sensitivity cohort (03a_build_design.R, FD_EXCLUDE_MM).

Input:   cohort_matched.csv (01_cohort_assembly.R)
Outputs: PLOT_DIR/fd_hist.{png,svg}, PLOT_DIR/fd_by_group.{png,svg},
         PLOT_DIR/fd_summary.csv
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# =====================================================
# CONFIG  --- edit paths for your environment
# =====================================================
COHORT_CSV   = "/path/to/derivatives/cohort_matched.csv"
PLOT_DIR     = "/path/to/nbs/04_qc/qc_216_main/plots"
MOTION_COL   = "MeanFD_i2"      # UKB field 25741, instance 2
FD_THRESHOLD = 0.30             # mm; high-motion threshold of the sensitivity cohort

os.makedirs(PLOT_DIR, exist_ok=True)


def save_both(basename):
    plt.savefig(os.path.join(PLOT_DIR, basename + ".png"), dpi=200, bbox_inches="tight")
    plt.savefig(os.path.join(PLOT_DIR, basename + ".svg"), format="svg", bbox_inches="tight")


df = pd.read_csv(COHORT_CSV, dtype={"eid": str})
if MOTION_COL not in df.columns:
    raise ValueError(f"Motion column '{MOTION_COL}' not in {COHORT_CSV}")
df[MOTION_COL] = pd.to_numeric(df[MOTION_COL], errors="coerce")
df = df.dropna(subset=[MOTION_COL]).reset_index(drop=True)
fd = df[MOTION_COL].values
print(f"Main cohort N = {len(df)}")

g_af   = df.loc[df["af"] == 1, MOTION_COL]
g_ctrl = df.loc[df["af"] == 0, MOTION_COL]
u, p_u = mannwhitneyu(g_af, g_ctrl, alternative="two-sided")
n_above      = int((fd > FD_THRESHOLD).sum())
n_af_above   = int((g_af > FD_THRESHOLD).sum())
n_ctrl_above = int((g_ctrl > FD_THRESHOLD).sum())

print(f"\nMean rs-fMRI head motion (mm): AF median = {g_af.median():.3f}, "
      f"Control median = {g_ctrl.median():.3f}, Mann-Whitney P = {p_u:.3g}")
print(f"Above {FD_THRESHOLD:.2f} mm: {n_above} ({100 * n_above / len(fd):.1f}%); "
      f"AF = {n_af_above}, Control = {n_ctrl_above}")

pd.DataFrame([{
    "n": len(df), "threshold_mm": FD_THRESHOLD, "n_above": n_above,
    "n_af_above": n_af_above, "n_control_above": n_ctrl_above,
    "af_median": g_af.median(), "control_median": g_ctrl.median(),
    "mannwhitney_U": u, "mannwhitney_P": p_u,
}]).to_csv(os.path.join(PLOT_DIR, "fd_summary.csv"), index=False)

# (A) histogram with the threshold
plt.figure(figsize=(7, 4.5))
plt.hist(fd, bins=50, color="#4C72B0", edgecolor="white")
plt.axvline(FD_THRESHOLD, ls="--", color="red",
            label=f"{FD_THRESHOLD:.2f} mm threshold ({n_above} participants excluded in the sensitivity cohort)")
plt.xlabel("Mean rs-fMRI head motion (UK Biobank field 25741, mm)")
plt.ylabel("Number of participants")
plt.title(f"Mean rs-fMRI head motion, main cohort (N = {len(df)})")
plt.legend(fontsize=9)
plt.tight_layout()
save_both("fd_hist")
plt.close()

# (B) by group
plt.figure(figsize=(5, 4.5))
bp = plt.boxplot([g_ctrl, g_af], patch_artist=True, widths=0.55)
plt.xticks([1, 2], ["Control", "AF"])
for patch, color in zip(bp["boxes"], ["#B0C4DE", "#E68A8A"]):
    patch.set_facecolor(color)
plt.axhline(FD_THRESHOLD, ls="--", color="red", alpha=0.6, label=f"{FD_THRESHOLD:.2f} mm threshold")
plt.ylabel("Mean rs-fMRI head motion (mm)")
plt.title(f"Head motion by group (N = {len(df)})\nMann-Whitney P = {p_u:.3g}")
plt.legend(loc="upper left", fontsize=9)
plt.tight_layout()
save_both("fd_by_group")
plt.close()
print(f"\nSaved fd_hist and fd_by_group (PNG + SVG) and fd_summary.csv in {PLOT_DIR}")
