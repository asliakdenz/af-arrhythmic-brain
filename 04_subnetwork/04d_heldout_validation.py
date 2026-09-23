#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04d_heldout_validation.py
=========================
Within- versus outside-subnetwork connectivity: internal held-out validation
of the AF > Control functional subnetwork (Supplementary Fig. 7) and the
magnitude comparison in the primary cohort (Supplementary Table 4).

Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
          in atrial fibrillation" (Akdeniz et al.)
Methods: Internal held-out validation of the functional subnetwork;
         Results, "A hyperconnected functional subnetwork ..."

--mode heldout (default)
    Generalizability is assessed by frozen-mask validation in the held-out
    sample, which did not contribute to defining the subnetwork: mean Fisher-z
    connectivity across the subnetwork edges is compared between groups by
    OLS regression with HC3 robust standard errors and identical covariate
    adjustment to the primary analysis (age, sex, education, diastolic BP,
    log total WMH, mean framewise displacement). The same contrast is fitted
    on the edges outside the mask, and the mean per-edge group effect
    (t statistic) inside versus outside the mask is reported.

--mode primary
    Supplementary Table 4: in the primary (mask-defining) cohort, mean
    Fisher-z connectivity within and outside the subnetwork per participant,
    with unadjusted group differences (Cohen's d, Mann-Whitney U) and the
    covariate-adjusted OLS/HC3 group effect.

Two guards enforce the design in held-out mode: the held-out participants
must not overlap the mask-defining sample, and the edge list is checked
against the expected size of the primary solution.

Inputs:
    --edges        significant-edge CSV of the primary subnetwork
                   (04a_extract_edges.m; columns Region1, Region2)
    --labels       216-node atlas label file
    --sample       subject table of the sample analysed here
                   (cohort_heldout_matched.csv or cohort_matched.csv)
    --primary-ids  subject table of the mask-defining sample
                   (cohort_matched.csv); used for the overlap guard
    --matrix-dir   per-subject functional connectivity matrices of the sample
Outputs (--out-dir):
    per_subject_fc.csv, summary.csv and, in held-out mode,
    heldout_subnetwork_fc.{svg,pdf}
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu

N_NODES = 216
COVARS = ["age_std", "sex", "education_std", "dia_bp_std", "log_total_wmh_std", "mean_fd_std"]


def get_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--edges", required=True, type=Path)
    p.add_argument("--labels", required=True, type=Path)
    p.add_argument("--sample", required=True, type=Path)
    p.add_argument("--primary-ids", required=True, type=Path)
    p.add_argument("--matrix-dir", required=True, type=Path)
    p.add_argument("--matrix-pattern", default="conn_{sid}.txt")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--mode", choices=["heldout", "primary"], default="heldout")
    p.add_argument("--expected-edges", type=int, default=264,
                   help="Expected size of the primary subnetwork (264 edges at t = 4.5); "
                        "a mismatch stops the run in held-out mode.")
    p.add_argument("--id-col", default="eid")
    p.add_argument("--group-col", default="af")
    return p.parse_args()


def load_labels(path: Path):
    labels = [l.strip() for l in path.read_text().splitlines() if l.strip()]
    if len(labels) != N_NODES:
        raise ValueError(f"expected {N_NODES} labels, found {len(labels)}")
    return labels, {name: i for i, name in enumerate(labels)}


def load_mask(edges_path: Path, name_to_idx: dict):
    edges = pd.read_csv(edges_path)
    if "Region1" not in edges.columns:
        edges = pd.read_csv(edges_path, sep="\t")
    unknown = (set(edges["Region1"]) | set(edges["Region2"])) - set(name_to_idx)
    if unknown:
        raise SystemExit(f"region names absent from the atlas: {sorted(unknown)[:5]}")
    return [(min(name_to_idx[a], name_to_idx[b]), max(name_to_idx[a], name_to_idx[b]))
            for a, b in zip(edges["Region1"], edges["Region2"])]


def network_of(label):
    parts = label.split("_")
    return parts[2] if label.startswith("7Networks") and len(parts) > 2 else "Subcortex"


def cohens_d(a, c):
    pooled = np.sqrt(((a.var(ddof=1) * (len(a) - 1)) + (c.var(ddof=1) * (len(c) - 1)))
                     / (len(a) + len(c) - 2))
    return (a.mean() - c.mean()) / pooled


def main():
    args = get_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    group = args.group_col

    labels, name_to_idx = load_labels(args.labels)
    node_net = np.array([network_of(l) for l in labels])
    edge_pairs = load_mask(args.edges, name_to_idx)
    print(f"[mask] {len(edge_pairs)} edges")
    if args.mode == "heldout" and len(edge_pairs) != args.expected_edges:
        raise SystemExit(f"expected the primary {args.expected_edges}-edge solution, found "
                         f"{len(edge_pairs)}; check that --edges points at the correct run "
                         f"or adjust --expected-edges.")

    sample = pd.read_csv(args.sample, dtype={args.id_col: str})
    primary_ids = set(pd.read_csv(args.primary_ids, dtype={args.id_col: str})[args.id_col])
    sample = sample.dropna(subset=[group] + COVARS).reset_index(drop=True)
    overlap = set(sample[args.id_col]) & primary_ids
    if args.mode == "heldout" and overlap:
        raise SystemExit(f"{len(overlap)} held-out participants also appear in the mask-defining "
                         f"sample; the validation would be circular.")
    if args.mode == "primary" and set(sample[args.id_col]) != primary_ids:
        print("[warn] --sample differs from --primary-ids; primary mode expects the mask-defining cohort")
    print(f"[sample] n = {len(sample)}  cases = {int((sample[group] == 1).sum())}  "
          f"controls = {int((sample[group] == 0).sum())}")

    # --- per-subject connectivity -----------------------------------------
    iu = np.triu_indices(N_NODES, 1)
    mask_set = set(edge_pairs)
    inside = np.array([k for k in range(len(iu[0])) if (iu[0][k], iu[1][k]) in mask_set])
    outside = np.array([k for k in range(len(iu[0])) if (iu[0][k], iu[1][k]) not in mask_set])

    edge_mat = np.full((len(sample), len(iu[0])), np.nan, dtype=np.float32)
    fc_in, fc_out, missing = [], [], 0
    for k, sid in enumerate(sample[args.id_col].values):
        path = args.matrix_dir / args.matrix_pattern.format(sid=sid)
        if not path.exists():
            fc_in.append(np.nan); fc_out.append(np.nan); missing += 1
            continue
        M = np.loadtxt(path)
        if M.shape != (N_NODES, N_NODES):
            raise ValueError(f"{path.name}: expected {N_NODES} x {N_NODES}, found {M.shape}")
        v = M[iu]
        edge_mat[k] = v
        fc_in.append(float(v[inside].mean()))
        fc_out.append(float(v[outside].mean()))
    sample["fc_subnetwork"], sample["fc_outside"] = fc_in, fc_out
    if missing:
        print(f"[warn] {missing} matrices missing")
    keep = sample["fc_subnetwork"].notna().values
    edge_mat, sample = edge_mat[keep], sample.loc[keep].reset_index(drop=True)

    covs = " + ".join(f"C({c})" if c == "sex" else c for c in COVARS)
    rows = []

    def adjusted(col):
        m = smf.ols(f"{col} ~ {group} + {covs}", data=sample).fit(cov_type="HC3")
        return dict(beta=m.params[group], se=m.bse[group], t=m.tvalues[group],
                    p=m.pvalues[group], df_resid=int(m.df_resid))

    a_in  = sample.loc[sample[group] == 1, "fc_subnetwork"]; c_in  = sample.loc[sample[group] == 0, "fc_subnetwork"]
    a_out = sample.loc[sample[group] == 1, "fc_outside"];    c_out = sample.loc[sample[group] == 0, "fc_outside"]

    if args.mode == "primary":
        # ---- Supplementary Table 4 -------------------------------------------
        print("\n" + "=" * 74 + f"\n  PRIMARY COHORT: WITHIN VERSUS OUTSIDE THE SUBNETWORK (n = {len(sample)})\n" + "=" * 74)
        for label, a, c, col, n_edges in [("Within subnetwork", a_in, c_in, "fc_subnetwork", len(inside)),
                                          ("Outside subnetwork", a_out, c_out, "fc_outside", len(outside))]:
            u, p_mw = mannwhitneyu(a, c, alternative="two-sided")
            adj = adjusted(col)
            row = dict(edge_set=label, n_edges=n_edges, AF_median=a.median(), Control_median=c.median(),
                       cohens_d=cohens_d(a, c), mannwhitney_P=p_mw, **{f"adj_{k}": v for k, v in adj.items()})
            rows.append(row)
            print(f"  {label:<20s} edges = {n_edges:>6d}  AF median = {a.median():+.4f}  "
                  f"Control median = {c.median():+.4f}  d = {cohens_d(a, c):+.2f}  MW P = {p_mw:.2e}  "
                  f"adj beta = {adj['beta']:+.4f}  t = {adj['t']:+.2f}  P = {adj['p']:.2e}")
    else:
        # ---- held-out validation ----------------------------------------------
        print("\n" + "=" * 74 + f"\n  HELD-OUT VALIDATION  (n = {len(sample)})\n" + "=" * 74)
        m_in = smf.ols(f"fc_subnetwork ~ {group} + {covs}", data=sample).fit(cov_type="HC3")
        for label, col, a, c in [("(a) within-subnetwork FC", "fc_subnetwork", a_in, c_in),
                                 ("(b) outside-subnetwork FC", "fc_outside", a_out, c_out)]:
            adj = adjusted(col)
            rows.append(dict(model=label, **adj, d=cohens_d(a, c)))
            print(f"  {label:<30s} b = {adj['beta']:+.4f}  SE = {adj['se']:.4f}  t({adj['df_resid']}) = "
                  f"{adj['t']:+.2f}  P = {adj['p']:.2e}  d = {cohens_d(a, c):+.2f}")

        # (c) per-edge group effect inside versus outside the mask
        design = np.column_stack([np.ones(len(sample)), sample[group].values.astype(float)]
                                 + [sample[c].values.astype(float) for c in COVARS])
        XtXi = np.linalg.pinv(design.T @ design)
        B = XtXi @ design.T @ edge_mat
        resid = edge_mat - design @ B
        dof = design.shape[0] - np.linalg.matrix_rank(design)
        sigma2 = np.einsum("ne,ne->e", resid, resid) / dof
        cvec = np.zeros(design.shape[1]); cvec[1] = 1.0
        se_edge = np.sqrt(sigma2 * (cvec @ XtXi @ cvec)); se_edge[se_edge <= 0] = np.inf
        t_edge = (cvec @ B) / se_edge
        t_in, t_out = float(t_edge[inside].mean()), float(t_edge[outside].mean())
        print(f"\n  (c) mean per-edge t: inside = {t_in:+.2f}, outside = {t_out:+.2f}")
        rows.append(dict(model="(c) mean per-edge t inside mask", beta=t_in))
        rows.append(dict(model="(c) mean per-edge t outside mask", beta=t_out))

        is_between = node_net[iu[0]] != node_net[iu[1]]
        print(f"  mask composition: between-network {100 * is_between[inside].mean():.1f}% "
              f"vs {100 * is_between.mean():.1f}% connectome-wide")

        # ---- figure (Supplementary Fig. 7) ------------------------------------
        import matplotlib.pyplot as plt
        import seaborn as sns
        plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "svg.fonttype": "none",
                             "pdf.fonttype": 42})
        plot_df = sample.copy()
        plot_df["label"] = pd.Categorical(plot_df[group].map({0: "Control", 1: "AF"}),
                                          ["Control", "AF"], ordered=True)
        palette = {"Control": "#6A7A8F", "AF": "#E07060"}
        fig, ax = plt.subplots(figsize=(3.0, 3.6), dpi=300)
        sns.violinplot(data=plot_df, x="label", y="fc_subnetwork", hue="label", palette=palette,
                       inner=None, cut=0, linewidth=0.4, width=0.85, ax=ax, legend=False)
        sns.boxplot(data=plot_df, x="label", y="fc_subnetwork", width=0.10, showfliers=False,
                    linewidth=0.4, boxprops=dict(facecolor="white", edgecolor="#333333"), ax=ax)
        ax.axhline(0, color="#8c8c8c", linestyle="--", linewidth=0.35)
        p_val = m_in.pvalues[group]
        stars = "***" if p_val <= 1e-3 else "**" if p_val <= 0.01 else "*" if p_val <= 0.05 else "n.s."
        y = plot_df["fc_subnetwork"]; span = y.max() - y.min()
        bar, top = y.max() + span * 0.06, y.max() + span * 0.08
        ax.plot([0, 0, 1, 1], [bar, top, top, bar], color="#333333", linewidth=0.35)
        ax.text(0.5, top + span * 0.005, stars, ha="center", va="bottom", fontsize=11)
        ax.text(0.02, 0.98, f"d = {cohens_d(a_in, c_in):+.2f}\nt = {m_in.tvalues[group]:+.2f}",
                transform=ax.transAxes, fontsize=8, va="top", ha="left")
        ax.set_xlabel(""); ax.set_ylabel("Mean Fisher z connectivity\nacross subnetwork edges")
        ax.set_title(f"Held-out validation (n = {len(sample)})", fontweight="normal")
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        plt.tight_layout()
        for ext in ("svg", "pdf"):
            fig.savefig(args.out_dir / f"heldout_subnetwork_fc.{ext}", format=ext, bbox_inches="tight")

    sample[[args.id_col, group, "fc_subnetwork", "fc_outside"]].to_csv(
        args.out_dir / "per_subject_fc.csv", index=False)
    summary = pd.DataFrame(rows)
    summary["n_total"] = len(sample)
    summary["n_case"] = int((sample[group] == 1).sum())
    summary["n_control"] = int((sample[group] == 0).sum())
    summary["n_edges_mask"] = len(edge_pairs)
    summary.to_csv(args.out_dir / "summary.csv", index=False)
    print(f"\n[saved] results written to {args.out_dir}/")


if __name__ == "__main__":
    main()
