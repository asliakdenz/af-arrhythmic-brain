#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07d_nested_cross_validation.py
==============================
Fully nested cross-validation of the functional deviation score against
processing speed.

Paper:   "The arrhythmic brain: Interoceptive overload and cognitive slowing
          in atrial fibrillation" (Akdeniz et al.)
Methods: Brain-behaviour analyses -> Hierarchical regression and nested
         cross-validation

Because the subnetwork mask, the normative model and the outcome model were
all fitted in the same sample, in-sample increments in explained variance are
optimistic. The entire pipeline is therefore refitted within each training
fold of a fully nested 10-fold cross-validation (100 repeats, stratified by
group), with no test data used at any stage.

Within every training fold, and never using test data:

    1. Edge-wise GLM -> t-map -> primary threshold -> largest connected
       component. This reproduces the extent-based selection of the
       network-based statistics; permutations are not required here because
       the mask itself is needed, not its family-wise error rate.
    2. Per-edge normative model fitted on training controls only, yielding
       coefficients B_k and residual standard deviation s_k.
    3. Deviation score computed with the frozen fold objects (mask, B_k, s_k).
    4. Baseline and augmented outcome models fitted on training data only.

The frozen fold objects are then applied to the test fold and predictions
pooled out of fold.

The design, contrast and score construction replicate the primary pipeline:

    design   : [case, control, log total WMH, age, sex, education,
                diastolic BP, mean framewise displacement]
    contrast : [1, -1, 0, 0, 0, 0, 0, 0]
    score    : mean over mask edges of
               (|z_edge| - control prediction) / control residual SD

Absolute edge weights are used so the score reflects connectivity magnitude
irrespective of direction, matching the functional deviation score.

Reported quantities are the median out-of-fold change in explained variance
within cases and for the interaction model in the full sample, with
significance assessed by within-group permutation of the outcome.

Usage:
    python 07d_nested_cross_validation.py --subjects derivatives/cohort_matched.csv \
        --cognition derivatives/cohort_matched_cognition.csv \
        --matrix-dir <with_gsr/primary> --out output/nested_cv.csv
    python 07d_nested_cross_validation.py ... --quick     # smoke test
"""

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

# ==============================================================================
# FIXED ANALYSIS PARAMETERS
# ==============================================================================
N_NODES = 216
T_THRESH = 4.5          # primary network-based statistics threshold
N_FOLD = 10
N_REPEAT = 100
N_PERM = 1000
SEED = 20260717

# Design matrix for edge selection; the contrast is [1, -1, 0, ...].
NBS_COVS = ["log_total_wmh_std", "age_std", "sex",
            "education_std", "dia_bp_std", "mean_fd_std"]

# Normative model covariates (core covariates plus mean framewise displacement).
NORM_COVS = ["age_std", "sex", "education_std",
             "log_total_wmh_std", "dia_bp_std", "mean_fd_std"]

# Baseline outcome model (core covariates plus mean framewise displacement).
BASE_COVS = ["age_std", "sex", "education_std",
             "log_total_wmh_std", "dia_bp_std", "mean_fd_std"]

IU = np.triu_indices(N_NODES, 1)          # canonical edge order


# ==============================================================================
# CONFIGURATION
# ==============================================================================
def get_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--subjects", required=True, type=Path,
                   help="Subject table (cohort_matched.csv) with the identifier, "
                        "group and standardised covariate columns.")
    p.add_argument("--cognition", type=Path, default=None,
                   help="Cognition table (cohort_matched_cognition.csv from 07a) "
                        "holding the outcome column; merged on the identifier. "
                        "Omit if the outcome is already in --subjects.")
    p.add_argument("--matrix-dir", required=True, type=Path,
                   help="Directory of per-subject functional connectivity "
                        "matrices, one plain-text 216 x 216 file per subject.")
    p.add_argument("--matrix-pattern", default="conn_{sid}.txt",
                   help="Filename pattern for a subject matrix; '{sid}' is "
                        "replaced by the identifier. Default: 'conn_{sid}.txt'")
    p.add_argument("--out", required=True, type=Path,
                   help="Destination CSV for per-repeat results. The "
                        "permutation null is written alongside it.")
    p.add_argument("--cache", type=Path, default=None,
                   help="Optional path for the edge matrix cache, built once "
                        "and reused across runs.")
    p.add_argument("--id-col", default="eid", help="Identifier column name.")
    p.add_argument("--group-col", default="af",
                   help="Binary group column name (1 = case, 0 = control).")
    p.add_argument("--outcome", default="domain_processing_std",
                   help="Outcome column name.")
    p.add_argument("--quick", action="store_true",
                   help="5 repeats and 50 permutations, for a smoke test.")
    return p.parse_args()


# ==============================================================================
# DATA
# ==============================================================================
def load_edges(subjects, args):
    """Return (n_subjects, n_edges) signed Fisher z upper triangle."""
    if args.cache and args.cache.exists():
        Y = np.load(args.cache)
        if Y.shape[0] == len(subjects):
            print(f"  cache hit: {Y.shape}")
            return Y
        print("  cache size mismatch, rebuilding")

    Y = np.full((len(subjects), len(IU[0])), np.nan, dtype=np.float32)
    missing, t0 = [], time.time()

    for i, sid in enumerate(subjects[args.id_col].values):
        path = args.matrix_dir / args.matrix_pattern.format(sid=sid)
        if not path.exists():
            missing.append(sid)
            continue
        Y[i] = np.loadtxt(path)[IU]           # signed, as used for selection
        if (i + 1) % 100 == 0:
            print(f"    {i + 1}/{len(subjects)}  ({time.time() - t0:.0f}s)")

    if missing:
        print(f"  [warn] {len(missing)} matrices missing")
    if args.cache:
        args.cache.parent.mkdir(parents=True, exist_ok=True)
        np.save(args.cache, Y)
    print(f"  built {Y.shape} in {time.time() - t0:.0f}s")
    return Y


# ==============================================================================
# FOLD-INTERNAL STEPS
# ==============================================================================
def edgewise_t(Y, X, c):
    """Vectorised edge-wise GLM t-map for contrast c. Y:(n,e) X:(n,p) c:(p,)."""
    XtXi = np.linalg.pinv(X.T @ X)
    beta = XtXi @ X.T @ Y
    resid = Y - X @ beta
    dof = X.shape[0] - np.linalg.matrix_rank(X)
    sigma2 = np.einsum("ne,ne->e", resid, resid) / dof
    se = np.sqrt(sigma2 * (c @ XtXi @ c))
    se[se <= 0] = np.inf
    return (c @ beta) / se


def largest_component(tmap, thresh):
    """Extent-based selection: largest connected component of suprathreshold edges."""
    supra = np.where(tmap >= thresh)[0]
    if supra.size == 0:
        return np.array([], dtype=int)

    rows, cols = IU[0][supra], IU[1][supra]
    A = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(N_NODES, N_NODES))
    A = A + A.T
    n_comp, labels = connected_components(A, directed=False)
    if n_comp == N_NODES:
        return np.array([], dtype=int)

    edge_labels = labels[rows]              # endpoints share a component label
    best = max(np.unique(edge_labels), key=lambda L: (edge_labels == L).sum())
    return supra[edge_labels == best]


def fit_normative(Y_abs_train_control, X_train_control):
    """Per-edge OLS on training controls only. Returns (B, s)."""
    B = (np.linalg.pinv(X_train_control.T @ X_train_control)
         @ X_train_control.T @ Y_abs_train_control)
    resid = Y_abs_train_control - X_train_control @ B
    s = resid.std(axis=0, ddof=0)
    s[s <= 0] = np.inf
    return B, s


def apply_score(Y_abs, X, B, s):
    """Frozen normative model -> mean z across mask edges."""
    return ((Y_abs - X @ B) / s).mean(axis=1)


def ols_fit(X, y):
    return np.linalg.pinv(X.T @ X) @ X.T @ y


def r2_oof(y, yhat):
    if y.size == 0:
        return np.nan
    return 1.0 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)


def dice(a, b):
    sa, sb = set(a.tolist()), set(b.tolist())
    return np.nan if not sa and not sb else 2 * len(sa & sb) / (len(sa) + len(sb))


# ==============================================================================
# ONE FULL NESTED PASS
# ==============================================================================
def nested_pass(Y, Y_abs, df, y, seed, group_col):
    """Returns out-of-fold change in explained variance for each contrast."""
    rng = np.random.default_rng(seed)
    n = len(df)
    case = df[group_col].values.astype(float)

    fold = np.empty(n, dtype=int)
    for g in (0, 1):                        # stratified by group
        idx = np.where(case == g)[0]
        rng.shuffle(idx)
        fold[idx] = np.arange(len(idx)) % N_FOLD

    X_nbs = np.column_stack([case, 1 - case] + [df[c].values for c in NBS_COVS])
    X_norm = np.column_stack([np.ones(n)] + [df[c].values for c in NORM_COVS])
    X_base = np.column_stack([np.ones(n)] + [df[c].values for c in BASE_COVS])
    contrast = np.array([1.0, -1.0] + [0.0] * len(NBS_COVS))

    oof = {k: np.full(n, np.nan)
           for k in ["base_case", "aug_case", "base_int", "aug_int"]}
    masks = []

    for k in range(N_FOLD):
        test = fold == k
        train = ~test

        # 1. edge selection on training data only
        mask = largest_component(edgewise_t(Y[train], X_nbs[train], contrast),
                                 T_THRESH)
        if mask.size < 5:                   # degenerate fold
            continue
        masks.append(mask)

        # 2. normative model on training controls only
        train_control = train & (case == 0)
        B, s = fit_normative(Y_abs[train_control][:, mask], X_norm[train_control])

        # 3. score everyone with the frozen fold objects
        score = np.full(n, np.nan)
        score[train] = apply_score(Y_abs[train][:, mask], X_norm[train], B, s)
        score[test] = apply_score(Y_abs[test][:, mask], X_norm[test], B, s)

        # 4a. within-case contrast
        train_case, test_case = train & (case == 1), test & (case == 1)
        if test_case.sum() > 0 and train_case.sum() > 10:
            Xb_tr, Xb_te = X_base[train_case], X_base[test_case]
            oof["base_case"][test_case] = Xb_te @ ols_fit(Xb_tr, y[train_case])
            Xa_tr = np.column_stack([Xb_tr, score[train_case]])
            Xa_te = np.column_stack([Xb_te, score[test_case]])
            oof["aug_case"][test_case] = Xa_te @ ols_fit(Xa_tr, y[train_case])

        # 4b. full-sample interaction contrast
        Xb_tr = np.column_stack([X_base[train], case[train]])
        Xb_te = np.column_stack([X_base[test], case[test]])
        oof["base_int"][test] = Xb_te @ ols_fit(Xb_tr, y[train])
        Xa_tr = np.column_stack([Xb_tr, score[train], score[train] * case[train]])
        Xa_te = np.column_stack([Xb_te, score[test], score[test] * case[test]])
        oof["aug_int"][test] = Xa_te @ ols_fit(Xa_tr, y[train])

    res = {}
    valid = np.isfinite(oof["base_case"]) & np.isfinite(oof["aug_case"])
    if valid.sum() > 50:
        res["r2_base_case"] = r2_oof(y[valid], oof["base_case"][valid])
        res["r2_aug_case"] = r2_oof(y[valid], oof["aug_case"][valid])
        res["dr2_case"] = res["r2_aug_case"] - res["r2_base_case"]

    valid = np.isfinite(oof["base_int"]) & np.isfinite(oof["aug_int"])
    res["r2_base_int"] = r2_oof(y[valid], oof["base_int"][valid])
    res["r2_aug_int"] = r2_oof(y[valid], oof["aug_int"][valid])
    res["dr2_int"] = res["r2_aug_int"] - res["r2_base_int"]
    res["mask_size_mean"] = float(np.mean([len(m) for m in masks])) if masks else 0.0
    return res, masks


# ==============================================================================
# MAIN
# ==============================================================================
def main():
    args = get_args()
    n_repeat = 5 if args.quick else N_REPEAT
    n_perm = 50 if args.quick else N_PERM
    rng = np.random.default_rng(SEED)

    print("=" * 74)
    print("FULLY NESTED CROSS-VALIDATION")
    print("=" * 74)

    df = pd.read_csv(args.subjects, dtype={args.id_col: str})
    if args.cognition is not None:
        cog = pd.read_csv(args.cognition, dtype={args.id_col: str})
        df = df.merge(cog[[args.id_col, args.outcome]], on=args.id_col, how="inner")
    needed = [args.outcome, args.group_col] + \
        sorted(set(NBS_COVS + NORM_COVS + BASE_COVS))
    df = df.dropna(subset=needed).reset_index(drop=True)
    print(f"n = {len(df)}  cases = {int((df[args.group_col] == 1).sum())}  "
          f"controls = {int((df[args.group_col] == 0).sum())}")

    print("\nloading connectivity matrices...")
    Y = load_edges(df, args)
    complete = np.isfinite(Y).all(axis=1)
    if not complete.all():
        print(f"  dropping {int((~complete).sum())} subjects with incomplete matrices")
        df, Y = df[complete].reset_index(drop=True), Y[complete]

    Y_abs = np.abs(Y)
    y = df[args.outcome].values.astype(float)
    case = df[args.group_col].values.astype(float)

    # Sanity check: the full-sample mask should reproduce the primary solution
    # reported in the manuscript (264 edges, 114 nodes at t = 4.5).
    X_nbs = np.column_stack([case, 1 - case] + [df[c].values for c in NBS_COVS])
    contrast = np.array([1.0, -1.0] + [0.0] * len(NBS_COVS))
    full_mask = largest_component(edgewise_t(Y, X_nbs, contrast), T_THRESH)
    n_nodes = len(np.unique(np.concatenate([IU[0][full_mask], IU[1][full_mask]])))
    print(f"\nfull-sample mask at t >= {T_THRESH}: "
          f"{len(full_mask)} edges, {n_nodes} nodes")

    print(f"\nrunning {n_repeat} repeats x {N_FOLD} folds...")
    rows, first_masks, t0 = [], [], time.time()
    for r in range(n_repeat):
        res, masks = nested_pass(Y, Y_abs, df, y, SEED + r, args.group_col)
        res["repeat"] = r
        rows.append(res)
        if r == 0:
            first_masks = masks
        if (r + 1) % 10 == 0:
            print(f"  {r + 1}/{n_repeat}  ({time.time() - t0:.0f}s)")

    results = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out, index=False)

    print("\n" + "=" * 74)
    print("OUT-OF-FOLD RESULTS  (median [2.5th, 97.5th percentile] across repeats)")
    print("=" * 74)
    for tag, label in [("case", "within cases: baseline vs baseline + score"),
                       ("int", "full sample: + score + score x group")]:
        if f"dr2_{tag}" not in results:
            continue
        d = results[f"dr2_{tag}"]
        lo, hi = np.percentile(d, [2.5, 97.5])
        print(f"\n  {label}")
        print(f"    R2 baseline  = {results[f'r2_base_{tag}'].median():+.4f}")
        print(f"    R2 augmented = {results[f'r2_aug_{tag}'].median():+.4f}")
        print(f"    dR2_cv       = {d.median():+.4f}  [{lo:+.4f}, {hi:+.4f}]")
        print(f"    excludes 0   = {bool(lo > 0)}")
    print(f"\n  mean mask size across folds: "
          f"{results['mask_size_mean'].mean():.1f} edges")

    if len(first_masks) > 1:
        pairs = [dice(first_masks[i], first_masks[j])
                 for i in range(len(first_masks))
                 for j in range(i + 1, len(first_masks))]
        print(f"  mask stability (pairwise Dice, first repeat): "
              f"median = {np.median(pairs):.3f}  min = {np.min(pairs):.3f}")

    observed = {t: results[f"dr2_{t}"].median()
                for t in ["case", "int"] if f"dr2_{t}" in results}
    observed = {t: v for t, v in observed.items() if np.isfinite(v)}
    if not observed:
        print("\nno fold produced a suprathreshold component of at least 5 edges; "
              "the out-of-fold increments are undefined and the permutation test is skipped.")
        return
    print(f"\nrunning {n_perm} permutations (outcome shuffled within group)...")
    null = {t: [] for t in observed}
    t0 = time.time()
    for p in range(n_perm):
        y_perm = y.copy()
        for g in (0, 1):
            idx = np.where(case == g)[0]
            y_perm[idx] = rng.permutation(y_perm[idx])
        res, _ = nested_pass(Y, Y_abs, df, y_perm, SEED + 10000 + p,
                             args.group_col)
        for t in observed:
            if f"dr2_{t}" in res:
                null[t].append(res[f"dr2_{t}"])
        if (p + 1) % 100 == 0:
            print(f"  {p + 1}/{n_perm}  ({time.time() - t0:.0f}s)")

    print("\n" + "=" * 74)
    print("PERMUTATION TEST")
    print("=" * 74)
    for t in observed:
        draws = np.array(null[t])
        if draws.size == 0:
            continue
        pval = (np.sum(draws >= observed[t]) + 1) / (len(draws) + 1)
        print(f"  {t:5s}: observed dR2_cv = {observed[t]:+.4f}  "
              f"null median = {np.median(draws):+.4f}  P = {pval:.4f}")

    null_path = args.out.with_name(args.out.stem + "_null.csv")
    pd.DataFrame(null).to_csv(null_path, index=False)
    print(f"\nwrote {args.out.name} and {null_path.name}")


if __name__ == "__main__":
    main()
