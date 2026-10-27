"""Recruitment: the second regime the linear low->high map cannot reach.

The baseline PLS model ``Y_hat = X B`` predicts held-out high-concentration
responses at Q^2 ~ 0.27, and a large block of strong responders is predicted
near zero. Those are *recruited* pairs -- a glomerulus silent at low
concentration but driven at high. A linear scaling map has no input to scale,
so recruitment bounds what the linear rule can do.

This module owns the **computation**: fitting the baseline model, classifying
every (odorant, glomerulus) pair, asking whether recruitment is predictable,
comparing linear against nonlinear regressors, and estimating a population
noise ceiling. Drawing lives in :mod:`highconc.panels.recruitment`.

    from highconc import recruitment as rc
    m = rc.fit()                     # cached; a few minutes the first time
    stats = rc.analyze(m)
    stats.frac_recruited

Definitions, per pair, in dF/F and ReLU-floored at 0:

===============  ==========================================================
responding       ``high > threshold``
scaled           responding and ``low > 0``   (the linear rule has an input)
recruited        responding and ``low == 0``  (silent at low, active at high)
de-recruited     ``high == 0`` and ``low > 0``
===============  ==========================================================

`fit()` is slow (a 4-fold CV PLS plus, in `compare_models`, a random forest and
a kernel ridge over 376 outputs), so it is cached per parameter set. Do not
call it in a loop.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import clusters as _clusters
from . import pls as _pls

# The baseline model -- variant, log1p on both sides, mean-centring inside the
# pipeline, family-stratified 4-fold CV, 8 components -- is defined once in
# highconc.pls and re-exported here so recruitment code reads standalone.
DEFAULT_VARIANT = _pls.DEFAULT_VARIANT
N_COMPONENTS = _pls.N_COMPONENTS
CV_SPLITS = _pls.CV_SPLITS
CV_RANDOM_STATE = _pls.CV_RANDOM_STATE
SCALING = _pls.SCALING


#: The baseline model object, owned by :mod:`highconc.pls` -- every PLS-side
#: figure evaluates the same fit, so there is exactly one of it. Re-exported
#: here under its recruitment-facing name.
RecruitmentModel = _pls.PLSModel

@dataclass
class Recruitment:
    """Everything `analyze` measures about one model, at one threshold."""

    threshold: float
    resp: np.ndarray                 # (odorants, gloms) boolean masks
    scaled: np.ndarray
    recruited: np.ndarray
    derecruited: np.ndarray
    n_resp: int
    n_scaled: int
    n_recruited: int
    n_derecruited: int
    frac_recruited: float
    per_odorant_frac: np.ndarray
    per_glom_prop: np.ndarray        # NaN where a glom never responds at high
    class_frac: dict                 # chemical-class letter -> fraction
    cluster_frac: dict               # Joel's 1-7 -> fraction
    gain_median: float
    gain_supralinear_frac: float
    recruited_glom_low_rate_median: float
    pop_sparse_low: np.ndarray
    pop_sparse_high: np.ndarray


def fit(variant: str = DEFAULT_VARIANT, n_components: int = N_COMPONENTS,
        cv_splits: int = CV_SPLITS, random_state: int = CV_RANDOM_STATE,
        scaling: str = SCALING) -> RecruitmentModel:
    """The baseline model. Thin alias for :func:`highconc.pls.fit`.

    One fit is shared across every PLS-side figure, so a notebook that draws
    recruitment and performance pays for it once. Cached there, not here.
    """
    return _pls.fit(variant, n_components, cv_splits, random_state, scaling)


def analyze(m: RecruitmentModel, threshold: float = 0.0) -> Recruitment:
    """Classify every pair and summarise recruitment.

    `threshold` is the dF/F above which a high-concentration response counts as
    responding. 0 is the manuscript's rule (any non-zero response); the
    fraction recruited is stable as it rises -- see `threshold_sweep`.
    """
    low, high = m.low, m.high
    resp = high > threshold
    scaled = resp & (low > 0)
    recruited = resp & (low == 0)
    derecruited = (high == 0) & (low > 0)

    with np.errstate(invalid="ignore"):
        per_odorant_frac = recruited.sum(1) / resp.sum(1)

    rj, nj = recruited.sum(0), resp.sum(0)
    per_glom_prop = np.full(nj.shape, np.nan)
    per_glom_prop[nj > 0] = rj[nj > 0] / nj[nj > 0]

    class_frac = {}
    for c in sorted(set(m.chem)):
        sel = m.chem == c
        rr, nn = recruited[sel].sum(), resp[sel].sum()
        class_frac[c] = rr / nn if nn else np.nan

    # Keyed by Joel's 1-7, not by our component index: this dict is what the
    # panels read, and everything downstream of highconc.clusters is 1-7.
    cluster_frac = {}
    for n in range(1, _clusters.N_CLUSTERS + 1):
        sel = m.cluster_number == n
        rr, nn = recruited[:, sel].sum(), resp[:, sel].sum()
        cluster_frac[n] = rr / nn if nn else np.nan

    gains = high[scaled] / low[scaled]
    low_resp_rate = (low > 0).mean(0)          # per glom, frac odorants at low
    _, jj = np.where(recruited)

    return Recruitment(
        threshold=float(threshold),
        resp=resp, scaled=scaled, recruited=recruited, derecruited=derecruited,
        n_resp=int(resp.sum()), n_scaled=int(scaled.sum()),
        n_recruited=int(recruited.sum()), n_derecruited=int(derecruited.sum()),
        frac_recruited=float(recruited.sum() / resp.sum()),
        per_odorant_frac=per_odorant_frac, per_glom_prop=per_glom_prop,
        class_frac=class_frac, cluster_frac=cluster_frac,
        gain_median=float(np.median(gains)),
        gain_supralinear_frac=float((gains > 1).mean()),
        recruited_glom_low_rate_median=float(np.median(low_resp_rate[jj])),
        pop_sparse_low=(low > 0).mean(1), pop_sparse_high=(high > 0).mean(1),
    )


def baseline_q2(m: RecruitmentModel) -> float:
    """Global Q^2 of the cross-validated baseline model, over all pairs."""
    return _pls.baseline_q2(m)


def low_response_rate(m: RecruitmentModel) -> np.ndarray:
    """Per glomerulus, the fraction of odorants it responds to at LOW."""
    return (m.low > 0).mean(0)


def predictability(m: RecruitmentModel, r: Recruitment,
                   n_splits: int = 5, random_state: int = 0) -> dict:
    """Can recruitment be predicted from what the linear map ignores?

    Among responding pairs, predict recruited vs scaled from chemistry and
    anatomy with logistic regression, held out by odorant (GroupKFold), so the
    score is generalisation to *new odorants*. Returns CV AUC for the full
    feature set, chemistry only, anatomy only, a shuffled control, and the base
    rate.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import GroupKFold, cross_val_predict
    from sklearn.metrics import roc_auc_score

    chem, gc = m.chem, m.cluster_index
    glom_low_rate = low_response_rate(m)          # anatomy-ish
    odor_low_drive = m.Xlog.mean(1)               # chemistry-ish

    ii, jj = np.where(r.resp)
    label = (m.low[ii, jj] == 0).astype(int)      # recruited
    groups = ii                                   # odorant id

    chem_oh = np.array([[chem[i] == c for c in sorted(set(chem))] for i in ii],
                       float)
    clus_oh = np.eye(_clusters.N_CLUSTERS)[gc[jj]]
    feats = {
        "chemistry": np.column_stack([chem_oh, odor_low_drive[ii]]),
        "anatomy": np.column_stack([clus_oh, m.x[jj], m.y[jj],
                                    glom_low_rate[jj]]),
    }
    feats["full"] = np.column_stack([feats["chemistry"], feats["anatomy"]])

    gkf = GroupKFold(n_splits=n_splits)
    aucs = {}
    for name, Xf in feats.items():
        Xf = StandardScaler().fit_transform(Xf)
        clf = LogisticRegression(max_iter=2000, class_weight="balanced")
        proba = cross_val_predict(clf, Xf, label, cv=gkf, groups=groups,
                                  method="predict_proba")[:, 1]
        aucs[name] = roc_auc_score(label, proba)
    rng = np.random.default_rng(random_state)
    aucs["shuffled"] = roc_auc_score(
        rng.permutation(label),
        StandardScaler().fit_transform(feats["full"])[:, 0])
    aucs["base_rate"] = float(label.mean())
    return aucs


def compare_models(m: RecruitmentModel, n_components: int = N_COMPONENTS,
                   scaling: str = SCALING, random_state: int = 0) -> dict:
    """Linear vs nonlinear low->high prediction, on identical folds.

    Returns ``{name: (Q^2 over all pairs, Pearson r over responding pairs)}``.
    Slow: the random forest fits 300 trees per fold.
    """
    from pls_explorer.cv import cross_val_predict_pls
    from pls_explorer.metrics import q2_global
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.kernel_ridge import KernelRidge

    Xlog, Ylog, folds = m.Xlog, m.Ylog, m.folds
    resp = (m.high > 0).ravel()

    def score(yh):
        return (q2_global(Ylog, yh),
                pearson(Ylog.ravel()[resp], yh.ravel()[resp]))

    out = {}
    out["PLS (linear, 8 LV)"] = score(cross_val_predict_pls(
        Xlog, Ylog, folds, n_components=n_components, scaling=scaling,
        sparse=False, relu_clip=True))
    out["Random forest"] = score(_cv_predict(
        lambda: RandomForestRegressor(n_estimators=300,
                                      random_state=random_state, n_jobs=-1),
        Xlog, Ylog, folds, scale=False))
    out["Kernel ridge (rbf)"] = score(_cv_predict(
        lambda: KernelRidge(kernel="rbf", alpha=1.0, gamma=1.0 / Xlog.shape[1]),
        Xlog, Ylog, folds))
    return out


def noise_ceiling(m: RecruitmentModel) -> dict:
    """Population-level reliability, as a proxy ceiling on the model.

    The processed data hold no trial replicates, so reliability is estimated
    *across subjects*: per (odorant, cluster) take each subject's mean high
    response, then correlate subjects pairwise. The ceiling is the median of
    those correlations. Compare it against the model's own cluster-level r --
    if they are close, the population signal is captured and the per-pair
    shortfall is recruitment, not unmodelled population structure.
    """
    Ylog, Yhat = m.Ylog, m.Yhat
    gc, subj = m.cluster_number, m.subject
    ks = range(1, _clusters.N_CLUSTERS + 1)

    per_subject = []
    for s in sorted(set(subj)):
        cols = []
        for k in ks:
            sel = (subj == s) & (gc == k)
            cols.append(np.nanmean(Ylog[:, sel], axis=1) if sel.any()
                        else np.full(Ylog.shape[0], np.nan))
        per_subject.append(np.array(cols).T)

    rs = []
    for a in range(len(per_subject)):
        for b in range(a + 1, len(per_subject)):
            u, v = per_subject[a].ravel(), per_subject[b].ravel()
            ok = np.isfinite(u) & np.isfinite(v)
            rs.append(pearson(u[ok], v[ok]))

    obs = np.column_stack([Ylog[:, gc == k].mean(1) for k in ks]).ravel()
    pred = np.column_stack([Yhat[:, gc == k].mean(1) for k in ks]).ravel()
    return dict(ceiling=float(np.median(rs)), inter_subject_rs=rs,
                model_cluster_r=pearson(obs, pred))


def threshold_sweep(m: RecruitmentModel,
                    thresholds=(0, 1, 2, 5, 10, 20, 50)):
    """Recruited fraction among responders as the responding threshold rises."""
    thresholds = list(thresholds)
    frac = []
    for t in thresholds:
        resp = m.high > t
        rec = resp & (m.low == 0)
        frac.append(rec.sum() / resp.sum() if resp.sum() else np.nan)
    return thresholds, frac


def pearson(a, b) -> float:
    return float(np.corrcoef(a, b)[0, 1])


def _cv_predict(make_est, X, Y, folds, scale=True, relu=True):
    """Cross-validated prediction for any sklearn regressor, on given folds."""
    from sklearn.preprocessing import StandardScaler

    Yhat = np.zeros_like(Y)
    for tr, te in folds:
        Xtr, Xte = X[tr], X[te]
        if scale:
            sc = StandardScaler().fit(Xtr)
            Xtr, Xte = sc.transform(Xtr), sc.transform(Xte)
        pred = make_est().fit(Xtr, Y[tr]).predict(Xte)
        Yhat[te] = np.clip(pred, 0, None) if relu else pred
    return Yhat


# ---------------------------------------------------------------------------
# The recruitment matrix: which odorants recruit whose glomeruli
# ---------------------------------------------------------------------------
#
# `analyze` above asks how *much* recruitment there is. This section asks how it
# is organised: given the glomeruli a low-concentration odorant drives, which
# other odorants recruit them at high concentration -- and does the linear model
# predict that structure from held-out data?
#
# Ported directly from the lab's original MATLAB (2026-09-24), not the
# probability/fraction reimplementation this module used before. Rows are the
# low-concentration ("source") odorant, columns the recruiting ("target")
# odorant -- both blocked by chemical class, so the diagonal blocks are
# within-class. A cell is NOT a probability here: each row is independently
# rescaled to its own [0, 1] range (see step 5), so a value is a target's
# recruitment strength *relative to that source's own strongest target*, not
# an absolute fraction of glomeruli recruited.
#
# The algorithm, exactly:
#
#   1. Peak-normalise every glomerulus to [0, 1], separately per
#      concentration (each glomerulus's own maximum response, over all 59
#      odorants, in that concentration, is what 1 means).
#   2. `recruited[odor, glom]` -- this (odorant, glomerulus) pair counts as
#      recruited if the glomerulus's normalised LOW response to that odorant
#      is below `lo_threshold` (0.2: not really active) AND its normalised
#      HIGH response is at or above `frac` (0.5, `ROW_RULE_FRAC` --
#      `threshold_effodor` in the original) -- i.e. newly, clearly active.
#   3. `effective_low[odor, glom]` -- this glomerulus counts as tuned to this
#      odorant at low concentration if its normalised LOW response exceeds
#      `frac` (the *same* 0.5, on the *low* side this time -- it is one
#      shared threshold in the original, used for two different questions).
#   4. Tally: for every glomerulus, for every odorant it is effectively tuned
#      to at low concentration, add that glomerulus's whole recruited-by
#      vector into the matching source row, and force that row's own
#      diagonal entry to 1 for every such contribution -- an odorant's
#      relationship to itself is defined, not measured, the same way a
#      correlation matrix's diagonal is 1 by construction. Sum across all
#      376 glomeruli.
#   5. Rescale each row of the summed total to its own [0, 1] range (its own
#      min and max, not a fixed reference) -- MATLAB's
#      ``normalize(recruitmentbyroi, 2, 'range')``.
#
# The predicted side has no reference spec -- `recruitment_matrix_tally`
# documents this module's own choice (apply the identical pipeline to
# `m.Yhat` in place of `m.high`) rather than the old "rate_matched" cutoff
# search, which no longer applies now that observed and predicted go through
# the same rule.

#: Default fraction for `recruited_mask`'s high side and `effective_low_mask`
#: -- `threshold_effodor` in the reference algorithm, one shared value
#: answering two different questions.
ROW_RULE = "halfmax"      # historical label, kept as RecruitmentMatrix metadata
ROW_RULE_FRAC = 0.2

#: The tally algorithm's low-concentration "not really active" cutoff --
#: `lo_threshold` above. Independent of `ROW_RULE_FRAC`: this asks "is the
#: glomerulus basically silent to this odorant at low," not "is it tuned to
#: it" (that's `ROW_RULE_FRAC`, at the opposite, higher, sense).
LO_SILENCE_FRAC = 0.2

#: The supplement's sweep of the assignment fraction. `ROW_RULE_FRAC` (0.50)
#: is the default and the sweep's own middle point; the sweep exists to show
#: the matrix does not depend on that choice. `LO_SILENCE_FRAC` is held fixed
#: across the sweep -- see `matrix_sweep`.
SWEEP_FRACS = (0.20, 0.35, 0.50, 0.65, 0.80)


@dataclass
class RecruitmentMatrix:
    """One observed/predicted recruitment-matrix pair, at one setting."""

    observed: np.ndarray             # (n, n) probabilities, NaN where undefined
    predicted: np.ndarray
    labels: list                     # odorant names, or class display names
    classes: list                    # class letter per row/column
    boundaries: list                 # indices where a new class block starts
    class_ticks: list                # (position, class name) for axis labels
    level: str                       # "odorant" | "class"
    row_rule: str
    frac: float
    cutoff: float                    # the prediction cutoff actually used
    n_rows: np.ndarray               # glomeruli behind each row
    n_empty: int                     # rows with no glomeruli at all
    r: float                         # observed vs predicted, all defined cells
    r_offdiag: float                 # the same, excluding the diagonal


def _peak_normalized(x: np.ndarray) -> np.ndarray:
    """Each column (ROI -- odorant or glomerulus, whichever axis is columns
    here) scaled by its own peak, into [0, 1]. A column whose peak is 0
    (never fires) stays all-zero rather than dividing by zero, so it never
    crosses a positive fraction threshold."""
    peak = x.max(0, keepdims=True)
    peak = np.where(peak > 0, peak, np.inf)
    return x / peak


def predicted_high(m: RecruitmentModel) -> np.ndarray:
    """`m.Yhat` -- a held-out cross-validated prediction of `m.Ylog` (log1p
    dF/F) -- converted back to raw dF/F units and ReLU-floored at 0, so it is
    on the same scale as `m.high`. The conversion `recruitment_matrix_tally`
    applies to its predicted side; exposed here so a caller doesn't have to
    duplicate it to look at the raw numbers behind the predicted matrix.
    """
    return np.clip(np.expm1(m.Yhat), 0, None)


def recruited_mask(low: np.ndarray, high: np.ndarray, *,
                   lo_threshold: float = LO_SILENCE_FRAC,
                   frac: float = ROW_RULE_FRAC) -> np.ndarray:
    """Which (odorant, ROI) pairs count as recruited: not really active at
    low concentration, clearly active at high -- both relative to that ROI's
    own peak in each concentration, never an absolute dF/F cutoff.

    Takes plain arrays (not a `RecruitmentModel`) so it works the same way
    on `m.low`/`m.high` (the observed side) and on `m.Xlog`... no -- on
    whatever the caller's two concentration arrays are (observed or
    predicted), same shape, same axis convention (odorants x ROIs).

    Matches the reference algorithm exactly:
    ``recruited = (tmpdatalo < lo_threshold) & (tmpdatahi >= threshold_effodor)``
    """
    lo_norm = _peak_normalized(low)
    hi_norm = _peak_normalized(high)
    return (lo_norm < lo_threshold) & (hi_norm >= frac)


def effective_low_mask(low: np.ndarray, *, frac: float = ROW_RULE_FRAC
                       ) -> np.ndarray:
    """Which (odorant, ROI) pairs count as the ROI's 'effective' odorants at
    low concentration -- active above `frac` of that ROI's own low peak.

    Same `frac` as `recruited_mask`'s high-side cutoff by design
    (`threshold_effodor` in the reference algorithm is one shared value,
    asked of two different things), not a coincidence.
    """
    return _peak_normalized(low) > frac


def _tally(recruited: np.ndarray, effective_low: np.ndarray) -> np.ndarray:
    """Steps 4-5 of the reference algorithm, shared by observed and
    predicted: tally, force the diagonal, rescale each row to [0, 1].

    `recruited`, `effective_low` : (n_odorants, n_rois) bool, same axis
    convention as `m.low`/`m.high` -- odorants down, ROIs (glomeruli) across.

    For every ROI (column), for every odorant it's effectively tuned to at
    low concentration, add that ROI's whole recruited-by row into the
    matching source row of the running total, forcing that row's own
    diagonal entry to 1 for every such contribution -- an odorant's
    relationship to itself is defined, not measured. Then rescale each row
    of the total to its own [0, 1] range (its own min and max).

    NaN, not 0, for a row with no contributions at all (nothing was ever
    this row's odorant effectively) -- undefined, matching how every other
    empty-row case in this module is treated, not a measured floor.
    """
    n_odorants, n_rois = recruited.shape
    total = np.zeros((n_odorants, n_odorants))
    for roi in range(n_rois):
        eff_odors = np.where(effective_low[:, roi])[0]
        if eff_odors.size == 0:
            continue
        row = recruited[:, roi].astype(float)
        for o in eff_odors:
            contribution = row.copy()
            contribution[o] = 1.0        # forced, matching the reference exactly
            total[o] += contribution

    row_min = total.min(axis=1, keepdims=True)
    row_max = total.max(axis=1, keepdims=True)
    span = row_max - row_min
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(span > 0, (total - row_min) / span, np.nan)


def recruitment_matrix_tally(m: RecruitmentModel, *,
                             lo_threshold: float = LO_SILENCE_FRAC,
                             frac: float = ROW_RULE_FRAC,
                             predicted: str = "same_rule") -> np.ndarray:
    """The odor x odor recruitment matrix, by direct tally -- ported from the
    lab's original MATLAB, replacing this module's earlier probability
    reimplementation. See the module note above for the full algorithm.

    A cell here is NOT a probability: each row is independently rescaled to
    its own [0, 1] range, so a value is a target odorant's recruitment
    strength *relative to that source odorant's own strongest target*, not
    an absolute fraction of glomeruli recruited.

    `predicted` -- the reference algorithm only specifies the observed side;
    this is my extrapolation for the model comparison, not something to
    take on faith:
        "same_rule" (default): apply the identical `recruited_mask` /
        `effective_low_mask` / tally / row-rescale pipeline to `m.Yhat` in
        place of `m.high`, so observed and predicted end up on the same kind
        of scale (each their own row-relative strength) and stay comparable.
        `m.Yhat` is a prediction of `m.Ylog` (log1p dF/F), not of `m.high`
        directly, so it is converted back with `expm1` and ReLU-floored at
        0 first -- otherwise `m.low` (raw dF/F) and `m.Yhat` (log1p dF/F)
        would be peak-normalised in two different domains, and "50% of
        peak" would not mean the same thing on both sides of `recruited`.
        `None`: skip the predicted side; returns only the observed matrix.

    Returns
    -------
    observed : ndarray, (59, 59)
    predicted : ndarray, (59, 59), or None if `predicted=None`
    """
    obs_recruited = recruited_mask(m.low, m.high, lo_threshold=lo_threshold,
                                   frac=frac)
    low_eff = effective_low_mask(m.low, frac=frac)
    observed = _tally(obs_recruited, low_eff)

    if predicted is None:
        return observed, None
    if predicted != "same_rule":
        raise ValueError(f"predicted must be 'same_rule' or None, got {predicted!r}")

    high_hat = predicted_high(m)
    pred_recruited = recruited_mask(m.low, high_hat, lo_threshold=lo_threshold,
                                    frac=frac)
    predicted_arr = _tally(pred_recruited, low_eff)
    return observed, predicted_arr


def recruitment_matrix(m: RecruitmentModel, *, level: str = "odorant",
                       lo_threshold: float = LO_SILENCE_FRAC,
                       frac: float = ROW_RULE_FRAC,
                       predicted: str = "same_rule",
                       order="matt", within="table",
                       codes=None) -> RecruitmentMatrix:
    """Observed and held-out-predicted recruitment matrices, class-blocked.

    The tally algorithm (see the module note above): each row is
    independently rescaled to its own [0, 1] range, so a cell is a target
    odorant's recruitment strength *relative to that source odorant's own
    strongest target* -- not an absolute fraction of glomeruli recruited,
    which is what this function returned before this algorithm replaced it.

    `level` is ``"odorant"`` only -- the reference algorithm has no
    class-level aggregate, and no current figure asks for one.

    `lo_threshold`, `frac` -- see `recruited_mask`/`effective_low_mask` and
    the module note above. `predicted` -- see `recruitment_matrix_tally`;
    it is this function's own extrapolation for the model-comparison side,
    not part of the reference algorithm.

    `order` is the chemical-class order along both axes -- see
    `highconc.clusters.class_order`; ``"matt"`` is the default because it is
    the order the manuscript blocks are laid out in.

    `within` is the odorant order *inside* each class block -- the same
    mechanism `highconc.similarity.class_sorted` uses, so Fig 3H/Fig 6D can
    match Fig 3B/C/D's order exactly. ``"table"`` (default) makes no change
    from the class-block-only order every existing export uses. A
    non-trivial order (``"matt"``, ``"burton"``) needs `codes`: an array of
    odorant CODEs (the leading token of `Dataset.odorants["Odorant"]`)
    aligned with `m.chem`/`m.names` -- `RecruitmentModel` doesn't carry
    codes itself, so the caller supplies them, e.g.
    ``codes=np.asarray([s.split()[0] for s in d.odorants["Odorant"].values])``.
    """
    if level != "odorant":
        raise NotImplementedError(
            f"level={level!r} is not supported by the tally algorithm -- "
            "the reference spec only covers odorant level, and no current "
            "figure asks for a class-level aggregate.")

    observed, predicted_arr = recruitment_matrix_tally(
        m, lo_threshold=lo_threshold, frac=frac, predicted=predicted)

    order_letters = _clusters.class_order(order, present=list(m.chem))
    if isinstance(within, str):
        try:
            within_map = _clusters.WITHIN_BLOCK_ORDERS[within.lower()]
        except KeyError:
            raise ValueError(
                f"unknown within-block order {within!r}; expected one of "
                f"{sorted(_clusters.WITHIN_BLOCK_ORDERS)} or a "
                "{class: codes} mapping") from None
    else:
        within_map = within

    idx_chunks = []
    for c in order_letters:
        rows = np.where(m.chem == c)[0]
        want = within_map.get(c)
        if want is not None:
            if codes is None:
                raise ValueError(
                    f"within-block order for class {c!r} needs `codes` "
                    "(an array of odorant codes aligned with m.chem)")
            pos = {code: i for i, code in enumerate(np.asarray(codes)[rows])}
            if sorted(pos) != sorted(want):
                raise ValueError(
                    f"within-block order for class {c!r} lists "
                    f"{sorted(want)} but the data has {sorted(pos)}")
            rows = rows[[pos[code] for code in want]]
        idx_chunks.append(rows)
    idx = np.concatenate(idx_chunks)

    obs = observed[np.ix_(idx, idx)]
    labels = list(m.names[idx])
    classes = list(m.chem[idx])
    boundaries = _class_boundaries(classes)

    n_rows_native = effective_low_mask(m.low, frac=frac).sum(axis=1)
    n_rows = n_rows_native[idx]

    if predicted_arr is not None:
        pred = predicted_arr[np.ix_(idx, idx)]
        ok = np.isfinite(obs) & np.isfinite(pred)
        off = ok & ~np.eye(len(labels), dtype=bool)
        r = pearson(obs[ok], pred[ok]) if ok.sum() > 1 else np.nan
        r_offdiag = pearson(obs[off], pred[off]) if off.sum() > 1 else np.nan
    else:
        pred = np.full_like(obs, np.nan)
        r = r_offdiag = np.nan

    return RecruitmentMatrix(
        observed=obs, predicted=pred, labels=labels, classes=classes,
        boundaries=boundaries, class_ticks=_class_ticks(classes),
        level=level, row_rule=ROW_RULE, frac=float(frac),
        cutoff=np.nan,  # no single cutoff under the tally algorithm
        n_rows=np.asarray(n_rows), n_empty=int((np.asarray(n_rows) == 0).sum()),
        r=r, r_offdiag=r_offdiag,
    )


def matrix_sweep(m: RecruitmentModel, fracs=SWEEP_FRACS, **kwargs):
    """The matrix at a range of assignment fractions -- the supplement.

    Returns a list of `RecruitmentMatrix`, one per entry of `fracs`. Everything
    else is held fixed, so the sweep isolates the one choice the row rule makes.
    """
    return [recruitment_matrix(m, frac=f, **kwargs) for f in fracs]


def _class_boundaries(classes) -> list:
    return [i for i in range(1, len(classes)) if classes[i] != classes[i - 1]]


def _class_ticks(classes) -> list:
    """(centre position, display name) for one tick per contiguous class run."""
    out, start = [], 0
    for i in range(1, len(classes) + 1):
        if i == len(classes) or classes[i] != classes[start]:
            out.append(((start + i - 1) / 2, _clusters.class_name(classes[start])))
            start = i
    return out



