
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
# Cell (a, b) is a recruitment probability in [0, 1]: among the glomeruli
# assigned to odorant `a` at low concentration, the fraction that odorant `b`
# recruits (silent at low for b, responding at high for b). Rows and columns are
# both blocked by chemical class, so the diagonal blocks are within-class.
#
# Two spec-dependent choices are arguments, not assumptions:
#
# `row_rule` -- how a glomerulus is assigned to a low-concentration odorant.
#     "halfmax" (default, `frac=0.5`): a glomerulus counts for *every* odorant
#     driving it to at least `frac` of its own maximum raw dF/F. Rows are built
#     from ~22 glomeruli and 3 odorants end up with an empty row.
#     "argmax": each glomerulus counts for its single best odorant only. Rows
#     are built from ~4 glomeruli and 11 odorants end up empty -- which is why
#     halfmax is the default.
#
# `prediction_threshold` -- how the model's continuous held-out prediction is
#     turned into a binary "recruited" call. Cross-validated predictions are
#     shrunk toward the mean, so a cutoff chosen on the observed scale
#     under-produces recruitment by roughly 2x.
#     "rate_matched" (default): the cutoff that makes the predicted matrix carry
#     the same total number of recruited pairs as the observed one, so the two
#     panels differ in *pattern* and not in overall rate -- which is the
#     comparison the figure is making.
#     "halfmax": per glomerulus, `frac` of its own maximum predicted response.
#     a float: an absolute cutoff on the model's log1p scale.

#: Default assignment rule and its fraction. See `row_membership`.
ROW_RULE = "halfmax"
ROW_RULE_FRAC = 0.5

#: Default rule for binarising the held-out prediction. See `predicted_mask`.
PREDICTION_THRESHOLD = "rate_matched"

#: The supplement's sweep of the assignment fraction. 0.50 is the default, and
#: the sweep exists to show the matrix does not depend on that choice.
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


def row_membership(m: RecruitmentModel, rule: str = ROW_RULE,
                   frac: float = ROW_RULE_FRAC) -> np.ndarray:
    """Assign glomeruli to low-concentration odorants. (odorants, gloms) bool.

    See the module note above for `rule`. A glomerulus that is silent at low
    concentration everywhere belongs to no odorant under either rule.
    """
    low = m.low
    if rule == "halfmax":
        peak = low.max(0, keepdims=True)
        member = (low >= frac * peak) & (low > 0)
    elif rule == "argmax":
        member = np.zeros(low.shape, bool)
        member[low.argmax(0), np.arange(low.shape[1])] = True
        member &= low > 0
    else:
        raise ValueError(f"unknown row_rule {rule!r}; expected "
                         "'halfmax' or 'argmax'")
    return member


def predicted_mask(m: RecruitmentModel, observed: np.ndarray,
                   threshold=PREDICTION_THRESHOLD,
                   frac: float = ROW_RULE_FRAC):
    """Binarise the held-out prediction into recruited pairs.

    Returns ``(mask, cutoff)``. `observed` is the observed recruited mask, used
    only by the ``"rate_matched"`` rule to fix the total. See the module note.
    """
    silent = m.low == 0                       # recruitment lives here or nowhere
    Yhat = m.Yhat
    if threshold == "halfmax":
        peak = Yhat.max(0, keepdims=True)
        return silent & (Yhat >= frac * peak) & (peak > 0), float(frac)
    if threshold == "rate_matched":
        pool = Yhat[silent]
        k = int(observed.sum())
        if k == 0 or k >= pool.size:
            cutoff = -np.inf if k else np.inf
        else:
            cutoff = float(np.partition(pool, -k)[-k])
        return silent & (Yhat >= cutoff), float(cutoff)
    cutoff = float(threshold)
    return silent & (Yhat > cutoff), cutoff


def recruitment_matrix(m: RecruitmentModel, *, level: str = "odorant",
                       row_rule: str = ROW_RULE, frac: float = ROW_RULE_FRAC,
                       prediction_threshold=PREDICTION_THRESHOLD,
                       threshold: float = 0.0, order="matt", within="table",
                       codes=None,
                       aggregate: str = "pooled",
                       exclude_self: bool = True, recruited: str = "exact",
                       lo_frac: float = 0.2,
                       hi_frac: float = 0.2) -> RecruitmentMatrix:
    """Observed and held-out-predicted recruitment matrices, class-blocked.

    `level` is ``"odorant"`` (59 x 59, the main panel) or ``"class"`` (10 x 10,
    the aggregated version). `order` is the chemical-class order along both
    axes -- see `highconc.clusters.class_order`; ``"matt"`` is the default
    because it is the order the manuscript blocks are laid out in.

    `within` is the odorant order *inside* each class block, at ``level=
    "odorant"`` only -- the same mechanism `highconc.similarity.class_sorted`
    uses, so Fig 3H/Fig 6D can match Fig 3B/C/D's order exactly. ``"table"``
    (default) makes no change from the class-block-only order every existing
    export uses. A non-trivial order (``"matt"``, ``"burton"``) needs `codes`:
    an array of odorant CODEs (the leading token of
    `Dataset.odorants["Odorant"]`) aligned with `m.chem`/`m.names` --
    `RecruitmentModel` doesn't carry codes itself, so the caller supplies
    them, e.g. ``codes=np.asarray([s.split()[0] for s in
    d.odorants["Odorant"].values])``.

    `exclude_self` drops the a == b cells, which are structurally zero: a
    glomerulus assigned to odorant `a` responds to `a` at low concentration, so
    `a` cannot recruit it. At odorant level they become the masked diagonal; at
    class level they are dropped from the pooled counts rather than diluting the
    within-class blocks with zeros.

    `recruited` is how a pair counts as recruited. ``"exact"`` (default) is
    silent at low (``low == 0``) and responding at high (``high > threshold``),
    with the prediction binarised by `prediction_threshold`. ``"relative"``
    judges both concentrations against each glomerulus's own peak at that
    concentration: ``low / low_peak < lo_frac`` and ``high / high_peak >
    hi_frac``. The held-out prediction (``expm1(m.Yhat)``, back on the
    dF/F scale) goes through the same test, so `prediction_threshold` is
    ignored. Row membership (`row_rule`, `frac`) is independent of both.

    `aggregate` only applies at class level. ``"pooled"`` (default) sums the
    glomerulus counts across the block and divides once, so the cell is a
    probability over every glomerulus involved and odorants with more
    glomeruli behind them weigh more. ``"mean"`` averages the odorant-level
    cells unweighted, giving a rare odorant with four glomeruli the same say as
    a common one with forty. The two agree closely -- median cell difference
    0.008, worst case 0.10 -- so the choice is presentational; it ships as an
    argument because there is no way to tell which one the June export used.
    """
    member = row_membership(m, row_rule, frac)
    if recruited == "exact":
        rec_obs = (m.low == 0) & (m.high > threshold)
        rec_pred, cutoff = predicted_mask(m, rec_obs, prediction_threshold, frac)
    elif recruited == "relative":
        rec_obs = _relative_recruited(m.low, m.high, lo_frac, hi_frac)
        rec_pred = _relative_recruited(m.low, np.expm1(m.Yhat), lo_frac, hi_frac)
        cutoff = float(hi_frac)
    else:
        raise ValueError(f"unknown recruited {recruited!r}; expected "
                         "'exact' or 'relative'")

    mem = member.astype(float)
    # Numerator: glomeruli of row a that column b recruits. Denominator: the
    # glomeruli of row a, broadcast -- so the ratio is a probability per row.
    num_obs = mem @ rec_obs.T.astype(float)
    num_pred = mem @ rec_pred.T.astype(float)
    den = np.repeat(mem.sum(1)[:, None], mem.shape[0], axis=1)

    if exclude_self:
        eye = np.eye(den.shape[0], dtype=bool)
        num_obs, num_pred, den = num_obs.copy(), num_pred.copy(), den.copy()
        num_obs[eye] = num_pred[eye] = den[eye] = 0.0

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

    if level == "odorant":
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
        num_obs, num_pred = num_obs[np.ix_(idx, idx)], num_pred[np.ix_(idx, idx)]
        den = den[np.ix_(idx, idx)]
        labels = list(m.names[idx])
        classes = list(m.chem[idx])
    elif level == "class":
        blocks = [m.chem == c for c in order_letters]
        if aggregate == "pooled":
            num_obs = _block_sum(num_obs, blocks)
            num_pred = _block_sum(num_pred, blocks)
            den = _block_sum(den, blocks)
        elif aggregate == "mean":
            num_obs = _block_mean(num_obs, den, blocks)
            num_pred = _block_mean(num_pred, den, blocks)
            den = np.ones_like(num_obs)
        else:
            raise ValueError(f"unknown aggregate {aggregate!r}; expected "
                             "'pooled' or 'mean'")
        labels = [_clusters.class_name(c) for c in order_letters]
        classes = list(order_letters)
    else:
        raise ValueError(f"unknown level {level!r}; expected "
                         "'odorant' or 'class'")

    with np.errstate(invalid="ignore", divide="ignore"):
        obs = np.where(den > 0, num_obs / den, np.nan)
        pred = np.where(den > 0, num_pred / den, np.nan)

    boundaries = _class_boundaries(classes)
    ok = np.isfinite(obs) & np.isfinite(pred)
    off = ok & ~np.eye(len(labels), dtype=bool)
    n_rows = member.sum(1) if level == "odorant" else np.array(
        [member[b].any(0).sum() for b in [m.chem == c for c in order_letters]])

    return RecruitmentMatrix(
        observed=obs, predicted=pred, labels=labels, classes=classes,
        boundaries=boundaries, class_ticks=_class_ticks(classes),
        level=level, row_rule=row_rule, frac=float(frac), cutoff=cutoff,
        n_rows=np.asarray(n_rows), n_empty=int((np.asarray(n_rows) == 0).sum()),
        r=pearson(obs[ok], pred[ok]) if ok.sum() > 1 else np.nan,
        r_offdiag=pearson(obs[off], pred[off]) if off.sum() > 1 else np.nan,
    )


def _relative_recruited(low: np.ndarray, high: np.ndarray, lo_frac: float,
                        hi_frac: float) -> np.ndarray:
    """Weak at low and clear at high, each against the glomerulus's own peak.

    A glomerulus with a zero peak gives NaN fractions and is never recruited.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        lo = low / low.max(0, keepdims=True)
        hi = high / high.max(0, keepdims=True)
    return (lo < lo_frac) & (hi > hi_frac)


def matrix_sweep(m: RecruitmentModel, fracs=SWEEP_FRACS, **kwargs):
    """The matrix at a range of assignment fractions -- the supplement.

    Returns a list of `RecruitmentMatrix`, one per entry of `fracs`. Everything
    else is held fixed, so the sweep isolates the one choice the row rule makes.
    """
    return [recruitment_matrix(m, frac=f, **kwargs) for f in fracs]


def _block_sum(M: np.ndarray, blocks) -> np.ndarray:
    """Sum M over each (block x block) pair. Counts pool; ratios do not."""
    rows = np.array([M[b].sum(0) for b in blocks])
    return np.array([rows[:, b].sum(1) for b in blocks]).T


def _block_mean(num: np.ndarray, den: np.ndarray, blocks) -> np.ndarray:
    """Unweighted mean of the odorant-level ratios over each block pair."""
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(den > 0, num / den, np.nan)
    return np.array([[np.nanmean(ratio[np.ix_(r, c)]) if
                      np.isfinite(ratio[np.ix_(r, c)]).any() else np.nan
                      for c in blocks] for r in blocks])


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
