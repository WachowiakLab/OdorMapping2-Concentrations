"""The baseline low->high PLS model, and the quantities the PLS figures read.

One fit, shared. Every PLS-side figure evaluates the same manuscript baseline --
glom_explorer ``nonresponders_dropped``, ``log1p`` on both sides, mean-centred
inside the pipeline, family-stratified 4-fold CV, 8 components -- so it is
fitted once here, cached, and handed to the panels as a `PLSModel`.

    from highconc import pls
    m = pls.fit()                    # cached; ~2 s
    obs, pred = pls.cluster_by_odorant(m)

Per-pair Q^2 is modest (0.27) because the signal lives in pattern and
population structure rather than per-pixel magnitude, and because recruitment
(see `highconc.recruitment`) puts a hard ceiling on it. The functions here are
the honest re-evaluations at the grain where the model is actually good:
cluster x odorant means, class x cluster "rules", per-odorant pattern
correlation, odorant decoding, representational geometry, calibration, and
top-k retrieval.

Cluster axes come out in **Joel's 1-7 order**, not raw component order.
Chemical-class axes take an `order` argument -- see
`highconc.clusters.class_order`; the default `"table"` is alphabetical by
letter, which is what these figures have always shown.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from . import clusters as _clusters

#: The clustering variant the manuscript uses.
DEFAULT_VARIANT = "nonresponders_dropped"

#: Baseline model hyperparameters.
N_COMPONENTS = 8
CV_SPLITS = 4
CV_RANDOM_STATE = 42
SCALING = "mean_center"

#: The canonical run's result files, vendored into this bundle.
#:
#: Three functions here read files the pls-explorer run wrote to disk rather
#: than refitting them (`load_weights`, `component_selection`, `lv_loadings`).
#: `pls_explorer.config.RESULTS_DIR` is ``<repo>/results`` inside a *source
#: checkout* of pls-explorer -- but the pinned dependency installs as a wheel
#: into site-packages, where no ``results/`` sits beside the package. On any
#: fresh clone that path simply does not exist, which is how the S5
#: fresh-machine test found this.
#:
#: So the five files the figures actually read -- ``component_selection.json``,
#: ``real/loadings.npz``, and ``B.npy`` for real/null1/null2, 3.3 MB in all --
#: are copied into ``data/pls_results/`` here and ship with the clone.
BUNDLED_RESULTS = Path(__file__).resolve().parent.parent / "data" / "pls_results"

#: Environment variable that overrides where results are read from.
RESULTS_ENV_VAR = "HIGHCONC_PLS_RESULTS"


def results_root() -> Path:
    """Where to read the canonical run's result files from.

    First hit wins:

    1. ``$HIGHCONC_PLS_RESULTS`` -- the escape hatch. Point it at your own
       results tree to render the figures against a different run.
    2. :data:`BUNDLED_RESULTS` -- the default, and what a clone gets.
    3. ``pls_explorer.config.RESULTS_DIR`` -- a source checkout of
       pls-explorer, which is what a developer with the repo beside this one
       has. Reached only if the bundled copy has been deleted.

    Returns the *first path that exists*; if none do, returns the bundled path
    so the caller's error message names the expected location.
    """
    override = os.environ.get(RESULTS_ENV_VAR)
    if override:
        return Path(override).expanduser()
    if BUNDLED_RESULTS.is_dir():
        return BUNDLED_RESULTS
    try:
        from pls_explorer.config import RESULTS_DIR
    except ImportError:
        return BUNDLED_RESULTS
    return Path(RESULTS_DIR) if Path(RESULTS_DIR).is_dir() else BUNDLED_RESULTS


@dataclass(frozen=True)
class PLSModel:
    """Paired responses plus the cross-validated baseline prediction.

    Attributes
    ----------
    low, high : ndarray, (59 odorants, 376 gloms)
        Raw dF/F, ReLU-floored at 0.
    Xlog, Ylog : ndarray
        ``log1p`` of the above -- the space the model is fit and scored in.
    Yhat : ndarray
        Cross-validated prediction of ``Ylog``.
    folds : list of (train, test)
        The family-stratified fold assignment. Reused by every comparison
        model, so alternatives are scored on identical splits.
    chem : ndarray, (59,)
        Chemical-class letter per odorant.
    chem_pooled : ndarray, (59,)
        The same with classes below min_size pooled -- what the folds use.
    names : ndarray, (59,)
        Odorant names, in row order.
    cluster_index : ndarray, (376,)
        Our NMF component index 0-6 per glomerulus.
    cluster_number : ndarray, (376,)
        The same glomeruli in Joel's 1-7 numbering. Prefer this for display.
    xy : ndarray, (376, 2)
        Dorsal-bulb coordinates.
    subject : ndarray, (376,)
    odorants : DataFrame
    variant : str
    """

    low: np.ndarray
    high: np.ndarray
    Xlog: np.ndarray
    Ylog: np.ndarray
    Yhat: np.ndarray
    folds: list
    chem: np.ndarray
    chem_pooled: np.ndarray
    names: np.ndarray
    cluster_index: np.ndarray
    cluster_number: np.ndarray
    xy: np.ndarray
    subject: np.ndarray
    odorants: object
    variant: str

    @property
    def x(self) -> np.ndarray:
        return self.xy[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.xy[:, 1]

    @property
    def n_odorants(self) -> int:
        return self.Ylog.shape[0]

    @property
    def n_gloms(self) -> int:
        return self.Ylog.shape[1]

    @property
    def chem_classes(self) -> list:
        """The chemical-class letters present, sorted."""
        return sorted(set(self.chem))

    def __repr__(self) -> str:
        return (f"PLSModel(variant={self.variant!r}, "
                f"{self.n_odorants} odorants x {self.n_gloms} gloms, "
                f"{N_COMPONENTS} LV)")


@lru_cache(maxsize=2)
def fit(variant: str = DEFAULT_VARIANT, n_components: int = N_COMPONENTS,
        cv_splits: int = CV_SPLITS, random_state: int = CV_RANDOM_STATE,
        scaling: str = SCALING) -> PLSModel:
    """Load the paired data and cross-validate the baseline model.

    Cached, so every figure in a notebook session shares one fit. The first
    call takes a couple of seconds.
    """
    from pls_explorer.data import load_paired, pool_small_groups
    from pls_explorer.cv import make_folds, cross_val_predict_pls
    from pls_explorer.preprocess import log1p_positive

    p = load_paired(variant)
    low, high = p.X.to_numpy(float), p.Y.to_numpy(float)
    Xlog, Ylog = log1p_positive(low), log1p_positive(high)
    cgroup = pool_small_groups(p.chemical_group, min_size=3)
    folds = make_folds(n_samples=Xlog.shape[0], scheme="family_stratified",
                       n_splits=cv_splits, chemical_group=cgroup,
                       random_state=random_state)
    Yhat = cross_val_predict_pls(Xlog, Ylog, folds, n_components=n_components,
                                 scaling=scaling, sparse=False, relu_clip=True)

    cols = p.X.columns
    odor_idx = p.X.index.get_level_values("odorant").astype(int).to_numpy()
    cluster_index = np.asarray(p.glom_clusters.astype(int))
    return PLSModel(
        low=low, high=high, Xlog=Xlog, Ylog=Ylog, Yhat=Yhat, folds=folds,
        chem=np.asarray(p.chemical_group), chem_pooled=cgroup,
        names=p.odorants["Name"].iloc[odor_idx].to_numpy(),
        cluster_index=cluster_index,
        cluster_number=np.array([_clusters.to_joel(int(i))
                                 for i in cluster_index], dtype=int),
        xy=np.column_stack([cols.get_level_values("x").to_numpy(float),
                            cols.get_level_values("y").to_numpy(float)]),
        subject=cols.get_level_values("subject").to_numpy(),
        odorants=p.odorants, variant=variant,
    )


def baseline_q2(m: PLSModel) -> float:
    """Global Q^2 of the cross-validated baseline, over all pairs."""
    from pls_explorer.metrics import q2_global

    return float(q2_global(m.Ylog, m.Yhat))


def baseline_rmse(m: PLSModel) -> float:
    """Cross-validated RMSE of the baseline, over all pairs."""
    from pls_explorer.metrics import rmse

    return float(rmse(m.Ylog, m.Yhat))


# ---------------------------------------------------------------------------
# Population grain: clusters and class x cluster rules
# ---------------------------------------------------------------------------
def cluster_by_odorant(m: PLSModel):
    """Mean response per (odorant, cluster), observed and predicted.

    Returns two ``(59, 7)`` arrays with columns in Joel's 1-7 order.
    """
    gc = m.cluster_number
    ks = range(1, _clusters.N_CLUSTERS + 1)
    obs = np.column_stack([m.Ylog[:, gc == k].mean(1) for k in ks])
    pred = np.column_stack([m.Yhat[:, gc == k].mean(1) for k in ks])
    return obs, pred


def class_by_cluster(m: PLSModel, M: np.ndarray, order="table"):
    """Mean of `M` in each (chemical class x cluster) block.

    Returns ``(matrix, classes)``; rows follow `order` (see
    `highconc.clusters.class_order`), columns Joel's 1-7.
    """
    classes = _clusters.class_order(order, present=m.chem_classes)
    gc = m.cluster_number
    out = np.zeros((len(classes), _clusters.N_CLUSTERS))
    for ci, c in enumerate(classes):
        for k in range(1, _clusters.N_CLUSTERS + 1):
            blk = M[np.ix_(m.chem == c, gc == k)]
            out[ci, k - 1] = blk.mean() if blk.size else 0.0
    return out, classes


def rules_matrices(m: PLSModel, order="table", normalize=True):
    """The observed and predicted class x cluster "rules", on a shared scale.

    `normalize=True` divides each cluster column by its own maximum, so the
    question the panel asks is "which classes drive this cluster", not "which
    cluster is loudest".
    """
    obs, classes = class_by_cluster(m, m.Ylog, order=order)
    pred, _ = class_by_cluster(m, m.Yhat, order=order)
    if normalize:
        return _colnorm(obs), _colnorm(pred), classes
    return obs, pred, classes


def _colnorm(M):
    cm = M.max(0, keepdims=True)
    cm[cm == 0] = 1
    return M / cm


# ---------------------------------------------------------------------------
# Pattern grain: per-odorant maps, decoding, geometry
# ---------------------------------------------------------------------------
def zscore_rows(M: np.ndarray) -> np.ndarray:
    """Z-score each row (odorant), so magnitude drops out and pattern remains."""
    M = np.asarray(M, float)
    mu = M.mean(1, keepdims=True)
    sd = M.std(1, keepdims=True)
    sd = np.where(sd == 0, 1.0, sd)
    return (M - mu) / sd


def similarity_matrices(m: PLSModel):
    """``(S, corr_obs, corr_pred)``: predicted-vs-observed map correlations and
    the two within-modality odorant x odorant correlation matrices."""
    Zo, Zp = zscore_rows(m.Ylog), zscore_rows(m.Yhat)
    n = Zo.shape[1]
    return Zp @ Zo.T / n, Zo @ Zo.T / n, Zp @ Zp.T / n


def pattern_correlations(m: PLSModel):
    """Matched (diagonal) vs mismatched (off-diagonal) map correlations."""
    S, _, _ = similarity_matrices(m)
    n = S.shape[0]
    return np.diag(S).copy(), S[~np.eye(n, dtype=bool)]


def decoding(m: PLSModel) -> dict:
    """Can a predicted map be matched back to its own odorant?

    The shared across-odorant mean map is removed from both sides first: it is
    stimulus-independent, and leaving it in would let a constant map "decode".
    Scoring is by RANK of the true odorant, not top-1 alone -- with 59 odorants
    and correlated chemistry, rank is the honest statistic.
    """
    Yo = m.Ylog - m.Ylog.mean(0, keepdims=True)
    Yp = m.Yhat - m.Yhat.mean(0, keepdims=True)
    S = zscore_rows(Yp) @ zscore_rows(Yo).T / Yo.shape[1]
    rank = (S >= np.diag(S)[:, None]).sum(1)       # 1 = true is the best match
    n = S.shape[0]
    return dict(S=S, rank=rank, n=n,
                top1=float((rank <= 1).mean()),
                top5=float((rank <= 5).mean()),
                top10=float((rank <= 10).mean()),
                median_rank=float(np.median(rank)),
                chance=1.0 / n)


def rsa(m: PLSModel) -> dict:
    """Representational similarity: do predicted maps keep the odour-space
    geometry? Spearman over the upper triangles of the two correlation-distance
    matrices."""
    from scipy.stats import spearmanr

    _, corr_obs, corr_pred = similarity_matrices(m)
    iu = np.triu_indices(corr_obs.shape[0], 1)
    obs_d, pred_d = (1 - corr_obs)[iu], (1 - corr_pred)[iu]
    return dict(obs=obs_d, pred=pred_d,
                rho=float(spearmanr(obs_d, pred_d).correlation))


# ---------------------------------------------------------------------------
# Calibration and retrieval
# ---------------------------------------------------------------------------
def scaled_pairs(m: PLSModel):
    """Observed and predicted values for pairs the linear rule can act on --
    responding at high AND non-zero at low."""
    mask = (m.Ylog.ravel() > 0) & (m.low.ravel() > 0)
    return m.Ylog.ravel()[mask], m.Yhat.ravel()[mask]


def responding_pairs(m: PLSModel):
    """Observed and predicted values for every responding pair."""
    mask = m.Ylog.ravel() > 0
    return m.Ylog.ravel()[mask], m.Yhat.ravel()[mask]


def calibration_deciles(m: PLSModel, n_bins: int = 10):
    """Mean observed response per predicted-value bin, with SEM.

    Bins are equal-COUNT (rank-based), not equal-width: the prediction has a
    large tie-block at zero, and equal-width bins would put most of the data in
    one bar.
    """
    from scipy.stats import spearmanr

    obs, pred = m.Ylog.ravel(), m.Yhat.ravel()
    groups = np.array_split(np.argsort(pred, kind="stable"), n_bins)
    mean = np.array([obs[g].mean() for g in groups])
    sem = np.array([obs[g].std() / np.sqrt(len(g)) for g in groups])
    return mean, sem, float(spearmanr(pred, obs).correlation)


def topk_retrieval(m: PLSModel, ks=(5, 10, 20, 40)):
    """Precision@k for each odorant's strongest glomeruli, against chance."""
    ks = list(ks)
    prec = []
    for k in ks:
        overlap = [len(set(np.argsort(m.Ylog[i])[::-1][:k]) &
                       set(np.argsort(m.Yhat[i])[::-1][:k])) / k
                   for i in range(m.n_odorants)]
        prec.append(float(np.mean(overlap)))
    return ks, prec, [k / m.n_gloms for k in ks]


def per_odorant_r(m: PLSModel) -> np.ndarray:
    """Pattern correlation between predicted and observed map, per odorant."""
    return (zscore_rows(m.Ylog) * zscore_rows(m.Yhat)).mean(1)


def best_predicted(m: PLSModel, n: int = 4) -> np.ndarray:
    """Indices of the `n` odorants the model reproduces best, best first.

    These are examples, chosen by a stated rule rather than by eye. Pass a
    different `n`, or index `per_odorant_r` yourself, to show the other end.
    """
    return np.argsort(per_odorant_r(m))[::-1][:n]


def pearson(a, b) -> float:
    return float(np.corrcoef(a, b)[0, 1])


# ---------------------------------------------------------------------------
# Predicted vs experimental spectra -- Figure 6B/C/D
# ---------------------------------------------------------------------------
# Matt's question on the Figure 6 page: "Is it possible to show the same
# response spectra from Fig 3E, F, G and compare predicted high to experimental
# high?" So these are the SAME three glomeruli as 3E/F/G -- the odorant order,
# the class blocks and the glomerulus indices all come from
# `highconc.similarity`, and nothing here re-derives any of them.

#: Which space the observed/predicted comparison is drawn in.
#:
#: - ``"raw"`` undoes the model's ``log1p`` so the bars are dF/F, the same units
#:   as 3E/F/G. Honest about magnitude, but one catastrophic overshoot can own
#:   the axis: on glomerulus 211 the model predicts 193.6 against an observed
#:   66.2, and on a shared scale that single bar flattens everything else.
#: - ``"log"`` is the space the model is fitted and scored in. Correlations are
#:   higher here because this is where the loss lives -- which is exactly why it
#:   flatters. ``log1p`` expands the bottom of the range, so the many small
#:   non-zero predictions a CV regression always makes become visible bars where
#:   the truth is silent: on these three glomeruli the model "responds" above
#:   0.4 in 36 / 27 / 53 odorants against 10 / 11 / 29 observed.
#: - ``"sqrt"`` is the middle ground: it compresses the overshoot without
#:   inflating the floor, because sqrt has a finite derivative at 0 where log's
#:   is unbounded. Small spurious predictions stay small and the 193.6 bar
#:   stays on the page. It reads best of the three and is NOT the default.
#:
#: **The default is ``"log"``**: it is the space the model is fitted and scored
#: in, so the residuals drawn are the residuals the model minimised. ``raw``
#: and ``sqrt`` are transforms applied afterwards for readability, and neither
#: has a modelling meaning. Log does expand the bottom of the range, as the
#: counts above show; the notebook prints the observed and predicted support
#: counts under the figure on every run so that cost stays visible.
#:
#: `r` is reported in whichever space is drawn, and it genuinely differs
#: between them -- 0.30/0.53/0.55 raw, 0.48/0.70/0.60 sqrt, 0.53/0.78/0.64 log.
#: Log gives the highest numbers of the three, which is precisely why the
#: caption must state the space.
SPECTRUM_SPACES = ("raw", "log", "sqrt")


@dataclass(frozen=True)
class PredictionSpectrum:
    """One glomerulus's experimental vs held-out predicted high spectrum.

    The sibling of `highconc.similarity.Spectrum`, which carries low vs high
    for the same glomerulus in the same odorant order. Here the top half is
    the measured high-concentration response and the bottom half is the
    cross-validated prediction of it, so 6B/C/D read exactly like 3E/F/G with
    one series swapped.

    Attributes
    ----------
    observed, predicted : ndarray, (59,)
        High-concentration response and its held-out prediction, in plotting
        order. `predicted` is a fold-out prediction: the odorant it belongs to
        was never in the training set of the model that produced it.
    space : str
        ``"raw"`` (dF/F) or ``"log"`` (``log1p`` dF/F). See `SPECTRUM_SPACES`.
    r : float
        Pearson correlation of the two, in whichever space they are in.
    """

    glom: int
    observed: np.ndarray
    predicted: np.ndarray
    classes: np.ndarray
    names: list
    blocks: tuple
    boundaries: tuple
    cluster: int
    space: str
    r: float

    @property
    def limit(self) -> float:
        """Max over both series -- for a shared y-range across a row."""
        return float(max(self.observed.max(), self.predicted.max()))

    @property
    def responded(self) -> np.ndarray:
        """Odorants the glomerulus actually responded to at high. The
        prediction is dense by construction (a cross-validated regression
        almost never lands exactly on zero), so a support comparison is only
        meaningful against a threshold -- see `support`."""
        return self.observed > 0

    def support(self, threshold: float = 0.0) -> tuple:
        """``(n_observed, n_predicted, n_shared)`` above `threshold`."""
        o, p = self.observed > threshold, self.predicted > threshold
        return int(o.sum()), int(p.sum()), int((o & p).sum())


def _to_space(x: np.ndarray, space: str) -> np.ndarray:
    """Convert a log1p dF/F array into `space` -- the same transform
    `predicted_spectrum` applies to one glomerulus, factored out so
    `predicted_high_floor` applies it identically across the whole dataset.
    """
    if space in ("raw", "sqrt"):
        x = np.expm1(x)
    x = np.clip(x, 0, None)
    if space == "sqrt":
        x = np.sqrt(x)
    return x


def predicted_spectrum(m: PLSModel, d, glom: int, *, order="matt",
                       within="table", space="raw") -> PredictionSpectrum:
    """Fig 6B/C/D: one glomerulus, experimental high vs held-out predicted high.

    `d` supplies the odorant table the class ordering is built from, and `m`
    supplies the values -- they are two views of the same 59 odorants, and the
    row orders are checked here rather than assumed.
    """
    from . import similarity as _sim

    if space not in SPECTRUM_SPACES:
        raise ValueError(f"unknown space {space!r}; expected one of "
                         f"{list(SPECTRUM_SPACES)}")
    names = np.asarray(d.odorants["Name"].values)
    if not np.array_equal(np.asarray(m.names), names):
        raise ValueError("the model and the dataset disagree on odorant order")

    perm, classes, blocks, boundaries = _sim.class_sorted(d, order, within)
    g = int(glom)
    obs, pred = _to_space(m.Ylog[perm, g], space), _to_space(m.Yhat[perm, g], space)
    return PredictionSpectrum(
        glom=g, observed=obs, predicted=pred, classes=classes,
        names=[str(x) for x in names[perm]], blocks=blocks,
        boundaries=boundaries, cluster=int(d.cluster_number[g]),
        space=space, r=pearson(obs, pred))


def predicted_high_floor(m: PLSModel, space: str = "log") -> float:
    """Minimum predicted high-concentration value across the WHOLE dataset --
    every odorant x every glomerulus, not just one row -- in whichever space
    Fig 6B/C/D draw. Fig 6A/B/C mark this as a horizontal reference line on
    the predicted (bottom) half of each spectrum, so the three panels can be
    read against the same floor rather than each other's own minimum.
    """
    if space not in SPECTRUM_SPACES:
        raise ValueError(f"unknown space {space!r}; expected one of "
                         f"{list(SPECTRUM_SPACES)}")
    return float(_to_space(m.Yhat, space).min())


def observed_high_floor(m: PLSModel, space: str = "log") -> float:
    """Smallest NONZERO observed high-concentration response over the WHOLE
    dataset (every odorant x every glomerulus), in whichever space Fig 6C draws.

    `_to_space` clips at zero, so the plain minimum is always 0 (a glomerulus
    that did not respond) and says nothing; the smallest value that is an
    actual response is the lowest response the experiment ever measured. Fig 6C
    draws it as a horizontal line so the predicted bars can be read against it.
    """
    if space not in SPECTRUM_SPACES:
        raise ValueError(f"unknown space {space!r}; expected one of "
                         f"{list(SPECTRUM_SPACES)}")
    x = _to_space(m.Ylog, space)
    return float(x[x > 0].min())


def eg_predicted_spectra(m: PLSModel, d, gloms=None, **kw) -> list:
    """The three Fig 3E/F/G glomeruli, predicted against experimental.

    `gloms` defaults to `highconc.similarity.EG_GLOMERULI` -- imported rather
    than repeated, because 6B/C/D showing a different three than 3E/F/G is the
    one way this panel can be wrong without looking wrong.
    """
    from . import similarity as _sim

    gloms = _sim.EG_GLOMERULI if gloms is None else gloms
    return [predicted_spectrum(m, d, g, **kw) for g in gloms]


# ---------------------------------------------------------------------------
# The weighting matrix W
# ---------------------------------------------------------------------------
# W (called B in pls-explorer) is the 376x376 concentration-remapping rule the
# model learns: what a low-concentration glomerulus contributes to each
# high-concentration one. It is NOT refitted here -- it comes from the
# canonical run pls-explorer has already written to disk, alongside its two
# reference nulls, so every figure that shows W shows the same W.

#: The runs a canonical result directory holds.
WEIGHT_RUNS = ("real", "null1", "null2")

#: What each null controls for. Kept next to the loader because a figure that
#: shows the three side by side has to label them, and mislabelling which null
#: is which inverts the argument.
WEIGHT_RUN_LABELS = {
    "real": "real",
    "null1": "homogeneous-gain null",     # per-glomerulus scalar low->high fit
    "null2": "random-retuning null",      # each column copied from a random one
}


@dataclass(frozen=True)
class Weights:
    """The learned W and its nulls, with the cluster label of every glomerulus.

    Attributes
    ----------
    matrices : dict
        ``{"real": (376, 376), "null1": ..., "null2": ...}``.
    cluster_number : ndarray, (376,)
        Joel's 1-7 per glomerulus, aligned to W's rows and columns. Everything
        derived from a `Weights` therefore comes out in 1-7 order.
    run : str
        Which result directory this came from.
    """

    matrices: dict
    cluster_number: np.ndarray
    run: str

    def __getitem__(self, key: str) -> np.ndarray:
        return self.matrices[key]

    @property
    def real(self) -> np.ndarray:
        return self.matrices["real"]

    @property
    def n(self) -> int:
        return self.real.shape[0]

    def __repr__(self) -> str:
        return (f"Weights(run={self.run!r}, {self.n}x{self.n}, "
                f"{sorted(self.matrices)})")


@lru_cache(maxsize=2)
def load_weights(run: str = "baseline",
                 variant: str = DEFAULT_VARIANT) -> Weights:
    """Load W and its nulls from pls-explorer's canonical results directory.

    Raises if the run is missing rather than refitting one. A W computed on the
    fly would silently differ from the one the manuscript's numbers came from,
    and the failure would be invisible in the figure -- the same reasoning as
    `highconc.dimsweep.load_cached`.
    """
    from pls_explorer.data import load_paired

    root = results_root() / run
    if not root.is_dir():
        raise FileNotFoundError(
            f"no pls-explorer results at {root}. W is not refitted here -- "
            f"generate the '{run}' run in pls-explorer, or point "
            f"${RESULTS_ENV_VAR} at a results tree that has it."
        )
    matrices = {}
    for key in WEIGHT_RUNS:
        path = root / key / "B.npy"
        if not path.exists():
            raise FileNotFoundError(
                f"{path} is missing; the '{run}' run is incomplete "
                f"(need {', '.join(WEIGHT_RUNS)})"
            )
        matrices[key] = np.load(path)

    index = np.asarray(load_paired(variant).glom_clusters).astype(int)
    cluster_number = np.array([_clusters.to_joel(int(i)) for i in index],
                              dtype=int)
    n = matrices["real"].shape[0]
    if cluster_number.size != n:
        raise ValueError(
            f"W is {n}x{n} but variant {variant!r} has {cluster_number.size} "
            f"glomeruli -- the run and the variant disagree"
        )
    return Weights(matrices=matrices, cluster_number=cluster_number, run=run)


def diagonality(W: np.ndarray) -> float:
    """How much of W sits on its diagonal -- i.e. how close it is to pure
    per-glomerulus gain. Invariant to a shared row+column permutation, so it
    can be quoted next to a matrix drawn in any order."""
    from pls_explorer.metrics import diagonality_index

    return float(diagonality_index(W))


def top_weight_gloms(W: np.ndarray, k: int = 85) -> np.ndarray:
    """The `k` glomeruli carrying the most row+column energy in W.

    W is mostly near-dead rows and columns; drawn whole it reads as empty
    whatever its structure. Keeping the top `k` densifies the picture without
    reordering it. `k` is a display choice -- say it in the caption.
    """
    energy = (W ** 2).sum(0) + (W ** 2).sum(1)
    return np.argsort(energy)[::-1][:k]


def reorder_by_cluster(W: np.ndarray, cluster_number: np.ndarray):
    """Group W's rows and columns by cluster, in Joel's 1-7 order.

    Returns ``(W_reordered, order, boundaries, centers)`` -- the last two being
    where to draw the block lines and where to put the cluster tick labels.
    """
    from pls_explorer.block_structure import reorder_by_clusters

    W_reord, order, boundaries = reorder_by_clusters(W, cluster_number)
    sorted_c = np.asarray(cluster_number)[order]
    edges = np.concatenate([[0], boundaries, [sorted_c.size]])
    centers = (edges[:-1] + edges[1:]) / 2.0
    return W_reord, order, boundaries, centers


def block_energy(w: Weights, W: np.ndarray | None = None) -> dict:
    """Per-cluster-block energy of W, in Joel's 1-7 order.

    Returns pls-explorer's decomposition plus `block_signed`: per block, the
    positive-weight energy minus the negative-weight energy, which is what says
    whether a cluster pair is net excitatory or net suppressive. `block_norm`
    (a Frobenius norm) cannot say that -- it is non-negative by construction.

    Computed on the FULL W, never on a densified subset: this is the
    quantitative claim, not the picture.
    """
    from pls_explorer.block_structure import decompose_block_energy

    W = w.real if W is None else W
    out = dict(decompose_block_energy(W, w.cluster_number))
    numbers = np.asarray(out["unique_row"])
    signed = np.sign(W) * W ** 2
    out["block_signed"] = np.array(
        [[signed[np.ix_(w.cluster_number == r, w.cluster_number == c)].sum()
          for c in numbers] for r in numbers]
    )
    out["numbers"] = numbers
    return out


def top_off_diagonal(decomp: dict, k: int = 3) -> list:
    """The `k` strongest cross-cluster blocks, as ``(from, to, norm)`` with
    cluster numbers in Joel's 1-7."""
    from pls_explorer.block_structure import top_off_diagonal_blocks

    return [(int(a), int(b), float(v)) for a, b, v in
            top_off_diagonal_blocks(decomp["block_norm"], decomp["numbers"], k=k)]


# ---------------------------------------------------------------------------
# Component selection, LV loadings, and the model zoo
#
# All three are READ from the canonical pls-explorer run rather than recomputed,
# for the same reason as `load_weights`: a sweep or a zoo refitted on the fly
# would silently differ from the numbers the manuscript quotes, and nothing in
# the figure would show it. Missing inputs raise.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ComponentSelection:
    """The n-components CV sweep and the two rules read off it.

    `one_sigma` is the smallest n within one standard error of the best Q^2 --
    the conservative rule, and the one the manuscript follows. `van_der_voet`
    is the randomization test, which here picks a larger n; they disagree, and
    `chosen` records which was taken.
    """

    grid: np.ndarray
    q2: np.ndarray
    rmsecv: np.ndarray
    rmsec: np.ndarray
    chosen: int
    one_sigma: int
    van_der_voet: int
    method: str
    run: str

    @property
    def overfit(self) -> np.ndarray:
        """RMSECV / RMSEC -- 1.0 is no optimism, and it climbs with n."""
        return self.rmsecv / np.where(self.rmsec > 0, self.rmsec, np.nan)

    def at(self, n: int) -> dict:
        """Every quantity at one n, for annotating a chosen point."""
        i = int(np.where(self.grid == n)[0][0])
        return dict(n=n, q2=float(self.q2[i]), rmsecv=float(self.rmsecv[i]),
                    rmsec=float(self.rmsec[i]), overfit=float(self.overfit[i]))


@lru_cache(maxsize=2)
def component_selection(run: str = "baseline") -> ComponentSelection:
    """Load the component-selection sweep written by the canonical run."""
    import json

    path = results_root() / run / "component_selection.json"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing -- the component sweep is not run here; "
            f"generate the '{run}' run in pls-explorer, or point "
            f"${RESULTS_ENV_VAR} at a results tree that has it."
        )
    raw = json.loads(path.read_text())
    s = raw["selection"]
    return ComponentSelection(
        grid=np.asarray(raw["grid"]), q2=np.asarray(raw["q2"], float),
        rmsecv=np.asarray(raw["rmsecv"], float),
        rmsec=np.asarray(raw["rmsec"], float),
        chosen=int(s["chosen"]), one_sigma=int(s["one_sigma"]),
        van_der_voet=int(s["van_der_voet"]), method=str(s["method"]), run=run,
    )


@dataclass(frozen=True)
class Loadings:
    """One latent component's X- and Y-loadings, per glomerulus.

    `cluster_number` is Joel's 1-7, aligned to `x` and `y`, so anything grouped
    or coloured from a `Loadings` comes out in 1-7 order.
    """

    x: np.ndarray
    y: np.ndarray
    cluster_number: np.ndarray
    component: int
    run: str

    @property
    def outer(self) -> np.ndarray:
        """The rank-1 x ⊗ y map this component contributes to W."""
        return np.outer(self.x, self.y)

    def take(self, index) -> "Loadings":
        """The same component restricted (and reordered) to `index`."""
        idx = np.asarray(index, int)
        return Loadings(x=self.x[idx], y=self.y[idx],
                        cluster_number=self.cluster_number[idx],
                        component=self.component, run=self.run)


@lru_cache(maxsize=4)
def lv_loadings(run: str = "baseline", component: int = 0,
                variant: str = DEFAULT_VARIANT) -> Loadings:
    """Load one latent component's loadings from the canonical run.

    `component` is 0-based: LV1 is ``component=0``.
    """
    from pls_explorer.data import load_paired

    path = results_root() / run / "real" / "loadings.npz"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing -- loadings are not refitted here; "
            f"generate the '{run}' run in pls-explorer, or point "
            f"${RESULTS_ENV_VAR} at a results tree that has it."
        )
    ld = np.load(path)
    index = np.asarray(load_paired(variant).glom_clusters).astype(int)
    cluster_number = np.array([_clusters.to_joel(int(i)) for i in index],
                              dtype=int)
    x = np.asarray(ld["x_loadings"][:, component], float)
    if cluster_number.size != x.size:
        raise ValueError(
            f"loadings have {x.size} glomeruli but variant {variant!r} has "
            f"{cluster_number.size} -- the run and the variant disagree"
        )
    return Loadings(x=x, y=np.asarray(ld["y_loadings"][:, component], float),
                    cluster_number=cluster_number, component=int(component),
                    run=run)


def group_by_cluster(cluster_number):
    """Sort glomeruli into Joel's 1-7 blocks.

    Returns ``(order, boundaries, centers, numbers)`` -- the permutation, where
    to draw the block lines, where to put the tick labels, and which cluster
    each block is. The 1-D counterpart of `reorder_by_cluster`; a marginal bar
    plot and the heatmap it sits beside must use the SAME order, so they take
    it from one call.
    """
    c = np.asarray(cluster_number, int)
    order = np.argsort(c, kind="stable")
    sorted_c = c[order]
    boundaries = np.where(np.diff(sorted_c) != 0)[0] + 1
    edges = np.concatenate([[0], boundaries, [sorted_c.size]])
    centers = (edges[:-1] + edges[1:]) / 2.0
    return order, boundaries, centers, np.unique(sorted_c)


#: Held-out Q^2 for the model zoo, family-stratified 4-fold CV, global and on
#: log-transformed Y. These are CACHED headline values from pls-explorer's
#: notebooks/03_model_zoo.ipynb, not refitted here -- the zoo takes minutes and
#: the figure's whole point is the ranking, which is fixed once the run is.
#: Keep sorted descending; the panel draws them top to bottom as given.
MODEL_ZOO: tuple[tuple[str, float], ...] = (
    ("RidgeCV", +0.343),
    ("KernelRidge Poly2", +0.336),
    ("PLS-R k=15", +0.304),
    ("RandomForest", +0.276),
    ("PLS-R k=8 (canonical)", +0.261),
    ("PCR k=8", +0.231),
    ("ExtraTrees", +0.206),
    ("KernelRidge RBF", +0.189),
    ("IDENTITY (Y=X)", +0.147),
    ("MLP (64)", +0.138),
    ("kNN k=3", +0.087),
    ("MEAN baseline", -0.028),
)

#: Which zoo entry is the model the manuscript uses. Everything else in the
#: panel is coloured relative to it.
ZOO_CANONICAL = "PLS-R k=8 (canonical)"


def zoo_baseline(zoo=MODEL_ZOO, canonical: str = ZOO_CANONICAL) -> float:
    """The canonical model's Q^2 -- the line the rest of the zoo is read against."""
    for name, q2 in zoo:
        if name == canonical:
            return float(q2)
    raise KeyError(f"{canonical!r} is not in the zoo: {[n for n, _ in zoo]}")
