"""Concentration flow in glomerular tuning space.

The dynamic reading of the cluster-aligned weighting matrix W (see
`highconc.pls.block_energy`): each glomerulus is a point in tuning space -- its
response vector across the 59 odorants -- and its low- and high-concentration
tuning vectors are two points, the displacement between them its concentration
trajectory. The claim is that similarly tuned glomeruli move in COHERENT
directions.

Two independent measures of that, both here:

- `coherence` is the all-pairs version: the mean cosine between displacement
  vectors, within cluster vs. between, against a cluster-label shuffle null.
- `centroid_coherence` compares each glomerulus's direction to the seven
  cluster-MEAN directions instead. The centroid is a denoised reference, so
  this is the cleaner statement, and it also yields the 7x7 cross-coherence
  matrix: how well cluster a's glomeruli align with cluster b's mean.

Everything comes out in Joel's 1-7. The embedding is PCA fit on LOW-concentration
tuning with high projected in -- fitting on both, or separately, would let the
axes move with concentration and the displacement would no longer mean what it
says.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from . import clusters as _clusters

#: Same paired variant the PLS side uses, so the two analyses describe the same
#: glomeruli. See `highconc.pls.DEFAULT_VARIANT`.
DEFAULT_VARIANT = "nonresponders_dropped"

#: Displacement vectors shorter than this percentile of ||low->high|| are
#: dropped before measuring alignment. A near-zero displacement has an
#: essentially random direction, so keeping them adds noise to every cosine
#: without adding a single glomerulus's worth of signal.
MIN_DISPLACEMENT_PCT = 10

#: Default PCs for the embedding. Three is what the 3-D figure draws; the
#: alignment result holds in the full 59-D space too, but concentration of
#: measure pins every angle near 90 degrees there, so a low-dimensional
#: embedding is what can actually be shown.
N_PCS = 3

#: Shuffle counts for the two nulls. The all-pairs null is cheap; the centroid
#: null refits seven centroids per draw, so it runs fewer.
N_NULL_PAIRS, N_NULL_CENTROID = 1000, 500
NULL_SEED = 0

#: Display order for the 7x7 coherence matrix, in Joel's 1-7, with the split
#: after the third entry. NOT the numeric order and NOT arbitrary: hierarchical
#: clustering of the coherence matrix puts {1, 3, 2} and {6, 5, 4, 7} into two
#: super-groups, and the matrix only reads as blocked when drawn that way.
COHERENCE_ORDER = (1, 3, 2, 6, 5, 4, 7)
COHERENCE_SPLIT = 3


@dataclass(frozen=True)
class Tuning:
    """Paired low/high tuning matrices, log1p, with Joel's 1-7 per glomerulus.

    `low` and `high` are (n_odorants, n_gloms) -- odorants down the rows, which
    is the orientation a tuning VECTOR is a column of.
    """

    low: np.ndarray
    high: np.ndarray
    cluster_number: np.ndarray
    variant: str

    @property
    def displacement(self) -> np.ndarray:
        """high - low, the concentration trajectory of every glomerulus."""
        return self.high - self.low

    @property
    def n_gloms(self) -> int:
        return self.low.shape[1]

    @property
    def n_odorants(self) -> int:
        return self.low.shape[0]


@lru_cache(maxsize=2)
def load(variant: str = DEFAULT_VARIANT) -> Tuning:
    """Load the paired tuning matrices and attach Joel's cluster numbers."""
    from pls_explorer.data import load_paired
    from pls_explorer.preprocess import log1p_positive

    p = load_paired(variant)
    index = np.asarray(p.glom_clusters).astype(int)
    return Tuning(
        low=log1p_positive(p.X.to_numpy(float)),
        high=log1p_positive(p.Y.to_numpy(float)),
        cluster_number=np.array([_clusters.to_joel(int(i)) for i in index],
                                dtype=int),
        variant=variant,
    )


def keep_moving(t: Tuning, percentile: float = MIN_DISPLACEMENT_PCT) -> np.ndarray:
    """Boolean mask of glomeruli whose displacement is long enough to have a
    meaningful direction. One definition, shared by every measure below, so the
    figure's panels are all computed on the same set."""
    norm = np.linalg.norm(t.displacement, axis=0)
    return norm > np.percentile(norm, percentile)


@dataclass(frozen=True)
class Coherence:
    """All-pairs displacement alignment, within cluster vs. between."""

    within: float
    between: float
    sigma: float
    p: float
    keep: np.ndarray

    @property
    def n_kept(self) -> int:
        return int(self.keep.sum())

    @property
    def gap(self) -> float:
        return self.within - self.between


def coherence(t: Tuning, *, percentile: float = MIN_DISPLACEMENT_PCT,
              n_null: int = N_NULL_PAIRS, seed: int = NULL_SEED) -> Coherence:
    """Mean cosine between displacement directions, within vs. between cluster.

    `sigma` and `p` are against a cluster-label shuffle, which holds the
    displacement vectors fixed and only re-assigns who is in which cluster --
    so it tests the cluster structure, not the fact that the vectors are
    correlated at all. Computed over ALL SEVEN clusters, whatever subset a
    figure chooses to draw.
    """
    keep = keep_moving(t, percentile)
    D = t.displacement[:, keep]
    Dn = D / np.linalg.norm(D, axis=0)
    clk = t.cluster_number[keep]

    C = Dn.T @ Dn
    iu = np.triu_indices(C.shape[0], 1)
    same = (clk[:, None] == clk[None, :])[iu]
    cos = C[iu]
    w, b = float(cos[same].mean()), float(cos[~same].mean())

    rng = np.random.default_rng(seed)
    null = np.empty(n_null)
    for k in range(n_null):
        cs = rng.permutation(clk)
        sm = (cs[:, None] == cs[None, :])[iu]
        null[k] = cos[sm].mean() - cos[~sm].mean()
    return Coherence(within=w, between=b,
                     sigma=float(((w - b) - null.mean()) / null.std()),
                     p=float(np.mean(null >= (w - b))), keep=keep)


@dataclass(frozen=True)
class Embedding:
    """Tuning space in `n` PCs, fit on low concentration, high projected in."""

    low: np.ndarray
    high: np.ndarray
    explained: np.ndarray
    pca: object

    @property
    def displacement(self) -> np.ndarray:
        return self.high - self.low


def embed(t: Tuning, n_components: int = N_PCS) -> Embedding:
    """PCA on low-concentration tuning (glomeruli as samples); project high in.

    One basis for both concentrations. Fitting the two separately would rotate
    the axes between them and the "displacement" would be partly a change of
    coordinates.
    """
    from sklearn.decomposition import PCA

    low_t, high_t = t.low.T, t.high.T            # gloms x odorants
    mu = low_t.mean(0)
    pca = PCA(n_components=n_components).fit(low_t - mu)
    return Embedding(low=pca.transform(low_t - mu),
                     high=pca.transform(high_t - mu),
                     explained=pca.explained_variance_ratio_, pca=pca)


@dataclass(frozen=True)
class CentroidCoherence:
    """Each glomerulus's direction against the seven cluster-mean directions.

    `matrix` is 7x7 in Joel's 1-7: row a, column b is the mean cosine between
    cluster a's glomeruli and cluster b's mean direction. Its diagonal is
    within-cluster coherence, its off-diagonal cross-coherence. `across` is the
    chance level the matrix should be read against, and the natural pivot for a
    diverging colour scale.
    """

    within: float
    across: float
    within_angle: np.ndarray
    across_angle: np.ndarray
    matrix: np.ndarray
    numbers: np.ndarray
    sigma: float
    keep: np.ndarray
    cluster_number: np.ndarray = field(repr=False)

    def angles_by_cluster(self) -> list:
        """Per-cluster own-centroid angles, in Joel order -- one box each."""
        return [self.within_angle[self.cluster_number == n] for n in self.numbers]

    def reorder(self, order=COHERENCE_ORDER):
        """`(matrix, labels)` in a chosen cluster order -- default the
        coherence super-groups, see `COHERENCE_ORDER`."""
        pos = {int(n): i for i, n in enumerate(self.numbers)}
        idx = [pos[int(n)] for n in order]
        return self.matrix[np.ix_(idx, idx)], [str(int(n)) for n in order]


def centroid_coherence(t: Tuning, emb: "Embedding | None" = None, *,
                       percentile: float = MIN_DISPLACEMENT_PCT,
                       n_components: int = N_PCS,
                       n_null: int = N_NULL_CENTROID,
                       seed: int = NULL_SEED) -> CentroidCoherence:
    """Alignment of every glomerulus to each cluster-mean direction.

    Measured in the embedded space the 3-D figure draws, not the full 59-D one:
    the direction of the result is identical, but in 59-D concentration of
    measure pins every angle near 90 degrees and nothing is legible.

    `percentile` is applied once, globally, over all 376 glomeruli's
    displacement norms -- not per cluster. Raises if that happens to drop an
    entire cluster (`clk == a` empty for some cluster `a`): a mean over zero
    glomeruli is not a coherence measurement, and letting it through as a
    silent NaN turns panel D's row/column for that cluster blank rather than
    surfacing the problem.
    """
    emb = embed(t, n_components) if emb is None else emb
    keep = keep_moving(t, percentile)
    numbers = np.arange(1, _clusters.N_CLUSTERS + 1, dtype=int)

    dv = emb.displacement[keep]
    Dn = dv / (np.linalg.norm(dv, axis=1, keepdims=True) + 1e-12)
    clk = t.cluster_number[keep]
    n = clk.size
    col = np.array([int(np.where(numbers == c)[0][0]) for c in clk])

    empty = [int(a) for a in numbers if not (clk == a).any()]
    if empty:
        raise ValueError(
            f"cluster(s) {empty} have no glomeruli left after the "
            f"{percentile}th-percentile displacement filter ({keep.sum()}/"
            f"{keep.size} kept overall) -- lower `percentile` "
            "(MIN_DISPLACEMENT_PCT) so every cluster keeps at least one "
            "glomerulus, rather than silently NaN-ing that cluster's row/"
            "column in the coherence matrix")

    cent = np.array([dv[clk == num].mean(0) for num in numbers])
    cent = cent / (np.linalg.norm(cent, axis=1, keepdims=True) + 1e-12)
    cos_gc = Dn @ cent.T                            # n x 7

    other = np.ones((n, numbers.size), bool)
    other[np.arange(n), col] = False
    within_cos, across_cos = cos_gc[np.arange(n), col], cos_gc[other]
    w, b = float(within_cos.mean()), float(across_cos.mean())

    rng = np.random.default_rng(seed)
    null = np.empty(n_null)
    for k in range(n_null):
        cs = rng.permutation(col)
        cu = np.array([dv[cs == j].mean(0) for j in range(numbers.size)])
        cu = cu / (np.linalg.norm(cu, axis=1, keepdims=True) + 1e-12)
        cg = Dn @ cu.T
        om = np.ones((n, numbers.size), bool)
        om[np.arange(n), cs] = False
        null[k] = cg[np.arange(n), cs].mean() - cg[om].mean()

    M = np.array([[cos_gc[clk == a, j].mean() for j in range(numbers.size)]
                  for a in numbers])
    return CentroidCoherence(
        within=w, across=b,
        within_angle=_degrees(within_cos), across_angle=_degrees(across_cos),
        matrix=M, numbers=numbers,
        sigma=float(((w - b) - null.mean()) / null.std()), keep=keep,
        cluster_number=clk)


def _degrees(cos) -> np.ndarray:
    return np.degrees(np.arccos(np.clip(cos, -1, 1)))


def cluster_flow(emb: Embedding, cluster_number, numbers=None,
                 spread_pct: float = 68, max_half_angle: float = 55):
    """Per-cluster mean trajectory and the spread around it.

    Yields one dict per cluster with at least three members: `number`, the
    starting centroid `origin`, the mean displacement `mean`, its unit
    direction `unit`, the `half_angle` containing `spread_pct` of the cluster's
    directions (capped at `max_half_angle`, or a wide cluster's cone swallows
    the figure), and the median displacement `length`. That cone is a spread,
    not a confidence interval -- it says how tightly the cluster moves
    together, not how well its mean is estimated.
    """
    cluster_number = np.asarray(cluster_number, int)
    numbers = np.unique(cluster_number) if numbers is None else numbers
    dv = emb.displacement
    out = []
    for num in numbers:
        m = cluster_number == num
        if m.sum() < 3:
            continue
        mean = dv[m].mean(0)
        norm = np.linalg.norm(mean)
        if norm < 1e-9:
            continue
        unit = mean / norm
        units = dv[m] / (np.linalg.norm(dv[m], axis=1, keepdims=True) + 1e-12)
        half = min(np.percentile(_degrees(units @ unit), spread_pct),
                   max_half_angle)
        out.append(dict(number=int(num), origin=emb.low[m].mean(0), mean=mean,
                        unit=unit, half_angle=float(half),
                        length=float(np.median(np.linalg.norm(dv[m], axis=1)))))
    return out
