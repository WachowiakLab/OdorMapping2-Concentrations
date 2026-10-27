"""Motif transfer between glomerular groups — the analysis behind Figure 7.

The question the figure asks is directional and pairwise. For an ordered pair
of glomerular groups ``F -> T``: does knowing group **F**'s mean *low*-
concentration response improve held-out prediction of group **T**'s glomeruli
at *high*, beyond what T's own low-concentration mean already predicts?

    intra    = predict T's high from  [T's low mean]
    intra+F  = predict T's high from  [T's low mean, F's low mean]
    M[F, T]  = Q2(intra+F) - Q2(intra)          # out-of-fold, pooled over T

A positive ``M[F, T]`` says a tuning motif carried by F at low concentration
reappears inside T at high — that is what "motif transfer" means here. The
matrix is deliberately *not* symmetric: F helping T does not imply T helps F,
and the one-way edges are the result.

Two groupings, the same test:

* ``"clusters"``   — the seven NMF clusters, in Joel's 1-7 order.
* ``"chem_class"`` — each glomerulus assigned to the chemical class evoking its
  largest mean low-concentration response. A control: if the transfers survive
  a grouping that never saw the NMF factorisation, they are not an artefact of
  it.

    >>> from highconc import pls, transfer
    >>> t = transfer.pairwise(pls.fit())
    >>> t.M.shape, len(t.edges())
    ((7, 7), 8)

**Significance is permutation-based and the nulls are cached.** Each null draw
recomputes the whole matrix under a *broken pairing* — the odorant order of Y
is shuffled against X, which preserves the marginal distribution of high
responses and all correlation structure within X and within Y while destroying
the low-to-high correspondence. Ten thousand draws take 11 minutes for the
clusters and 16 for the classes, so both live in ``data/pairwise/*.npz`` and
:func:`pairwise` reads them. Nothing in a notebook run recomputes a null;
:func:`compute_null` exists to rebuild a cache deliberately, and says so.

The broken-pairing null is the one the manuscript publishes, and it is a
choice: a size-matched control run on 2026-08-20 leaves only three of the eight
cluster edges standing. That caveat is recorded in the methods audit; it is not
what this figure shows. See the ``null`` parameter.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import clusters as _clusters

#: Ridge penalty on the group-mean predictors. The predictors are one or two
#: columns of 59 samples, so this is regularisation against collinearity
#: between two group means, not against dimensionality.
ALPHA = 1.0

#: Permutation count behind the shipped caches.
N_PERM = 10_000

#: Seed for the permutation sequence. Fixed, and reproducing the shipped
#: cluster cache bit-for-bit depends on it.
SEED = 42

#: FDR level. Edges are called significant at BH q < ALPHA_FDR across all
#: ordered off-diagonal pairs — 42 of them for clusters, 72 for classes.
ALPHA_FDR = 0.05

#: Minimum glomeruli for a chemical class to be tested. Nine of the ten classes
#: clear 15 (372 of 376 glomeruli); the ester-preferring group has 4 and is
#: excluded, because a group mean over four glomeruli is not a population.
CHEM_MIN_SIZE = 15

#: The two groupings, and the cache each one reads.
GROUPINGS = ("clusters", "chem_class")

_CACHE_DIR = Path(__file__).resolve().parent.parent / "data" / "pairwise"
_CACHE_FILES = {
    "clusters": _CACHE_DIR / "pairwise_pvalues.npz",
    "chem_class": _CACHE_DIR / "pairwise_chem_class.npz",
}

#: Null models. Only ``"broken_pairing"`` ships a cache; the name is explicit
#: so that a future size-matched run has somewhere to go without the shipped
#: results silently changing meaning.
NULLS = ("broken_pairing",)


# ---------------------------------------------------------------------------
# Groupings
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Grouping:
    """Which glomeruli form which group, already in display order.

    Attributes
    ----------
    id : str
        ``"clusters"`` or ``"chem_class"``.
    columns : list of ndarray
        Glomerulus column indices per group, in display order.
    labels : tuple of str
        Axis labels, in the same order — ``"1".."7"`` or class names.
    keys : tuple
        The underlying group key per position (cluster number, class letter).
    colors : tuple of str
        One colour per group, from the canonical cluster/class palettes.
    dropped : int
        Glomeruli excluded because their group was below `CHEM_MIN_SIZE`.
    """

    id: str
    columns: list
    labels: tuple
    keys: tuple
    colors: tuple
    dropped: int = 0

    @property
    def n_groups(self) -> int:
        return len(self.columns)

    @property
    def sizes(self) -> tuple:
        return tuple(len(c) for c in self.columns)

    @property
    def n_pairs(self) -> int:
        """Ordered off-diagonal pairs — what BH corrects across."""
        return self.n_groups * (self.n_groups - 1)

    def __repr__(self) -> str:
        return (f"Grouping({self.id!r}, {self.n_groups} groups, "
                f"sizes={self.sizes}, dropped={self.dropped})")


def cluster_grouping(m) -> Grouping:
    """The seven NMF clusters, **in Joel's 1-7 order**.

    The stored analysis ran in our raw component-index order (historically
    displayed A-G); this is where the renumbering happens, once, so no panel
    and no cache ever has to think about it again.
    """
    idx = np.asarray(m.cluster_index)
    keys = tuple(range(1, _clusters.N_CLUSTERS + 1))
    columns = [np.where(idx == _clusters.to_index(n))[0] for n in keys]
    return Grouping(
        id="clusters", columns=columns,
        labels=tuple(str(n) for n in keys), keys=keys,
        colors=tuple(_rgba_hex(_clusters_color(n)) for n in keys),
    )


def chem_class_grouping(m, min_size: int = CHEM_MIN_SIZE) -> Grouping:
    """Each glomerulus grouped by the class evoking its largest low mean.

    Assignment is over `m.Xlog` — the same log space the test is run in — and
    is a hard argmax, so every glomerulus lands in exactly one class. Classes
    below `min_size` are dropped entirely rather than pooled: pooling would
    invent a group ("other") that has no chemical meaning and would then be
    tested as though it did.

    Class order is `highconc.clusters.CLASS_ORDER_TABLE`, restricted to the
    classes that survive.
    """
    chem = np.asarray(m.chem)
    present = [c for c in _clusters.CLASS_ORDER_TABLE if (chem == c).any()]
    means = np.stack([m.Xlog[chem == c].mean(0) for c in present])
    preferred = np.asarray(present)[means.argmax(0)]

    keys, columns = [], []
    dropped = 0
    for c in present:
        col = np.where(preferred == c)[0]
        if len(col) < min_size:
            dropped += len(col)
            continue
        keys.append(c)
        columns.append(col)
    return Grouping(
        id="chem_class", columns=columns,
        labels=tuple(_clusters.CHEM_CLASS_NAMES[c] for c in keys),
        keys=tuple(keys),
        colors=tuple(_clusters.CHEM_CLASS_COLORS[c] for c in keys),
        dropped=dropped,
    )


def grouping(m, which: str = "clusters", **kw) -> Grouping:
    """Dispatch to :func:`cluster_grouping` or :func:`chem_class_grouping`."""
    if which == "clusters":
        return cluster_grouping(m)
    if which == "chem_class":
        return chem_class_grouping(m, **kw)
    raise ValueError(f"grouping must be one of {GROUPINGS}, got {which!r}")


def _clusters_color(n):
    from . import style as hs
    return hs.cluster_color(n)


def _rgba_hex(c):
    from matplotlib.colors import to_hex
    return to_hex(c)


# ---------------------------------------------------------------------------
# The test itself
# ---------------------------------------------------------------------------
def group_means(m, g: Grouping) -> np.ndarray:
    """``(59, K)`` mean low-concentration response per group — the predictors."""
    return np.column_stack([m.Xlog[:, col].mean(1) for col in g.columns])


def increments(m, g: Grouping, Ymat: np.ndarray | None = None,
               *, alpha: float = ALPHA) -> np.ndarray:
    """``M[F, T]`` = out-of-fold Q2 gain from adding source F to target T.

    Rows are sources, columns are targets, diagonal is NaN. Pass `Ymat` to
    score a permuted response matrix — that is how the null is built.

    Both models are fitted inside `m.folds`, the same family-stratified splits
    every other model in this bundle is scored on, and Q2 is pooled over all of
    T's glomeruli at once (`q2_global`), not averaged per glomerulus: a group's
    prediction is a population claim.
    """
    from sklearn.linear_model import Ridge
    from pls_explorer.metrics import q2_global

    Y = m.Ylog if Ymat is None else Ymat
    K = g.n_groups
    CM = group_means(m, g)
    n = Y.shape[0]

    M = np.full((K, K), np.nan)
    for T in range(K):
        gl = g.columns[T]
        intra = np.zeros((n, len(gl)))
        plus = {F: np.zeros((n, len(gl))) for F in range(K) if F != T}
        for tr, te in m.folds:
            dT = CM[:, [T]]
            intra[te] = Ridge(alpha).fit(dT[tr], Y[np.ix_(tr, gl)]).predict(dT[te])
            for F in plus:
                dTF = CM[:, [T, F]]
                plus[F][te] = Ridge(alpha).fit(dTF[tr],
                                               Y[np.ix_(tr, gl)]).predict(dTF[te])
        q0 = q2_global(Y[:, gl], intra)
        for F in plus:
            M[F, T] = q2_global(Y[:, gl], plus[F]) - q0
    return M


def intra_baseline(m, g: Grouping, *, alpha: float = ALPHA) -> np.ndarray:
    """``(K,)`` out-of-fold Q2 of each group predicted from its own low mean.

    The denominator for both percentage readings of `M`. Fast — no
    permutations — so it is recomputed rather than cached.
    """
    from sklearn.linear_model import Ridge
    from pls_explorer.metrics import q2_global

    CM = group_means(m, g)
    q2 = np.zeros(g.n_groups)
    for T, gl in enumerate(g.columns):
        intra = np.zeros((m.Ylog.shape[0], len(gl)))
        for tr, te in m.folds:
            dT = CM[:, [T]]
            intra[te] = Ridge(alpha).fit(dT[tr],
                                         m.Ylog[np.ix_(tr, gl)]).predict(dT[te])
        q2[T] = q2_global(m.Ylog[:, gl], intra)
    return q2


def bh_fdr(pvals) -> np.ndarray:
    """Benjamini-Hochberg q-values for a flat array of p-values."""
    p = np.asarray(pvals, float)
    mlen = len(p)
    order = np.argsort(p)
    q = np.empty(mlen)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        i = mlen - rank
        prev = min(prev, p[idx] * mlen / i)
        q[idx] = prev
    return q


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Transfer:
    """One motif-transfer test, ready to draw.

    All matrices are ``(K, K)`` with **rows = source F, columns = target T**
    and a NaN diagonal, already in the grouping's display order.

    Attributes
    ----------
    M : ndarray
        Raw increment, ``Q2(intra+F) - Q2(intra)``.
    P : ndarray
        One-sided permutation p, ``(1 + #{null >= observed}) / (1 + n_perm)``.
    Q : ndarray
        BH q-values across the off-diagonal pairs.
    q2_intra : ndarray, (K,)
        Each target's own-mean baseline.
    """

    grouping: Grouping
    M: np.ndarray
    P: np.ndarray
    Q: np.ndarray
    q2_intra: np.ndarray
    n_perm: int = N_PERM
    null: str = "broken_pairing"

    @property
    def labels(self) -> tuple:
        return self.grouping.labels

    @property
    def off(self) -> np.ndarray:
        """Boolean mask of the ordered off-diagonal pairs."""
        return ~np.eye(self.grouping.n_groups, dtype=bool)

    @property
    def p_floor(self) -> float:
        """Smallest p the permutation count can resolve."""
        return 1.0 / (self.n_perm + 1)

    @property
    def pct_var(self) -> np.ndarray:
        """`M` as a percentage of the target's own baseline Q2.

        "How much more of T's variance does F explain." Blows up when a target
        has a small baseline, which is exactly when it is least meaningful.
        """
        return 100.0 * self.M / self.q2_intra[None, :]

    @property
    def pct_err(self) -> np.ndarray:
        """`M` as a percentage of the target's *remaining* error.

        "How much of what T's own mean could not explain does F explain."
        Bounded by construction and the default the figure shows — see the
        `metric` parameter on `panels.transfer.improvement_matrix`.
        """
        return 100.0 * self.M / (1.0 - self.q2_intra[None, :])

    def significant(self, alpha: float = ALPHA_FDR) -> np.ndarray:
        return (self.Q < alpha) & self.off

    def edges(self, alpha: float = ALPHA_FDR) -> list:
        """Significant transfers, strongest increment first.

        Each entry is a dict with ``source``, ``target`` (labels), ``dq2``,
        ``pct_var``, ``pct_err``, ``p``, ``q`` and ``reciprocal`` — the last
        True when the reverse edge is significant too, which is the difference
        between a sink and a two-way exchange and is a result in its own right.
        """
        sig = self.significant(alpha)
        out = []
        for F, T in zip(*np.where(sig)):
            out.append(dict(
                source=self.labels[F], target=self.labels[T],
                source_index=int(F), target_index=int(T),
                dq2=float(self.M[F, T]),
                pct_var=float(self.pct_var[F, T]),
                pct_err=float(self.pct_err[F, T]),
                p=float(self.P[F, T]), q=float(self.Q[F, T]),
                reciprocal=bool(sig[T, F]),
            ))
        return sorted(out, key=lambda e: -e["dq2"])

    def summary(self, alpha: float = ALPHA_FDR) -> dict:
        e = self.edges(alpha)
        return dict(
            grouping=self.grouping.id, n_groups=self.grouping.n_groups,
            n_pairs=self.grouping.n_pairs, n_perm=self.n_perm,
            n_significant=len(e),
            n_one_way=sum(not x["reciprocal"] for x in e),
            p_floor=self.p_floor,
        )

    def __repr__(self) -> str:
        return (f"Transfer({self.grouping.id!r}, {self.grouping.n_groups} groups, "
                f"{len(self.edges())} of {self.grouping.n_pairs} significant, "
                f"N={self.n_perm:,})")


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------
def cache_path(which: str = "clusters") -> Path:
    if which not in _CACHE_FILES:
        raise ValueError(f"grouping must be one of {GROUPINGS}, got {which!r}")
    return _CACHE_FILES[which]


def pairwise(m, which: str = "clusters", *, alpha_fdr: float = ALPHA_FDR,
             null: str = "broken_pairing") -> Transfer:
    """Load the cached permutation result for one grouping.

    Reads ``M_real``, ``P`` and the raw ``null`` draws from disk and recomputes
    only the intra baselines, which are seconds. **Never runs permutations** —
    if the cache is missing this raises and tells you which
    :func:`compute_null` call would build it, because a silent 16-minute
    recompute inside a notebook is worse than an error.

    The cluster cache was written in our raw component-index order; the stored
    ``labels`` say which, and the matrices are permuted into the grouping's
    display order (Joel 1-7) on the way out.
    """
    if null not in NULLS:
        raise ValueError(f"null must be one of {NULLS}, got {null!r}")
    path = cache_path(which)
    if not path.exists():
        raise FileNotFoundError(
            f"No cached permutation null at {path}.\n"
            f"Build it with:  highconc.transfer.compute_null(pls.fit(), "
            f"{which!r})   (~{11 if which == 'clusters' else 16} minutes)"
        )
    g = grouping(m, which)
    d = np.load(path, allow_pickle=True)
    stored = [str(s) for s in d["labels"]]
    perm = _reindex(stored, g)

    M = np.asarray(d["M_real"], float)[np.ix_(perm, perm)]
    P = np.asarray(d["P"], float)[np.ix_(perm, perm)]
    n_perm = int(d["n_perm"]) if "n_perm" in d.files else int(d["null"].shape[0])

    off = ~np.eye(g.n_groups, dtype=bool)
    Q = np.full((g.n_groups, g.n_groups), np.nan)
    Q[off] = bh_fdr(P[off])
    return Transfer(grouping=g, M=M, P=P, Q=Q,
                    q2_intra=intra_baseline(m, g), n_perm=n_perm, null=null)


def _reindex(stored_labels, g: Grouping) -> list:
    """Row/column permutation taking the cache's order into display order.

    Cluster caches store the legacy A-G labels (our component index); class
    caches store class letters. Anything else is a mismatch worth failing on.
    """
    if g.id == "clusters":
        legacy = list(_clusters.LEGACY_LABELS)
        if stored_labels == legacy:
            return [_clusters.to_index(n) for n in g.keys]
        if stored_labels == list(g.labels):
            return list(range(g.n_groups))
    else:
        if stored_labels == list(g.keys):
            return list(range(g.n_groups))
        if stored_labels == list(g.labels):
            return list(range(g.n_groups))
    raise ValueError(
        f"cached labels {stored_labels} do not match grouping {g.id!r} "
        f"(keys {list(g.keys)}). The cache was built for a different grouping."
    )


def compute_null(m, which: str = "clusters", *, n_perm: int = N_PERM,
                 seed: int = SEED, path: Path | None = None,
                 progress: bool = True, n_jobs: int = 1) -> Path:
    """Run the broken-pairing permutation null and write the cache.

    **Slow on purpose to call, not to run**: 10,000 draws is 11 minutes for the
    clusters, 16 for the classes. Nothing calls this automatically.

    The permutations are drawn up-front from ``default_rng(seed)`` in the same
    sequence the original 2026-06-04 script drew them one at a time, so a
    ``n_jobs > 1`` rebuild is bit-for-bit identical to the shipped serial one.
    """
    from joblib import Parallel, delayed

    g = grouping(m, which)
    path = Path(path) if path is not None else cache_path(which)
    rng = np.random.default_rng(seed)
    n = m.Ylog.shape[0]
    perms = [rng.permutation(n) for _ in range(n_perm)]

    M_real = increments(m, g)
    if n_jobs == 1:
        null = np.empty((n_perm, g.n_groups, g.n_groups))
        for b, pi in enumerate(perms):
            null[b] = increments(m, g, m.Ylog[pi])
            if progress and (b + 1) % 500 == 0:
                print(f"   {b + 1}/{n_perm}", flush=True)
    else:
        out = Parallel(n_jobs=n_jobs, verbose=5 if progress else 0)(
            delayed(increments)(m, g, m.Ylog[pi]) for pi in perms)
        null = np.stack(out)

    K = g.n_groups
    off = ~np.eye(K, dtype=bool)
    P = np.full((K, K), np.nan)
    for F in range(K):
        for T in range(K):
            if F != T:
                P[F, T] = (1 + np.sum(null[:, F, T] >= M_real[F, T])) / (1 + n_perm)
    Q = np.full((K, K), np.nan)
    Q[off] = bh_fdr(P[off])

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, M_real=M_real, P=P, Q=Q, null=null,
                        labels=np.array([str(k) for k in g.keys]),
                        n_perm=n_perm, seed=seed, grouping=g.id,
                        null_model="broken_pairing")
    return path


# ---------------------------------------------------------------------------
# Cross-grouping comparison — the control's actual claim
# ---------------------------------------------------------------------------
def matching_edges(cluster_t: Transfer, class_t: Transfer, m,
                   *, alpha: float = ALPHA_FDR) -> list:
    """Pair each significant cluster edge with its chemical-class counterpart.

    A cluster edge ``i -> j`` has a counterpart when some significant class
    edge ``a -> b`` has ``a`` as cluster ``i``'s dominant class and ``b`` as
    cluster ``j``'s. "Dominant class" is the class most of the cluster's
    glomeruli prefer — the same preferred-class assignment
    :func:`chem_class_grouping` uses, so the two analyses are related by one
    definition rather than by eye.

    Returns one entry per cluster edge, ``counterpart`` None when there is
    none. Seven of the eight cluster edges are expected to match.
    """
    dom = dominant_class(m, cluster_t.grouping, class_t.grouping)
    class_edges = {(e["source"], e["target"]): e for e in class_t.edges(alpha)}
    out = []
    for e in cluster_t.edges(alpha):
        key = (dom.get(e["source"]), dom.get(e["target"]))
        out.append(dict(cluster_edge=e, counterpart=class_edges.get(key),
                        via=key))
    return out


def dominant_class(m, cluster_g: Grouping, class_g: Grouping) -> dict:
    """``{cluster label: class label}`` — each cluster's most common preference.

    Only classes that survived `CHEM_MIN_SIZE` can be returned, since those are
    the only ones the class-level test has an edge for.
    """
    member = {}
    for pos, col in enumerate(class_g.columns):
        for gl in col:
            member[int(gl)] = class_g.labels[pos]
    out = {}
    for pos, col in enumerate(cluster_g.columns):
        votes = [member[int(gl)] for gl in col if int(gl) in member]
        if votes:
            names, counts = np.unique(votes, return_counts=True)
            out[cluster_g.labels[pos]] = str(names[counts.argmax()])
    return out
