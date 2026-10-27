"""Odorant x odorant similarity — the top row of Figure 3.

Four matrices over the same 59 odorants, in the same order, on the same
0-1 scale:

* **3A** response correlation, 'most sensitive' -- `most_sensitive_correlation(d)`
* **3B** response correlation at low (1x)      -- `response_correlation(d, "low")`
* **3C** response correlation at high (25x)    -- `response_correlation(d, "high")`
* **3D** MACCS chemical similarity             -- `chemical_similarity(d)`

They are here together on purpose. The figure's whole claim is that the
*response* structure changes with concentration while the *chemical* structure
does not, and that comparison only holds if all four matrices carry the same
row order and the same class blocks. Computing the ordering once, in
`class_sorted()`, is what makes that structural rather than a convention four
call sites have to remember.

    >>> from highconc.data import load
    >>> from highconc import similarity as sim
    >>> s = sim.response_correlation(load(), "low")
    >>> s.M.shape, s.blocks[:2]
    ((59, 59), ('C', 'D'))

Every matrix comes back as an `OdorantSimilarity`, already permuted into
plotting order, with the diagonal set to NaN so a panel draws it as "not a
measurement" rather than as a stripe of 1.0 that dominates the colour scale.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd

from . import clusters as _clusters

#: data/burton_et_al/ -- the "most-sensitive" condition matrix and the
#: odorant order it (and now every other within-block order) follows.
DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "burton_et_al"

#: Response transform applied before correlating. Raw is the default: the
#: matrices it produces are the ones in Matt's 2026-08-23 PDF (verified panel
#: by panel against page 3), and dF/F is already the physiological quantity.
#: "log" is offered because response amplitudes span two orders of magnitude
#: and a log correlation is the more conservative statistic; it changes the
#: picture very little, which is itself worth being able to show.
TRANSFORMS = ("raw", "log")

#: Within-block odorant orders ("table", "matt", "burton") and the aromatics
#: subclass taxonomy live in `highconc.clusters` now — both this module and
#: `recruitment.py` need them, and `clusters` is the taxonomy home both
#: already import.
WITHIN_BLOCK_ORDERS = _clusters.WITHIN_BLOCK_ORDERS

#: The three glomeruli drawn in Matt's Fig 3E, 3F and 3G — **recovered, not
#: chosen**. His panels were reduced to bar geometry straight out of the PDF's
#: vector content and scored against all 376 glomeruli; these three reproduce
#: every one of his 60 bars under a single per-panel scale factor, with a
#: maximum residual of 0.008 pt. The scale that fit panel E independently
#: predicts its own axis tick spacing (11.98 pt for a 20-unit tick against
#: 11.980 measured), so the amplitudes are right and not just the shape.
#:
#: Fig 6B/C/D must draw the SAME three. That is the entire reason this is one
#: named constant instead of three literals in two notebook sections.
EG_GLOMERULI: tuple[int, int, int] = (173, 211, 97)

#: Which cluster each of them belongs to, recorded so a mismatch after any
#: re-clustering is loud rather than silent. E -> 5, F -> 2, G -> 6.
EG_GLOMERULI_CLUSTERS: tuple[int, int, int] = (5, 2, 6)


@dataclass(frozen=True)
class Spectrum:
    """One glomerulus's response across the 59 odorants — Fig 3E/F/G, 6B/C/D.

    `low` and `high` are in plotting order (same `class_sorted` layout as the
    matrices), so a spectrum and a matrix always agree on which column is which
    odorant.
    """

    glom: int
    low: np.ndarray
    high: np.ndarray
    classes: np.ndarray
    names: list[str]
    blocks: tuple[str, ...]
    boundaries: tuple[int, ...]
    cluster: int

    @property
    def limit(self) -> float:
        """Symmetric-friendly max, for setting a shared y-range across panels."""
        return float(max(self.low.max(), self.high.max()))

    @property
    def responded(self) -> np.ndarray:
        """Boolean mask of odorants with a nonzero response at either
        concentration — the support that identified these three."""
        return (self.low > 0) | (self.high > 0)


def spectrum(d, glom: int, *, order="matt", within="table") -> Spectrum:
    """One glomerulus's low and high response spectra, in plotting order."""
    perm, classes, blocks, boundaries = class_sorted(d, order, within)
    g = int(glom)
    return Spectrum(glom=g,
                    low=d.low.values[perm, g].astype(float),
                    high=d.high.values[perm, g].astype(float),
                    classes=classes,
                    names=[str(x) for x in d.odorants["Name"].values[perm]],
                    blocks=blocks, boundaries=boundaries,
                    cluster=int(d.cluster_number[g]))


def eg_spectra(d, gloms=EG_GLOMERULI, **kw) -> list[Spectrum]:
    """The three example spectra, in Matt's E, F, G order. Pass `gloms` to
    override — the whole point of `EG_GLOMERULI` being a parameter."""
    return [spectrum(d, g, **kw) for g in gloms]


@dataclass(frozen=True)
class OdorantSimilarity:
    """A 59 x 59 odorant similarity matrix, already in plotting order.

    Attributes
    ----------
    M : ndarray, (59, 59)
        The matrix, permuted into class-block order, diagonal set to NaN.
    order : ndarray, (59,)
        Indices into the dataset's original odorant order. `M` is
        ``raw[np.ix_(order, order)]``.
    classes : ndarray of str, (59,)
        Chemical-class letter per row, in plotting order. Grouped by
        construction.
    names : list of str, (59,)
        Odorant display names, in plotting order.
    blocks : tuple of str
        The class letters in block order, e.g. ``('C', 'D', 'E', ...)``.
    boundaries : tuple of int
        Row indices where a new class block starts, excluding 0. Feed straight
        to `style.heatmap(row_boundaries=..., col_boundaries=...)`.
    what : str
        Human-readable description, used as a panel title and for provenance.
    """

    M: np.ndarray
    order: np.ndarray
    classes: np.ndarray
    names: list[str]
    blocks: tuple[str, ...]
    boundaries: tuple[int, ...]
    what: str

    @property
    def centers(self) -> tuple[float, ...]:
        """Mid-row of each class block, for tick placement."""
        edges = (0,) + self.boundaries + (len(self.classes),)
        return tuple((a + b - 1) / 2 for a, b in zip(edges, edges[1:]))

    @property
    def offdiagonal(self) -> np.ndarray:
        """The off-diagonal values as a flat array, NaNs dropped."""
        v = self.M[~np.eye(len(self.M), dtype=bool)]
        return v[~np.isnan(v)]

    def block_means(self) -> pd.DataFrame:
        """Mean similarity within and between each pair of class blocks.

        The numeric counterpart of looking at the picture: it is what tells you
        whether the bright square you think you see on the diagonal is really
        brighter than its row.
        """
        n = len(self.blocks)
        out = np.full((n, n), np.nan)
        for i, a in enumerate(self.blocks):
            ra = self.classes == a
            for j, b in enumerate(self.blocks):
                rb = self.classes == b
                sub = self.M[np.ix_(ra, rb)]
                out[i, j] = np.nanmean(sub) if np.isfinite(sub).any() else np.nan
        names = [_clusters.class_name(c) for c in self.blocks]
        return pd.DataFrame(out, index=names, columns=names)


@dataclass(frozen=True)
class OdorantGlomMatrix:
    """Figure 2: an odorant x glomerulus response matrix, already in
    plotting order on the odorant (row) axis.

    Unlike `OdorantSimilarity` this is not square and has no diagonal --
    it's a response matrix, not a similarity matrix -- so only the row axis
    carries chemical-class structure. The column axis is whatever glomeruli
    the source file ships (Burton's own 113, not this dataset's 376) with no
    ordering or identity information, so it is left exactly as given.

    Attributes
    ----------
    M : ndarray, (59, n_gloms)
        The matrix, rows permuted into class-block order. Columns untouched.
    classes : ndarray of str, (59,)
        Chemical-class letter per row, in plotting order.
    names : list of str, (59,)
        Odorant display names, in plotting order.
    blocks : tuple of str
        The class letters in block order.
    boundaries : tuple of int
        Row indices where a new class block starts, excluding 0. Feed to
        `style.heatmap(row_boundaries=...)`.
    what : str
        Human-readable description, used as a panel title.
    """

    M: np.ndarray
    classes: np.ndarray
    names: list[str]
    blocks: tuple[str, ...]
    boundaries: tuple[int, ...]
    what: str


def class_sorted(d, order="matt", within="table"
                 ) -> tuple[np.ndarray, np.ndarray, tuple, tuple]:
    """``(perm, classes, blocks, boundaries)`` for a dataset's odorants.

    `order` is ``"table"``, ``"matt"``, ``"matt_esters_first"`` or an explicit
    sequence of class letters — see `highconc.clusters.class_order`. It sets
    the order of the BLOCKS.

    `within` sets the order *inside* each block: ``"table"`` (default, the
    odorant table's own order) or ``"matt"``, which additionally applies
    `WITHIN_BLOCK_MATT` — a ten-odorant reshuffle of the pyrazine block, the
    only place his figure and our table disagree. It is not the default
    because it changes nothing structural and every existing export uses the
    table order; it is here so a panel can be made cell-for-cell identical to
    his when that is what someone wants.

    Default block order is ``"matt"``: that is the order actually drawn in the
    2026-08-23 figure PDF, confirmed against page 3 rather than assumed.
    """
    grp = np.asarray(d.odorants["Group"].values)
    blocks = tuple(_clusters.class_order(order, present=grp))
    if isinstance(within, str):
        try:
            within = WITHIN_BLOCK_ORDERS[within.lower()]
        except KeyError:
            raise ValueError(f"unknown within-block order {within!r}; expected "
                             f"one of {sorted(WITHIN_BLOCK_ORDERS)} or a "
                             "{class: codes} mapping") from None
    codes = np.asarray([str(x).split()[0] for x in d.odorants["Odorant"].values])
    chunks = []
    for g in blocks:
        rows = np.where(grp == g)[0]
        want = within.get(g)
        if want is not None:
            pos = {c: i for i, c in enumerate(codes[rows])}
            if sorted(pos) != sorted(want):
                raise ValueError(
                    f"within-block order for class {g!r} lists {sorted(want)} "
                    f"but the data has {sorted(pos)}")
            rows = rows[[pos[c] for c in want]]
        chunks.append(rows)
    perm = np.concatenate(chunks)
    if len(perm) != len(grp):
        missing = sorted(set(range(len(grp))) - set(perm.tolist()))
        raise ValueError(f"{len(missing)} odorants dropped by class order "
                         f"{order!r}: rows {missing}")
    classes = grp[perm]
    sizes = np.cumsum([int((grp == g).sum()) for g in blocks])
    return perm, classes, blocks, tuple(int(b) for b in sizes[:-1])


def _package(M, d, perm, classes, blocks, boundaries, what) -> OdorantSimilarity:
    M = np.asarray(M, float)[np.ix_(perm, perm)].copy()
    np.fill_diagonal(M, np.nan)     # self-similarity is 1 by construction
    names = [str(x) for x in d.odorants["Name"].values[perm]]
    return OdorantSimilarity(M=M, order=perm, classes=classes, names=names,
                             blocks=blocks, boundaries=boundaries, what=what)


def response_correlation(d, conc="low", *, order="matt", within="table",
                         transform="raw") -> OdorantSimilarity:
    """Figure 3B/3C: odorant x odorant Pearson correlation of response maps.

    Each odorant is a 376-glomerulus vector; the matrix is the correlation
    between every pair of those vectors at one concentration. `conc` is
    ``"low"`` (1x) or ``"high"`` (25x).

    Responses are non-negative, so negative correlations are rare and small —
    which is why the figure's 0-1 scale is honest rather than a crop. Anything
    below zero is left in the matrix, not clipped; a panel that pins vmin=0
    renders it at the floor, and `offdiagonal.min()` still reports it.
    """
    if conc not in ("low", "high"):
        raise ValueError(f"conc must be 'low' or 'high', got {conc!r}")
    if transform not in TRANSFORMS:
        raise ValueError(f"transform must be one of {TRANSFORMS}, "
                         f"got {transform!r}")
    X = getattr(d, conc).values.astype(float)
    if transform == "log":
        X = np.log10(X + 1.0)
    perm, classes, blocks, boundaries = class_sorted(d, order, within)
    label = {"low": "'low' (1x)", "high": "'high' (25x)"}[conc]
    what = f"response correlation, {label}"
    if transform != "raw":
        what += f" [{transform}]"
    return _package(np.corrcoef(X), d, perm, classes, blocks, boundaries, what)


def _load_odorant_order() -> list[str]:
    """The 59 `Odorant` strings in data/burton_et_al/odorant_order.txt, in
    file order. Each line is a quoted, comma-terminated literal (not valid
    Python on its own without wrapping in brackets), so this just strips the
    quoting rather than eval-ing it."""
    lines = (DATA_DIR / "odorant_order.txt").read_text().splitlines()
    return [ln.strip().rstrip(",").strip("'\"") for ln in lines if ln.strip()]


def most_sensitive_correlation(d, *, order="matt", within="table"
                               ) -> OdorantSimilarity:
    """Figure 3A: odorant x odorant response correlation at the
    'most-sensitive' condition -- a third concentration this dataset does not
    otherwise have.

    Loaded from `data/burton_et_al/most_sensitive.csv`, which ships as a
    61x61 matrix in `odorant_order.txt`'s row order rather than this
    dataset's native order: the first 59 rows/columns match the 59 named
    odorants exactly (its four all-NaN rows land precisely on D15, E1, P1 and
    P4's positions); rows/columns 60-61 are two unnamed extras and are
    dropped. The remaining NaNs -- those same four odorants, which have no
    data at this condition -- are filled with 0.0 rather than left as gaps:
    the panel pins vmin=0, so 0.0 renders at the colour floor rather than as
    a hole. `_package` re-masks the true diagonal to NaN afterward regardless
    of what was filled there, so self-similarity still reads as "not a
    measurement" exactly like 3B/3C/3D.

    Reindexed onto the dataset's own odorant order before permuting into
    `class_sorted`'s plotting order, so it composes with `split_block`
    identically to `response_correlation`/`chemical_similarity`.
    """
    names = _load_odorant_order()
    n = len(names)
    raw = pd.read_csv(DATA_DIR / "most_sensitive.csv", header=None).values
    if raw.shape[0] < n or raw.shape[1] < n:
        raise ValueError(f"most_sensitive.csv is {raw.shape}, too small for "
                         f"odorant_order.txt's {n} odorants")
    raw = raw[:n, :n]
    raw = np.nan_to_num(raw.astype(float), nan=0.0)

    have = set(d.odorants["Odorant"])
    missing = sorted(set(names) - have)
    if missing:
        raise ValueError("odorant_order.txt names not found in the dataset: "
                         f"{missing}")

    wide = pd.DataFrame(raw, index=names, columns=names)
    wide = wide.reindex(index=d.odorants["Odorant"].values,
                        columns=d.odorants["Odorant"].values)
    perm, classes, blocks, boundaries = class_sorted(d, order, within)
    return _package(wide.values, d, perm, classes, blocks, boundaries,
                    "response correlation, 'most sensitive'")


def burton_response_matrix(d, conc="high", *, order="matt", within="burton"
                           ) -> OdorantGlomMatrix:
    """Figure 2: odorant x glomerulus response matrix from Burton's own data,
    at ``"high"`` (`data/burton_et_al/sortedmatrix.csv`, 113 glomeruli),
    ``"low"`` (`sortedmatrixlo.csv`, 113 glomeruli) or ``"most_sensitive"``
    (`sortedmatrix_omp114L_omp8x.csv`, 76 glomeruli -- a different
    subject/preparation, not a 113-column slice, hence the different width).
    None of these are this dataset's own 376 glomeruli.

    None of the files have a header or odorant labels, so rows are matched by
    position against `data/burton_et_al/odorant_order.txt`, the same
    convention `most_sensitive_correlation` uses. Reindexed onto this
    dataset's own odorant order before permuting into `class_sorted`'s
    plotting order, so the row axis is directly comparable to Figure 3's.

    Columns are Burton's own glomeruli, in whatever order each file ships --
    there is no identity or ordering information for them, so they are left
    untouched and unlabeled.
    """
    fnames = {"high": "sortedmatrix.csv", "low": "sortedmatrixlo.csv",
             "most_sensitive": "sortedmatrix_omp114L_omp8x.csv"}
    if conc not in fnames:
        raise ValueError(f"conc must be one of {sorted(fnames)}, got {conc!r}")
    fname = fnames[conc]

    names = _load_odorant_order()
    n = len(names)
    raw = pd.read_csv(DATA_DIR / fname, header=None).values.astype(float)
    if raw.shape[0] != n:
        raise ValueError(f"{fname} has {raw.shape[0]} rows, expected {n} to "
                         "match odorant_order.txt")

    have = set(d.odorants["Odorant"])
    missing = sorted(set(names) - have)
    if missing:
        raise ValueError("odorant_order.txt names not found in the dataset: "
                         f"{missing}")

    wide = pd.DataFrame(raw, index=names)
    wide = wide.reindex(index=d.odorants["Odorant"].values)
    perm, classes, blocks, boundaries = class_sorted(d, order, within)
    row_names = [str(x) for x in d.odorants["Name"].values[perm]]
    label = {"low": "'low'", "high": "'high'",
             "most_sensitive": "'most sensitive'"}[conc]
    return OdorantGlomMatrix(M=wide.values[perm, :], classes=classes,
                             names=row_names, blocks=blocks,
                             boundaries=boundaries,
                             what=f"odorant x glomerulus response, {label}")


def split_block(sim: OdorantSimilarity, d, letter: str,
                subclass_by_code: dict[str, str]) -> OdorantSimilarity:
    """Re-label one class block's rows into finer display subclasses.

    Built for Figure 3's aromatics breakdown: `letter="R"` with
    `highconc.clusters.AROMATIC_SUBCLASS_BY_CODE` splits the ten aromatic
    odorants into Ar. Ketone / Ar. Ketone/Ester / Phenol/Methoxy for the axis
    blocking and colours, without touching `Dataset.chem_class` or any other
    panel -- `sim.classes` is relabelled and `blocks`/`boundaries` recomputed
    from the result; everything else (`M`, `order`, `names`) is untouched.

    The split only reads cleanly when `subclass_by_code`'s values are
    contiguous in `sim`'s row order, which they are here because `sim` was
    built with `within="burton"` -- that order already lays the block out in
    the same 1/7/2 sequence the subclasses follow.
    """
    codes = np.asarray([str(x).split()[0]
                        for x in d.odorants["Odorant"].values])[sim.order]
    classes = sim.classes.astype(object).copy()
    mask = classes == letter
    missing = sorted(set(codes[mask]) - set(subclass_by_code))
    if missing:
        raise ValueError(f"no subclass given for {letter!r} odorant(s) "
                         f"{missing}")
    classes[mask] = [subclass_by_code[c] for c in codes[mask]]

    blocks, boundaries = [], []
    for i, c in enumerate(classes):
        if not blocks or blocks[-1] != c:
            blocks.append(c)
            if i:
                boundaries.append(i)
    return replace(sim, classes=classes, blocks=tuple(blocks),
                   boundaries=tuple(boundaries))


def chemical_similarity(d, *, order="matt", within="table", **kw
                        ) -> OdorantSimilarity:
    """Figure 3D: MACCS-key cosine similarity between the 59 odorants.

    Delegates the chemistry to `highconc.chem.maccs_cosine`, which reads a
    cached matrix and only needs RDKit when that cache has to be rebuilt. The
    ordering comes from the same `class_sorted()` call as 3B and 3C, so the
    three panels cannot drift apart.
    """
    from . import chem as _chem

    M = _chem.maccs_cosine(d, **kw)
    perm, classes, blocks, boundaries = class_sorted(d, order, within)
    return _package(M, d, perm, classes, blocks, boundaries,
                    "MACCS chemical similarity (cosine)")


def concentration_shift(low: OdorantSimilarity, high: OdorantSimilarity
                        ) -> dict:
    """How much did the correlational structure move between 1x and 25x?

    Returns the two mean off-diagonal correlations, their difference, and the
    Spearman correlation between the two matrices' upper triangles — the
    figure's claim in three numbers, so a caption never has to be eyeballed
    off the picture.
    """
    from scipy.stats import spearmanr

    if low.order.tolist() != high.order.tolist():
        raise ValueError("the two matrices are in different odorant orders; "
                         "build both with the same `order` argument")
    iu = np.triu_indices(len(low.M), k=1)
    a, b = low.M[iu], high.M[iu]
    ok = np.isfinite(a) & np.isfinite(b)
    rho, p = spearmanr(a[ok], b[ok])
    return dict(mean_low=float(np.nanmean(a)), mean_high=float(np.nanmean(b)),
                delta=float(np.nanmean(b) - np.nanmean(a)),
                spearman=float(rho), p=float(p), n_pairs=int(ok.sum()))


# ---------------------------------------------------------------------------
# Figure 4A — the same idea one axis over: glomerulus x glomerulus
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GlomerulusSimilarity:
    """A 376 x 376 glomerulus similarity matrix, already in cluster order.

    The glomerular counterpart of `OdorantSimilarity`. It blocks by NMF
    cluster rather than by chemical class, which is why it is a separate type
    rather than the same one with different labels: the blocks here are the
    figure's own clustering, so a reader is entitled to ask whether the bright
    diagonal squares are the clusters, and `block_means()` answers that.

    Attributes
    ----------
    M : ndarray, (376, 376)
        Permuted into cluster order, diagonal set to NaN.
    order : ndarray, (376,)
        Indices into the dataset's glomerulus order; `M` is
        ``raw[np.ix_(order, order)]``.
    numbers : ndarray of int, (376,)
        Joel 1-7 cluster number per row, in plotting order.
    blocks : tuple of int
        The cluster numbers in block order, normally ``(1, ..., 7)``.
    boundaries : tuple of int
        Row indices where a new cluster block starts, excluding 0.
    what : str
        Human-readable description, used as a panel title.
    """

    M: np.ndarray
    order: np.ndarray
    numbers: np.ndarray
    blocks: tuple[int, ...]
    boundaries: tuple[int, ...]
    what: str

    @property
    def centers(self) -> tuple[float, ...]:
        """Mid-row of each cluster block, for tick placement."""
        edges = (0,) + self.boundaries + (len(self.numbers),)
        return tuple((a + b - 1) / 2 for a, b in zip(edges, edges[1:]))

    @property
    def offdiagonal(self) -> np.ndarray:
        """The off-diagonal values as a flat array, NaNs dropped."""
        v = self.M[~np.eye(len(self.M), dtype=bool)]
        return v[~np.isnan(v)]

    def block_means(self) -> pd.DataFrame:
        """Mean correlation within and between each pair of clusters, 7 x 7."""
        n = len(self.blocks)
        out = np.full((n, n), np.nan)
        for i, a in enumerate(self.blocks):
            ra = self.numbers == a
            for j, b in enumerate(self.blocks):
                rb = self.numbers == b
                sub = self.M[np.ix_(ra, rb)]
                out[i, j] = np.nanmean(sub) if np.isfinite(sub).any() else np.nan
        labels = [str(b) for b in self.blocks]
        return pd.DataFrame(out, index=labels, columns=labels)


def glomerulus_correlation(d, conc="low", *, order="joel", transform="raw"
                           ) -> GlomerulusSimilarity:
    """Figure 4A: glomerulus x glomerulus Pearson correlation of tuning curves.

    Each glomerulus is a 59-odorant vector; the matrix is the correlation
    between every pair of those vectors at one concentration, permuted so the
    NMF clusters run down the diagonal. Matt's 4A is the ``"low"`` one.

    Parameters
    ----------
    conc : {"low", "high"}
    order : {"joel", "component"}
        Which cluster ordering the blocks run in — the same argument, with the
        same meaning, as `panels.clustering.bicluster_order`. ``"joel"`` is the
        default so 4A's blocks read 1-7 left to right like every other panel of
        the figure.
    transform : {"raw", "log"}
        As in `response_correlation`.
    """
    if conc not in ("low", "high"):
        raise ValueError(f"conc must be 'low' or 'high', got {conc!r}")
    if transform not in TRANSFORMS:
        raise ValueError(f"transform must be one of {TRANSFORMS}, "
                         f"got {transform!r}")
    if order == "joel":
        key = np.asarray(d.cluster_number, dtype=int)
    elif order == "component":
        key = np.asarray(d.cluster_index, dtype=int)
    else:
        raise ValueError(f"order must be 'joel' or 'component', got {order!r}")

    X = getattr(d, conc).values.astype(float)          # 59 odorants x 376 gloms
    if transform == "log":
        X = np.log10(X + 1.0)

    perm = np.argsort(key, kind="stable")
    sorted_key = key[perm]
    boundaries = tuple(int(b) for b in np.where(np.diff(sorted_key) != 0)[0] + 1)
    blocks = tuple(int(k) for k in dict.fromkeys(sorted_key.tolist()))

    M = np.corrcoef(X, rowvar=False)                   # glom x glom
    M = M[np.ix_(perm, perm)].copy()
    np.fill_diagonal(M, np.nan)                        # 1 by construction

    label = {"low": "'low' (1x)", "high": "'high' (25x)"}[conc]
    what = f"glomerulus correlation, {label}"
    if transform != "raw":
        what += f" [{transform}]"
    if order != "joel":
        what += f" [{order} order]"
    return GlomerulusSimilarity(M=M, order=perm, numbers=sorted_key,
                                blocks=blocks, boundaries=boundaries, what=what)
