"""Clustering panels — the bi-clustered response heatmaps and their companions.

Extracted from ``figures/nmf_clusters_main/generate.py``. Each function draws
into an ``ax`` you pass in and returns the artist a colorbar would need, so the
notebook owns the layout and the colorbars.

    from highconc import load, style as hs
    from highconc.panels import clustering

    d = load()
    fig, (ax_lo, ax_hi) = plt.subplots(2, 1)
    clustering.biclustered(d, ax_lo, conc="low", title="low")
    mesh = clustering.biclustered(d, ax_hi, conc="high", title="high")
    fig.colorbar(mesh, ax=[ax_lo, ax_hi], label="dF/F")
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patheffects as patheffects
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

from .. import clusters as _clusters
from ..data import Dataset

__all__ = [
    "BiclusterOrder", "bicluster_order", "biclustered",
    "glomerulus_correlation_matrix",
    "class_cluster_glom_composition", "cluster_summary",
    "CHEM_CLASS_NAMES",
    "topography", "cluster_contrast", "NONRESPONDER_MARKERS",
    "bulb_limits", "per_subject_heatmap", "block_means", "per_subject_blocks",
    "per_subject_topography",
]

# Defaults carried over from the pre-refactor script. Every one is exposed as a
# keyword argument; these are the values the manuscript figure was drawn with.
VMIN, VMAX = 1.0, 140.0
CMAP_NAME = "magma"
EMPTY_COLOR = "gainsboro"       # cells below vmin read as "no response"
BLOCK_COLOR = "#c4c4c4"         # diagonal block highlight, lighter than empty
BLOCK_EDGE_COLOR = "#888888"
BLOCK_EDGE_LW = 1.0
BOUNDARY_COLOR = "blue"
BOUNDARY_LW = 1.5
BOUNDARY_ALPHA = 0.5


@dataclass(frozen=True)
class BiclusterOrder:
    """The row/column sort shared by both concentrations.

    Attributes
    ----------
    rows, cols : ndarray
        Index arrays that sort odorants by their bi-cluster and glomeruli by
        their NMF component. Applying them to `low` and `high` alike is what
        makes a glomerulus mean the same thing in both heatmaps.
    col_boundaries, row_boundaries : ndarray
        Positions where the sorted cluster label changes — where to draw lines.
    col_edges, row_edges : ndarray
        The same boundaries with 0 and N appended, so ``edges[k]:edges[k+1]``
        is the span of cluster k.
    row_class_boundaries : ndarray
        Positions where the sorted chemical-class letter changes, within the
        odorant row order actually produced (finer-grained than
        `row_boundaries`, and independent of it -- a class can change without
        a cluster changing). Meaningful as a clean, blocky set of lines only
        when `bicluster_order`'s `chem_order` was set; with `chem_order=None`
        this still reflects wherever native order happens to cross a class,
        which will look scattered rather than blocked.
    """

    rows: np.ndarray
    cols: np.ndarray
    col_boundaries: np.ndarray
    row_boundaries: np.ndarray
    col_edges: np.ndarray
    row_edges: np.ndarray
    row_class_boundaries: np.ndarray

    @property
    def n_blocks(self) -> int:
        """Diagonal blocks: one per glom cluster that has an odorant cluster."""
        return min(len(self.col_edges), len(self.row_edges)) - 1


def bicluster_order(
    dataset: Dataset, *, order: str = "joel", chem_order: str | None = None,
    within_order: str = "burton", row_order: str = "cluster"
) -> BiclusterOrder:
    """Sort glomeruli by their cluster and odorants by their bi-cluster.

    Glom cluster k pairs with odorant cluster k, so applying the *same* cluster
    ordering to rows and columns is what places each pair's highlighted block
    on the diagonal. Any ordering does that, as long as both axes use it.

    Parameters
    ----------
    order : {"joel", "component"}
        Which cluster ordering the blocks run in. ``"joel"`` is the canonical
        1-7 numbering, and is the default so that this panel's blocks line up
        with the 1-7 columns of `cluster_summary` beside it — left-to-right
        block k really is cluster k. ``"component"`` is the raw NMF component
        order 0-6 the pre-refactor script used; it is here so an old render can
        be reproduced, not because anything should be laid out that way now.
    chem_order : str, optional
        Secondary sort key for odorant rows, within each cluster block --
        glomeruli have no chemical class, so this only ever touches `rows`,
        never `cols`. `None` (default) leaves odorants in their native order
        within a block, same as before this was added. Pass ``"matt"`` (or
        any `highconc.clusters.class_order` name) to additionally sort
        odorants by chemical class inside each block, at the fine-grained
        per-odorant resolution `within_order` sets -- not just grouped by
        letter -- so the result matches Figure 3/6's own odorant order,
        block for block. The block boundaries (and therefore the diagonal
        alignment with the glomerulus columns) are unchanged either way;
        only the order *within* a block moves. Required (not `None`) when
        `row_order="odorant"` -- there is no other order to fall back to.
    within_order : str
        Fine order *inside* each chemical class, only used when `chem_order`
        is set -- see `highconc.similarity.class_sorted`'s `within`.
        ``"burton"`` (default) is the full 59-odorant canonical order in
        `data/burton_et_al/odorant_order.txt`, the same one §3/§6 use via
        `WITHIN_CLASS_ORDER`.
    row_order : {"cluster", "odorant"}
        ``"cluster"`` (default): odorant rows are primarily sorted by their
        bi-cluster (`chem_order`, if set, only breaks ties *within* a
        cluster's block) -- this is what makes the diagonal block structure
        show up, and is what every figure before this option existed drew.
        ``"odorant"``: odorant rows ignore cluster membership entirely and
        follow `chem_order`/`within_order` directly (Fig 3/6's own odorant
        order, if that's what's passed) -- a cluster's rows are then
        scattered rather than contiguous, so the diagonal block highlight in
        `biclustered` (`diagonal_blocks=True`) no longer means anything and
        should be turned off by the caller. Glomerulus columns are
        untouched either way -- they always sort by cluster, so the
        vertical `col_boundaries` lines stay meaningful in both modes.
    """
    if order == "joel":
        col_key = dataset.cluster_number
        # An odorant cluster carries a component index, so map it the same way.
        row_key = np.array(
            [_clusters.to_joel(k) for k in dataset.odorant_cluster], dtype=int
        )
    elif order == "component":
        col_key = dataset.cluster_index
        row_key = dataset.odorant_cluster
    else:
        raise ValueError(f"order must be 'joel' or 'component', got {order!r}")

    cols = np.argsort(col_key, kind="stable")
    if row_order == "odorant":
        if chem_order is None:
            raise ValueError(
                "row_order='odorant' needs chem_order set -- there is no "
                "cluster fallback to sort rows by in this mode")
        from .. import similarity as _similarity

        rows, *_ = _similarity.class_sorted(dataset, order=chem_order,
                                            within=within_order)
    elif row_order != "cluster":
        raise ValueError(f"row_order must be 'cluster' or 'odorant', got "
                         f"{row_order!r}")
    elif chem_order is None:
        rows = np.argsort(row_key, kind="stable")
    else:
        from .. import similarity as _similarity

        # The exact canonical odorant permutation Fig 3/6 draw in, inverted
        # into a per-odorant rank -- so the secondary sort matches them block
        # for block, not just letter for letter.
        perm, *_ = _similarity.class_sorted(dataset, order=chem_order,
                                            within=within_order)
        chem_key = np.empty(len(perm), dtype=int)
        chem_key[perm] = np.arange(len(perm))
        # lexsort's last key is primary: cluster first, chemical class breaks
        # ties within a cluster's block.
        rows = np.lexsort((chem_key, row_key))

    sorted_glom = col_key[cols]
    sorted_odor = row_key[rows]
    sorted_class = np.asarray(dataset.chem_class)[rows]

    col_boundaries = np.where(np.diff(sorted_glom) != 0)[0] + 1
    row_boundaries = np.where(np.diff(sorted_odor) != 0)[0] + 1
    row_class_boundaries = np.where(sorted_class[:-1] != sorted_class[1:])[0] + 1

    return BiclusterOrder(
        rows=rows,
        cols=cols,
        col_boundaries=col_boundaries,
        row_boundaries=row_boundaries,
        col_edges=np.concatenate([[0], col_boundaries, [len(sorted_glom)]]),
        row_edges=np.concatenate([[0], row_boundaries, [len(sorted_odor)]]),
        row_class_boundaries=row_class_boundaries,
    )


def biclustered(
    dataset: Dataset,
    ax,
    *,
    conc: str = "low",
    order: BiclusterOrder | None = None,
    col_normalize: bool = False,
    raw_floor: float | None = None,
    vmin: float = VMIN,
    vmax: float = VMAX,
    cmap: str = CMAP_NAME,
    empty_color: str = EMPTY_COLOR,
    block_color: str = BLOCK_COLOR,
    block_edge_color: str = BLOCK_EDGE_COLOR,
    block_edge_lw: float = BLOCK_EDGE_LW,
    diagonal_blocks: bool = True,
    boundaries: bool = True,
    boundary_color: str = BOUNDARY_COLOR,
    boundary_lw: float = BOUNDARY_LW,
    boundary_alpha: float = BOUNDARY_ALPHA,
    band_colors: tuple[str, str] | None = None,
    class_names: bool = False,
    class_name_fontsize: float = 9,
    title: str | None = None,
    xlabel: str | None = "sorted glomeruli",
    ylabel: str | None = "sorted odorants",
    title_fontsize: float = 15,
    label_fontsize: float = 12,
):
    """Draw one concentration's bi-clustered response heatmap into `ax`.

    Odorants (rows) and glomeruli (columns) are sorted into their clusters, the
    same way for both concentrations, and each glom cluster's paired odorant
    block is highlighted a shade lighter than the background so the diagonal
    structure is visible where responses are sparse. Responses are on a log
    colour scale (matching `per_subject_heatmap`'s); anything under `vmin`
    falls through to the background.

    Parameters
    ----------
    dataset : Dataset
    ax : matplotlib Axes
    conc : {"low", "high"}
        Which concentration to draw. Both use the same sort.
    order : BiclusterOrder, optional
        Precomputed sort. Pass the same object to both concentrations to avoid
        recomputing it; omitted, it is derived from `dataset`.
    col_normalize : bool
        Peak-normalise each column (glomerulus) to [0, 1] by its own maximum
        response across all odorants, on the raw (linear) dF/F values --
        before any log transform, and instead of one. A column that never
        fires stays all-zero rather than dividing by zero. Off by default
        (plots raw dF/F on a log scale, as below). On, this draws on a
        LINEAR scale instead of log -- the two are meant to be set together,
        since a [0, 1] normalised quantity has no absolute floor to put a log
        scale's dynamic range to use. Pass `vmin=0.0, vmax=1.0` alongside it;
        the defaults below (`VMIN`/`VMAX`, tuned for raw dF/F) are wrong for
        this mode and are not adjusted for you.
    raw_floor : float, optional
        Only meaningful with `col_normalize=True`. Normalising to a
        column's own peak can make a cell that was noise-floor-weak in
        absolute dF/F terms (a fraction of 1) look just as bright as a
        genuine response, if its glomerulus's peak also happens to be small
        -- normalising erases the very floor `vmin` enforces when
        `col_normalize=False`. `raw_floor` restores it: any cell whose
        *pre-normalisation* raw value is below `raw_floor` is drawn as
        `empty_color`, regardless of what its normalised value would have
        been. `None` (default) applies no such floor -- every normalised
        value draws on the ramp, including a glomerulus's weakest response
        to its own best odorant, however small that response was in
        absolute terms.
    vmin, vmax : float
        Colour-scale limits, in dF/F (log scale, `col_normalize=False`) or in
        the normalised [0, 1] units (linear scale, `col_normalize=True`).
        Cells below `vmin` are drawn as `empty_color` rather than as the
        colormap's darkest end.
    diagonal_blocks : bool
        Draw the lighter diagonal-block highlight (and its edge lines). On by
        default -- this is what makes sense when `order`'s rows are
        cluster-ordered (`row_order="cluster"`). Turn off when `order` was
        built with `row_order="odorant"`: a cluster's rows are no longer
        contiguous there, so the highlight would pick out the wrong cells
        entirely rather than just looking redundant.
    boundaries : bool
        Draw the glomerulus-cluster (vertical) and chemical-class (horizontal)
        boundary lines. On by default; turn off for a plain heatmap with no
        rules. The horizontal *cluster* lines stay dropped deliberately (see
        below) even when this is on -- only the finer class lines draw.
    band_colors : pair of colors, optional
        Alternating per-cluster background bands instead of the uniform
        background. Off by default; it only reads well when the empty cells are
        transparent.
    class_names : bool
        Mark each contiguous chemical-class row-run with a coloured block and
        its written-out name ("aldehydes"), via `highconc.style.class_strip`
        (`side="left"`, same mechanism Figure 2's `response_matrix` uses on
        the right) -- colours from `highconc.clusters.class_color`, names
        from `cluster_summary` (panel D)'s own `CHEM_CLASS_NAMES`. Off by
        default. Reads `order.row_class_boundaries` to find the runs, so it
        marks however many contiguous class blocks the current row order
        actually produces -- exactly 10 when `order` was built with
        `row_order="odorant"` (one clean block per class); with
        `row_order="cluster"` a class split across multiple cluster blocks
        gets one block per run, so its colour and name can repeat.
    class_name_fontsize : float
        Font size for those name labels.
    title, xlabel, ylabel : str, optional
        Pass None to leave off — the usual case for the upper of two stacked
        heatmaps, which shares the lower one's x axis.

    Returns
    -------
    QuadMesh
        The response mesh, for handing to `fig.colorbar`.
    """
    if conc not in ("low", "high"):
        raise ValueError(f"conc must be 'low' or 'high', got {conc!r}")
    if order is None:
        order = bicluster_order(dataset)

    data = getattr(dataset, conc).values
    data = data[np.ix_(order.rows, order.cols)]

    if col_normalize:
        peak = data.max(0, keepdims=True)
        peak = np.where(peak > 0, peak, np.inf)
        normed = data / peak
        if raw_floor is not None:
            # NaN, not 0 -- 0 would draw as the ramp's own floor colour
            # (a measurement), where this cell was actually below the noise
            # floor before normalising (mesh_cmap.set_bad below sends NaN to
            # the background instead).
            normed = np.where(data >= raw_floor, normed, np.nan)
        data = normed

    if band_colors is not None:
        for i in range(len(order.col_edges) - 1):
            ax.axvspan(
                order.col_edges[i], order.col_edges[i + 1],
                facecolor=band_colors[i % 2], zorder=0,
            )

    # Transparent empty/masked cells so the background and the block
    # highlights, both drawn as explicit patches behind the data, show
    # through -- a drawn rectangle rather than the axes facecolor, so it
    # renders the same way regardless of backend.
    mesh_cmap = plt.get_cmap(cmap).copy()
    mesh_cmap.set_under("none")
    mesh_cmap.set_bad("none")           # LogNorm masks zeros; `raw_floor` uses this too

    n_rows, n_cols = data.shape
    ax.add_patch(Rectangle((0, 0), n_cols, n_rows,
                           facecolor=empty_color, edgecolor="none", zorder=0))

    norm = (plt.Normalize(vmin=vmin, vmax=vmax) if col_normalize
           else LogNorm(vmin=vmin, vmax=vmax))
    sns.heatmap(
        data, cmap=mesh_cmap, norm=norm, ax=ax,
        xticklabels=[], yticklabels=[],
        cbar=False,                     # the notebook owns the colorbar
    )
    ax.collections[0].set_zorder(2)     # response mesh above the block patches

    if class_names:
        # Coloured block + name per contiguous chemical-class row-run, the
        # same `highconc.style.class_strip` mechanism
        # `panels.similarity.response_matrix` uses for Figure 2's own
        # odorant x glomerulus matrix -- `side="left"` here, away from this
        # panel's colourbar on the right.
        from .. import style as _style

        sorted_class = np.asarray(dataset.chem_class)[order.rows]
        present = sorted(set(sorted_class))
        display = {c: CHEM_CLASS_NAMES.get(c, c) for c in present}
        strip_colors = {display[c]: _clusters.class_color(c) for c in present}
        _style.class_strip(ax, [display[c] for c in sorted_class],
                           colors=strip_colors, orientation="y", side="left", offset=0,
                           labels=True, rotation=0,
                           fontsize=class_name_fontsize)

    if diagonal_blocks:
        for k in range(order.n_blocks):
            c0, c1 = order.col_edges[k], order.col_edges[k + 1]
            r0, r1 = order.row_edges[k], order.row_edges[k + 1]
            ax.add_patch(Rectangle((c0, r0), c1 - c0, r1 - r0,
                                   facecolor=block_color, edgecolor="none",
                                   zorder=1))
            # Rules along each block's top and bottom, spanning only its
            # columns. Skip the two that coincide with the edge of the matrix.
            if k != 0:
                ax.plot([c0, c1], [r0, r0], color=block_edge_color, ls="-",
                        lw=block_edge_lw, zorder=4)
            if k != order.n_blocks - 1:
                ax.plot([c0, c1], [r1, r1], color=block_edge_color, ls="-",
                        lw=block_edge_lw, zorder=4)

    if title is not None:
        ax.set_title(title, fontsize=title_fontsize)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
    # No set_aspect: the heatmap fills its grid cell, so it stays aligned with
    # whatever sits beside it and with the other concentration above or below.

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_edgecolor("black")

    # Column (glomerulus-cluster) boundaries. The horizontal odorant-*cluster*
    # lines were dropped deliberately — with both, the panel reads as a
    # checkerboard rather than as a set of blocks. The finer chemical-class
    # lines are a different, thinner annotation (only meaningful as tidy
    # blocks when `chem_order` was set) so they draw regardless.
    if boundaries:
        for b in order.col_boundaries:
            ax.axvline(b, color=boundary_color, linewidth=boundary_lw,
                       alpha=boundary_alpha, zorder=3)
        for b in order.row_class_boundaries:
            ax.axhline(b, color="white", linewidth=boundary_lw * 0.6,
                       alpha=boundary_alpha, zorder=3)

    return ax.collections[0]


# ---------------------------------------------------------------------------
# Fig 4A — glomerulus x glomerulus correlation
# ---------------------------------------------------------------------------

def glomerulus_correlation_matrix(
    ax,
    gsim,
    *,
    vmin: float = 0.0,
    vmax: float = 1.0,
    cmap=None,
    boundaries: bool = True,
    boundary_color: str = BOUNDARY_COLOR,
    boundary_lw: float = 0.7,
    boundary_alpha: float = BOUNDARY_ALPHA,
    cluster_boxes: bool = False,
    cluster_box_color: str = "green",
    cluster_box_lw: float = 1.2,
    cluster_labels: bool = True,
    label_fontsize: float = 9,
    title: str | None = None,
    xlabel: str | None = "glomeruli",
    ylabel: str | None = "glomeruli",
    title_fontsize: float = 15,
    axis_label_fontsize: float = 12,
):
    """Draw a `similarity.GlomerulusSimilarity` matrix into `ax`.

    Deliberately close to `panels.similarity.correlation_matrix`, and for the
    same reason: 4A is 3B's picture one axis over, so it gets the same 0-1
    scale, the same treatment of the NaN diagonal, and the same ramp. What
    differs is only the blocking — clusters, in Joel's numbering, rather than
    chemical classes — and the boundary colour, which matches the blue rules
    the rest of Figure 4 uses.

    Parameters
    ----------
    gsim : GlomerulusSimilarity
        From `highconc.similarity.glomerulus_correlation` — already permuted,
        diagonal already NaN.
    vmin, vmax : float
        Pinned to 0-1 like Figure 3's matrices. A few correlations are slightly
        negative; those are measurements and render at the floor rather than
        grey (see `cmap`).
    cluster_labels : bool
        Number each block along the top, in that cluster's colour.
    boundaries : bool
        Draw the cluster boundary lines -- full-width/height rules at every
        block edge, across the whole matrix. On by default.
    cluster_boxes : bool
        Draw an unfilled rectangle around each of the seven **diagonal**
        blocks only (the within-cluster squares), rather than rules spanning
        the full matrix. Off by default; independent of `boundaries` -- pass
        either, neither, or both.

    Returns the image, for the colourbar.
    """
    from .. import style as hs

    if cmap is None:
        # Same reasoning as Figure 3's matrices: sub-vmin cells are real,
        # slightly-negative correlations, so they belong at the bottom of the
        # ramp. Grey stays reserved for the NaN diagonal, which is not a
        # measurement.
        cmap = hs.sequential_cmap(under=hs.SEQ(0.0))

    im = hs.heatmap(ax, gsim.M, vmin=vmin, vmax=vmax, cmap=cmap,
                    aspect="equal")

    n = len(gsim.numbers)
    if boundaries:
        for b in gsim.boundaries:
            ax.axvline(b - 0.5, color=boundary_color, lw=boundary_lw,
                       alpha=boundary_alpha)
            ax.axhline(b - 0.5, color=boundary_color, lw=boundary_lw,
                       alpha=boundary_alpha)

    if cluster_boxes:
        edges = (0,) + gsim.boundaries + (n,)
        for e0, e1 in zip(edges, edges[1:]):
            ax.add_patch(Rectangle(
                (e0 - 0.5, e0 - 0.5), e1 - e0, e1 - e0,
                facecolor="none", edgecolor=cluster_box_color,
                lw=cluster_box_lw, zorder=5))

    ax.set_xticks([])
    ax.set_yticks([])
    if cluster_labels:
        from .. import style as _style
        for k, c in zip(gsim.blocks, gsim.centers):
            ax.text(c, -0.015 * n, str(k), ha="center", va="bottom",
                    fontsize=label_fontsize, color=_style.cluster_color(k))

    if title is not None:
        # Padded clear of the cluster numbers along the top.
        ax.set_title(title, fontsize=title_fontsize,
                     pad=18 if cluster_labels else 6)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=axis_label_fontsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=axis_label_fontsize)
    return im


# ---------------------------------------------------------------------------
# Fig 4D — chemical class x cluster summary
# ---------------------------------------------------------------------------

#: Written-out names for the chemical-class letters, in Group_num order.
CHEM_CLASS_NAMES = {
    "C": "acids", "D": "aldehydes", "E": "esters", "F": "furans",
    "K": "ketones", "L": "alcohols", "N": "amines", "P": "pyrazines",
    "R": "aromatics", "S": "thiazoles",
}


def _responding_glom_counts(dataset: Dataset, conc: str, order: str = "table"
                            ) -> tuple[np.ndarray, list[str]]:
    """For each (chemical class, cluster) cell, the count of distinct
    glomeruli in that cluster with a nonzero response to at least one
    odorant of that class, at concentration `conc`. Not normalised. A
    glomerulus that responds to several odorants of a class still counts
    once.

    `order` is the row (class) order -- see `highconc.clusters.class_order`.
    ``"table"`` (default) is alphabetical by letter; pass ``"matt"`` to match
    Figure 3/6's broad chemical-class block order instead.
    """
    nonzero = getattr(dataset, conc).values != 0

    letters = _clusters.class_order(order, present=list(dataset.chem_class))
    counts = np.zeros((len(letters), _clusters.N_CLUSTERS))
    for i, letter in enumerate(letters):
        odors = dataset.chem_class == letter
        for n in range(1, _clusters.N_CLUSTERS + 1):
            gloms = dataset.cluster_mask(n)
            if not (odors.any() and gloms.any()):
                continue
            block = nonzero[np.ix_(odors, gloms)]
            counts[i, n - 1] = block.any(axis=0).sum()
    return counts, letters


def class_cluster_glom_composition(dataset: Dataset, conc: str = "low",
                                   order: str = "table"
                                   ) -> tuple[np.ndarray, list[str]]:
    """Fig 4D's metric: each cluster's chemical-class composition by
    responding *glomerulus*, scaled so the class with the most responding
    glomeruli is 1.

    For each (chemical class, cluster) cell, counts the distinct glomeruli
    in that cluster with a nonzero response to at least one odorant of that
    class, then divides each column (cluster) by its own maximum: the
    chemical class reaching the most responding glomeruli in a cluster is 1,
    and the rest read as a fraction of it. Columns do **not** sum to 1 --
    this is a relative-magnitude scale, not a share of a whole. A glomerulus
    that responds to five odorants of a class still counts once (the count
    is of glomeruli, not (odorant, glomerulus) pairs); a glomerulus
    responding to more than one class's odorants is counted in each such
    class.

    Parameters
    ----------
    dataset : Dataset
    conc : {"low", "high"}
    order : str
        Row (class) order -- see `highconc.clusters.class_order`. ``"table"``
        (default) is alphabetical; ``"matt"`` matches Figure 3/6's broad
        chemical-class block order.

    Returns
    -------
    matrix : ndarray, shape (n_classes, 7)
        Each column's maximum is 1 (the class with the most responding
        glomeruli in that cluster); other cells are that class's count as a
        fraction of the max. A cluster with zero responding glomeruli
        anywhere (none in this dataset) is left all-zero rather than
        divided by zero.
    class_letters : list of str
        Row labels, the chemical-class letters in `order`.
    """
    if conc not in ("low", "high"):
        raise ValueError(f"conc must be 'low' or 'high', got {conc!r}")

    counts, letters = _responding_glom_counts(dataset, conc, order)

    maxima = counts.max(axis=0, keepdims=True)
    matrix = np.divide(counts, maxima, out=np.zeros_like(counts), where=maxima > 0)
    return matrix, letters


def cluster_summary(
    dataset: Dataset,
    ax,
    *,
    conc: str = "low",
    class_order: str = "table",
    cmap: str = CMAP_NAME,
    vmin: float = 0.0,
    vmax: float | None = None,
    empty_color: str = EMPTY_COLOR,
    class_names: bool = True,
    boundaries: bool = True,
    boundary_color: str = BOUNDARY_COLOR,
    boundary_lw: float = 0.8,
    boundary_alpha: float = BOUNDARY_ALPHA,
    xlabel: str | None = "glom cluster",
    ylabel: str | None = None,
    label_fontsize: float = 12,
    tick_fontsize: float = 10,
):
    """Fig 4D — each cluster's chemical-class preference profile, from
    `class_cluster_glom_composition`: for each (chemical class, cluster)
    cell, the distinct glomeruli in that cluster responding to that class,
    scaled so the cluster's top class is 1. Columns do not sum to 1 -- a
    relative magnitude within the cluster, not a share of a whole.

    Reading across a row gives one class's reach across the seven clusters.

    Columns are in the canonical **1-7** numbering (`highconc.clusters`),
    not the retired A-G letters.

    Parameters
    ----------
    class_order : str
        Row order -- see `highconc.clusters.class_order`. ``"table"``
        (default) is alphabetical; ``"matt"`` matches Figure 3/6's broad
        chemical-class block order.
    vmax : float, optional
        Colour ceiling. Left as None it is this call's own matrix maximum
        (by construction, 1 unless a cluster has zero responding
        glomeruli), so a low and a high call each autoscale independently.
        Pass a number to pin both to a shared ceiling instead.
    class_names : bool
        Label rows with written-out names ("aldehydes") rather than the raw
        Group letters.
    boundaries : bool
        Draw the vertical cluster-boundary rules. On by default.

    Returns
    -------
    QuadMesh
        For handing to `fig.colorbar`.
    """
    matrix, letters = class_cluster_glom_composition(dataset, conc, class_order)

    if vmax is None:
        vmax = float(matrix.max()) if matrix.max() > 0 else 1.0

    panel_cmap = plt.get_cmap(cmap).copy()
    panel_cmap.set_under(empty_color)
    panel_cmap.set_bad(empty_color)

    yticklabels = (
        [CHEM_CLASS_NAMES.get(c, c) for c in letters] if class_names else letters
    )
    sns.heatmap(
        matrix, cmap=panel_cmap, vmin=vmin, vmax=vmax, ax=ax,
        xticklabels=list(_clusters.LABELS_BY_JOEL),
        yticklabels=yticklabels,
        cbar=False,                     # the notebook owns the colorbar
    )
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    ax.set_ylabel("" if ylabel is None else ylabel, fontsize=label_fontsize)
    ax.tick_params(axis="x", labelsize=tick_fontsize, rotation=0)
    ax.tick_params(axis="y", labelsize=tick_fontsize, rotation=0)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_edgecolor("black")

    if boundaries:
        for x in range(1, _clusters.N_CLUSTERS):
            ax.axvline(x, color=boundary_color, linewidth=boundary_lw,
                       alpha=boundary_alpha)

    return ax.collections[0]


# ---------------------------------------------------------------------------
# Per-cluster glomerular maps
# ---------------------------------------------------------------------------

#: Colour ceiling for the topography maps, in mean dF/F.
#:
#: **This is a LOW-concentration ceiling.** At low concentration, mean dF/F
#: runs median 1.19 / p95 4.85 / max 13.88, so a ceiling of 15.0 sits *above
#: the data maximum*: nothing reaches full scale and the median glomerulus
#: renders at 8% of the ramp, which makes every map look uniformly near-zero.
#: At 5.0 the median sits at 24% of scale and 18 of 376 glomeruli (4.8%)
#: saturate.
#:
#: **Do not reuse it for `conc="high"`.** High runs median 5.30 / p95 19.82, so
#: 51% of glomeruli would saturate. Use `TOPO_VMAX_HIGH`, which clips the same
#: ~5% at that concentration. A map drawing low and high on one shared scale
#: needs a single ceiling that is neither of these.
TOPO_VMIN, TOPO_VMAX = 0.0, 5.0
TOPO_VMAX_HIGH = 20.0
CONTEXT_COLOR = "darkgray"

#: How the glomeruli that are *not* in the cluster being drawn are marked.
#:
#: ``"context"`` keeps them small and semi-transparent behind the cluster:
#: anatomical context and nothing more. ``"requested"`` is Figure 5A's — larger
#: than the cluster's own points, grey, with no border, and pale and
#: transparent (``0.85`` at ``alpha=0.5``). Drawing them larger than the
#: cluster's points inverts the usual figure-ground relation, so the cluster
#: reads as the sparse thing it is against a dense grey bulb.
#:
#: ``alpha=0.35`` was tried and rejected: the bulb's outline starts to
#: dissolve, and the map stops reading as a bulb.
#:
#: A dict of `ax.scatter` keywords may be passed instead of a name.
NONRESPONDER_MARKERS: dict[str, dict] = {
    "context":   dict(s=16, c=CONTEXT_COLOR, alpha=0.35, linewidths=0),
    "requested": dict(s=46, c="0.85", alpha=0.5, linewidths=0),
    "hollow":    dict(s=16, facecolors="none", edgecolors=CONTEXT_COLOR,
                      linewidths=0.4, alpha=0.5),
}

def _nonresponder_kw(spec) -> dict:
    """Resolve `nonresponder_marker` to a dict of scatter keywords."""
    if isinstance(spec, str):
        try:
            return dict(NONRESPONDER_MARKERS[spec])
        except KeyError:
            raise ValueError(
                f"unknown nonresponder_marker {spec!r}; expected one of "
                f"{sorted(NONRESPONDER_MARKERS)} or a dict of scatter kwargs"
            ) from None
    return dict(spec)


def topography(
    dataset: Dataset,
    ax,
    cluster: int,
    *,
    conc: str = "low",
    vmin: float = TOPO_VMIN,
    vmax: float = TOPO_VMAX,
    cmap: str = CMAP_NAME,
    size: float = 30,
    edgecolor: str | None = None,
    linewidth: float = 0.0,
    context: bool = True,
    color: str | None = None,
    nonresponder_marker="context",
    context_color: str | None = None,
    context_size: float | None = None,
    context_alpha: float | None = None,
    context_linewidth: float | None = None,
    title: str | bool = True,
    ylabel: str | None = None,
    title_fontsize: float = 13,
    label_fontsize: float = 13,
):
    """One cluster's glomeruli on the dorsal bulb, coloured by mean response.

    Every other glomerulus is drawn in grey behind it as anatomical context, so
    each map shows where the cluster sits rather than just its shape. All four
    subjects are superimposed in raw pixel coordinates, which is meaningful
    here because the field of view is consistent between them.

    Parameters
    ----------
    cluster : int
        Which cluster to highlight, in **Joel's 1-7** numbering.
    conc : {"low", "high"}
        Which concentration's mean response colours the points.
    vmin, vmax : float
        Colour limits in mean dF/F. The default ceiling sits below the data
        maximum on purpose (see `TOPO_VMAX`); keep it the same for both
        concentrations or the broadening stops being visible.
    edgecolor, linewidth : optional
        Outline on the cluster's own points. Off by default. Worth turning on
        (white, ~0.6) whenever the non-cluster glomeruli are drawn solid —
        the bottom of the sequential ramp is nearly white, so a weakly
        responding cluster member would otherwise read as a hole in the grey
        rather than as a point.
    context : bool
        Draw the non-cluster glomeruli. Off gives the cluster alone.
    color : str, optional
        Draw every cluster member as one solid colour, so the map shows
        *membership* and nothing else -- no response value is encoded, and
        `conc`, `vmin`, `vmax` and `cmap` are ignored. Off by default (points
        coloured by mean dF/F as above). The returned collection then has no
        colour scale, so there is nothing to hand to `fig.colorbar`.
    nonresponder_marker : str or dict
        How those non-cluster glomeruli are drawn — a key of
        `NONRESPONDER_MARKERS` or an explicit dict of scatter keywords.
        ``"context"`` (default) keeps them small and faint behind the cluster;
        ``"requested"`` is Figure 5A's treatment, larger and solid grey with no
        border. `context_color`, `context_size`, `context_alpha` and
        `context_linewidth` override individual fields of whichever preset
        is in force.
    title : str or bool
        True titles the map with the cluster number; False or None leaves it
        off; a string is used as given.
    ylabel : str, optional
        Row label — used on the leftmost map of a row ("low", "high").

    Returns
    -------
    PathCollection
        The cluster's points, for handing to `fig.colorbar`.
    """
    if conc not in ("low", "high"):
        raise ValueError(f"conc must be 'low' or 'high', got {conc!r}")

    values = getattr(dataset, conc).values.mean(axis=0)
    mask = dataset.cluster_mask(cluster)
    x, y = dataset.xy[:, 0], dataset.xy[:, 1]

    if context:
        ctx = _nonresponder_kw(nonresponder_marker)
        if context_color is not None:
            if "facecolors" in ctx:
                ctx["edgecolors"] = context_color
            else:
                ctx["c"] = context_color
        if context_size is not None:
            ctx["s"] = context_size
        if context_alpha is not None:
            ctx["alpha"] = context_alpha
        if context_linewidth is not None:
            ctx["linewidths"] = context_linewidth
        ax.scatter(x[~mask], y[~mask], zorder=1, **ctx)

    if color is not None:
        points = ax.scatter(x[mask], y[mask], color=color, s=size,
                            edgecolors=edgecolor, linewidths=linewidth,
                            zorder=2)
    else:
        points = ax.scatter(x[mask], y[mask], c=values[mask], cmap=cmap,
                            vmin=vmin, vmax=vmax, s=size, edgecolors=edgecolor,
                            linewidths=linewidth, zorder=2)

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():     # no bounding box on a bulb map
        spine.set_visible(False)

    if title is True:
        ax.set_title(str(cluster), fontsize=title_fontsize)
    elif title:
        ax.set_title(title, fontsize=title_fontsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=label_fontsize)

    return points


# ---------------------------------------------------------------------------
# Two-cluster spatial contrast (the A-vs-E breakout)
# ---------------------------------------------------------------------------

def cluster_contrast(
    dataset: Dataset,
    ax,
    first: int = 1,
    second: int = 2,
    *,
    conc: str = "high",
    cmap=None,
    size: float = 99,
    edgecolor: str = "0.4",
    edge_lw: float = 0.4,
    context: bool = True,
    context_color: str = CONTEXT_COLOR,
    context_alpha: float = 0.30,
    territory: bool = True,
    territory_level: float = 0.35,
    territory_alpha: float = 0.22,
    bw_adjust: float = 1.1,
    labels: bool = True,
    label_fontsize: float = 13,
):
    """Do two clusters occupy different territory on the bulb?

    Draws the union of two clusters' glomeruli, each coloured by which of the
    two clusters' preferred odorants it responds to more::

        index = resp_first / (resp_first + resp_second)

    where a "first-cluster odorant" is one the first cluster's canonical
    spectrum ranks above the second's, and vice versa. 0 is wholly
    second-leaning, 1 wholly first-leaning. Optional KDE shading outlines each
    cluster's territory as a guide to the eye — it is a density contour, not a
    statistical claim about a boundary.

    Parameters
    ----------
    first, second : int
        The two clusters to contrast, in **Joel's 1-7** numbering. The default
        pair is the one the manuscript figure breaks out. It was written as
        A-vs-E before the renumbering: legacy A is component 0, which is
        cluster 1, and legacy E is component 4, which is cluster **2** — not
        cluster 5, which is what reading the letter E as the fifth cluster
        would wrongly give.
    conc : {"low", "high"}
        Which concentration defines both the canonical spectra and the
        responses. High by default: at low concentration too few glomeruli are
        active for the index to be stable.
    cmap : Colormap, optional
        Diverging map, second-leaning to first-leaning. Defaults to
        `highconc.style.DIV`, the house teal-to-red diverging map. Independent
        of `highconc.style.SEQ` (the sequential ramp) — DIV is for signed
        quantities and is not changed by picking a different sequential map.
    territory : bool
        Draw the KDE territory shading behind the points.

    Returns
    -------
    PathCollection
        The coloured points, for handing to `fig.colorbar`.
    """
    if conc not in ("low", "high"):
        raise ValueError(f"conc must be 'low' or 'high', got {conc!r}")
    if cmap is None:
        from .. import style as _style
        cmap = _style.DIV

    data = getattr(dataset, conc).values
    mask_first = dataset.cluster_mask(first)
    mask_second = dataset.cluster_mask(second)

    # Canonical spectrum of each cluster: its mean response per odorant.
    spec_first = data[:, mask_first].mean(axis=1)
    spec_second = data[:, mask_second].mean(axis=1)
    odors_first = spec_first > spec_second
    odors_second = ~odors_first

    resp_first = data[odors_first, :].mean(axis=0)
    resp_second = data[odors_second, :].mean(axis=0)
    index = resp_first / (resp_first + resp_second + 1e-9)

    both = mask_first | mask_second
    x, y = dataset.xy[:, 0], dataset.xy[:, 1]

    if context:
        ax.scatter(x[~both], y[~both], c=context_color, s=size,
                   alpha=context_alpha, linewidths=0, zorder=1)
    points = ax.scatter(x[both], y[both], c=index[both], cmap=cmap,
                        vmin=0, vmax=1, s=size, linewidths=edge_lw,
                        edgecolors=edgecolor, zorder=2)

    if territory:
        # Shade after the scatter so the axis limits are already set by the
        # data, then restore them — a KDE will otherwise expand them.
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        for mask, color in ((mask_first, cmap(0.95)), (mask_second, cmap(0.05))):
            base = dict(x=x[mask], y=y[mask], ax=ax, bw_adjust=bw_adjust)
            sns.kdeplot(**base, fill=True, color=color, alpha=territory_alpha,
                        levels=5, thresh=territory_level, zorder=0.5)
            sns.kdeplot(**base, fill=False, color=color, linewidths=1.5,
                        levels=[territory_level], zorder=0.6)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    if labels:
        # The words themselves carry the colour, so no marker is needed. A
        # white halo keeps them readable where the glomeruli reach the corner —
        # they do, at this panel's size.
        halo = [patheffects.withStroke(linewidth=2.5, foreground="white")]
        for i, cluster in enumerate((first, second)):
            ax.text(0.02, 0.98 - 0.08 * i, f"cluster {cluster}",
                    transform=ax.transAxes, ha="left", va="top",
                    color=cmap(1.0 - i), fontsize=label_fontsize,
                    zorder=5, path_effects=halo)

    return points


# ---------------------------------------------------------------------------
# Per-subject reproducibility — is the block structure an artifact of pooling?
# ---------------------------------------------------------------------------

#: Padding on the shared spatial frame, as a fraction of the data's extent. All
#: four subjects' maps are drawn in the same frame, so they are comparable by
#: position and not merely by shape.
BULB_PAD = 0.04
SUBJECT_CONTEXT_COLOR = "0.84"


def bulb_limits(dataset: Dataset, pad: float = BULB_PAD):
    """(xlim, ylim) covering every glomerulus, padded. One frame for all maps."""
    x, y = dataset.xy[:, 0], dataset.xy[:, 1]
    (x0, x1), (y0, y1) = (x.min(), x.max()), (y.min(), y.max())
    return ((x0 - pad * (x1 - x0), x1 + pad * (x1 - x0)),
            (y0 - pad * (y1 - y0), y1 + pad * (y1 - y0)))


def per_subject_heatmap(
    dataset: Dataset,
    ax,
    subject,
    *,
    cluster: int | None = None,
    conc: str = "low",
    order: BiclusterOrder | None = None,
    vmin: float = VMIN,
    vmax: float = VMAX,
    cmap: str = CMAP_NAME,
    empty_color: str = "gainsboro",
    col_normalize: bool = False,
    raw_floor: float | None = None,
    boundaries: bool = True,
    boundary_color: str = "#1f3fff",
    boundary_lw: float = 0.7,
    boundary_alpha: float = 0.45,
    cluster_labels: bool = True,
    label_fontsize: float = 7.5,
    title: str | bool = True,
    title_fontsize: float = 12,
    ylabel: str | None = None,
    ylabel_fontsize: float = 10,
):
    """One mouse's own glomeruli, sorted into clusters, against all 59 odorants.

    The point of the panel: the diagonal block structure shows up *within* each
    subject, so it is not something pooling four fields of view manufactured.
    Only the columns are subject-specific — the odorant rows are the shared
    sort, so the four panels line up row for row.

    Parameters
    ----------
    subject : str
        Which mouse, as it appears in `dataset.subjects`.
    order : BiclusterOrder, optional
        Supplies the odorant row sort and its boundaries. Pass the shared one
        so all four subjects (and the pooled figures) agree; omitted, it is
        derived from `dataset`.
    cluster : int, optional
        Restrict the columns to a single cluster, in Joel's 1-7 numbering.
        This is what Figure 4C is: one cluster's glomeruli, four subjects side
        by side, so that the cluster's tuning can be seen to repeat across
        mice. Omitted, the panel shows every cluster (Supp. Fig. B's use).
    col_normalize, raw_floor
        Same meaning as in `biclustered`: peak-normalise each glomerulus to
        [0, 1] by its own maximum over all odorants, on the linear raw dF/F
        values, and draw on a LINEAR scale (pass ``vmin=0.0, vmax=1.0``);
        `raw_floor` greys out cells whose raw value was below it. Because this
        panel always carries all 59 odorant rows, a glomerulus normalises to
        exactly the same values here as in `biclustered`. Off by default (raw
        dF/F on a log scale), so the other figures using this panel are
        unchanged.
    cluster_labels : bool
        Write each block's cluster number above it, in that cluster's colour.
        With `cluster` set there is only one block, so this writes one number.
    boundaries : bool
        Draw the cluster-boundary rules (both the vertical glomerulus-cluster
        lines and the lighter horizontal odorant-cluster ones). On by default.

    Returns
    -------
    AxesImage or None
        For handing to `fig.colorbar`. `None` when the subject has no
        glomeruli in the requested cluster — an empty frame is drawn instead
        of raising, because the absence is itself a result.
    """
    if conc not in ("low", "high"):
        raise ValueError(f"conc must be 'low' or 'high', got {conc!r}")
    if order is None:
        order = bicluster_order(dataset)

    mask = dataset.subject_mask(subject)
    if cluster is not None:
        mask = mask & dataset.cluster_mask(cluster)
    cols = np.where(mask)[0]
    cols = cols[np.argsort(dataset.cluster_number[cols], kind="stable")]
    data = getattr(dataset, conc).values[np.ix_(order.rows, cols)]
    numbers = dataset.cluster_number[cols]

    if len(cols) == 0:
        # A real result, not an error: cluster 5 is present in only two of the
        # four mice. Draw the empty frame rather than dropping the panel, so a
        # row of subjects stays a row of subjects and the absence is visible.
        _empty_subject_panel(ax, order, empty_color, boundary_color,
                             boundary_lw, boundary_alpha, cluster=cluster,
                             boundaries=boundaries)
        _finish_subject_panel(ax, subject, title, title_fontsize,
                              ylabel, ylabel_fontsize)
        return None

    if col_normalize:
        peak = data.max(0, keepdims=True)
        peak = np.where(peak > 0, peak, np.inf)
        normed = data / peak
        if raw_floor is not None:
            normed = np.where(data >= raw_floor, normed, np.nan)
        data = normed

    heat_cmap = plt.get_cmap(cmap).copy()
    heat_cmap.set_under(empty_color)
    heat_cmap.set_bad(empty_color)      # LogNorm masks zeros; raw_floor uses this too
    norm = (plt.Normalize(vmin=vmin, vmax=vmax) if col_normalize
            else LogNorm(vmin=vmin, vmax=vmax))
    im = ax.imshow(data, aspect="auto", cmap=heat_cmap, norm=norm,
                   interpolation="nearest")
    ax.set_facecolor(empty_color)

    edges = np.concatenate(
        [[0], np.where(np.diff(numbers) != 0)[0] + 1, [len(numbers)]]
    )
    if boundaries:
        for b in edges[1:-1]:
            ax.axvline(b - 0.5, color=boundary_color, lw=boundary_lw,
                       alpha=boundary_alpha)
    if cluster_labels:
        from .. import style as _style
        for a, b in zip(edges[:-1], edges[1:]):
            n = int(numbers[a])
            ax.text((a + b - 1) / 2, -1.5, str(n), ha="center", va="bottom",
                    fontsize=label_fontsize, color=_style.cluster_color(n))

    # Odorant-cluster rules. Lighter than the vertical ones: they are a reading
    # aid for the block structure, not the panel's subject.
    if boundaries:
        for b in order.row_boundaries:
            ax.axhline(b - 0.5, color=boundary_color, lw=boundary_lw * 0.7,
                       alpha=boundary_alpha * 0.7)

    _finish_subject_panel(ax, subject, title, title_fontsize,
                          ylabel, ylabel_fontsize)
    return im


def _empty_subject_panel(ax, order, empty_color, boundary_color,
                         boundary_lw, boundary_alpha, cluster=None,
                         boundaries=True):
    """The frame a subject panel keeps when that subject has no such glomeruli.

    Deliberately NOT `empty_color`. In a populated panel grey means "below the
    colour floor" — a glomerulus that was measured and did not respond. If an
    absent subject were also grey, the two would be indistinguishable and the
    absence would read as a failed render. So this panel is white, hatched, and
    says in words what is missing: structural absence has to look structurally
    different from a floored measurement.
    """
    n = len(order.rows)
    ax.set_facecolor("white")
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    ax.add_patch(Rectangle(
        (-0.5, -0.5), 1.0, n, transform=ax.transData, facecolor="white",
        edgecolor="0.82", hatch="////", lw=0.8, zorder=0))
    if boundaries:
        for b in order.row_boundaries:
            ax.axhline(b - 0.5, color=boundary_color, lw=boundary_lw * 0.7,
                       alpha=boundary_alpha * 0.35, zorder=1)
    label = ("no glomeruli\nin this subject" if cluster is None
             else f"no cluster-{cluster}\nglomeruli\nin this subject")
    ax.text(0.0, (n - 1) / 2, label, ha="center", va="center", rotation=90,
            fontsize=7.5, color="0.35", zorder=2)


def _finish_subject_panel(ax, subject, title, title_fontsize,
                          ylabel, ylabel_fontsize):
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if title is True:
        ax.set_title(f"S{subject}", fontsize=title_fontsize, pad=10)
    elif title:
        ax.set_title(title, fontsize=title_fontsize, pad=10)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=ylabel_fontsize)


def block_means(
    dataset: Dataset,
    subject=None,
    *,
    conc: str = "low",
    normalize: str = "column",
) -> np.ma.MaskedArray:
    """Mean response per (odorant cluster, glom cluster) block, 7 x 7.

    Rows are odorant clusters and columns glomerular clusters, both in Joel's
    1-7 order, so the reproducible signal sits on the diagonal.

    Parameters
    ----------
    subject : str, optional
        Restrict to one mouse's glomeruli. None pools all four.
    normalize : {"column", "none"}
        ``"column"`` divides each glom cluster's column by its own maximum, so
        every cluster's strongest odorant group reads as 1.0 and the diagonal
        is comparable between clusters of very different size and brightness.
        That is what makes this a reproducibility panel rather than a
        brightness panel. ``"none"`` leaves raw mean dF/F.

    Returns
    -------
    MaskedArray, shape (7, 7)
        Blocks with no glomeruli — subjects 3 and 4 have no cluster-5
        glomeruli — are masked, and draw as the colormap's `bad` colour.
    """
    if conc not in ("low", "high"):
        raise ValueError(f"conc must be 'low' or 'high', got {conc!r}")
    if normalize not in ("column", "none"):
        raise ValueError(f"normalize must be 'column' or 'none', got {normalize!r}")

    data = getattr(dataset, conc).values
    in_subject = (
        np.ones(dataset.n_gloms, dtype=bool) if subject is None
        else dataset.subject_mask(subject)
    )
    odorant_number = np.array(
        [_clusters.to_joel(k) for k in dataset.odorant_cluster], dtype=int
    )

    M = np.full((_clusters.N_CLUSTERS, _clusters.N_CLUSTERS), np.nan)
    for col, glom_n in enumerate(range(1, _clusters.N_CLUSTERS + 1)):
        gloms = dataset.cluster_mask(glom_n) & in_subject
        if not gloms.any():
            continue                    # this mouse has none of this cluster
        for row, odor_n in enumerate(range(1, _clusters.N_CLUSTERS + 1)):
            odors = odorant_number == odor_n
            if odors.any():
                M[row, col] = data[np.ix_(odors, gloms)].mean()

    if normalize == "column":
        with np.errstate(invalid="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                ceiling = np.nanmax(M, axis=0, keepdims=True)
        M = M / np.where(ceiling > 0, ceiling, np.nan)
    return np.ma.masked_invalid(M)


def per_subject_blocks(
    dataset: Dataset,
    ax,
    subject=None,
    *,
    conc: str = "low",
    normalize: str = "column",
    gamma: float = 2.0,
    vmin: float = 0.02,
    vmax: float = 1.0,
    cmap: str = CMAP_NAME,
    empty_color: str = "gainsboro",
    boundaries: bool = True,
    grid_color: str = "white",
    grid_lw: float = 0.6,
    title: str | bool = True,
    title_fontsize: float = 11,
    xlabel: str | None = "glom cluster",
    ylabel: str | None = None,
    label_fontsize: float = 9,
    tick_fontsize: float = 7,
):
    """One mouse's 7 x 7 block-average matrix — the diagonal, subject by subject.

    Same claim as `per_subject_heatmap`, coarsened: if cluster k's glomeruli
    really do prefer odorant cluster k, the bright cell is on the diagonal in
    every mouse separately.

    Parameters
    ----------
    gamma : float
        Exponent of the colour scale. The pre-refactor default of 2.0 pushes
        the mid-range down so the diagonal separates from its neighbours; 1.0
        is a linear scale, which is the honest choice if the panel is ever read
        quantitatively rather than as a pattern. It is a display gamma either
        way — `block_means` returns the numbers.
    normalize : {"column", "none"}
        See `block_means`.
    boundaries : bool
        Draw the white cell-separator grid lines. On by default.

    Returns
    -------
    AxesImage
        For handing to `fig.colorbar`.
    """
    from .. import style as _style

    M = block_means(dataset, subject, conc=conc, normalize=normalize)
    block_cmap = plt.get_cmap(cmap).copy()
    block_cmap.set_under(empty_color)
    block_cmap.set_bad(empty_color)

    from matplotlib.colors import PowerNorm
    im = ax.imshow(M, cmap=block_cmap,
                   norm=PowerNorm(gamma, vmin=vmin, vmax=vmax),
                   interpolation="nearest")
    ax.set_facecolor(empty_color)

    labels = list(_clusters.LABELS_BY_JOEL)
    ax.set_xticks(range(_clusters.N_CLUSTERS))
    ax.set_yticks(range(_clusters.N_CLUSTERS))
    ax.set_xticklabels(labels, fontsize=tick_fontsize)
    ax.set_yticklabels(labels, fontsize=tick_fontsize)
    for n, tick in zip(range(1, _clusters.N_CLUSTERS + 1), ax.get_xticklabels()):
        tick.set_color(_style.cluster_color(n))
    ax.tick_params(length=0)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_edgecolor("0.4")
    if boundaries:
        for g in range(1, _clusters.N_CLUSTERS):
            ax.axvline(g - 0.5, color=grid_color, lw=grid_lw)
            ax.axhline(g - 0.5, color=grid_color, lw=grid_lw)

    if title is True:
        ax.set_title("pooled" if subject is None else f"S{subject}",
                     fontsize=title_fontsize, pad=4)
    elif title:
        ax.set_title(title, fontsize=title_fontsize, pad=4)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=label_fontsize)

    return im


def per_subject_topography(
    dataset: Dataset,
    ax,
    subject,
    cluster: int,
    *,
    color=None,
    size: float = 54,
    edgecolor: str = "black",
    linewidth: float = 0.5,
    context: bool = True,
    context_color: str = SUBJECT_CONTEXT_COLOR,
    context_size: float = 14,
    context_linewidth: float = 0.4,
    context_alpha: float | None = None,
    limits=None,
    empty_marker: str = "—",
    empty_fontsize: float = 10,
    title: str | bool = False,
    title_fontsize: float = 12,
    ylabel: str | None = None,
    label_fontsize: float = 12,
):
    """Where one mouse's members of one cluster sit on the bulb.

    Unlike `topography`, the points carry no response value — they are one flat
    colour, the cluster's own. The question is spatial: does cluster k occupy
    the same territory in each mouse? Every glomerulus from *every* subject is
    drawn grey behind, so all 28 maps share one visual backdrop and can be
    compared by position.

    Parameters
    ----------
    cluster : int
        Which cluster, in **Joel's 1-7** numbering.
    color : optional
        Point colour. Defaults to the cluster's own colour, which is indexed by
        component so it matches the exports and the other panels' block labels.
    limits : ((x0, x1), (y0, y1)), optional
        The shared frame. Omitted, it is `bulb_limits(dataset)` — pass the
        precomputed value when drawing a grid of these, which is the usual case.
    empty_marker : str
        Drawn when this mouse has no glomeruli in this cluster, which is a real
        result and not a failure: subjects 3 and 4 have no cluster-5 glomeruli.
        Pass None to leave the cell blank.

    Returns
    -------
    PathCollection or None
        The cluster's points, or None when the mouse has none.
    """
    from .. import style as _style

    if color is None:
        color = _style.cluster_color(cluster)
    if limits is None:
        limits = bulb_limits(dataset)

    x, y = dataset.xy[:, 0], dataset.xy[:, 1]
    if context:
        ax.scatter(x, y, facecolors="none", edgecolors=context_color,
                   s=context_size, linewidths=context_linewidth,
                   alpha=context_alpha, zorder=1)

    mask = dataset.subject_mask(subject) & dataset.cluster_mask(cluster)
    points = None
    if mask.any():
        points = ax.scatter(x[mask], y[mask], color=color, s=size,
                            edgecolors=edgecolor, linewidths=linewidth, zorder=2)
    elif empty_marker:
        ax.text(0.5, 0.5, empty_marker, transform=ax.transAxes,
                ha="center", va="center", color="0.6",
                fontsize=empty_fontsize)

    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    if title is True:
        ax.set_title(str(cluster), fontsize=title_fontsize,
                     color=color, pad=4)
    elif title:
        ax.set_title(title, fontsize=title_fontsize, color=color, pad=4)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=label_fontsize, rotation=0,
                      ha="right", va="center", labelpad=10)

    return points
