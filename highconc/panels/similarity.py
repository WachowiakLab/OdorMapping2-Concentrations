"""Panels for Figure 3's top row — the odorant x odorant similarity matrices.

One panel function, used three times: 3B (response correlation at low), 3C (at
high) and 3D (MACCS chemical similarity) are the same picture of different
matrices, and drawing them with one function is what guarantees they share a
colour scale, a block order and a class palette.

The same module also holds the mirrored bar spectra — 3E/F/G (low vs high) and
6B/C/D (experimental vs predicted high). They are one function, `mirror_spectrum`,
called two ways, for the same reason the matrices are: the manuscript reads
6B/C/D against 3E/F/G, and that comparison only holds if they are literally the
same picture with one series swapped.

Composition — the row of three, the shared colourbar, the panel letters — stays
in the notebook, as everywhere else in this package.
"""

from __future__ import annotations

import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

from .. import clusters as cl
from .. import style as hs

__all__ = ["correlation_matrix", "response_matrix", "class_colors",
           "similarity_colorbar", "mirror_spectrum", "response_spectrum",
           "prediction_spectrum", "shared_spectrum_limit", "SPECTRUM_YLABELS"]


def class_colors(sim) -> dict:
    """``{class letter: colour}`` for one similarity matrix's blocks.

    Keyed by letter, from `highconc.clusters.CHEM_CLASS_COLORS`, so a class
    keeps its colour no matter which block order the panel is drawn in.
    """
    return {c: cl.class_color(c) for c in sim.blocks}


def correlation_matrix(ax, sim, *, vmin=0.0, vmax=1.0, cmap=None,
                       title=None, names=False, boundaries=True,
                       strips=("x",), strip_labels=True, strip_fontsize=None,
                       rotation=90, tick_step=3,
                       box_boundaries=None, box_color="green", box_lw=1.2):
    """Draw one odorant x odorant similarity matrix.

    Parameters
    ----------
    sim : OdorantSimilarity
        From `highconc.similarity` — already permuted, diagonal already NaN.
    vmin, vmax : float
        Pinned to 0-1 by default, and they should stay pinned: 3B, 3C and 3D
        are only comparable if they share one scale, and autoscaling 3B (which
        is sparse) against 3C (which is not) would draw away the figure's
        entire claim.
    cmap : Colormap or str, optional
        Defaults to the house sequential ramp. Matt's 2026-08-23 draft used
        ``"inferno"``; pass it to reproduce that look exactly.
    names : bool
        Label every row and column with the odorant name. Legible only on a
        large panel — off by default, which is what the figure uses.
    strips : sequence
        Which axes get the coloured chemical-class strip: any of ``"x"``,
        ``"y"``. Default is the x-axis only, matching 3B and 3C; 3D in Matt's
        draft carries the strip on both.
    tick_step : int or None
        Numeric tick spacing, drawn **only on axes with no class strip** —
        a strip and a row of numbers in the same margin collide, and the strip
        is the more informative of the two. `None` turns numbers off entirely.
    strip_fontsize : float, optional
        Class-name label size on the strip. `None` uses `class_strip`'s own
        default (`FS_TICK - 1`); pass a number to override, e.g. for a small
        fixed-size panel where the house tick size doesn't fit.
    box_boundaries : sequence of int, optional
        Interior cut points (same convention as `sim.boundaries`) at which to
        draw unfilled diagonal boxes, same mechanism as Figure 4A's
        `cluster_boxes`: `(0,) + box_boundaries + (n,)` gives the edges, and
        one square is drawn between each consecutive pair. Unlike `sim`'s own
        `boundaries` (every chemical-class edge, drawn as full-width/height
        rules), this is an independent, caller-chosen grouping — pass a
        coarser or finer set of cuts than the class blocks. `None` (default)
        draws no boxes.
    box_color, box_lw : str, float
        Box edge colour and line width. Defaults match `cluster_boxes`'
        (`"green"`, `1.2`).

    Returns the image, which you need for the colourbar.
    """
    M = np.asarray(sim.M, float)
    if cmap is None:
        # The house ramp, but with the sub-floor colour set to the ramp's own
        # zero instead of grey. A handful of correlations are very slightly
        # negative (min ~ -0.2); those are MEASUREMENTS and belong at the
        # bottom of the scale. Leaving `under` grey would paint a third of the
        # panel the colour that means "not measured". The NaN diagonal keeps
        # the grey, which is right — self-similarity is not a measurement.
        cmap = hs.sequential_cmap(under=hs.SEQ(0.0))
    bounds = sim.boundaries if boundaries else ()
    im = hs.heatmap(ax, M, vmin=vmin, vmax=vmax, cmap=cmap,
                    row_boundaries=bounds, col_boundaries=bounds,
                    aspect="equal")

    if box_boundaries is not None:
        n = len(sim.classes)
        edges = (0,) + tuple(box_boundaries) + (n,)
        for e0, e1 in zip(edges, edges[1:]):
            ax.add_patch(Rectangle(
                (e0 - 0.5, e0 - 0.5), e1 - e0, e1 - e0,
                facecolor="none", edgecolor=box_color, lw=box_lw, zorder=5))

    if names:
        ax.set_xticks(range(len(sim.names)))
        ax.set_yticks(range(len(sim.names)))
        ax.set_xticklabels(sim.names, rotation=90, fontsize=hs.FS_TICK - 3)
        ax.set_yticklabels(sim.names, fontsize=hs.FS_TICK - 3)
    else:
        t = list(range(0, len(sim.names), int(tick_step))) if tick_step else []
        ax.set_xticks([] if "x" in strips else t)
        ax.set_yticks([] if "y" in strips else t)
        ax.tick_params(labelsize=hs.FS_TICK - 1)

    colors = class_colors(sim)
    labels = {c: cl.class_name(c) for c in sim.blocks}
    for orientation in strips:
        hs.class_strip(ax, [labels[c] for c in sim.classes],
                       colors={labels[c]: v for c, v in colors.items()},
                       orientation=orientation, labels=strip_labels,
                       fontsize=strip_fontsize,
                       rotation=rotation if orientation == "x" else 0)

    ax.set_title(title if title is not None else sim.what,
                 fontsize=hs.FS_LABEL)
    return im


def response_matrix(ax, mat, *, vmin=0.0, vmax=1.0, cmap=None, title=None,
                    boundaries=True, strip_labels=True, strip_fontsize=None):
    """Figure 2: one odorant x glomerulus response matrix.

    Parameters
    ----------
    mat : OdorantGlomMatrix
        From `highconc.similarity.burton_response_matrix` -- rows already
        permuted into class-block order, columns untouched.
    vmin, vmax : float
        Pinned to 0-1 by default, matching Figure 3's correlation panels --
        panels A and B should share one scale for the same reason 3B/3C do.

    Unlike `correlation_matrix`, this is not square: it's a response matrix,
    not a similarity matrix, so there's no diagonal to grey out and only the
    row (odorant) axis carries chemical-class structure -- the strip goes on
    `y`, and the glomerulus columns are left unlabeled (Burton's own 113
    glomeruli, with no ordering or identity information supplied).

    Returns the image, which you need for the colourbar.
    """
    M = np.asarray(mat.M, float)
    if cmap is None:
        # Same reasoning as `correlation_matrix`: this is a sparse,
        # non-negative response matrix, so zero is a measurement and belongs
        # at the ramp's own floor, not at the "not measured" grey.
        cmap = hs.sequential_cmap(under=hs.SEQ(0.0))
    bounds = mat.boundaries if boundaries else ()
    im = hs.heatmap(ax, M, vmin=vmin, vmax=vmax, cmap=cmap,
                    row_boundaries=bounds, aspect="auto")

    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(labelsize=hs.FS_TICK - 1)

    colors = class_colors(mat)
    labels = {c: cl.class_name(c) for c in mat.blocks}
    hs.class_strip(ax, [labels[c] for c in mat.classes],
                   colors={labels[c]: v for c, v in colors.items()},
                   orientation="y", labels=strip_labels, rotation=0,
                   fontsize=strip_fontsize)

    ax.set_title(title if title is not None else mat.what,
                 fontsize=hs.FS_LABEL)
    return im


def similarity_colorbar(fig, im, cax, *, label=r"similarity $\rho$",
                        ticks=(0.0, 0.5, 1.0)):
    """The one colourbar the three matrices share.

    Takes the image from whichever panel was drawn last — they are all on the
    same 0-1 scale, so which one is arbitrary, and that is the point.
    """
    cb = fig.colorbar(im, cax=cax, ticks=list(ticks))
    cb.set_label(label, fontsize=hs.FS_LABEL)
    cb.ax.tick_params(labelsize=hs.FS_TICK)
    cb.outline.set_visible(False)
    return cb


def shared_spectrum_limit(spectra, pad=1.05) -> float:
    """One y-limit for a row of spectra.

    Matt's 3E/F/G each autoscale, and so his three panels are on three
    different axes — 30, 60 and 100 units of the same quantity, drawn the same
    height. Pass this to `response_spectrum(limit=...)` to put them on one
    scale instead. Which is right depends on whether the figure is making a
    point about SHAPE (autoscale) or about MAGNITUDE (shared), so it is a call
    at the composition site, not a default buried in the panel.
    """
    return float(max(s.limit for s in spectra)) * pad


def mirror_spectrum(ax, up, down, classes, blocks, *, limit=None, width=0.8,
                    colors=None, strip=True, strip_labels=True, rotation=90,
                    label_up="low", label_down="high", annotate=True,
                    ylabel="Resp. Amplitude", linewidth=0.0):
    """Two series over the 59 odorants, mirrored about zero on one axis.

    The shape both spectrum panels are: `up` above the line, `down` below it,
    bars coloured by chemical class. Drawing `down` as a negative and then
    formatting the tick labels back to positive is what stops the lower half
    reading as a decrease.

    Fig 3E/F/G pass low and high (`response_spectrum`); Fig 6B/C/D pass
    experimental and predicted high (`prediction_spectrum`). They share this
    function so the two rows of the manuscript are the same picture with one
    series swapped, which is the comparison the figures are meant to support.
    """
    up, down = np.asarray(up, float), np.asarray(down, float)
    n = len(up)
    x = np.arange(n)
    if colors is None:
        colors = {c: cl.class_color(c) for c in blocks}
    bar_colors = [colors[c] for c in classes]
    kw = dict(width=width, linewidth=linewidth)
    ax.bar(x, up, color=bar_colors, **kw)
    ax.bar(x, -down, color=bar_colors, **kw)
    ax.axhline(0, color=hs.IDENTITY, lw=hs.BOUNDARY_LW, zorder=3)

    lim = float(limit) if limit else max(up.max(), down.max()) * 1.08
    ax.set_ylim(-lim, lim)
    ax.set_xlim(-0.5, n - 0.5)
    # the lower half is a magnitude, not a negative number
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, _pos: f"{abs(v):g}"))
    ax.set_ylabel(ylabel, fontsize=hs.FS_LABEL)
    ax.set_xticks([])
    ax.tick_params(labelsize=hs.FS_TICK)
    hs.despine(ax, keep=("left",))

    if annotate:
        ax.text(0.015, 0.97, label_up, transform=ax.transAxes,
                ha="left", va="top", fontsize=hs.FS_LABEL)
        ax.text(0.015, 0.03, label_down, transform=ax.transAxes,
                ha="left", va="bottom", fontsize=hs.FS_LABEL)
    if strip:
        labels = {c: cl.class_name(c) for c in blocks}
        hs.class_strip(ax, [labels[c] for c in classes],
                       colors={labels[c]: v for c, v in colors.items()},
                       orientation="x", labels=strip_labels, rotation=rotation)
    return ax


def response_spectrum(ax, spec, *, label_low="low", label_high="high", **kw):
    """Fig 3E/F/G: one glomerulus's spectrum, low up and high down.

    The two concentrations share one axis, mirrored about zero, which is what
    lets the eye read "the same odorant got bigger" as a single shape rather
    than as two charts to compare.

    Parameters
    ----------
    spec : Spectrum
        From `highconc.similarity.spectrum` / `eg_spectra`.
    limit : float, optional
        Symmetric y-limit. `None` autoscales this panel alone, which is what
        Matt's draft does; pass `shared_spectrum_limit(...)` to put a row of
        panels on one scale.
    colors : dict, optional
        Class letter -> colour. Defaults to the shared chemical-class palette,
        so a bar is the same colour as its class strip and as the same class in
        3B/3C/3D.

    Everything else is `mirror_spectrum`'s.
    """
    return mirror_spectrum(ax, spec.low, spec.high, spec.classes, spec.blocks,
                           label_up=label_low, label_down=label_high, **kw)


#: Y-axis label per prediction space, so a panel cannot claim dF/F while
#: drawing log1p. Keyed by `highconc.pls.PredictionSpectrum.space`.
SPECTRUM_YLABELS = {"raw": "Resp. Amplitude", "log": "log(1 + ΔF/F)",
                    "sqrt": r"$\sqrt{\Delta F/F}$"}


def prediction_spectrum(ax, spec, *, label_observed="experimental high",
                        label_predicted="predicted high", annotate_r=True,
                        r_fontsize=None, predicted_floor=None,
                        observed_floor=None, floor_color="0.4", **kw):
    """Fig 6A/B/C: experimental high up, held-out predicted high down.

    Deliberately the same picture as `response_spectrum` -- same odorant
    order, same class colours, same mirrored axis -- because 6A/B/C are read
    against 3E/F/G. Takes a `highconc.pls.PredictionSpectrum`.

    `annotate_r` prints the panel's own observed-vs-predicted correlation in
    the top right. It is a per-glomerulus number and belongs on the panel: the
    three glomeruli are predicted very differently well, and a row of three
    without it invites reading the best one as typical.

    `predicted_floor`, if given, draws a dashed horizontal line on the
    predicted (bottom) half at that value -- meant for
    `highconc.pls.predicted_high_floor`, the minimum predicted value over the
    *whole* dataset, so all three panels are read against the same floor
    rather than each other's own minimum. `None` (default) draws nothing.

    `observed_floor`, if given, draws a solid horizontal line at that value on
    the same (bottom) half -- meant for `highconc.pls.observed_high_floor`, the
    smallest nonzero response the experiment ever measured, so a predicted bar
    can be read against the lowest real response. It is drawn at the mirrored
    (negative) position, like `predicted_floor`.
    """
    kw.setdefault("ylabel", SPECTRUM_YLABELS.get(spec.space, "Resp. Amplitude"))

    ax = mirror_spectrum(ax, spec.observed, spec.predicted, spec.classes,
                         spec.blocks, label_up=label_observed,
                         label_down=label_predicted, **kw)
    if predicted_floor is not None:
        # Predicted is drawn as `-down` in mirror_spectrum, so its floor sits
        # at the mirrored (negative) position too.
        ax.axhline(-predicted_floor, color=floor_color, ls="--", lw=1.0,
                  zorder=2)
    if observed_floor is not None:
        ax.axhline(-observed_floor, color="#cdcdcd", lw=1.0, zorder=2)
    if annotate_r:
        ax.text(0.985, 0.97, f"r = {spec.r:.2f}", transform=ax.transAxes,
                ha="right", va="top", color="0.3",
                fontsize=hs.FS_LABEL if r_fontsize is None else r_fontsize)
    return ax
