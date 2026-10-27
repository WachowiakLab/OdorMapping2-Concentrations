"""Panels for Figure 7 — the motif-transfer matrices.

Two panels over the same `highconc.transfer.Transfer`, and they must be read
together: the left says *whether* a transfer is there, the right says *how big*
it is. Neither alone is the result.

* `pvalue_matrix`      — permutation p per ordered pair, clipped at 0.1.
* `improvement_matrix` — the increment as a percentage of the target's own
  baseline, on the diverging ramp with zero at grey.

Both are drawn as **one vector rectangle per cell**, not `imshow`. A 7x7 raster
is a few hundred pixels that Illustrator is free to drop or resample on import
(STYLE_GUIDE 6.3); rectangles survive, and at this size there is no cost. The
same reason `highconc.style.vector_matrix` exists — these panels do their own
grid because they additionally need per-cell text and a significance box.

Orientation is fixed and stated on both axes, because it is the one thing a
reader can get backwards: **rows are the source, columns are the target.**
Row F, column T is "adding F's low-concentration mean to the prediction of T's
high". The matrix is not symmetric and is not meant to be.
"""

from __future__ import annotations

import matplotlib as mpl
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

from .. import style as hs

__all__ = ["pvalue_matrix", "improvement_matrix", "cell_grid",
           "significance_legend", "edge_table", "METRICS", "METRIC_LABELS"]

#: Which reading of the increment the right-hand panel shows.
#:
#: ``"pct_err"`` is the default and the manuscript's headline: the increment as
#: a fraction of the error the target's own mean left behind, ``dQ2 / (1 -
#: Q2_intra)``. It is bounded and comparable across targets.
#:
#: ``"pct_var"`` (``dQ2 / Q2_intra``) reads as "how much more of T's variance",
#: which is the more intuitive sentence but divides by a baseline that is 0.046
#: for cluster 5 and 0.438 for cluster 2 — so the same absolute increment looks
#: ten times larger on the weakly-predicted target. ``"dq2"`` shows the raw
#: increment and leaves the comparison to the reader.
METRICS = ("pct_err", "pct_var", "dq2")

METRIC_LABELS = {
    "pct_err": "% reduction in prediction error\n($\\Delta Q^2$ / (1 $-$ $Q^2_{intra}$))",
    "pct_var": "% increase in explained variance\n($\\Delta Q^2$ / $Q^2_{intra}$)",
    "dq2": "$\\Delta Q^2$ (held-out increment)",
}

#: Axis wording. Spelled out rather than "F" and "T" because the direction is
#: the result and an abbreviation is where a reader loses it.
XLABEL = "target group (high predicted)"
YLABEL = "source group (low added)"


def cell_grid(ax, t, *, facecolor, text=None, textcolor=None, alpha=0.05,
              box_lw=2.0, fontsize=None, labels=None, rotation=0,
              xlabel=XLABEL, ylabel=YLABEL, diagonal=hs.GAINSBORO,
              box_color=hs.BOUNDARY, text_significant_only=False):
    """Draw one K x K matrix of vector cells, with the significance boxes.

    `facecolor`, `text` and `textcolor` are callables taking ``(F, T)`` — the
    shared machinery behind both panels, so the two can never drift apart in
    cell size, tick handling or which pairs get a box.

    `box_color` is the colour of the significance box. `text_significant_only`
    prints a cell's number only where the pair is significant (BH q < `alpha`),
    leaving every other cell blank -- the panel then reads as the claim alone.

    The diagonal is drawn as an explicit grey square rather than left blank: a
    group predicting itself is not a missing measurement, it is the baseline
    the rest of the row is measured against, and an empty cell reads as data we
    failed to collect.
    """
    K = t.grouping.n_groups
    sig = t.significant(alpha)
    labels = list(t.labels) if labels is None else list(labels)
    fontsize = hs.FS_TICK - 1.5 if fontsize is None else fontsize

    for F in range(K):
        for T in range(K):
            if F == T:
                ax.add_patch(Rectangle((T - .5, F - .5), 1, 1,
                                       facecolor=diagonal, edgecolor="white",
                                       lw=0.6, zorder=0))
                continue
            ax.add_patch(Rectangle((T - .5, F - .5), 1, 1,
                                   facecolor=facecolor(F, T),
                                   edgecolor="white", lw=0.6, zorder=0))
            if text is not None and (sig[F, T] or not text_significant_only):
                ax.text(T, F, text(F, T), ha="center", va="center",
                        fontsize=fontsize, zorder=4,
                        color="black" if textcolor is None else textcolor(F, T))
            if sig[F, T]:
                ax.add_patch(Rectangle((T - .5, F - .5), 1, 1, fill=False,
                                       edgecolor=box_color, lw=box_lw,
                                       zorder=5))

    ax.set_xlim(-.5, K - .5)
    ax.set_ylim(K - .5, -.5)
    ax.set_aspect("equal")
    ax.set_xticks(range(K))
    ax.set_xticklabels(labels, rotation=rotation,
                       ha="right" if rotation else "center",
                       rotation_mode="anchor" if rotation else None)
    ax.set_yticks(range(K))
    ax.set_yticklabels(labels)
    ax.tick_params(length=0, labelsize=hs.FS_TICK)
    ax.set_xlabel(xlabel, fontsize=hs.FS_LABEL)
    ax.set_ylabel(ylabel, fontsize=hs.FS_LABEL)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


def _p_formatter(t):
    floor = t.p_floor

    def fmt(pv):
        if not np.isfinite(pv):
            return ""
        if pv <= floor + 1e-12:
            return f"<{floor:.0e}".replace("e-0", "e-")
        if pv < 0.001:
            return f"{pv:.0e}".replace("e-0", "e-")
        return f"{pv:.3f}"

    return fmt


def _pct_fmt(v):
    """Percentages to the nearest whole number, with a signed zero collapsed.

    ``f"{-0.4:+.0f}"`` is ``"-0"``, which reads as a negative quantity that is
    not there. A cell that rounds to zero should say ``0``."""
    if not np.isfinite(v):
        return ""
    s = f"{v:+.0f}"
    return "0" if s in ("+0", "-0") else s


def _dq2_fmt(v):
    return "" if not np.isfinite(v) else f"{v:+.3f}"


def pvalue_matrix(ax, t, *, alpha=0.05, clip=0.1, values=True, cbar=True,
                  cax=None, title=None, cmap=None, **kw):
    """Permutation p per ordered pair, dark = small.

    Colour is ``-log10 p`` on the sequential ramp, **clipped at `clip`**: every
    pair with p above 0.1 is the palest cell there is, because the panel's job
    is to separate the transfers that happened from the ones that did not, not
    to rank the 34 that did not. The printed number is the unclipped p, so
    nothing is hidden — and a p at the resolution floor prints as ``<1e-4``
    rather than as a number the 10,000 draws cannot support.

    `cmap` is a Colormap or its name (default `highconc.style.SEQ`).

    The box, not the colour, is the claim: it marks BH q < `alpha` across
    all ordered off-diagonal pairs.
    """
    nlp = -np.log10(np.where(np.isfinite(t.P), t.P, np.nan))
    lo, hi = -np.log10(clip), -np.log10(t.p_floor)
    norm = Normalize(lo, hi)
    cmap = hs.SEQ if cmap is None else cmap
    if isinstance(cmap, str):            # a name such as "OrRd", not just a Colormap
        cmap = mpl.colormaps[cmap]
    sig = t.significant(alpha)
    fmt = _p_formatter(t)

    cell_grid(ax, t,
              facecolor=lambda F, T: cmap(norm(np.clip(nlp[F, T], lo, hi))),
              text=(lambda F, T: fmt(t.P[F, T])) if values else None,
              textcolor=lambda F, T: "white" if sig[F, T] else "black",
              alpha=alpha, **kw)
    ax.set_title(title if title is not None
                 else f"permutation $p$  (N = {t.n_perm:,})", fontsize=hs.FS_TITLE)

    if cbar:
        ticks = np.arange(int(np.ceil(lo)), int(np.floor(hi)) + 1)
        cb = _colorbar(ax, cmap, norm, cax, ticks=ticks)
        cb.ax.set_yticklabels([f"{10.0 ** -k:g}" for k in ticks])
        cb.set_label(f"$p$ (clipped at {clip:g})", fontsize=hs.FS_LABEL)
        return cb
    return None


def improvement_matrix(ax, t, *, metric="pct_err", alpha=0.05, values=True,
                       vmax=None, cbar=True, cax=None, title=None,
                       fmt=None, **kw):
    """The increment itself, on the diverging ramp with zero at grey.

    Symmetric limits around zero by default (`vmax` overrides), so grey means
    "adding this source changed nothing" and the two directions of failure and
    success are equally visible. Most cells *are* grey, and that is the honest
    picture: 34 of 42 ordered pairs carry no transfer.

    Negative cells are real and not an error — adding an uninformative
    predictor costs held-out variance.
    """
    if metric not in METRICS:
        raise ValueError(f"metric must be one of {METRICS}, got {metric!r}")
    R = {"pct_err": t.pct_err, "pct_var": t.pct_var, "dq2": t.M}[metric]
    if fmt is None:
        fmt = _dq2_fmt if metric == "dq2" else _pct_fmt

    vmax = float(np.nanmax(np.abs(R[t.off]))) if vmax is None else float(vmax)
    norm = Normalize(-vmax, vmax)
    cmap = hs.DIV

    cell_grid(ax, t,
              facecolor=lambda F, T: cmap(norm(R[F, T])),
              text=(lambda F, T: fmt(R[F, T])) if values else None,
              textcolor=lambda F, T: ("white" if abs(R[F, T]) > 0.7 * vmax
                                      else "0.15"),
              alpha=alpha, **kw)
    ax.set_title(title if title is not None
                 else "improvement over intra-group baseline",
                 fontsize=hs.FS_TITLE)

    if cbar:
        cb = _colorbar(ax, cmap, norm, cax)
        cb.set_label(METRIC_LABELS[metric], fontsize=hs.FS_LABEL - 1)
        return cb
    return None


def _colorbar(ax, cmap, norm, cax=None, ticks=None):
    import matplotlib.pyplot as plt

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig = ax.get_figure()
    kw = dict(cax=cax) if cax is not None else dict(ax=ax, fraction=0.046,
                                                    pad=0.04)
    cb = fig.colorbar(sm, ticks=ticks, **kw)
    cb.ax.tick_params(labelsize=hs.FS_TICK)
    cb.outline.set_visible(False)
    return cb


def significance_legend(ax, t, *, alpha=0.05, loc="upper right",
                        box_color=hs.BOUNDARY):
    """One proxy handle explaining the significance box. Cheap, and it removes
    the only piece of the panel a reader has to be told about in prose."""
    handle = Rectangle((0, 0), 1, 1, fill=False, edgecolor=box_color, lw=2.0)
    leg = ax.legend([handle], [f"FDR $q$ < {alpha:g}"], loc=loc,
                    fontsize=hs.FS_TICK, frameon=False,
                    handlelength=1.1, handleheight=1.1, borderpad=0.2)
    return leg


def edge_table(t, *, alpha=0.05, metric="pct_err") -> str:
    """The significant edges as fixed-width text — for the notebook, not the
    figure. The matrix shows the pattern; this is what gets read into a
    caption or checked against the methods table."""
    key = {"pct_err": "pct_err", "pct_var": "pct_var", "dq2": "dq2"}[metric]
    fmt = _p_formatter(t)
    rows = [f"{'source':>12} -> {'target':<12} {'dQ2':>8} {metric:>8} "
            f"{'p':>8} {'q':>8}  way"]
    for e in t.edges(alpha):
        rows.append(f"{e['source']:>12} -> {e['target']:<12} "
                    f"{e['dq2']:+8.4f} {e[key]:+8.1f} {fmt(e['p']):>8} "
                    f"{e['q']:8.4f}  {'both' if e['reciprocal'] else 'one-way'}")
    return "\n".join(rows)
