"""Panels for the Y = X · W schematic.

These draw `highconc.schematic.Cartoon` data -- **synthetic and illustrative,
never measured**. They are deliberately not in `highconc.panels.pls`, whose
panels all read a real `highconc.pls.PLSModel`: keeping the cartoon in its own
module means a figure that makes a claim cannot reach these by accident.

Each takes an `ax` and draws one frameless, tickless glyph. The equation's
layout -- where the maps, the strips and the matrix sit, and the arrows and
operators between them -- is composition and stays in the calling script.
"""

from __future__ import annotations

import numpy as np

from .. import style as hs

#: Edge colour for the cartoon glomeruli and cell borders.
EDGE = "0.35"

#: Cap the weight scale below the largest gain, so the diagonal reads uniformly
#: "strong" and the weak signed off-diagonal mixing terms stay visible. Without
#: it the diagonal takes the whole range and the cross terms vanish -- which is
#: the opposite of what the panel is for.
W_VMAX = 1.3

#: Colour range for the Fig 6F retuning cartoons. Their entries are built on a
#: 0-1 scale (`highconc.schematic.RETUNING_GAINS`), and all three must share
#: one range or "structured" would look stronger than "unstructured" merely by
#: autoscaling.
RETUNING_VMAX = 1.0


def cartoon_map(ax, px, py, values, *, vmax, cmap=None, size=430,
                edgecolors=EDGE, linewidths=0.8, pad=0.18, title=None,
                title_fontsize=11):
    """A cartoon glomerular map: one disc per glomerulus, coloured by activation.

    `vmax` is required and has no default on purpose -- the low- and
    high-concentration maps must share one scale, or the figure's whole claim
    (that Y is bigger than X) is drawn away by autoscaling.
    """
    ax.scatter(px, py, c=values, cmap=hs.SEQ if cmap is None else cmap,
               vmin=0, vmax=vmax, s=size, edgecolors=edgecolors,
               linewidths=linewidths, zorder=2)
    ax.set_xlim(-pad, 1 + pad)
    ax.set_ylim(-pad, 1 + pad)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(False)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    return ax


def vector_strip(ax, values, *, vmax, cmap=None, edgecolors="white",
                 linewidth=1.3):
    """The same activations as a 1 x n row of cells -- X or Y as a vector.

    Shares `vmax` with `cartoon_map`; the strip and the map above it are the
    same numbers drawn two ways, and must be read on one scale.
    """
    ax.pcolormesh(np.asarray(values, float)[None, :],
                  cmap=hs.SEQ if cmap is None else cmap, vmin=0, vmax=vmax,
                  edgecolors=edgecolors, linewidth=linewidth)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("auto")
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


def weight_glyph(ax, W, *, vmax=W_VMAX, cmap=None, edgecolors="white",
                 linewidth=1.2):
    """The cartoon W as a gridded matrix, signed, on the diverging map.

    Row-major with the origin at top left (`imshow` convention), so a cell's
    position matches how W is written in the equation. For the real 376x376 W
    use `highconc.panels.pls.weight_matrix` instead -- this one draws cell
    borders, which only make sense at cartoon size.
    """
    ax.pcolormesh(np.asarray(W, float), cmap=hs.DIV if cmap is None else cmap,
                  vmin=-vmax, vmax=vmax, edgecolors=edgecolors,
                  linewidth=linewidth)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return ax


def retuning_glyph(ax, cartoon, *, vmax=RETUNING_VMAX, cmap=None,
                   edgecolors="white", linewidth=1.0, frame="0.35",
                   frame_width=1.0, title=True, title_fontsize=9,
                   title_color="0.2"):
    """One Fig 6F cartoon: a 7 x 7 hypothesis about what W could look like.

    Takes a `highconc.schematic.Retuning`. Same glyph as `weight_glyph` — the
    same diverging map, the same cell borders — plus the thin frame and the
    caption that make it a labelled hypothesis rather than a matrix. `vmax` is
    shared across the three by default, which is the point: they differ in
    where the weight is, not in how much of it there is.

    Pass ``title=False`` to place the captions yourself, or a string to
    override this cartoon's own name.
    """
    weight_glyph(ax, cartoon.W, vmax=vmax, cmap=cmap, edgecolors=edgecolors,
                 linewidth=linewidth)
    if frame:
        for s in ax.spines.values():
            s.set_visible(True)
            s.set_color(frame)
            s.set_linewidth(frame_width)
    if title:
        ax.set_title(cartoon.name if title is True else title,
                     fontsize=title_fontsize, color=title_color, loc="left")
    return ax


def flow_field(ax, groups, *, colors=None, point_size=22, arrow_width=0.006,
               arrow_alpha=0.85, recruit_size=40, recruit_width=0.004,
               recruit_alpha=0.5, label_fontsize=11, xlim=(-5.2, 5.2),
               ylim=(-5.0, 5.2), xlabel="representational (tuning) space",
               label_fontsize_axis=11, origin=True, origin_size=70,
               origin_note="silent at low conc.\n(recruited gloms start here)"):
    """The concentration-flow cartoon: four groups each moving coherently.

    Takes `highconc.schematic.make_flow()`. Each group gets its filled points
    with solid step arrows, plus one OPEN marker at the origin with a lighter
    arrow -- a glomerulus that is silent at low concentration and is recruited
    into its cluster's motion. Open vs. filled is the distinction that carries
    the recruitment half of the story, so it survives greyscale.

    `colors` defaults to the tab10 entries the groups name. These are cartoon
    clusters lettered A-D, not the real 1-7 -- see `highconc.schematic`.
    """
    colors = ([hs.CLUSTER_COLORS[g.palette] for g in groups]
              if colors is None else list(colors))
    for g, color in zip(groups, colors):
        ax.quiver(g.points[:, 0], g.points[:, 1], g.steps[:, 0], g.steps[:, 1],
                  color=color, angles="xy", scale_units="xy", scale=1,
                  width=arrow_width, alpha=arrow_alpha, zorder=2)
        ax.scatter(g.points[:, 0], g.points[:, 1], s=point_size, color=color,
                   edgecolors="white", linewidths=0.6, zorder=3)
        ax.text(*g.label_at, g.name, color=color, fontsize=label_fontsize,
                weight="bold", ha="center")
        ax.scatter(*g.recruit_start, s=recruit_size, facecolors="none",
                   edgecolors=color, linewidths=1.3, zorder=3)
        ax.quiver(*g.recruit_start, *g.recruit_step, color=color, angles="xy",
                  scale_units="xy", scale=1, width=recruit_width,
                  alpha=recruit_alpha, linestyle="--", zorder=2)

    if origin:
        ax.scatter([0], [0], s=origin_size, marker="x", color="0.5", zorder=4)
        if origin_note:
            ax.text(0, -0.45, origin_note, fontsize=8.5, color="0.45",
                    ha="center", va="top")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=label_fontsize_axis)
    for s in ax.spines.values():
        s.set_visible(False)
    return ax
