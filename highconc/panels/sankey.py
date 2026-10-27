"""Panel for Figure 4E -- the flow from the original clusters to the folded-in clusters.

One function over `highconc.folding.FoldedFlows`. Drawn with matplotlib patches (a rectangle per
cluster, a Bezier ribbon per flow) so it is vector, editable, and inherits the house type scale
like every other panel; no plotting library is needed beyond matplotlib.

Cluster numbers sit at the sides: the original clusters 1-7 to the left of the left column, the
new clusters 1-7 to the right of the right column. A ribbon is coloured by where it came from, so
a cluster that carries over reads as one band of its own colour; glomeruli that were not in the
original clustering come in grey.
"""

from __future__ import annotations

import numpy as np
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path

from .. import style as hs

__all__ = ["flow_diagram"]


def _ribbon(x0, x1, y0a, y0b, y1a, y1b):
    """Closed Bezier band from the left span [y0a, y0b] at x0 to the right span [y1a, y1b] at x1."""
    xm = (x0 + x1) / 2
    verts = [(x0, y0a), (xm, y0a), (xm, y1a), (x1, y1a),
             (x1, y1b), (xm, y1b), (xm, y0b), (x0, y0b), (x0, y0a)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    return Path(verts, codes)


def flow_diagram(ax, flows, *, gap=0.022, node_width=0.035, ribbon_alpha=0.5,
                 unmatched_color=None, left_title="original clustering",
                 right_title="high folded in", unmatched_label="not in\noriginal",
                 number_fontsize=None, title_fontsize=None):
    """Draw the flow from the original clusters (left) to the folded-in clusters (right).

    Parameters
    ----------
    flows : highconc.folding.FoldedFlows
        ``counts[i, j]``: glomeruli from original cluster i+1 (last row: not in the original
        clustering) now in new cluster j+1.
    gap : float
        Vertical gap between nodes, as a fraction of the panel height.
    node_width : float
        Width of each node column, as a fraction of the panel width.
    unmatched_color : optional
        Colour of the not-in-original node and its ribbons. Defaults to `highconc.style.GAINSBORO`.
    left_title, right_title, unmatched_label : str
        Column headings and the label of the grey node.

    Returns
    -------
    dict
        ``{"left": [Rectangle, ...], "right": [Rectangle, ...]}`` -- the node patches.
    """
    counts = np.asarray(flows.counts)
    n_src, n_tgt = counts.shape                      # 8 sources (7 + not-in-original), 7 targets
    k = n_tgt
    unmatched_color = hs.GAINSBORO if unmatched_color is None else unmatched_color
    number_fontsize = hs.FS_TICK if number_fontsize is None else number_fontsize
    title_fontsize = hs.FS_LABEL if title_fontsize is None else title_fontsize

    def colour(i):                                   # source i: cluster i+1, or the grey row
        return unmatched_color if i == k else hs.cluster_color(i + 1)

    total = counts.sum()
    scale = (1.0 - gap * (n_src - 1)) / total        # one scale for both columns
    src_sizes, tgt_sizes = counts.sum(axis=1), counts.sum(axis=0)

    def tops(sizes, extra=0.0):
        t, y = [], extra
        for s in sizes:
            t.append(y)
            y += s * scale + gap
        return np.array(t)

    src_top = tops(src_sizes)
    tgt_top = tops(tgt_sizes, extra=gap * (n_src - n_tgt) / 2)   # centre the shorter column

    x_l1, x_r0 = node_width, 1.0 - node_width
    src_fill = src_top.copy()                        # running cursor inside each source node
    tgt_fill = tgt_top.copy()
    for i in range(n_src):                           # ribbons first, so the nodes sit on top
        for j in range(n_tgt):
            n = counts[i, j]
            if n == 0:
                continue
            h = n * scale
            path = _ribbon(x_l1, x_r0, src_fill[i], src_fill[i] + h, tgt_fill[j], tgt_fill[j] + h)
            ax.add_patch(PathPatch(path, facecolor=colour(i), edgecolor="none",
                                   alpha=ribbon_alpha, lw=0, zorder=1))
            src_fill[i] += h
            tgt_fill[j] += h

    left_nodes, right_nodes = [], []
    for i in range(n_src):
        r = Rectangle((0, src_top[i]), node_width, src_sizes[i] * scale,
                      facecolor=colour(i), edgecolor="none", zorder=2)
        ax.add_patch(r)
        left_nodes.append(r)
        mid = src_top[i] + src_sizes[i] * scale / 2
        if i == k:
            ax.text(-0.012, mid, unmatched_label, ha="right", va="center",
                    fontsize=number_fontsize - 1, color="0.35", linespacing=0.95)
        else:
            ax.text(-0.012, mid, str(i + 1), ha="right", va="center",
                    fontsize=number_fontsize, color=colour(i), weight="bold")
    for j in range(n_tgt):
        r = Rectangle((x_r0, tgt_top[j]), node_width, tgt_sizes[j] * scale,
                      facecolor=hs.cluster_color(j + 1), edgecolor="none", zorder=2)
        ax.add_patch(r)
        right_nodes.append(r)
        mid = tgt_top[j] + tgt_sizes[j] * scale / 2
        ax.text(1.0 + 0.012, mid, str(j + 1), ha="left", va="center",
                fontsize=number_fontsize, color=hs.cluster_color(j + 1), weight="bold")

    ax.text(node_width / 2, -0.035, left_title, ha="center", va="bottom", fontsize=title_fontsize)
    ax.text(1.0 - node_width / 2, -0.035, right_title, ha="center", va="bottom",
            fontsize=title_fontsize)
    ax.set_xlim(-0.12, 1.08)
    ax.set_ylim(1.02, -0.20)                         # cluster 1 at the top; room above for headings
    ax.set_axis_off()
    return {"left": left_nodes, "right": right_nodes}
