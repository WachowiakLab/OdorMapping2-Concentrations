"""Fig 5A -- tuning by cluster, as mirrored bars, drawn with plotly.

The one panel in this package that does not take a matplotlib ``ax``: plotly
figures are their own canvas. It follows the other rules -- it never loads
data, never saves, and takes finished numbers (`highconc.clusters.
fraction_responding`) -- and returns a ``plotly.graph_objects.Figure`` for the
notebook to export. Static export needs ``kaleido``.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from matplotlib.colors import to_hex
from plotly.subplots import make_subplots

__all__ = ["mirrored_tuning_figure"]

PX_PER_PT = 100 / 72      # the figure is laid out at 100 px per inch


def mirrored_tuning_figure(low, high, names, classes, *, class_names, class_colors,
                           cluster_colors, n_glom, width_in=7.0, height_in=5.5,
                           xlim=1.0, low_color="#0072b2", high_color="#d55e00",
                           class_gap=0.9, left_in=0.95, right_in=1.60,
                           top_in=0.62, bottom_in=0.40, font_pt=8,
                           name_font_pt=4.6, label_font_pt=9) -> go.Figure:
    """Seven columns of mirrored bars: low to the left of zero, high to the right.

    `low` and `high` are ``(n_odorants, n_clusters)`` fractions, rows already in
    display order (top to bottom); `names` and `classes` give each row's odorant
    name and chemical-class letter. `class_names` and `class_colors` map a class
    letter to its display name and colour; `cluster_colors` is one colour per
    cluster; `n_glom` the glomerulus count per cluster. Classes are separated by
    `class_gap` rows, with a rule across each cluster's frame between them. The margins are in inches, so the plot area
    keeps its size whatever the page is.
    """
    low, high = np.asarray(low, float), np.asarray(high, float)
    class_colors = {k: to_hex(v) for k, v in class_colors.items()}
    cluster_colors = [to_hex(c) for c in cluster_colors]
    low_color, high_color = to_hex(low_color), to_hex(high_color)
    classes = np.asarray(classes)
    n_odor, n_cl = low.shape
    px = lambda pt: pt * PX_PER_PT
    W, H = width_in * 100, height_in * 100

    shift = np.concatenate([[0], np.cumsum(classes[1:] != classes[:-1])]) * class_gap
    y = np.arange(n_odor) + shift
    blocks, start = [], 0
    for i in range(1, n_odor + 1):
        if i == n_odor or classes[i] != classes[start]:
            blocks.append((start, i - 1))
            start = i
    ylo, yhi = y[-1] + 0.7, y[0] - 0.7                     # reversed: first row on top

    plot_w = W - (left_in + right_in) * 100
    gap = 0.05 * 100
    fig = make_subplots(rows=1, cols=n_cl, shared_yaxes=True,
                        horizontal_spacing=gap / plot_w)
    for n in range(n_cl):
        fig.add_trace(go.Bar(y=y, x=-low[:, n], orientation="h", width=0.8,
                             marker=dict(color=low_color, line_width=0),
                             name="low", legendgroup="low", showlegend=False), row=1, col=n + 1)
        fig.add_trace(go.Bar(y=y, x=high[:, n], orientation="h", width=0.8,
                             marker=dict(color=high_color, line_width=0),
                             name="high", legendgroup="high", showlegend=False), row=1, col=n + 1)
        ax = "" if n == 0 else str(n + 1)
        fig.update_xaxes(range=[-xlim, xlim], tickvals=[-0.5, 0, 0.5],
                         ticktext=["0.5", "0", "0.5"], showline=False, zeroline=True,
                         zerolinecolor="#cccccc", zerolinewidth=1, showgrid=False,
                         ticks="outside", ticklen=3, tickfont=dict(size=px(font_pt - 1)),
                         row=1, col=n + 1)
        fig.update_yaxes(range=[ylo, yhi], showticklabels=False, showgrid=False,
                         zeroline=False, showline=False, ticks="", row=1, col=n + 1)
        fig.add_shape(type="rect", xref=f"x{ax}", yref=f"y{ax}", x0=-xlim, x1=xlim,
                      y0=y[0] - 0.6, y1=y[-1] + 0.6, line=dict(color="#595959", width=0.8),
                      fillcolor="rgba(0,0,0,0)", layer="above")
        for (_, b0), (a1, _) in zip(blocks[:-1], blocks[1:]):   # a rule between chemical classes
            ym = 0.5 * (y[b0] + y[a1])
            fig.add_shape(type="line", xref=f"x{ax}", yref=f"y{ax}", x0=-xlim, x1=xlim,
                          y0=ym, y1=ym, line=dict(color="#595959", width=0.8), layer="above")
        fig.add_annotation(xref=f"x{ax} domain", yref="paper", x=0.5, y=1.0,
                           text=str(n + 1), showarrow=False, yanchor="bottom",
                           font=dict(size=px(label_font_pt), color=cluster_colors[n]))
        fig.add_annotation(xref=f"x{ax} domain", yref=f"y{ax}", x=0.04, y=y[-1] + 0.5,
                           text=str(int(n_glom[n])), showarrow=False, xanchor="left",
                           yanchor="bottom", font=dict(size=px(font_pt - 1)))
    # low / high key above the first cluster, title and x label over the whole plot
    fig.add_annotation(xref="x domain", yref="paper", x=0.46, y=1.0 + 23 / (H - top_in * 100 - bottom_in * 100),
                       text="low", showarrow=False, xanchor="right", yanchor="bottom",
                       font=dict(size=px(font_pt), color=low_color))
    fig.add_annotation(xref="x domain", yref="paper", x=0.54, y=1.0 + 23 / (H - top_in * 100 - bottom_in * 100),
                       text="high", showarrow=False, xanchor="left", yanchor="bottom",
                       font=dict(size=px(font_pt), color=high_color))
    fig.add_annotation(xref="paper", yref="paper", x=0.5, y=1 + (top_in * 100 - 5) / (H - top_in * 100 - bottom_in * 100),
                       text="Tuning by cluster", showarrow=False, yanchor="top",
                       font=dict(size=px(label_font_pt)))
    fig.add_annotation(xref="paper", yref="paper", x=0.5, y=-(bottom_in * 100 - 4) / (H - top_in * 100 - bottom_in * 100),
                       text="Fraction responsive glomeruli per cluster", showarrow=False,
                       yanchor="bottom", font=dict(size=px(font_pt)))
    for a, b in blocks:                                    # class names left, odorant names right
        c = classes[a]
        fig.add_annotation(xref="paper", yref="y", x=-0.12 * 100 / plot_w,
                           y=0.5 * (y[a] + y[b]), text=class_names[c], showarrow=False,
                           xanchor="right", yanchor="middle",
                           font=dict(size=px(label_font_pt), color=class_colors[c]))
    for k, nm in enumerate(names):
        fig.add_annotation(xref="paper", yref=f"y{n_cl}", x=1.0 + 0.05 * 100 / plot_w, y=y[k],
                           text=nm, showarrow=False, xanchor="left", yanchor="middle",
                           font=dict(size=px(name_font_pt), color=class_colors[classes[k]]))
    fig.add_annotation(xref="paper", yref="paper", x=-(left_in * 100 - 4) / plot_w,
                       y=1 + (top_in * 100 - 4) / (H - top_in * 100 - bottom_in * 100),
                       text="<b>A</b>", showarrow=False, xanchor="left", yanchor="top",
                       font=dict(size=px(11), color="black"))
    fig.update_layout(
        width=W, height=H, barmode="overlay", bargap=0, plot_bgcolor="white",
        paper_bgcolor="white", showlegend=False, font=dict(family="Arial", size=px(font_pt), color="black"),
        margin=dict(l=left_in * 100, r=right_in * 100, t=top_in * 100, b=bottom_in * 100))
    return fig
