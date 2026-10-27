"""highconc.style — house figure style for the highconcmss manuscript.

The single source of truth for the manuscript's look. Every figure imports it;
nothing defines its own rcParams. Depends only on numpy + matplotlib, so
collaborators without the repo can still copy this one file next to their
script and `import style as hs`.

    from highconc import style as hs
    hs.apply()

    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    ax.scatter(x, y, **hs.CLOUD)
    hs.despine(ax)
    hs.save(fig, "my_panel")          # -> my_panel.png / .pdf / .svg

Everything exported here is the reference implementation of STYLE_GUIDE.md.
If a value in the guide and a value here disagree, this file is right and the
guide needs updating.

Self-test:  python -m highconc.style        (renders a demo of every helper)
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

try:                                    # in-package: cluster identity has an owner
    from . import clusters as _clusters
    _LABELS = list(_clusters.LABELS)
except ImportError:                     # copied out as a standalone file
    _LABELS = ["1", "6", "3", "5", "2", "4", "7"]

__all__ = [
    "apply", "SEQ", "DIV", "CLUSTER_COLORS", "CLUSTER_LABELS",
    "CLUSTER_LABELS_LEGACY", "cluster_color", "GAINSBORO",
    "IDENTITY", "BOUNDARY", "BOUNDARY_LW", "C_PRIMARY", "C_OVERFIT", "C_WIN",
    "C_LOSE", "C_SHUFFLED", "C_NMF", "C_PCA", "CLOUD", "EMPHASIS",
    "C_MOST_SENSITIVE", "C_LOW", "C_HIGH", "CONDITIONS", "CONDITION_COLORS",
    "condition_color",
    "FS_TICK", "FS_LABEL", "FS_TITLE", "FS_PANEL", "WIDTH_1COL", "WIDTH_2COL",
    "despine", "panel_label", "identity_line", "save", "check_export",
    "heatmap", "topomap", "parity", "vector_matrix", "sequential_cmap",
    "ecdf", "box", "rank_align", "class_strip",
]

# ---------------------------------------------------------------------------
# Type scale (points) and canvas widths (inches), final on-the-page sizes
# ---------------------------------------------------------------------------
FS_TICK, FS_LABEL, FS_TITLE, FS_PANEL = 8, 9, 10, 11
WIDTH_1COL, WIDTH_2COL = 3.5, 7.0

# ---------------------------------------------------------------------------
# Palette
# ---------------------------------------------------------------------------
SEQ = plt.cm.magma                      # sequential: activation, magnitudes
DIV = LinearSegmentedColormap.from_list  # (built below)
DIV = DIV("teal_gainsboro_orrd",
          ["#01665e", "#5ab4ac", "gainsboro", "#fc8d59", "#b30000"])

# Cluster identity is owned by highconc.clusters, not by the style module.
# Labels are Joel's 1-7, in our component-index order; A-G is legacy only.
CLUSTER_LABELS = _LABELS
CLUSTER_LABELS_LEGACY = list("ABCDEFG")
CLUSTER_COLORS = plt.cm.tab10(np.arange(7))   # fixed per component across the paper


def cluster_color(number: int):
    """Colour for one of Joel's clusters, 1-7.

    The palette is indexed by *component*, not by display number, so a cluster
    keeps the colour it has always had in the exports — the renumbering changed
    what a cluster is called, not what colour it is drawn in.
    """
    return CLUSTER_COLORS[_clusters.to_index(number)]


GAINSBORO = "#dcdcdc"    # below floor / not measured / no data
IDENTITY = "0.45"        # identity and unity lines  (== #737373)
BOUNDARY = "black"       # cluster / block boundaries
BOUNDARY_LW = 0.8

C_PRIMARY = "#0072b2"    # Okabe-Ito blue      - primary metric
C_OVERFIT = "#d55e00"    # Okabe-Ito vermilion - overfit / canonical comparison
C_WIN = "#009e73"        # Okabe-Ito green     - "wins"
C_LOSE = "#888888"       #                     - "loses" / neutral
C_SHUFFLED = "#bbbbbb"   # nulls and shuffles are ALWAYS grey
C_NMF = "#9b5bd0"
C_PCA = "#777777"

# --- the concentration conditions -----------------------------------------
# The paper's primary categorical scheme. Okabe-Ito amber / blue / vermilion:
# colourblind-safe, and 'low'/'high' reuse the accent blue and vermilion above.
# Ordered weakest -> strongest, which is also the order they appear in legends.
C_MOST_SENSITIVE = "#e69f00"   # Okabe-Ito amber
C_LOW = C_PRIMARY              # 1x
C_HIGH = C_OVERFIT             # 25x 'low'

CONDITIONS = ("most-sensitive", "low", "high")
CONDITION_COLORS = {
    "most-sensitive": C_MOST_SENSITIVE,
    "low": C_LOW,
    "high": C_HIGH,
}
# Tolerated spellings, so a column that says "1x" or "most_sensitive" still works
_CONDITION_ALIASES = {
    "most sensitive": "most-sensitive", "most_sensitive": "most-sensitive",
    "mostsensitive": "most-sensitive", "ms": "most-sensitive",
    "0.01x": "most-sensitive",
    "1x": "low", "lo": "low",
    "25x": "high", "hi": "high",
}


def condition_color(name):
    """Colour for a concentration condition, tolerant of spelling.

    >>> condition_color("25x") == C_HIGH
    True
    """
    key = str(name).strip().strip("'\"").lower()
    key = _CONDITION_ALIASES.get(key, key)
    if key not in CONDITION_COLORS:
        raise KeyError(
            f"unknown condition {name!r}; expected one of {CONDITIONS} "
            f"(or an alias: {sorted(_CONDITION_ALIASES)})")
    return CONDITION_COLORS[key]

# Ready-made scatter kwargs for the two marker regimes (STYLE_GUIDE §4.1)
CLOUD = dict(s=6.5, linewidths=0, alpha=0.55, zorder=2)
EMPHASIS = dict(s=55, edgecolors="white", linewidths=1.0, zorder=3)


def apply() -> None:
    """Pin the export settings and type scale. Call once, before plotting."""
    # The Illustrator contract (STYLE_GUIDE §5)
    mpl.rcParams["svg.fonttype"] = "none"   # SVG text stays editable <text>
    mpl.rcParams["pdf.fonttype"] = 42       # TrueType, not Type 3
    mpl.rcParams["ps.fonttype"] = 42
    # Arial first: it exists on the Macs where figures are assembled, so
    # Illustrator does not flag imported text as a missing font.
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

    mpl.rcParams["font.size"] = FS_LABEL
    mpl.rcParams["axes.labelsize"] = FS_LABEL
    mpl.rcParams["axes.titlesize"] = FS_TITLE
    mpl.rcParams["xtick.labelsize"] = FS_TICK
    mpl.rcParams["ytick.labelsize"] = FS_TICK
    mpl.rcParams["legend.fontsize"] = FS_TICK
    mpl.rcParams["legend.frameon"] = False
    mpl.rcParams["axes.grid"] = False
    mpl.rcParams["figure.dpi"] = 120        # on-screen only; save() sets output dpi


def resolved_font() -> str:
    """Which sans-serif actually resolved. 'DejaVu Sans' means Arial is missing."""
    return mpl.font_manager.FontProperties(
        family=mpl.rcParams["font.sans-serif"]).get_name()


# ---------------------------------------------------------------------------
# Axis helpers
# ---------------------------------------------------------------------------
def despine(ax, keep=("left", "bottom")) -> None:
    """Top and right spines off. Use on data axes; leave matrices framed."""
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)


def panel_label(ax, letter: str, x: float = -0.10, y: float = 1.05) -> None:
    """Bold capital outside the axes, upper-left. 'A', not '(a)'."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=FS_PANEL,
            fontweight="bold", va="bottom", ha="right")


def identity_line(ax, lim: float, lo: float = 0.0) -> None:
    """Dashed grey unity line, drawn under the data."""
    ax.plot([lo, lim], [lo, lim], ls="--", lw=1.2, color=IDENTITY, zorder=1)


def sequential_cmap(under=GAINSBORO, bad=GAINSBORO):
    """magma with opaque grey for missing (and, if you set one, sub-floor) values.

    Zero is not missing. A zero response is a measurement and belongs at the
    bottom of the colour ramp (black); grey means "not measured". Keep
    vmin=0 on sparse response matrices, or most of the panel turns grey and the
    figure claims the data are absent when they are merely zero.

    A floor above zero (vmin > 0) is a separate, deliberate choice — it makes
    weak-but-nonzero values recede so strong ones read. It hides data, so say so
    in the caption.

    Opaque matters: a transparent raster becomes a PDF soft mask, which
    Illustrator renders wrong (STYLE_GUIDE §6.2).
    """
    cm = SEQ.copy()
    cm.set_under(under)
    cm.set_bad(bad)
    return cm


# ---------------------------------------------------------------------------
# Figure-type helpers
# ---------------------------------------------------------------------------
def parity(ax, observed, predicted, *, color=None, dense=None, label_stat=True,
           xlabel="observed", ylabel="predicted"):
    """Square parity panel with unity line and an in-panel r.

    Marker regime follows n, because the same size cannot serve 50 points and
    50,000 (STYLE_GUIDE §4.1). Override with dense=True/False.
    """
    observed = np.asarray(observed, float).ravel()
    predicted = np.asarray(predicted, float).ravel()
    ok = np.isfinite(observed) & np.isfinite(predicted)
    observed, predicted = observed[ok], predicted[ok]
    n = observed.size

    if dense is None:
        dense = n > 200
    if dense:
        # let overplotting carry density: tiny, edgeless, translucent
        kw = dict(s=1 if n > 2000 else 4, linewidths=0,
                  alpha=0.20 if n > 2000 else 0.40, zorder=2,
                  color=color or "0.25")
    elif n > 60:
        kw = dict(s=14, edgecolors="white", linewidths=0.5, zorder=3,
                  color=color or C_PRIMARY)
    else:
        kw = dict(EMPHASIS, color=color or C_PRIMARY)

    ax.scatter(observed, predicted, **kw)
    lim = float(np.nanpercentile(np.concatenate([observed, predicted]), 99.8)) * 1.03
    identity_line(ax, lim)
    ax.set(xlim=(0, lim), ylim=(0, lim), xlabel=xlabel, ylabel=ylabel)
    ax.set_aspect("equal")
    if label_stat and observed.size > 1:
        r = float(np.corrcoef(observed, predicted)[0, 1])
        ax.text(0.04, 0.95, f"r = {r:.2f}", transform=ax.transAxes,
                fontsize=FS_TICK, va="top")
    despine(ax)
    return ax


def heatmap(ax, M, *, vmin=None, vmax=None, diverging=False, row_boundaries=(),
            col_boundaries=(), cmap=None, aspect="auto"):
    """Matrix panel: nearest-neighbour, opaque missing cells, black boundaries.

    `diverging=True` pins the colour scale symmetrically about zero — required
    whenever the quantity is signed (STYLE_GUIDE §1.2).

    On a sparse non-negative matrix leave `vmin=0`: zeros then render as the
    palest end of the ramp (they are measurements), and grey is reserved for
    genuinely missing cells. See `sequential_cmap`.
    """
    M = np.asarray(M, float)
    if cmap is None:
        cmap = DIV.copy() if diverging else sequential_cmap()
        if diverging:
            cmap.set_bad(GAINSBORO)
    if diverging:
        v = float(np.nanmax(np.abs(M))) if vmax is None else vmax
        vmin, vmax = -v, v
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation="nearest", aspect=aspect)
    for b in row_boundaries:
        ax.axhline(b - 0.5, color=BOUNDARY, lw=BOUNDARY_LW)
    for b in col_boundaries:
        ax.axvline(b - 0.5, color=BOUNDARY, lw=BOUNDARY_LW)
    ax.tick_params(length=0)
    return im          # keep it: you need it for the colorbar


def topomap(ax, x, y, values, *, vmin, vmax, size=34, title=None):
    """Spatial map. vmin/vmax are REQUIRED and must be shared across every map
    that will be compared — that is the whole point (STYLE_GUIDE §4.3)."""
    sc = ax.scatter(x, y, c=values, s=size, cmap=sequential_cmap(),
                    vmin=vmin, vmax=vmax, linewidths=0.3, edgecolors="white")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    despine(ax, keep=())
    if title:
        ax.set_title(title, fontsize=FS_LABEL)
    return sc


def vector_matrix(ax, M, *, cmap=None, vmin=None, vmax=None, labels=None):
    """Draw a small matrix as one vector rectangle per cell.

    Use for anything up to roughly 20x20. A tiny raster can be dropped entirely
    on Illustrator import (STYLE_GUIDE §6.3); rectangles always survive and look
    identical in the PNG.
    """
    M = np.asarray(M, float)
    cmap = cmap or sequential_cmap()
    vmin = np.nanmin(M) if vmin is None else vmin
    vmax = np.nanmax(M) if vmax is None else vmax
    span = (vmax - vmin) or 1.0
    nr, nc = M.shape
    for r in range(nr):
        for c in range(nc):
            v = M[r, c]
            fc = GAINSBORO if not np.isfinite(v) else cmap((v - vmin) / span)
            ax.add_patch(plt.Rectangle((c - .5, r - .5), 1, 1, facecolor=fc,
                                       edgecolor="none", zorder=0))
    ax.set_xlim(-.5, nc - .5)
    ax.set_ylim(nr - .5, -.5)
    ax.set_aspect("equal")
    if labels is not None:
        ax.set_xticks(range(nc))
        ax.set_xticklabels(labels[:nc])
        ax.set_yticks(range(nr))
        ax.set_yticklabels(labels[:nr])
    ax.tick_params(length=0)


def ecdf(ax, series, *, xlabel=None, ylabel="cumul. prob.", legend=True,
         lw=1.6):
    """Cumulative distributions, one line per condition.

    `series` is {condition_name: values}. Conditions are drawn in the canonical
    weakest-to-strongest order and coloured from CONDITION_COLORS, so the same
    condition is the same colour in every figure of the paper.
    """
    names = [c for c in CONDITIONS if c in series] + \
            [c for c in series if c not in CONDITIONS]
    for name in names:
        v = np.asarray(series[name], float)
        v = np.sort(v[np.isfinite(v)])
        if v.size == 0:
            continue
        y = np.arange(1, v.size + 1) / v.size
        try:
            color = condition_color(name)
        except KeyError:
            color = C_SHUFFLED           # anything unrecognised reads as a null
        ax.step(v, y, where="post", lw=lw, color=color, label=str(name))
    ax.set_ylim(0, 1.02)
    ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(xlabel)
    if legend:
        ax.legend(frameon=False, loc="upper left")
    despine(ax)
    return ax


def box(ax, groups, *, log=False, ylabel=None, widths=0.6, showfliers=True):
    """Box plot, one box per condition, filled with the condition colour.

    `groups` is {condition_name: values}. Zero baseline is meaningless on a log
    axis, so `log=True` is the right choice for delivered concentrations.
    """
    names = [c for c in CONDITIONS if c in groups] + \
            [c for c in groups if c not in CONDITIONS]
    data = [np.asarray(groups[n], float)[np.isfinite(groups[n])] for n in names]
    bp = ax.boxplot(data, patch_artist=True, widths=widths,
                    showfliers=showfliers,
                    flierprops=dict(marker="D", markersize=2.5,
                                    markerfacecolor="black",
                                    markeredgecolor="none"),
                    medianprops=dict(color="black", lw=1.0),
                    whiskerprops=dict(color="black", lw=0.8),
                    capprops=dict(color="black", lw=0.8),
                    boxprops=dict(edgecolor="black", lw=0.8))
    for patch, name in zip(bp["boxes"], names):
        try:
            patch.set_facecolor(condition_color(name))
        except KeyError:
            patch.set_facecolor(C_SHUFFLED)
    ax.set_xticks(range(1, len(names) + 1))
    ax.set_xticklabels(names)
    if log:
        ax.set_yscale("log")
    if ylabel:
        ax.set_ylabel(ylabel)
    despine(ax)
    return bp


def rank_align(M, half_width=30, *, axis=1):
    """Rank-order each row and centre it on its peak, padding with NaN.

    Turns a glomerulus x odorant response matrix into the "ranked, aligned"
    view (each row sorted descending, peak in the middle, +/- half_width either
    side). NaN padding renders as opaque grey via sequential_cmap(), so a row
    with fewer odorants than the window does not masquerade as zero response.

    Returns an array of shape (n_rows, 2 * half_width + 1).
    """
    M = np.asarray(M, float)
    if axis == 0:
        M = M.T
    n, w = M.shape[0], 2 * half_width + 1
    out = np.full((n, w), np.nan)
    for i, row in enumerate(M):
        v = np.sort(row[np.isfinite(row)])[::-1]      # descending
        if v.size == 0:
            continue
        take = min(v.size, half_width + 1)
        out[i, half_width:half_width + take] = v[:take]   # peak at centre
        out[i, half_width - take + 1:half_width + 1] = v[:take][::-1]
    return out


def class_strip(ax, classes, *, colors=None, orientation="x", labels=True,
                rotation=90, fontsize=None, pad=0.012, height=0.018,
                side="right", offset=-0.5):
    """Coloured category strip along a matrix axis, with rotated labels.

    `classes` is a per-column (or per-row) sequence of category names, already
    in plotting order. Draws one coloured segment per contiguous run — which
    means the categories must be grouped, and if they are not you will see it
    immediately, which is the point.

    `side` only matters for `orientation="y"`: `"right"` (default, what
    `panels.similarity.response_matrix` uses) draws just past the axes'
    right edge; `"left"` mirrors the same rectangle-then-text layout onto the
    left edge instead, for a panel that already has something else (a
    colourbar, say) occupying the right side.

    `offset` is where cell 0 starts in data coordinates. The default -0.5 fits
    `imshow`/`hs.heatmap` (cell i centred on i); pass `offset=0` for
    `sns.heatmap` or `pcolormesh` (cell i spans i..i+1).

    Returns {category: colour} so the same mapping can be reused on sibling
    panels. Pass it back in as `colors` to keep panels consistent.
    """
    classes = list(classes)
    if colors is None:
        uniq = list(dict.fromkeys(classes))
        # Chemical families are NOT clusters — they are a different categorical
        # dimension and get their own palette, wide enough that no two families
        # share a colour. Reusing CLUSTER_COLORS here would both collide (7
        # colours, ~10 families) and imply a relationship that does not exist.
        src = plt.cm.tab10 if len(uniq) <= 10 else plt.cm.tab20
        n = 10 if len(uniq) <= 10 else 20
        if len(uniq) > n:
            raise ValueError(
                f"{len(uniq)} categories is more than the {n}-colour palette can "
                "distinguish; group them further or pass an explicit `colors` map")
        colors = {c: src(i) for i, c in enumerate(uniq)}
    fontsize = fontsize or FS_TICK - 1

    runs, start = [], 0
    for i in range(1, len(classes) + 1):
        if i == len(classes) or classes[i] != classes[start]:
            runs.append((start, i, classes[start]))
            start = i

    trans = ax.get_xaxis_transform() if orientation == "x" else \
        ax.get_yaxis_transform()
    for a, b, name in runs:
        if orientation == "x":
            ax.add_patch(plt.Rectangle((a + offset, -pad - height), b - a, height,
                                       transform=trans, clip_on=False,
                                       facecolor=colors[name], edgecolor="none"))
            if labels:
                ax.text((a + b) / 2 + offset, -pad - height - .012, name,
                        transform=trans, rotation=rotation, ha="center",
                        va="top", fontsize=fontsize, color=colors[name])
        elif side == "right":
            ax.add_patch(plt.Rectangle((1 + pad, a + offset), height, b - a,
                                       transform=trans, clip_on=False,
                                       facecolor=colors[name], edgecolor="none"))
            if labels:
                ax.text(1 + pad + height + .008, (a + b) / 2 + offset, name,
                        transform=trans, ha="left", va="center",
                        fontsize=fontsize, color=colors[name])
        else:
            ax.add_patch(plt.Rectangle((-pad - height, a + offset), height, b - a,
                                       transform=trans, clip_on=False,
                                       facecolor=colors[name], edgecolor="none"))
            if labels:
                ax.text(-pad - height - .008, (a + b) / 2 + offset, name,
                        transform=trans, ha="right", va="center",
                        fontsize=fontsize, color=colors[name])
    return colors


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
def save(fig, stem, outdir=".", dpi=300, formats=("png", "pdf", "svg"),
         verbose=True):
    """Write the export trio: white PNG to look at, transparent PDF + SVG with
    live text to assemble in Illustrator (STYLE_GUIDE §5)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for ext in formats:
        p = outdir / f"{stem}.{ext}"
        fig.savefig(p, dpi=dpi, bbox_inches="tight",
                    facecolor="white" if ext == "png" else "none",
                    transparent=(ext != "png"))
        written.append(p)
    if verbose:
        print(f"saved: {', '.join(p.name for p in written)}")
    return written


def check_export(svg_path, verbose=True) -> bool:
    """Confirm an exported SVG still carries live, editable text.

    This is the single most common handover failure: text exported as outlines
    or pixels cannot be corrected at assembly. Returns True if it passed.
    """
    svg = Path(svg_path).read_text()
    n_text = len(re.findall(r"<text[^>]*>", svg))
    fonts = sorted(set(re.findall(r"font-family:\s*([^;\"]+)", svg)))
    ok = n_text > 0
    if verbose:
        print(f"  <text> elements : {n_text}")
        print(f"  fonts referenced: {fonts or ['none declared']}")
        print("  PASS - text is live and editable" if ok else
              "  FAIL - text was outlined or rasterised; check your export settings")
        if ok and not any("Arial" in f for f in fonts):
            print("  NOTE - Arial is not in the font list; see STYLE_GUIDE §2")
    return ok


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------
def _selftest(outdir="."):
    apply()
    print(f"resolved sans-serif font: {resolved_font()}")
    rng = np.random.default_rng(0)

    n = 400
    xs, ys = rng.random(n), rng.random(n)
    obs = rng.gamma(2, 1, n)
    prd = 0.7 * obs + 0.3 * obs.mean() + rng.normal(0, .5, n)
    M = rng.normal(0, 1, (7, 18))

    fig = plt.figure(figsize=(WIDTH_2COL, 4.6), layout="constrained")
    gs = fig.add_gridspec(2, 3)

    ax = fig.add_subplot(gs[0, 0])
    parity(ax, obs, np.clip(prd, 0, None))
    panel_label(ax, "A")

    ax = fig.add_subplot(gs[0, 1])
    sc = topomap(ax, xs, ys, obs, vmin=0.4, vmax=obs.max(), size=18,
                 title="spatial map")
    fig.colorbar(sc, ax=ax, fraction=.05, shrink=.85)
    panel_label(ax, "B", x=-0.02)

    ax = fig.add_subplot(gs[0, 2])
    vector_matrix(ax, np.abs(M[:, :7]), labels=CLUSTER_LABELS)
    ax.set_title("vector matrix", fontsize=FS_LABEL)
    panel_label(ax, "C")

    ax = fig.add_subplot(gs[1, :])
    im = heatmap(ax, M, diverging=True, row_boundaries=[3])
    cb = fig.colorbar(im, ax=ax, fraction=.03, shrink=.9)
    cb.ax.yaxis.set_label_position("left")
    cb.set_label("weight", fontsize=FS_TICK)
    cb.ax.tick_params(labelsize=FS_TICK)
    ax.set_yticks(range(7))
    ax.set_yticklabels(CLUSTER_LABELS)
    ax.set_xticks([])
    ax.set_xlabel("component")
    panel_label(ax, "D", x=-0.045)

    paths = save(fig, "highconc_style_selftest", outdir=outdir)
    plt.close(fig)
    print("export check:")
    ok = check_export([p for p in paths if p.suffix == ".svg"][0])
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if _selftest(*sys.argv[1:2]) else 1)
