"""Panels for the recruitment analysis (`highconc.recruitment`).

Twelve axes-level panels, grouped as the six composite figures the report is
built from:

===============================  ============================================
landscape                        `by_chem_class`, `by_cluster`,
                                 `propensity_hist`
topography                       `propensity_map`, `low_rate_map`
predictability                   `predictability`
model comparison                 `model_scores`
mechanism                        `gain_hist`, `recruited_low_rate_hist`,
                                 `sparseness_hist`
ceiling + robustness             `ceiling_bars`, `threshold_robustness`
===============================  ============================================

Each takes an existing `ax` and a result object from `highconc.recruitment`;
none creates a figure or loads anything.

Two orderings are keyword arguments rather than assumptions. `by_chem_class`
takes `order` (``"value"``, ``"table"``, ``"matt"``, or an explicit sequence),
and `by_cluster` lays clusters out in Joel's 1-7 order, which is the numbering
everything downstream of `highconc.clusters` uses -- the pre-refactor script
plotted them in raw component order under the retired A-G letters.
"""

from __future__ import annotations

import numpy as np

from .. import clusters as cl
from .. import style as hs
from ..recruitment import low_response_rate

# The recruitment palette. Fixed here rather than at each call site so the six
# composites agree: recruitment reads in the sequential ramp's mid-tone, and
# the categorical accents come from the style module's Okabe-Ito set.
C_RECRUIT = hs.SEQ(0.6)          # bars and histograms of recruitment itself
C_HIGHLIGHT = "#b30000"          # the quantity being argued about
C_SCALED = "#1b7837"             # scaled pairs / the linear backbone
C_CHEM = "#4575b4"               # chemistry features, and 'low' in comparisons
C_ANATOMY = "#91bfdb"            # anatomy features
C_NEUTRAL = "0.5"                # references and ceilings
C_SHUFFLED = "0.7"               # nulls are always grey (STYLE_GUIDE)

_MODEL_COLORS = (C_HIGHLIGHT, C_SCALED, C_CHEM)


# ---------------------------------------------------------------------------
# 1. Recruitment landscape
# ---------------------------------------------------------------------------
def by_chem_class(ax, r, *, order="value", names=True, color=C_RECRUIT):
    """Fraction of responding pairs that are recruited, per chemical class.

    `order` is ``"value"`` (descending, the pre-refactor default),
    ``"table"``, ``"matt"``, or an explicit sequence of class letters -- see
    `highconc.clusters.class_order` for the two named orderings and why both
    ship. `names=False` labels with the bare letters instead.
    """
    letters = list(r.class_frac)
    if order == "value":
        letters = [letters[i]
                   for i in np.argsort([r.class_frac[c] for c in letters])[::-1]]
    else:
        letters = cl.class_order(order, present=letters)
    vals = [r.class_frac[c] for c in letters]
    labels = [cl.class_name(c) for c in letters] if names else letters

    ax.bar(range(len(letters)), vals, color=color)
    ax.set_xticks(range(len(letters)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("fraction recruited")
    ax.set_ylim(0, 1)
    ax.set_title("recruitment by chemical class", fontsize=11, weight="bold")
    hs.despine(ax)
    return ax


def by_cluster(ax, r, *, color=C_RECRUIT, per_cluster_color=False):
    """Fraction recruited per glomerular cluster, in Joel's 1-7 order.

    `per_cluster_color=True` paints each bar in that cluster's identity colour
    (`highconc.style.cluster_color`) instead of the flat recruitment colour --
    use it when the panel sits beside a topography that is already colour-coded.
    """
    numbers = list(range(1, cl.N_CLUSTERS + 1))
    vals = [r.cluster_frac[n] for n in numbers]
    colors = ([hs.cluster_color(n) for n in numbers] if per_cluster_color
              else color)
    ax.bar(range(len(numbers)), vals, color=colors)
    ax.set_xticks(range(len(numbers)))
    ax.set_xticklabels([str(n) for n in numbers])
    ax.set_ylabel("fraction recruited")
    ax.set_ylim(0, 1)
    ax.set_xlabel("glomerular cluster")
    ax.set_title("recruitment by glom cluster", fontsize=11, weight="bold")
    hs.despine(ax)
    return ax


def propensity_hist(ax, r, *, bins=20, color=C_RECRUIT):
    """Distribution of per-glomerulus recruitment propensity, with the overall
    fraction marked. Glomeruli that never respond at high are excluded (their
    propensity is undefined, not zero)."""
    prop = r.per_glom_prop
    prop = prop[np.isfinite(prop)]
    ax.hist(prop, bins=bins, color=color, edgecolor="white")
    ax.axvline(r.frac_recruited, ls="--", color="0.4",
               label=f"overall = {r.frac_recruited:.2f}")
    ax.set_xlabel("per-glomerulus recruitment propensity")
    ax.set_ylabel("# glomeruli")
    ax.legend(fontsize=9)
    ax.set_title("per-glomerulus propensity", fontsize=11, weight="bold")
    hs.despine(ax)
    return ax


# ---------------------------------------------------------------------------
# 2. Topography
# ---------------------------------------------------------------------------
def propensity_map(ax, m, r, *, cmap="plasma", vmin=0.0, vmax=1.0, size=55,
                   missing=hs.GAINSBORO, title="recruitment propensity (all subjects)"):
    """Recruitment propensity over the bulb, all subjects superimposed.

    Glomeruli with no high-concentration response are drawn in `missing` grey:
    their propensity is undefined, and colouring them 0 would claim they are
    never recruited.
    """
    x, y = m.x, m.y
    ok = np.isfinite(r.per_glom_prop)
    ax.scatter(x[~ok], y[~ok], c=missing, s=45, linewidths=0)
    sc = ax.scatter(x[ok], y[ok], c=r.per_glom_prop[ok], cmap=cmap, s=size,
                    vmin=vmin, vmax=vmax, linewidths=0.3, edgecolors="0.4")
    _bulb_axes(ax, title)
    return sc


def low_rate_map(ax, m, *, cmap="viridis", size=55,
                 title="low-conc responsiveness"):
    """Each glomerulus's low-concentration response rate, same frame as
    `propensity_map` -- the contrast that shows recruitment filling the gaps."""
    sc = ax.scatter(m.x, m.y, c=low_response_rate(m), cmap=cmap, s=size,
                    linewidths=0.3, edgecolors="0.4")
    _bulb_axes(ax, title)
    return sc


def _bulb_axes(ax, title):
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    hs.despine(ax, keep=())
    if title:
        ax.set_title(title, fontsize=11, weight="bold")


# ---------------------------------------------------------------------------
# 3. Is recruitment predictable?
# ---------------------------------------------------------------------------
_AUC_KEYS = ("shuffled", "chemistry", "anatomy", "full")
_AUC_LABELS = ("shuffled", "chemistry\nonly", "anatomy\nonly", "full")
_AUC_COLORS = (C_SHUFFLED, C_CHEM, C_ANATOMY, C_HIGHLIGHT)


def predictability(ax, aucs, *, keys=_AUC_KEYS, ylim=(0.4, 1.0),
                   annotate=True):
    """CV AUC for predicting recruited vs scaled, feature set by feature set.

    The chance line at 0.5 is drawn because the interesting claim is distance
    from it, not the absolute height of the bars.
    """
    keys = list(keys)
    vals = [aucs[k] for k in keys]
    labels = [dict(zip(_AUC_KEYS, _AUC_LABELS)).get(k, k) for k in keys]
    colors = [dict(zip(_AUC_KEYS, _AUC_COLORS)).get(k, C_NEUTRAL) for k in keys]
    ax.bar(range(len(keys)), vals, color=colors)
    ax.axhline(0.5, ls="--", color="0.4", lw=1)
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("CV AUC  (recruited vs. scaled)")
    ax.set_ylim(*ylim)
    if annotate:
        for i, v in enumerate(vals):
            ax.text(i, v + 0.01, f"{v:.2f}", ha="center", fontsize=9)
    ax.set_title("Is recruitment predictable?\n(logistic, GroupKFold by odorant)",
                 fontsize=11, weight="bold")
    hs.despine(ax)
    return ax


# ---------------------------------------------------------------------------
# 4. Linear vs nonlinear
# ---------------------------------------------------------------------------
def model_scores(ax, mc, *, metric="q2", colors=_MODEL_COLORS, annotate=True):
    """One bar per model. `metric` is ``"q2"`` (all pairs) or ``"r"``
    (Pearson r over responding pairs) -- the two halves of the comparison
    figure, so the notebook can place them independently."""
    which = {"q2": 0, "r": 1}
    if metric not in which:
        raise ValueError(f"metric must be 'q2' or 'r', got {metric!r}")
    names = list(mc)
    vals = [mc[n][which[metric]] for n in names]
    ax.bar(range(len(names)), vals,
           color=[colors[i % len(colors)] for i in range(len(names))])
    if metric == "q2":
        ax.set_ylabel("$Q^2$ (all pairs)")
        ax.axhline(0, color="0.4", lw=0.8)
    else:
        ax.set_ylabel("Pearson r (responding pairs)")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8.5)
    if annotate:
        for i, v in enumerate(vals):
            ax.text(i, v + 0.005, f"{v:.2f}", ha="center", fontsize=8.5)
    hs.despine(ax)
    return ax


# ---------------------------------------------------------------------------
# 5. Mechanism
# ---------------------------------------------------------------------------
def gain_hist(ax, m, r, *, bins=40, clip=20, color=C_SCALED):
    """High/low gain among *scaled* pairs. Gains are clipped at `clip` for
    display only; the median line is computed on the unclipped values."""
    gains = m.high[r.scaled] / m.low[r.scaled]
    ax.hist(np.clip(gains, 0, clip), bins=bins, color=color, edgecolor="white")
    ax.axvline(1, ls="--", color="0.4", label="gain = 1")
    ax.axvline(r.gain_median, ls="-", color=C_HIGHLIGHT,
               label=f"median = {r.gain_median:.1f}")
    ax.set_xlabel("high / low gain (scaled pairs)")
    ax.set_ylabel("# pairs")
    ax.legend(fontsize=9)
    ax.set_title("supralinear scaling", fontsize=11, weight="bold")
    hs.despine(ax)
    return ax


def recruited_low_rate_hist(ax, m, r, *, bins=25, color=C_HIGHLIGHT):
    """For every recruited pair, how responsive that glomerulus is at low
    concentration to *other* odorants -- the panel that separates "weakly
    tuned, gated by concentration" from "dead glomerulus"."""
    rate = low_response_rate(m)
    _, jj = np.where(r.recruited)
    ax.hist(rate[jj], bins=bins, color=color, edgecolor="white")
    ax.axvline(np.median(rate[jj]), ls="-", color="0.2",
               label=f"median = {np.median(rate[jj]):.2f}")
    ax.set_xlabel("recruited glom's low-conc\nresponse rate to OTHER odorants")
    ax.set_ylabel("# recruited pairs")
    ax.legend(fontsize=9)
    ax.set_title("recruited gloms: weak, not silent, at low",
                 fontsize=11, weight="bold")
    hs.despine(ax)
    return ax


def sparseness_hist(ax, r, *, bins=20, low_color=C_CHEM,
                    high_color=C_HIGHLIGHT, alpha=0.6):
    """Population sparseness per odorant at low vs high -- tuning broadening."""
    ax.hist(r.pop_sparse_low, bins=bins, alpha=alpha, color=low_color,
            label=f"low (med {np.median(r.pop_sparse_low):.2f})")
    ax.hist(r.pop_sparse_high, bins=bins, alpha=alpha, color=high_color,
            label=f"high (med {np.median(r.pop_sparse_high):.2f})")
    ax.set_xlabel("fraction of glomeruli active (per odorant)")
    ax.set_ylabel("# odorants")
    ax.legend(fontsize=9)
    ax.set_title("tuning broadens with conc.", fontsize=11, weight="bold")
    hs.despine(ax)
    return ax


# ---------------------------------------------------------------------------
# 6. Noise ceiling + robustness
# ---------------------------------------------------------------------------
def ceiling_bars(ax, nc, *, annotate=True):
    """The model's cluster-level r against the inter-subject reliability proxy.
    Both are population-level; the gap, not the height, is the claim."""
    vals = [nc["ceiling"], nc["model_cluster_r"]]
    ax.bar([0, 1], vals, color=[C_NEUTRAL, C_HIGHLIGHT])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["inter-subject\nreliability\n(ceiling, proxy)",
                        "model\ncluster-level r"], fontsize=9)
    ax.set_ylabel("Pearson r (cluster level)")
    ax.set_ylim(0, 1)
    if annotate:
        for i, v in enumerate(vals):
            ax.text(i, v + 0.01, f"{v:.2f}", ha="center", fontsize=9)
    ax.set_title("noise ceiling (population level)", fontsize=11, weight="bold")
    hs.despine(ax)
    return ax


def threshold_robustness(ax, thresholds, fracs, *, color=C_HIGHLIGHT):
    """Recruited fraction as the responding threshold rises: the check that the
    headline number is not an artifact of calling any non-zero pixel a
    response."""
    ax.plot(thresholds, fracs, "o-", color=color)
    ax.set_xlabel("'responding' threshold (high dF/F >)")
    ax.set_ylabel("fraction recruited")
    ax.set_ylim(0, 1)
    ax.set_title("recruitment is threshold-robust", fontsize=11, weight="bold")
    hs.despine(ax)
    return ax


# ---------------------------------------------------------------------------
# 7. The recruitment matrix
# ---------------------------------------------------------------------------
def matrix(ax, rm, which="observed", *, cmap=None, vmax=None, title=None,
           xlabel=True, ylabel=True, tick_labels="class", boundaries=True):
    """One side of the recruitment matrix -- observed or held-out prediction.

    `rm` is a `highconc.recruitment.RecruitmentMatrix`. Rows are the
    low-concentration odorant (or class) a glomerulus is assigned to, columns
    the high-concentration odorant that recruits it, and the cell is a
    recruitment probability in [0, 1]. Undefined cells -- the structurally zero
    diagonal, and rows no glomerulus is assigned to -- render grey, the style
    module's missing-data colour, not as a probability of zero.

    Draw the pair on a shared scale. `vmax=None` (default) takes the 98th
    percentile over both matrices together -- shared by construction, and
    necessary because the cell distribution is heavily skewed: the median cell
    is 0.03 and a handful reach 1.0, so the full 0-1 scale leaves the whole
    panel in the palest two tenths of the ramp. Pass `vmax=1` for the absolute
    probability scale. Returns the image, which you need for the colorbar.

    `tick_labels` is ``"class"`` (one centred label per chemical-class block),
    ``"all"`` (every odorant, only legible at figure width), or ``None``.

    `boundaries` draws the chemical-class block rules. On by default.
    """
    if which not in ("observed", "predicted"):
        raise ValueError(f"unknown which {which!r}; expected "
                         "'observed' or 'predicted'")
    M = getattr(rm, which)
    if vmax is None:
        vmax = shared_vmax(rm)
    rb = rm.boundaries if boundaries else ()
    im = hs.heatmap(ax, M, vmin=0, vmax=vmax, cmap=cmap, aspect="equal",
                    row_boundaries=rb, col_boundaries=rb)
    if tick_labels == "class":
        pos = [p for p, _ in rm.class_ticks]
        names = [n for _, n in rm.class_ticks]
        ax.set_xticks(pos); ax.set_xticklabels(names, rotation=45, ha="right",
                                               fontsize=8)
        ax.set_yticks(pos); ax.set_yticklabels(names, fontsize=8)
    elif tick_labels == "all":
        ax.set_xticks(range(len(rm.labels)))
        ax.set_xticklabels(rm.labels, rotation=90, fontsize=5)
        ax.set_yticks(range(len(rm.labels)))
        ax.set_yticklabels(rm.labels, fontsize=5)
    else:
        ax.set_xticks([]); ax.set_yticks([])
    unit = "odorant" if rm.level == "odorant" else "class"
    if xlabel:
        ax.set_xlabel(f"recruited high-conc {unit}")
    if ylabel:
        ax.set_ylabel(f"driving low-conc {unit}")
    if title is None:
        title = ("observed" if which == "observed"
                 else "predicted (held-out PLS)")
    ax.set_title(title, fontsize=11, weight="bold")
    return im


def matrix_annotate(ax, rm, which="observed", *, vmax=None, fmt="{:.2f}",
                    fontsize=6):
    """Write each cell's value onto a `matrix` panel. Class level only.

    59 x 59 is far too dense to annotate; the 10 x 10 class-level matrix is the
    one worth reading off, which is why the older export printed its numbers.
    Pass the same `vmax` you gave `matrix`, so the light/dark text switch lands
    where the colour ramp actually turns dark.
    """
    M = getattr(rm, which)
    if vmax is None:
        vmax = shared_vmax(rm)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if not np.isfinite(v):
                continue
            ax.text(j, i, fmt.format(v), ha="center", va="center",
                    fontsize=fontsize,
                    color="white" if v > 0.6 * vmax else "black")
    return ax


def matrix_sweep_label(ax, rm, *, fontsize=9):
    """Label a sweep panel with the assignment fraction it was drawn at.

    The supplement is a row of matrices that differ only in `frac`, so the
    fraction, the empty-row count it produces, and the observed-vs-predicted
    correlation are the three things a reader needs beside each one.
    """
    ax.set_title(f"frac = {rm.frac:.2f}", fontsize=fontsize, weight="bold")
    ax.text(0.02, -0.14, f"{rm.n_empty} empty rows   r$_{{off}}$ = "
            f"{rm.r_offdiag:.2f}", transform=ax.transAxes, fontsize=fontsize - 1,
            va="top")
    return ax


def shared_vmax(*matrices, percentile: float = 98.0) -> float:
    """One colour ceiling for every panel that will be compared.

    Pass every `RecruitmentMatrix` that shares a colorbar -- the two halves of
    one figure, or all ten panels of the sweep -- and give the result to each
    `matrix` call. Comparing panels drawn on different scales is the mistake
    this exists to prevent (STYLE_GUIDE 4.3).
    """
    vals = []
    for rm in matrices:
        for M in (rm.observed, rm.predicted):
            vals.append(M[np.isfinite(M)])
    v = np.concatenate(vals)
    return float(np.percentile(v, percentile)) if v.size else 1.0
