"""Panels for the PLS-side figures, all reading one `highconc.pls.PLSModel`.

The claim these panels make together is that per-pair Q^2 understates the
model: evaluated at the population grain (`cluster_parity`, `rules_matrix`),
at the pattern grain (`pattern_correlation_hist`, `decoding_matrix`,
`rsa_scatter`, `pattern_norm_parity`) and on the scaled backbone
(`scaled_parity`, `decile_calibration`, `topk_retrieval`), the linear rule does
well. Each panel takes an `ax` and returns whatever the caller needs for a
colorbar.

Two conventions differ from the pre-refactor scripts, deliberately: clusters
are Joel's 1-7 rather than A-G in raw component order, and chemical-class rows
follow `highconc.clusters.class_order` (default ``"table"``, which is the
alphabetical order these figures already used).
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from .. import clusters as cl
from .. import pls as _pls
from .. import style as hs

#: In-panel statistic box. One definition, so every annotation matches.
NOTE_BOX = dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.7)

C_BLOCK = "#b30000"        # class x cluster block means
C_RETRIEVAL = "#1b7837"
C_RSA = "#4575b4"
C_CHANCE = "0.6"

#: Axis label for every panel in log(1+dF/F) space.
POINT_LABEL = "log(1 + ΔF/F)"

#: Example-map scale, log(1+dF/F). Responders sit around 3.5-4.5, so a low
#: ceiling saturates them and lets the shrunken predictions read; values below
#: the floor go flat grey so near-silent glomeruli recede to context.
MAP_VMIN, MAP_VMAX = 0.5, 3.5


def note(ax, text, x=0.04, y=0.96, va="top", ha="left", fontsize=10, **kw):
    """Statistic box in axes coordinates -- the r, rho or slope being claimed."""
    return ax.text(x, y, text, transform=ax.transAxes, va=va, ha=ha,
                   fontsize=fontsize, bbox=NOTE_BOX, **kw)


# ---------------------------------------------------------------------------
# Population grain
# ---------------------------------------------------------------------------
def cluster_parity(ax, obs, pred, *, size=26, alpha=0.75, legend=True,
                   edgecolors="none", linewidths=0, annotate=True,
                   note_pos=(0.04, 0.96), title="cluster × odorant parity",
                   title_fontsize=12):
    """Observed vs predicted cluster means, one colour per cluster.

    `obs`/`pred` are the ``(odorants, 7)`` arrays from
    `highconc.pls.cluster_by_odorant`, columns in Joel's 1-7 order; colours come
    from `highconc.style.cluster_color`, which indexes by component, so a
    cluster keeps the colour it has in every other figure.
    """
    lim = max(obs.max(), pred.max()) * 1.05
    hs.identity_line(ax, lim)
    for j in range(obs.shape[1]):
        number = j + 1
        ax.scatter(obs[:, j], pred[:, j], s=size, color=hs.cluster_color(number),
                   alpha=alpha, linewidths=linewidths, edgecolors=edgecolors,
                   label=str(number))
    r = _pls.pearson(obs.ravel(), pred.ravel())
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal", "box")
    ax.set_xlabel("observed cluster mean  log(1+ΔF/F)")
    ax.set_ylabel("predicted cluster mean")
    if legend:
        ax.legend(ncol=2, fontsize=8, frameon=False, title="cluster")
    if annotate:
        note(ax, f"r = {r:.2f}\n({obs.shape[1]} clusters × {obs.shape[0]} odorants)",
             *note_pos)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    hs.despine(ax)
    return r


def block_parity(ax, obs, pred, *, size=28, color=C_BLOCK, alpha=0.6,
                 edgecolors="none", linewidths=0, annotate=True,
                 title="class × cluster rules parity", title_fontsize=12):
    """Observed vs predicted class x cluster block means -- does the model
    reproduce the chemical "rules", not just the totals?"""
    lim = max(obs.max(), pred.max()) * 1.05
    hs.identity_line(ax, lim)
    ax.scatter(obs.ravel(), pred.ravel(), s=size, color=color, alpha=alpha,
               linewidths=linewidths, edgecolors=edgecolors)
    r = _pls.pearson(obs.ravel(), pred.ravel())
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal", "box")
    ax.set_xlabel("observed block mean  log(1+ΔF/F)")
    ax.set_ylabel("predicted block mean")
    if annotate:
        note(ax, f"r = {r:.2f}\n({obs.shape[0]} classes × {obs.shape[1]} clusters)")
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    hs.despine(ax)
    return r


def rules_matrix(ax, M, classes, *, vmin=0.0, vmax=1.0, names=True,
                 xticks=True, title=None, title_fontsize=12, label_fontsize=8):
    """One class x cluster rules matrix. Returns the image, for the colorbar.

    Observed and predicted must be drawn with the same `vmin`/`vmax` -- the
    panel's whole point is that the two look alike.
    """
    im = ax.imshow(M, cmap=hs.SEQ, vmin=vmin, vmax=vmax, aspect="auto")
    labels = [cl.class_name(c) for c in classes] if names else list(classes)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=label_fontsize)
    if xticks:
        ax.set_xticks(range(M.shape[1]))
        ax.set_xticklabels([str(k + 1) for k in range(M.shape[1])],
                           fontsize=label_fontsize)
    else:
        ax.set_xticks([])
    ax.tick_params(length=0)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    return im


# ---------------------------------------------------------------------------
# Pattern grain
# ---------------------------------------------------------------------------
def pattern_correlation_hist(ax, matched, mismatched, *, bins_matched=20,
                             bins_mismatched=40, color=C_BLOCK,
                             mismatch_color="0.7", alpha=0.8,
                             title="right pattern (per odorant)",
                             title_fontsize=12):
    """Matched vs mismatched per-odorant map correlations.

    The mismatched distribution is the null the matched one has to beat, and it
    is drawn from the same matrix, so no extra assumption enters.
    """
    ax.hist(mismatched, bins=bins_mismatched, color=mismatch_color, alpha=alpha,
            density=True, label=f"mismatched (med {np.median(mismatched):.2f})")
    ax.hist(matched, bins=bins_matched, color=color, alpha=alpha, density=True,
            label=f"matched (med {np.median(matched):.2f})")
    ax.axvline(np.median(matched), color=color, ls="-", lw=1.2)
    ax.set_xlabel("per-odorant map correlation (pred vs obs)")
    ax.set_ylabel("density")
    ax.legend(fontsize=8.5, frameon=False)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    hs.despine(ax)
    return ax


def decoding_matrix(ax, m, dec, *, cmap="magma", order_by="chem",
                    annotate=True, title="odorant decoding (rank)",
                    title_fontsize=12):
    """Predicted-vs-observed identity-map correlations, odorants sorted so the
    chemical block structure is visible. Returns the image, for the colorbar."""
    S = dec["S"]
    order = (np.argsort(m.chem, kind="stable") if order_by == "chem"
             else np.arange(S.shape[0]))
    im = ax.imshow(S[np.ix_(order, order)], cmap=cmap, aspect="equal",
                   vmin=np.percentile(S, 2), vmax=np.percentile(S, 99))
    ax.set_xlabel("observed map  (odorants, by class)")
    ax.set_ylabel("predicted map")
    ax.set_xticks([])
    ax.set_yticks([])
    if annotate:
        n = dec["n"]
        ax.text(0.03, 0.97,
                f"true odorant median rank {int(dec['median_rank'])}/{n} (top 10%)\n"
                f"top-5 = {dec['top5']*100:.0f}%  top-10 = {dec['top10']*100:.0f}%\n"
                f"chance top-5 = {5/n*100:.0f}%, top-10 = {10/n*100:.0f}%",
                transform=ax.transAxes, va="top", fontsize=8.5, color="white",
                bbox=dict(boxstyle="round,pad=0.3", fc="black", ec="none",
                          alpha=0.55))
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    return im


def rsa_scatter(ax, r, *, size=4, color=C_RSA, alpha=0.15,
                title="odor-space geometry (RSA)", title_fontsize=12):
    """Observed vs predicted odorant dissimilarity, one point per odorant pair."""
    ax.scatter(r["obs"], r["pred"], s=size, color=color, alpha=alpha,
               linewidths=0)
    ax.set_xlabel("observed dissimilarity  (1 − r)")
    ax.set_ylabel("predicted dissimilarity")
    note(ax, f"RSA  ρ = {r['rho']:.2f}\n(odor × odor geometry)")
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    hs.despine(ax)
    return ax


def pattern_norm_parity(ax, m, *, gridsize=45, percentile=99.5,
                        title="relative pattern, not magnitude",
                        title_fontsize=12):
    """Parity after z-scoring within each odorant: pattern with magnitude
    removed, which is the grain the model is actually good at."""
    zo, zp = _pls.zscore_rows(m.Ylog).ravel(), _pls.zscore_rows(m.Yhat).ravel()
    lim = np.percentile(np.abs(np.concatenate([zo, zp])), percentile)
    ax.plot([-lim, lim], [-lim, lim], ls="--", lw=1.2, color=hs.IDENTITY)
    ax.hexbin(zo, zp, gridsize=gridsize, cmap=hs.SEQ, mincnt=1,
              norm=mpl.colors.LogNorm(), extent=(-lim, lim, -lim, lim))
    r = _pls.pearson(zo, zp)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", "box")
    ax.set_xlabel("observed (z within odorant)")
    ax.set_ylabel("predicted (z within odorant)")
    note(ax, f"r = {r:.2f}\n(pattern-normalized)")
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    hs.despine(ax)
    return r


# ---------------------------------------------------------------------------
# Calibration, retrieval, and the raw parity
# ---------------------------------------------------------------------------
def scaled_parity(ax, m, *, gridsize=40, annotate=True,
                  title="calibrated where it can scale", title_fontsize=12):
    """Parity on the scaled backbone only (low > 0), where the linear rule has
    an input to act on. Returns ``(slope, r)``."""
    obs, pred = _pls.scaled_pairs(m)
    lim = max(obs.max(), pred.max()) * 1.03
    hs.identity_line(ax, lim)
    ax.hexbin(obs, pred, gridsize=gridsize, cmap=hs.SEQ, mincnt=1,
              norm=mpl.colors.LogNorm(), extent=(0, lim, 0, lim))
    slope = float(np.polyfit(obs, pred, 1)[0])
    r = _pls.pearson(obs, pred)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal", "box")
    ax.set_xlabel("observed  log(1+ΔF/F)")
    ax.set_ylabel("predicted")
    if annotate:
        note(ax, f"scaled pairs (low>0)\nslope = {slope:.2f}  r = {r:.2f}")
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    hs.despine(ax)
    return slope, r


def square_parity_axes(ax, lim, *, label=POINT_LABEL, label_fontsize=11,
                       tick_fontsize=9):
    """Square parity frame: unity line, equal aspect, both axes on [0, lim]."""
    hs.identity_line(ax, lim)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal", "box")
    ax.set_xlabel(f"observed  {label}", fontsize=label_fontsize)
    ax.set_ylabel(f"predicted  {label}", fontsize=label_fontsize)
    ax.tick_params(labelsize=tick_fontsize)
    hs.despine(ax)
    return ax


def responding_parity(ax, m, *, gridsize=34, extent=(1.0, 5.0), saturation=3.0,
                      linewidths=None, frame="pinned", stat="q2", annotate=True,
                      cbar=False, title="responding pairs", title_fontsize=11,
                      note_fontsize=9, cmap=None, floor=None, floor_color="0.35"):
    """Per-pair parity over responding pairs (observed > 0). Returns the hexbin.

    Two framings, because the two figures that use this panel ask different
    questions of it:

    ``frame="pinned"``  axes fixed to `extent` and the colour ceiling raised
        above the densest bin by `saturation`, so the cloud's shape reads
        instead of one saturated blob. Pairs outside the window are drawn but
        not scaled to.
    ``frame="square"``  a full-range square parity frame with the unity line
        and axis labels (`square_parity_axes`), scaled to the data.

    `stat` selects the annotation: ``"q2"`` (Q² and n) or ``"r+q2"`` (n, the
    correlation over responding pairs, and Q²). `linewidths=None` leaves the
    hexagon edges at matplotlib's default -- do not "tidy" that to 0, it
    changes how much ink is on the page.

    `cmap` is the hexbin colormap (default `highconc.style.SEQ`). `floor`, if
    given, draws a dotted horizontal line at that predicted value -- meant for
    `highconc.pls.observed_high_floor`, the smallest nonzero response ever
    measured, so a prediction under the line is one below anything observed.
    """
    obs, pred = _pls.responding_pairs(m)
    if frame == "square":
        lim = max(obs.max(), pred.max()) * 1.03
        extent = (0.0, lim)
        square_parity_axes(ax, lim)
    elif frame != "pinned":
        raise ValueError(f"frame must be 'pinned' or 'square', got {frame!r}")

    lo, hi = extent
    kw = {} if linewidths is None else {"linewidths": linewidths}
    hb = ax.hexbin(obs, pred, gridsize=gridsize,
                   cmap=hs.SEQ if cmap is None else cmap, mincnt=1,
                   norm=mpl.colors.LogNorm(), extent=(lo, hi, lo, hi),
                   zorder=2, **kw)
    if saturation:
        hb.set_clim(1, hb.get_array().max() * saturation)

    if frame == "pinned":
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ticks = list(range(int(lo), int(hi) + 1))
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xlabel(f"observed  {POINT_LABEL}")
        ax.set_ylabel(f"predicted  {POINT_LABEL}")
        hs.despine(ax)

    if floor is not None:
        ax.axhline(floor, color=floor_color, ls=":", lw=1.0, zorder=3)

    if annotate:
        q2 = _pls.baseline_q2(m)
        if stat == "r+q2":
            text = (f"responding pairs (n={obs.size:,})\n"
                    f"r = {_pls.pearson(obs, pred):.2f}"
                    f"   $Q^2_{{all}}$ = {q2:.2f}")
        else:
            text = f"Q² = {q2:.2f}  (all pairs)\nresponding (n={obs.size:,})"
        note(ax, text, y=0.96 if frame == "square" else 0.97,
             fontsize=note_fontsize)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    if cbar:
        cb = ax.figure.colorbar(hb, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("glom × odorant pairs", fontsize=9)
        cb.ax.tick_params(labelsize=8)
    return hb


def mean_parity(ax, obs, pred, *, size=42, color=C_BLOCK, alpha=0.55,
                unit="odorant", pad=1.08, title=None, title_fontsize=11.5,
                note_fontsize=9):
    """Parity after averaging one axis away -- one point per odorant, or per
    glomerulus.

    Note which is which when you read the r: averaging over odorants leaves
    each glomerulus's baseline excitability, which mean-centring reproduces
    almost by construction. The per-glomerulus panel is the easy axis and is
    included for completeness, not as evidence.
    """
    lim = max(obs.max(), pred.max()) * pad
    square_parity_axes(ax, lim)
    ax.scatter(obs, pred, s=size, c=color, alpha=alpha, linewidths=0, zorder=2)
    note(ax, f"per {unit} (n={obs.size})\nr = {_pls.pearson(obs, pred):.2f}",
         fontsize=note_fontsize)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    return _pls.pearson(obs, pred)


def residual_hist(ax, m, *, bins=70, percentile=99.5, color=C_BLOCK,
                  alpha=0.85, log=True, title="residuals", title_fontsize=11.5,
                  note_fontsize=9):
    """Observed minus predicted over all pairs, on a log count axis.

    The log axis is doing real work: the distribution is dominated by a huge
    near-zero peak, and on a linear axis the recruited tail -- the thing the
    figure is about -- is invisible.
    """
    resid = (m.Ylog - m.Yhat).ravel()
    lim = np.percentile(np.abs(resid), percentile)
    ax.hist(resid, bins=bins, range=(-lim, lim), color=color, alpha=alpha,
            edgecolor="white", linewidth=0.2)
    ax.axvline(0, ls="--", lw=1.2, color=hs.IDENTITY)
    if log:
        ax.set_yscale("log")
    ax.set_box_aspect(1)
    ax.set_xlabel(f"residual  (obs − pred)  {POINT_LABEL}", fontsize=11)
    ax.set_ylabel("glom × odorant pairs", fontsize=11)
    ax.tick_params(labelsize=9)
    hs.despine(ax)
    note(ax, f"all pairs (n={resid.size:,})\nmean = {resid.mean():+.2f}\n"
             f"median = {np.median(resid):+.2f}", fontsize=note_fontsize)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    return resid


def parity_by_low_input(ax, m, *, cmap="plasma", size=14, alpha=0.6,
                        cbar=False, title="colored by low-conc input",
                        title_fontsize=11.5, note_fontsize=9):
    """Responding pairs coloured by each pair's own low-concentration input.

    This is the panel that names the failure mode: the cluster along the x-axis
    (large observed, ~0 predicted) is uniformly dark, i.e. low input ~ 0. A
    linear map cannot scale a near-zero input up. Points are drawn in order of
    input so the high-input ones sit on top.
    """
    obs, pred = _pls.responding_pairs(m)
    low = m.Xlog.ravel()[m.Ylog.ravel() > 0]
    lim = max(obs.max(), pred.max()) * 1.03
    square_parity_axes(ax, lim)
    order = np.argsort(low)
    sc = ax.scatter(obs[order], pred[order], c=low[order], cmap=cmap, s=size,
                    alpha=alpha, linewidths=0, zorder=2)
    note(ax, "x-axis cluster = low input ≈ 0\n(recruited at high conc.)",
         fontsize=note_fontsize)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    if cbar:
        cb = ax.figure.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(f"low-conc input  {POINT_LABEL}", fontsize=9)
        cb.ax.tick_params(labelsize=8)
    return sc


def scaled_vs_recruited(ax, m, r, *, size=12, scaled_color=C_RETRIEVAL,
                        recruited_color=C_BLOCK, scaled_alpha=0.35,
                        recruited_alpha=0.30, fit=True, legend=True,
                        title="scaled vs. recruited", title_fontsize=11.5,
                        note_fontsize=9):
    """The same responding pairs, split by whether the linear rule had an input.

    `r` is a `highconc.recruitment.Recruitment` -- the masks come from there
    rather than being recomputed, so this panel and the recruitment figure
    cannot drift apart on what "recruited" means. The two fitted slopes are the
    quantitative claim: the scaled backbone is near-calibrated, the recruited
    population collapses toward the x-axis.
    """
    responding = m.Ylog.ravel() > 0
    obs, pred = _pls.responding_pairs(m)
    scaled = r.scaled.ravel()[responding]
    recruited = r.recruited.ravel()[responding]
    lim = max(obs.max(), pred.max()) * 1.03
    square_parity_axes(ax, lim)
    ax.scatter(obs[recruited], pred[recruited], s=size, c=recruited_color,
               alpha=recruited_alpha, linewidths=0, zorder=2)
    ax.scatter(obs[scaled], pred[scaled], s=size, c=scaled_color,
               alpha=scaled_alpha, linewidths=0, zorder=3)

    slopes = {}
    if fit:
        for name, mask, color in (("scaled", scaled, scaled_color),
                                  ("recruited", recruited, recruited_color)):
            slope, intercept = np.polyfit(obs[mask], pred[mask], 1)
            xs = np.array([0, lim])
            ax.plot(xs, slope * xs + intercept, color=color, lw=2.0, zorder=4)
            slopes[name] = float(slope)
        note(ax, f"scaled (low>0):  slope = {slopes['scaled']:.2f}\n"
                 f"recruited (low=0):  slope = {slopes['recruited']:.2f}\n"
                 f"{recruited.mean() * 100:.0f}% of responders recruited",
             fontsize=note_fontsize)
    if legend:
        ax.text(0.97, 0.11, "scaled", color=scaled_color, fontsize=11,
                weight="bold", ha="right", va="bottom", transform=ax.transAxes)
        ax.text(0.97, 0.03, "recruited", color=recruited_color, fontsize=11,
                weight="bold", ha="right", va="bottom", transform=ax.transAxes)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    return slopes


def decile_calibration(ax, m, *, n_bins=10, color=C_BLOCK, annotate=True,
                       title="monotone calibration (lift)", title_fontsize=12):
    """Mean observed response per predicted-value bin: is the model monotone
    even where it is not calibrated? Returns Spearman rho."""
    mean, sem, rho = _pls.calibration_deciles(m, n_bins=n_bins)
    ax.errorbar(range(1, n_bins + 1), mean, yerr=sem, fmt="o-", color=color,
                capsize=3)
    ax.set_xlabel(f"predicted-value decile (1=low → {n_bins}=high)")
    ax.set_ylabel("mean observed  log(1+ΔF/F)")
    if annotate:
        note(ax, f"Spearman ρ = {rho:.2f}\n(all pairs)")
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    hs.despine(ax)
    return rho


def topk_retrieval(ax, m, *, ks=(5, 10, 20, 40), color=C_RETRIEVAL,
                   chance_color=C_CHANCE,
                   title="retrieves the right strong gloms", title_fontsize=12):
    """Precision@k for each odorant's strongest glomeruli, against the chance
    line k/n_gloms. Returns ``(ks, precision)``."""
    ks, prec, chance = _pls.topk_retrieval(m, ks=ks)
    ax.plot(ks, prec, "o-", color=color, label="model")
    ax.plot(ks, chance, "s--", color=chance_color, label="chance")
    ax.set_xlabel("k (top glomeruli per odorant)")
    ax.set_ylabel("precision@k")
    ax.set_ylim(0, 1)
    ax.legend(fontsize=9, frameon=False)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    hs.despine(ax)
    return ks, prec


# ---------------------------------------------------------------------------
# Example maps
# ---------------------------------------------------------------------------
def glom_map(ax, m, values, *, vmin=MAP_VMIN, vmax=MAP_VMAX, size=42,
             under=hs.GAINSBORO, limits=None, title=None, title_fontsize=9):
    """One odorant's glomerular map, observed or predicted.

    Glomeruli below `vmin` are drawn flat grey *behind* the rest, so a
    near-silent field recedes to context and never occludes a responder. Pass
    the same `limits` (from `map_limits`) to every map that will be compared --
    that is what makes the grid one frame instead of twelve.
    """
    values = np.asarray(values, float)
    cmap = hs.SEQ.copy()
    cmap.set_under(under)
    below = values < vmin
    ax.scatter(m.x[below], m.y[below], color=under, s=size, linewidths=0,
               zorder=1)
    sc = ax.scatter(m.x[~below], m.y[~below], c=values[~below], cmap=cmap,
                    vmin=vmin, vmax=vmax, s=size, linewidths=0.25,
                    edgecolors="0.45", zorder=2)
    if limits is not None:
        (xlo, xhi), (ylo, yhi) = limits
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    hs.despine(ax, keep=())
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    return sc


def map_limits(m, pad=0.05):
    """``((xlo, xhi), (ylo, yhi))`` -- the shared frame for every example map."""
    xr = m.x.max() - m.x.min()
    yr = m.y.max() - m.y.min()
    return ((m.x.min() - pad * xr, m.x.max() + pad * xr),
            (m.y.min() - pad * yr, m.y.max() + pad * yr))


def map_aspect(m, pad=0.05):
    """Data height / width of the shared frame, for sizing map columns."""
    (xlo, xhi), (ylo, yhi) = map_limits(m, pad)
    return (yhi - ylo) / (xhi - xlo)


def map_colorbar(fig, cax, *, vmin=MAP_VMIN, vmax=MAP_VMAX, under=hs.GAINSBORO,
                 label="log(1 + ΔF/F)", label_fontsize=8, tick_fontsize=7):
    """Colorbar for a set of `glom_map` panels, with the below-floor extension
    that says the grey is a floor and not missing data."""
    cmap = hs.SEQ.copy()
    cmap.set_under(under)
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=mpl.colors.Normalize(vmin, vmax))
    cb = fig.colorbar(sm, cax=cax, extend="min")
    cb.set_label(label, fontsize=label_fontsize)
    cb.ax.tick_params(labelsize=tick_fontsize)
    return cb


# ---------------------------------------------------------------------------
# The weighting matrix W
# ---------------------------------------------------------------------------
def weight_matrix(ax, W, *, norm=None, percentile=99, cmap=None,
                  frame="0.6", title=None, title_fontsize=12, xlabel=None,
                  ylabel=None, label_fontsize=10):
    """One W-like matrix, signed, on the diverging map centred at zero.

    With no `norm`, the panel scales to its own `percentile` of |W| and shows
    ±1 -- which is what lets three matrices of very different magnitude be
    compared for STRUCTURE side by side. Say "per-panel scaled" in the caption
    when you do that, and read nothing into the absolute values.
    """
    W = np.asarray(W, float)
    cmap = hs.DIV if cmap is None else cmap
    if norm is None:
        scale = np.percentile(np.abs(W), percentile) or 1.0
        ax.imshow(W / scale, cmap=cmap, vmin=-1, vmax=1,
                  interpolation="nearest")
    else:
        ax.imshow(W, cmap=cmap, norm=norm, interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    if frame:
        for s in ax.spines.values():
            s.set_color(frame)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold", pad=6)
    return ax


def weight_blocks(ax, W_reord, boundaries, centers, *, norm=None,
                  cmap=None, line_color="0.15", line_width=0.8,
                  labels=None, label_fontsize=10,
                  xlabel="high-conc glomerulus (by cluster)",
                  ylabel="low-conc glomerulus (by cluster)",
                  title="W reordered by NMF cluster", title_fontsize=12):
    """W grouped by cluster, with the block boundaries drawn.

    Take `W_reord`, `boundaries` and `centers` from
    `highconc.pls.reorder_by_cluster`, which sorts into Joel's 1-7. A signed
    `norm` (typically a `TwoSlopeNorm` anchored at zero) is worth passing here:
    W's negative weights are smaller than its positive ones, so a symmetric
    range washes them out.
    """
    cmap = hs.DIV if cmap is None else cmap
    if norm is None:
        scale = np.percentile(np.abs(W_reord), 99) or 1.0
        ax.imshow(W_reord, cmap=cmap, vmin=-scale, vmax=scale,
                  interpolation="nearest")
    else:
        ax.imshow(W_reord, cmap=cmap, norm=norm, interpolation="nearest")
    for b in boundaries:
        ax.axvline(b - 0.5, color=line_color, lw=line_width)
        ax.axhline(b - 0.5, color=line_color, lw=line_width)
    labels = list(cl.LABELS_BY_JOEL) if labels is None else list(labels)
    ax.set_xticks(centers)
    ax.set_xticklabels(labels, fontsize=label_fontsize)
    ax.set_yticks(centers)
    ax.set_yticklabels(labels, fontsize=label_fontsize)
    ax.tick_params(length=0)
    ax.set_xlabel(xlabel, fontsize=label_fontsize)
    ax.set_ylabel(ylabel, fontsize=label_fontsize)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    return ax


def remixing_table(ax, M, *, norm=None, cmap=None, highlight=(),
                   labels=None, label_fontsize=10, highlight_width=2.0,
                   xlabel="high-conc cluster", ylabel="low-conc cluster",
                   title="cross-cluster remixing table\n(net signed block energy)",
                   title_fontsize=12):
    """The 7x7 cluster-to-cluster table, with chosen blocks outlined.

    `highlight` takes ``(from, to, value)`` triples in Joel's 1-7 -- pass
    `highconc.pls.top_off_diagonal(...)` straight in. Returns the image, for
    the colorbar.
    """
    cmap = hs.DIV if cmap is None else cmap
    im = ax.imshow(np.asarray(M, float), cmap=cmap, norm=norm, aspect="equal")
    labels = list(cl.LABELS_BY_JOEL) if labels is None else list(labels)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=label_fontsize)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=label_fontsize)
    ax.tick_params(length=0)
    ax.set_xlabel(xlabel, fontsize=label_fontsize)
    ax.set_ylabel(ylabel, fontsize=label_fontsize)
    for item in highlight:
        frm, to = int(item[0]), int(item[1])
        ax.add_patch(mpl.patches.Rectangle((to - 1.5, frm - 1.5), 1, 1,
                                           fill=False, edgecolor="black",
                                           lw=highlight_width))
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    return im


# ---------------------------------------------------------------------------
# Model overview: is there a linear rule, what does LV1 look like, and how
# does the linear rule compare to the field?
# ---------------------------------------------------------------------------
def component_selection(ax, sel, *, chosen=None, ms_q2=4, ms_overfit=3.5,
                        lw=1.8, mew=0.0, overfit_ref=2.0, tick_every=2,
                        title="Does a linear rule exist?  (component selection)",
                        title_fontsize=11):
    """Q^2 and the RMSECV/RMSEC overfit ratio against number of components.

    Two y-axes on purpose: the choice of n is a trade between the two curves,
    and separating them into two panels hides that. Takes a
    `highconc.pls.ComponentSelection`. Returns the twin axis so the caller can
    add to it.

    `chosen` defaults to the sweep's own `chosen` -- pass a number to mark a
    different n (e.g. to show what the van der Voet rule would have picked).
    The marker sizes are arguments because this panel is drawn both small, as
    one cell of the overview, and large, as a standalone figure.
    """
    n = sel.chosen if chosen is None else int(chosen)
    ax.plot(sel.grid, sel.q2, "-o", color=hs.C_PRIMARY, ms=ms_q2, lw=lw,
            mec="white", mew=mew, label="Q² (CV)")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.axvline(n, color="black", ls="--", lw=1.0)
    ax.text(n, sel.q2.max(), f" n={n}", fontsize=9, va="bottom", ha="left")
    ax.set_xlabel("n components")
    ax.set_ylabel("Q²  (CV)", color=hs.C_PRIMARY)
    ax.tick_params(axis="y", labelcolor=hs.C_PRIMARY)
    ax.set_xticks(sel.grid[::tick_every])
    ax.spines["top"].set_visible(False)

    ax2 = ax.twinx()
    ax2.plot(sel.grid, sel.overfit, ":s", color=hs.C_OVERFIT, ms=ms_overfit,
             lw=lw * 0.85, mec="white", mew=mew, label="RMSECV / RMSEC")
    if overfit_ref:
        ax2.axhline(overfit_ref, color="0.6", lw=0.6, ls=":")
    ax2.set_ylabel("RMSECV / RMSEC", color=hs.C_OVERFIT)
    ax2.tick_params(axis="y", labelcolor=hs.C_OVERFIT)
    ax2.spines["top"].set_visible(False)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    return ax2


def loading_outer(ax, outer, *, boundaries=(), centers=(), labels=None,
                  percentile=99, cmap=None, line_color="0.3", line_width=0.5,
                  xlabel="Y-loadings (high-conc, by cluster)", ylabel=None,
                  label_fontsize=9, tick_fontsize=8):
    """The rank-1 outer product x ⊗ y of one latent component, by cluster.

    Rows and columns must already be grouped -- pass the block `boundaries` and
    tick `centers` from `highconc.pls.group_by_cluster` so the heatmap and the
    two marginal skylines beside it cannot drift apart. Scaled symmetrically to
    its own `percentile`, on the diverging map centred at zero. Returns the
    image, for the colorbar.
    """
    outer = np.asarray(outer, float)
    vmax = np.percentile(np.abs(outer), percentile) or 1.0
    im = ax.imshow(outer, cmap=hs.DIV if cmap is None else cmap,
                   vmin=-vmax, vmax=vmax, interpolation="nearest")
    for b in boundaries:
        ax.axvline(b - 0.5, color=line_color, lw=line_width)
        ax.axhline(b - 0.5, color=line_color, lw=line_width)
    labels = list(cl.LABELS_BY_JOEL) if labels is None else list(labels)
    ax.set_xticks(list(centers))
    ax.set_xticklabels(labels, fontsize=tick_fontsize)
    ax.set_yticks(list(centers))
    ax.set_yticklabels(labels, fontsize=tick_fontsize)
    ax.tick_params(length=0)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
    return im


def loading_skyline(ax, values, cluster_number, *, orientation="left",
                    label=None, label_fontsize=9, title=None,
                    title_fontsize=11, title_pad=4):
    """One margin of `loading_outer`: per-glomerulus loadings as a bar skyline.

    Coloured by cluster through `highconc.style.cluster_color`, so a bar has
    the colour its cluster has everywhere else in the paper. `values` and
    `cluster_number` must already be in the heatmap's order.

    ``orientation="left"`` draws the X-loading margin -- horizontal bars
    growing leftward, y running top-down to match `imshow`. ``"top"`` draws the
    Y-loading margin as vertical bars. Frameless and tickless by design: it is
    a margin, and the numbers are read from the heatmap it annotates.
    """
    values = np.asarray(values, float)
    colors = [hs.cluster_color(int(c)) for c in np.asarray(cluster_number)]
    n = values.size
    if orientation == "left":
        ax.barh(np.arange(n), values, color=colors, height=1.0)
        ax.set_ylim(n - 0.5, -0.5)
        ax.axvline(0, color="0.6", lw=0.5)
        ax.invert_xaxis()
        if label:
            ax.set_ylabel(label, fontsize=label_fontsize)
    elif orientation == "top":
        ax.bar(np.arange(n), values, color=colors, width=1.0)
        ax.set_xlim(-0.5, n - 0.5)
        ax.axhline(0, color="0.6", lw=0.5)
        if label:
            ax.set_xlabel(label, fontsize=label_fontsize)
    else:
        raise ValueError(f"orientation must be 'left' or 'top', got {orientation!r}")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold", pad=title_pad)
    return ax


def model_zoo(ax, zoo=None, *, canonical=None, xlim=(-0.10, 0.42),
              value_fontsize=7.5, label_fontsize=9, annotate_baseline=True,
              xlabel="Q²  (global, log-transformed Y)",
              title="Model zoo — how much room above the linear baseline?\n"
                    "(family-stratified 4-fold CV)", title_fontsize=11):
    """Held-out Q^2 for the model zoo, ranked, against the canonical PLS-R.

    Colour carries the claim and nothing else: the canonical model in the
    comparison colour, models that beat it in the "wins" green, everything else
    neutral grey. Defaults to `highconc.pls.MODEL_ZOO`, whose values are cached
    from pls-explorer rather than refitted -- see that constant.
    """
    zoo = list(_pls.MODEL_ZOO if zoo is None else zoo)
    canonical = _pls.ZOO_CANONICAL if canonical is None else canonical
    base = _pls.zoo_baseline(zoo, canonical)
    labels = [name for name, _ in zoo]
    q2s = np.array([q for _, q in zoo], float)
    colors = [hs.C_OVERFIT if name == canonical
              else (hs.C_WIN if q > base + 1e-6 else hs.C_LOSE)
              for name, q in zoo]

    ypos = np.arange(len(labels))
    ax.barh(ypos, q2s, color=colors)
    ax.invert_yaxis()
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=label_fontsize)
    ax.axvline(base, color=hs.C_OVERFIT, ls="--", lw=1.2)
    ax.axvline(0, color="0.7", lw=1.0)
    for y, q in zip(ypos, q2s):
        ax.text(q + (0.006 if q >= 0 else -0.006), y, f"{q:+.3f}", va="center",
                ha="left" if q >= 0 else "right", fontsize=value_fontsize)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    if annotate_baseline:
        short = canonical.split(" (")[0]
        ax.text(base, -0.9, f"{short} baseline ({base:+.3f})",
                color=hs.C_OVERFIT, fontsize=8, ha="center", va="bottom")
    hs.despine(ax)
    if title:
        ax.set_title(title, fontsize=title_fontsize, weight="bold")
    return ax
