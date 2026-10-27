"""Dimensionality-choice panels — why k = 7 at 1x and k = 5 at 25x.

Extracted from ``figures/nmf_dim_choice_1x/generate.py`` and its 25x sibling,
which were near-identical: the same four panels over a different concentration's
cache. There is one implementation here, and concentration is an argument to
:func:`highconc.dimsweep.load_cached`, not a copy of the file.

    from highconc import dimsweep, style as hs
    from highconc.panels import dimchoice

    sweep, consensus = dimsweep.load_cached(cache_dir, conc="1x", ks=(7, 9))
    dimchoice.variance_explained(sweep, ax_a)
    dimchoice.variance_slope(sweep, ax_b)
    dimchoice.stability(sweep, ax_c)
    im = dimchoice.consensus(consensus[7], ax_d)

Every panel reads an already-loaded cache. None of them computes a sweep, and
none of them can: they take result objects, not a response matrix. That is the
guard against a notebook run quietly turning into a several-minute refit.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .. import style as hs
from ..dimsweep import ConsensusResult, DimSweepResult, N_NMF_SEEDS, N_SHUFFLES

__all__ = [
    "variance_explained", "variance_slope", "stability", "consensus",
    "CONSENSUS_VMIN", "CONSENSUS_VMAX", "consensus_cmap",
]

# Marker and line weights carried over from the pre-refactor scripts.
LW_CURVE = 1.6
LW_NMF = 1.4
MS_NMF = 6.5
BAND_ALPHA_NMF = 0.20
BAND_ALPHA_SHUFFLED = 0.25

#: Consensus colour range. Values below `vmin` are pairs that essentially never
#: co-cluster and read as background; 0.5-1.0 is where the instability gradient
#: that separates a good k from a bad one actually lives. The pre-refactor
#: notebook used (0.93, 1.01), tuned to a dataset whose variation sat just below
#: 1.0 — that range flattens these data to a single colour.
CONSENSUS_VMIN, CONSENSUS_VMAX = 0.5, 1.0


def _sem(std, n):
    """Standard error from a standard deviation. `n` is the count behind it."""
    return np.asarray(std) / np.sqrt(n)


def variance_explained(
    sweep: DimSweepResult,
    ax,
    *,
    n_nmf_seeds: int = N_NMF_SEEDS,
    n_shuffles: int = N_SHUFFLES,
    legend: bool = True,
):
    """Panel A — cumulative variance explained vs k, for PCA, shuffled PCA, NMF.

    The shuffled-PCA curve is the null: where NMF stops rising faster than
    chance is the argument for a k. Bands are standard errors, so they need the
    counts the sweep was run with — pass them if you recomputed the cache with
    something other than the defaults.
    """
    k = sweep.k_values
    ax.plot(k, sweep.pca_var_explained,
            color=hs.C_PCA, lw=LW_CURVE, label="PCA")
    ax.plot(k, sweep.pca_shuffled_var_mean,
            color=hs.C_SHUFFLED, lw=LW_CURVE, ls="--", label="PCA (shuffled)")
    shuf_sem = _sem(sweep.pca_shuffled_var_std, n_shuffles)
    ax.fill_between(k,
                    sweep.pca_shuffled_var_mean - shuf_sem,
                    sweep.pca_shuffled_var_mean + shuf_sem,
                    color=hs.C_SHUFFLED, alpha=BAND_ALPHA_SHUFFLED)
    ax.plot(k, sweep.nmf_var_explained_mean,
            color=hs.C_NMF, marker="o", ms=MS_NMF, lw=LW_NMF, label="NMF")
    nmf_sem = _sem(sweep.nmf_var_explained_std, n_nmf_seeds)
    ax.fill_between(k,
                    sweep.nmf_var_explained_mean - nmf_sem,
                    sweep.nmf_var_explained_mean + nmf_sem,
                    color=hs.C_NMF, alpha=BAND_ALPHA_NMF)

    ax.set_xlabel("Number of Components")
    ax.set_ylabel("Fraction of Variance Explained")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    if legend:
        ax.legend(fontsize=hs.FS_TICK, frameon=False, loc="lower right")
    return ax


def variance_slope(sweep: DimSweepResult, ax):
    """Panel B — the first difference of panel A's curves.

    Same three series, no error bands: the slope of a mean, not a mean of
    slopes, so the sweep's per-seed spread does not carry through.
    """
    k_mid = sweep.k_values[1:]
    ax.plot(k_mid, np.diff(sweep.pca_var_explained),
            color=hs.C_PCA, lw=LW_CURVE)
    ax.plot(k_mid, np.diff(sweep.pca_shuffled_var_mean),
            color=hs.C_SHUFFLED, lw=LW_CURVE, ls="--")
    ax.plot(k_mid, np.diff(sweep.nmf_var_explained_mean),
            color=hs.C_NMF, marker="o", ms=MS_NMF, lw=LW_NMF)

    ax.set_xlabel("Number of Components")
    ax.set_ylabel("Slope of Explained Variance")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    return ax


def stability(
    sweep: DimSweepResult,
    ax,
    *,
    n_nmf_seeds: int = N_NMF_SEEDS,
    ylim: tuple[float, float] | None = (0.75, 1.02),
):
    """Panel C — adjusted Rand index across NMF initialisations vs k.

    How reproducibly a given k recovers the same partition. `ylim` defaults to
    the pre-refactor window, which crops the k = 1..2 end; pass ``None`` to see
    the whole curve.
    """
    k = sweep.k_values
    ax.plot(k, sweep.nmf_ari_mean,
            color=hs.C_NMF, marker="o", ms=MS_NMF, lw=LW_CURVE)
    ari_sem = _sem(sweep.nmf_ari_std, n_nmf_seeds)
    ax.fill_between(k,
                    sweep.nmf_ari_mean - ari_sem,
                    sweep.nmf_ari_mean + ari_sem,
                    color=hs.C_NMF, alpha=BAND_ALPHA_NMF)

    ax.set_xlabel("number of components")
    ax.set_ylabel("adjusted Rand index")
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    return ax


def consensus_cmap(under: str = "whitesmoke"):
    """Reversed plasma: dark purple = always co-clusters, orange = sometimes.

    Reversed so the stable, uninteresting blocks are dark and the intermediate
    values — the ones that say a k is unstable — are the bright ones.
    """
    cm = plt.cm.plasma_r.copy()
    cm.set_under(under)
    return cm


def consensus(
    result: ConsensusResult,
    ax,
    *,
    vmin: float = CONSENSUS_VMIN,
    vmax: float = CONSENSUS_VMAX,
    cmap=None,
    colorbar: bool = True,
    ticks=(0.5, 0.75, 1.0),
    anchor: str = "E",
):
    """Panel D — one consensus matrix, square, with a flush colorbar.

    Clean diagonal blocks mean a k the data support; intermediate values
    scattered through the would-be blocks mean it is over-split.

    `colorbar` is on by default because the bar sits inside this panel's own
    space (an axes divider), not in the figure's layout — so unlike the other
    panels' colorbars, the notebook has nothing useful to do with it. Turn it
    off for the second matrix in a stack if one bar serves both. Returns the
    image, so a caller that wants its own bar still can.
    """
    im = ax.imshow(result.consensus,
                   cmap=consensus_cmap() if cmap is None else cmap,
                   vmin=vmin, vmax=vmax, aspect="equal")
    ax.set_anchor(anchor)   # right-align, so a wide column pads on the left
    ax.set_xticks([])
    ax.set_yticks([])
    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.03)
        cbar = ax.figure.colorbar(im, cax=cax)
        cbar.ax.tick_params(labelsize=hs.FS_TICK)
        if ticks is not None:
            cbar.set_ticks(list(ticks))
    return im
