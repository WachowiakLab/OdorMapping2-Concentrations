# %% [markdown]
# <!-- GENERATED FILE — do not edit. Source: notebooks/figures.py.
#      Rebuild with `python dev/build_notebooks.py`. -->
#
# *This is the **human** edition of the figure notebook — the friendly one — what each figure shows, and how it was made.*
# *Same code, same figures, different words: the other edition is `figures_LLM`.*
# *Edits go in `notebooks/figures.py`, never here — this file is overwritten.*

# %% [markdown]
# # High-concentration odor mapping — the figures
#
# This notebook draws the figures of the paper that can be built from the data, using the paper's own
# numbering and panel letters. Run it top to bottom and it rebuilds all of them into
# `notebooks/outputs/` as PNG, PDF and SVG.
#
# ### The basic story the figures tell
#
# - At three odorant concentrations — most-sensitive, low and high — the same 59 odorants drive very
#   different numbers of glomeruli, and more of them as the concentration rises (**Fig 2**).
# - Raising concentration changes which odorants evoke similar glomerular response patterns, not only
#   how strongly they respond (**Fig 3**).
# - Glomeruli group into seven clusters by shared odorant preference (**Fig 4**). The groups are stable
#   when glomeruli that only respond at high concentration are folded in: most keep their group, and the
#   newcomers spread across all seven (**Fig 4D**).
# - Each group favours particular chemical classes, several widen their reach at high concentration (**Fig 5**),
#   and the groups occupy different parts of the bulb (**Fig 6**).
# - A linear model fitted across all odorants predicts held-out high-concentration responses from the
#   low-concentration ones (**Fig 7**).
# - The broadening is selective and asymmetric: particular groups take on the tuning of particular others,
#   and not the reverse (**Fig 8**).
#
# ### What's in here
#
# Section numbers are the paper's figure numbers, so §3 draws Figure 3; the supplements follow, §9 to §12
# being Figures S1 to S4. A section draws only the panels of its figure that come from this data and code.
# Everything else in the paper — Figure 1, Figure 2B–E, Fig 6B, the schematics in Fig 7A and 7C, and
# the procedure diagram in Fig 8A — is drawn elsewhere, so a figure here can start at B or D.
#
# | § | Figure | In one line |
# |---|---|---|
# | 1 | setup | Loads the data and fits the model. Run this first. |
# | 2 | **Fig 2** | Odorant responses across Burton's own glomeruli, at three concentrations |
# | 3 | **Fig 3** | Odorant–odorant similarity, the recruitment matrix, and four example glomeruli |
# | 4 | **Fig 4** | The seven glomerular groups, and how they carry over when high-only glomeruli are folded in |
# | 5 | **Fig 5** | Which chemical classes each group responds to |
# | 6 | **Fig 6** | Where those groups sit on the bulb |
# | 7 | **Fig 7** | A linear model predicting the low → high change |
# | 8 | **Fig 8** | Directional transfer of tuning between groups |
# | 9 | **Fig S1** | How the number of groups was chosen |
# | 10 | **Fig S2** | The same structure within each mouse |
# | 11 | **Fig S3** | Where each mouse's groups sit on the bulb |
# | 12 | **Fig S4** | Predicted vs measured, over every responding pair |
#
# The paper has no Figure S3 from the old concentration-flow analysis, so that figure is no longer built.
#
# ### How to change something
#
# - **A panel's contents** (what's plotted, its colours): edit the function in `highconc/panels/`. Every
#   figure using it follows.
# - **A figure's layout** (grid, colourbars, panel letters): edit the cell that builds it, right here in the
#   notebook.
# - **Fonts, colours, export settings**: `highconc/style.py`, once, for everything. No cell sets its own.
#
# ### Two things that will bite you
#
# 1. **Run it in order.** §1 loads the data and fits the model that later sections share. Running a later
#    section on a fresh kernel refits it, which takes about a minute.
# 2. **Some results are cached on disk and are not recomputed** — the model fit, the dimensionality sweep in
#    §9, and the re-clustering with high-only glomeruli folded in (§4D, built once by
#    `dev/fold_nonresponders.py`). If a cache is missing the notebook stops with an error instead of quietly
#    recomputing, because a recomputed number might not match the one in the paper.
#
# Glomerular groups are numbered **1–7** throughout. Older figures and drafts label them A–G; the letters do
# not map onto the numbers in order, so don't assume cluster C is cluster 3. `highconc.clusters` holds the
# translation.

# %% [markdown]
# ## 1. Setup
#
# Loads the imaging data and fits the model that most of the later sections use.
# Run this before anything else. It takes about a minute the first time and is
# then cached, so re-running costs nothing.
#
# The cell after the imports also explains **how to export figures** — the short
# version is that a normal run already writes PNG, PDF and SVG of everything.

# %%
import dataclasses
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, PowerNorm, TwoSlopeNorm
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the 3d projection)

from highconc import (
    clusters as cl,
    dimsweep,
    flow,
    folding,
    load,
    pls,
    recruitment as rc,
    schematic,
    similarity,
    style as hs,
    transfer,
)
from highconc import panels as P

# Export contract (editable SVG text, TrueType in PDF), type scale, colours.
hs.apply()

# %%
# Where we are, and where renders go. `outputs/` is gitignored; the notebook can
# be re-run from a clean clone and will rebuild everything in it.
ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents]
            if (p / "pyproject.toml").exists())
FIGURES = ROOT / "figures"
OUT = ROOT / "notebooks" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

SAVE = True          # write png/pdf/svg to OUT as well as showing inline
FIG_DPI = 200        # raster preview density; vector formats are resolution-free
VARIANT = "nonresponders_dropped"

#: Every figure this notebook draws is scaled to span an exact physical
#: width: `MAIN_TARGET_W` for the numbered figures (2-7), `SUPP_TARGET_W`
#: for the supplements (S1-S4) -- smaller, and shared by all four rather than
#: each carrying its own inherited size. Each section solves its own
#: `*_SCALE` (or, for Fig 3, `FIG3_SCALE` applied to every margin/gap in the
#: section) so that section's own total width comes out to exactly the
#: target, then multiplies every other length in the section by that same
#: factor -- so raising or lowering either target here rescales every figure
#: in step, proportions held, no section needing a by-hand edit.
#:
#: This makes each figure's own physical size the invariant, at the cost of
#: an older one: before this, every section computed its scale from
#: `HEATMAP_IN` so that Fig 2 A, Fig 3 A-E, Fig 4 A-B, Fig 7 D, Fig 8 B,
#: S2's per-subject heatmaps and S3's coherence matrix all rendered at the
#: *same* physical size across figures. Two exact-width targets and one
#: shared heatmap size cannot both hold in general -- the figures' own
#: layouts (margins, panel counts, aspect ratios) differ too much -- so
#: pick one. `HEATMAP_PX`/`FIG_SCALE`/`HEATMAP_IN` below still set the
#: reference heatmap size a section's `*_SCALE` is *measured against*
#: (unchanged from before), but no longer guarantee it comes out identical
#: figure to figure once that section's scale is re-solved for its own
#: width target.
MAIN_TARGET_W = 7.0      # inches -- Figures 2-7
SUPP_TARGET_W = 5.5      # inches -- Figures S1-S4, smaller and shared
FIG_SCALE = 1.2
HEATMAP_PX = 260
HEATMAP_IN = HEATMAP_PX / FIG_DPI * FIG_SCALE

#: Chemical-class strip label size ("acids", "aldehydes", ...) on those same
#: panels -- small enough to fit next to a 0.6in heatmap without crowding it.
CLASS_STRIP_FONTSIZE = 5

#: Title size for the same panels. The house `FS_LABEL` (9pt) titles like
#: "'most sensitive'" are wider than a 0.6in/120px panel, so they overflow
#: into the next panel or the panel letter beside them; this is small enough
#: to mostly fit within one panel's width instead.
HEATMAP_TITLE_FONTSIZE = 7


EXPORT_FORMATS = ("png", "pdf", "svg")

#: Every figure this notebook has drawn, by stem, in the order drawn.
#: `save()` fills it; `export_figure()` reads it. Holding the figure objects
#: is what lets you re-export without recomputing anything.
DRAWN: dict[str, "plt.Figure"] = {}


def _write(fig, stem, outdir, formats, dpi):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(outdir / f"{stem}.{ext}", dpi=dpi,
                    bbox_inches="tight",
                    facecolor="white" if ext == "png" else "none",
                    transparent=(ext != "png"))
    return [outdir / f"{stem}.{ext}" for ext in formats]


def save(fig, stem, dpi=None):
    """Write one figure three ways, remember it, then show it.

    PNG on white for preview; PDF and SVG transparent with editable text, so
    either drops into Illustrator with the labels still real text objects.

    Called at the end of every figure cell. A plain top-to-bottom run therefore
    already writes all three formats of everything into `notebooks/outputs/` —
    you do not have to do anything to get them. `export_figure()` is for when
    you want them somewhere else, in fewer formats, or at a different DPI.
    """
    DRAWN[stem] = fig
    if SAVE:
        _write(fig, stem, OUT, EXPORT_FORMATS, dpi or FIG_DPI)
        print(f"saved: {stem}.{{png,pdf,svg}}")
    return fig


def export_figure(which=None, *, formats=EXPORT_FORMATS, outdir=OUT, dpi=None):
    """Re-export figures that have already been drawn. Nothing is recomputed.

    A normal run already writes `png`, `pdf` and `svg` of every figure into
    `notebooks/outputs/`. Use this when you want something different:

        export_figure()                        # the figure just drawn, again
        export_figure("fig7_pls")              # one, by stem
        export_figure("all")                   # every figure drawn so far
        export_figure("all", outdir="~/Desktop/highconc")
        export_figure("fig4_tuning_by_cluster", formats=("pdf",), dpi=600)

    Parameters
    ----------
    which : str, optional
        A figure's stem, or ``"all"``. Omitted, the most recently drawn figure.
    formats : sequence of str
        Any of ``"png"``, ``"pdf"``, ``"svg"``. PDF and SVG keep live text.
    outdir : path
        Created if missing. `~` is expanded.
    dpi : int, optional
        Raster density for the PNG. Vector formats ignore it.

    Notes
    -----
    Only figures drawn **in this kernel session** can be exported — the figure
    objects live in `DRAWN`. If you restart the kernel, run the notebook again
    first. `"all"` on a fresh kernel says so rather than silently doing nothing.
    """
    outdir = Path(outdir).expanduser()
    bad = [f for f in formats if f not in EXPORT_FORMATS]
    if bad:
        raise ValueError(f"unknown format(s) {bad}; expected any of "
                         f"{list(EXPORT_FORMATS)}")
    if not DRAWN:
        print("nothing drawn in this kernel yet — run the notebook first, "
              "then export.")
        return []
    if which == "all":
        stems = list(DRAWN)
    elif which is None:
        stems = [list(DRAWN)[-1]]
    elif which in DRAWN:
        stems = [which]
    else:
        raise KeyError(f"{which!r} has not been drawn in this session. "
                       f"Available: {sorted(DRAWN)}")

    written = []
    for stem in stems:
        written += _write(DRAWN[stem], stem, outdir, formats, dpi or FIG_DPI)
    print(f"exported {len(stems)} figure(s) x {len(formats)} format(s) "
          f"-> {outdir}")
    return written

# %%
# ---------------------------------------------------------------------------
# EXPORTING FIGURES
#
# You do not need to do anything: a top-to-bottom run already writes png, pdf
# and svg of every figure into `notebooks/outputs/`. The lines below are for
# when you want something else. Uncomment and run one AFTER the notebook has
# been run, so there is something to export.
#
# export_figure("all")                                # everything, again
# export_figure("all", outdir="~/Desktop/highconc")   # somewhere else
# export_figure("all", formats=("pdf",), dpi=600)     # just the vector PDFs
# ---------------------------------------------------------------------------

# %%
# The clustered dataset — 59 odorants x 376 glomeruli at two concentrations,
# with NMF cluster identity and bulb coordinates attached. Ships inside
# glom-explorer; there is no external data path to configure.
d = load(VARIANT)
print(f"{d.n_odorants} odorants x {d.n_gloms} gloms | "
      f"subjects {list(d.subjects)} | clusters 1-{cl.N_CLUSTERS}")

# %%
# The baseline PLS model, fitted once and shared by the figure sections below:
# log1p, mean-centered inside the pipeline, family-stratified 4-fold CV,
# 8 components. `pls.fit` is memoized, so the later sections are free.
m = pls.fit(VARIANT)
print(f"Q2_global = {pls.baseline_q2(m):+.4f} | RMSECV = {pls.baseline_rmse(m):.3f} "
      f"| {m.n_odorants} odorants x {m.n_gloms} gloms")

# %% [markdown]
# ## 2. Figure 2 — odorant responses across Burton's own glomeruli, at three concentrations
#
# Panel A of the paper's Figure 2: three response matrices over the same 59 odorants, stacked top to bottom
# at the "most sensitive", "low" and "high" concentrations. These are Burton et al.'s own glomeruli rather
# than this dataset's 376 — the "most sensitive" matrix has 76 glomeruli (a different, smaller preparation,
# `omp114L`/`omp8x`) and the other two have 113.
#
# - **Each cell** is one odorant–glomerulus pair; brighter means a bigger response. All three matrices share
#   one colour scale. Odorants run down the rows grouped by chemical class, with the coloured strip and the
#   class names on the right of each matrix — the same convention as Figure 3's strip, turned onto the row
#   axis since odorants are rows here. Columns (glomeruli) carry no identity or ordering information and are
#   left unlabelled.
# - **The matrices get denser as concentration rises**: 3.8% of cells are nonzero at "most sensitive", 4.9%
#   at low and 12.1% at high (the counts are printed under the figure).
# - **The top matrix is narrower** (76 columns against 113) only because its source preparation has fewer
#   glomeruli, not because anything was cropped.
#
# The rest of the paper's Figure 2 — the sorted response spectra, the histogram of glomeruli per odorant, and
# the lifetime-sparseness panels (B–E) — is drawn elsewhere.

# %%
# Own constants — this section inherits no bare name from anywhere, even
# though its legal values mirror §3's.
FIG2_ORDER = "matt"      # chemical-class block order -- highconc.clusters.class_order
FIG2_WITHIN = "burton"   # within-block odorant order -- highconc.clusters.WITHIN_BLOCK_ORDERS
FIG2_CONCS = ("most_sensitive", "low", "high")   # top to bottom, as in the paper's panel A

fig2_mats = [similarity.burton_response_matrix(d, c, order=FIG2_ORDER,
                                               within=FIG2_WITHIN)
            for c in FIG2_CONCS]
for label, mat in zip(FIG2_CONCS, fig2_mats):
    print(f"panel {label!r}: {mat.M.shape[1]} gloms, "
          f"{(mat.M != 0).mean():.1%} of cells nonzero")

# %%
# --- Figure 2 ---------------------------------------------------------------
# A over B over C rather than side by side: every panel is 59 odorants x
# however many gloms it has, so a wide, short box (rather than a squarish
# one) is what lets the glomerulus axis actually read as wider than the
# odorant axis is tall.
FIG2_W_MM, FIG2_H_MM = 200, 250                        # original target print size
# Scale the whole figure down by one factor rather than rebuild its layout,
# solved so the rendered width comes out to exactly MAIN_TARGET_W.
FIG2_SCALE = MAIN_TARGET_W * 25.4 / FIG2_W_MM
FIG2_W_MM, FIG2_H_MM = FIG2_W_MM * FIG2_SCALE, FIG2_H_MM * FIG2_SCALE
FIG2_W, FIG2_H = FIG2_W_MM / 25.4, FIG2_H_MM / 25.4    # -> inches for figsize

fig, axes = plt.subplots(3, 1, figsize=(FIG2_W, FIG2_H), dpi=200,
                         gridspec_kw=dict(left=0.16, right=0.68, top=0.95,
                                          bottom=0.05, hspace=0.30))
im = None
titles = ["'most sensitive'", "'low'", "'high'"]
for ax, mat, title in zip(axes, fig2_mats, titles):
    im = P.similarity.response_matrix(ax, mat, title=title, boundaries=False,
                                      strip_fontsize=CLASS_STRIP_FONTSIZE)
    ax.set_title(title, fontsize=HEATMAP_TITLE_FONTSIZE)
hs.panel_label(axes[0], "A", x=-0.16, y=1.14)   # one panel, A: the three matrices stacked
pos_b = axes[1].get_position()
cax = fig.add_axes([0.84, pos_b.y0, 0.03, pos_b.height])
P.similarity.similarity_colorbar(fig, im, cax, label="response")

save(fig, "fig2_response_matrix")
# plt.show()

# %% [markdown]
# ## 3. Figure 3 — odorant–odorant similarity, the recruitment matrix, and four example glomeruli
#
# Panel A is a chemical-structure similarity for the 59 odorants (MACCS fingerprints, cosine), which does not
# change with concentration. Panels B, C and D ask, for every pair of odorants, how similar the two patterns of
# glomerular activity are — at the "most sensitive" condition (B), at low concentration (C) and at high (D).
# Average similarity between different odorants is 0.03 at low and 0.10 at high, and the low and high matrices
# agree with each other at a rank correlation of 0.55. Chemical similarity tracks response similarity more
# closely at high concentration than at low (rank correlation 0.21 → 0.31).
#
# **Reading order is A, B, C, D, left to right** — the fixed chemical ruler first, then the three
# concentration-dependent conditions it is judged against. Panels E and F are the second row.
#
# **Reading the panels.**
#
# - **A–D** — every cell is one odorant pair; brighter means more alike. All four share one colour scale, so
#   differences between them are not a scaling artifact. Odorants are grouped along both edges by chemical class
#   (the coloured strip carries the class names under every panel), with the aromatics split into three finer
#   groups (Ar. Ketone, Ar. Ketone/Ester, Phenol/Methoxy). Green boxes mark the five broad groups (acids
#   through alcohols, aromatics and furans, amines, pyrazines, thiazoles). The grey diagonal is each odorant compared with
#   itself, which is 1 by definition; four odorants (D15, E1, P1, P4) have no data at the "most sensitive"
#   condition, so panel B reads them at the colour floor rather than as a gap.
# - **E** — recruitment. Each row is an odorant that drove a glomerulus at low concentration; each column an
#   odorant that newly drives it at high. The cell is the fraction of that row's glomeruli recruited. Grey is
#   undefined, not zero.
# - **F** — four individual glomeruli (173 and 211 on top, 97 and 155 below). Each bar is one odorant: the
#   low-concentration response points up, the high-concentration response points down. Bar colours are chemical
#   class — the plain ten-class version; the aromatics split is a panel A–D thing only. The downward bars are
#   magnitudes, not negative numbers.
#
# **Panel B's data.** `data/burton_et_al/most_sensitive.csv` and `odorant_order.txt`, supplied 2026-09-02 — the
# "most sensitive" condition this dataset didn't otherwise have.
#
# **The example glomeruli** are numbers 173, 211, 97 and 155. Nothing on record said which glomeruli the original
# panels used for the first three, so they were identified by measuring the bars in the draft PDF and matching them
# against all 376; they match to within the PDF's own rounding. 211 and 97 appear again as the predicted-versus-
# measured examples in Figure 7B. Glomerulus 155 was added later, by its number in the dataset.
#
# **Panel letters follow the manuscript figure PDF** (`allfigs_hiconcms_penultimate_102126.pdf`).

# %%
CLASS_ORDER_FIG3 = "matt"      # "table" | "matt" | "matt_esters_first"
WITHIN_CLASS_ORDER = "matt"    # "table" | "matt" (pyrazine only) | "burton" (full order)
FIG3_SPLIT_AROMATICS = True    # split "R" into 3 display subclasses on panels A-D
# Diagonal boxes on A-D, as in the paper: five groups under CLASS_ORDER_FIG3="matt" with
# FIG3_SPLIT_AROMATICS -- acids through alcohols (C+D+E+K+L, 23), aromatics+furans (R+F, 12),
# amines (N, 7), pyrazines (P, 10), thiazoles (S, 7). A different CLASS_ORDER_FIG3 or
# FIG3_SPLIT_AROMATICS setting changes the block sizes these cuts land on, so they are not
# derived automatically.
FIG3_BOX_BOUNDARIES = (23, 35, 42, 52)
FIG3_BOX_COLOR = "green"
SPECTRUM_SCALE = "per-panel"   # "per-panel" | "shared"
EG_GLOMERULI = similarity.EG_GLOMERULI     # panel F's first three example glomeruli
# Fourth example glomerulus (panel F), by dataset glomerulus id like EG_GLOMERULI.
# Kept out of EG_GLOMERULI, which other code imports.
EXTRA_GLOM = 155

# Odorant CODEs (the leading token of `Odorant`), positionally aligned with
# `d`'s native odorant order -- `recruitment_matrix`'s `within` needs these,
# (only this section draws a recruitment matrix now, so only this section uses it).
ODORANT_CODES = np.asarray([s.split()[0] for s in d.odorants["Odorant"].values])

mats = [similarity.most_sensitive_correlation(d, order=CLASS_ORDER_FIG3,
                                              within=WITHIN_CLASS_ORDER),
        similarity.response_correlation(d, "low", order=CLASS_ORDER_FIG3,
                                        within=WITHIN_CLASS_ORDER),
        similarity.response_correlation(d, "high", order=CLASS_ORDER_FIG3,
                                        within=WITHIN_CLASS_ORDER),
        similarity.chemical_similarity(d, order=CLASS_ORDER_FIG3,
                                       within=WITHIN_CLASS_ORDER)]
if FIG3_SPLIT_AROMATICS:
    mats = [similarity.split_block(s, d, "R", cl.AROMATIC_SUBCLASS_BY_CODE)
           for s in mats]
spectra = similarity.eg_spectra(d, [*EG_GLOMERULI, EXTRA_GLOM],
                                order=CLASS_ORDER_FIG3, within=WITHIN_CLASS_ORDER)

shift = similarity.concentration_shift(mats[1], mats[2])
print(f"mean off-diagonal correlation: {shift['mean_low']:.4f} (low) -> "
      f"{shift['mean_high']:.4f} (high), delta {shift['delta']:+.4f}")
print(f"the two matrices agree at Spearman {shift['spearman']:.3f} "
      f"over {shift['n_pairs']} odorant pairs")
for s in spectra:
    print(f"glom {s.glom:3d} (cluster {s.cluster}): "
          f"{int(s.responded.sum()):2d} odorants responding, max {s.limit:.1f}")

# %%
# Panel E's matrix is computed here rather than inherited from §7, so the two
# sections cannot drift and either can be run without the other.
# `recruitment.fit` is memoized, so the second call costs nothing.
# A glomerulus is assigned to a low-concentration odorant when it responds at 0.20 of its own
# low-concentration maximum ("halfmax" is the rule's name; the fraction is ROW_FRAC_FIG3), and is
# recruited by the 0.20/0.20 relative test below. The supplement below sweeps the fraction.
ROW_RULE_FIG3 = "halfmax"
ROW_FRAC_FIG3 = 0.20
PRED_THRESHOLD_FIG3 = "rate_matched"
RECRUITED_FIG3 = "relative"     # "exact" = low == 0 & high > 0 (the original)
LO_FRAC_FIG3 = 0.20             # recruited: low  < 0.20 x the glomerulus's low peak
HI_FRAC_FIG3 = 0.20             # recruited: high > 0.20 x the glomerulus's high peak

rmx3 = rc.recruitment_matrix(m, level="odorant", row_rule=ROW_RULE_FIG3,
                             frac=ROW_FRAC_FIG3, order=CLASS_ORDER_FIG3,
                             within=WITHIN_CLASS_ORDER, codes=ODORANT_CODES,
                             prediction_threshold=PRED_THRESHOLD_FIG3,
                             recruited=RECRUITED_FIG3, lo_frac=LO_FRAC_FIG3,
                             hi_frac=HI_FRAC_FIG3)
print(f"panel E: rows built from {rmx3.n_rows.mean():.1f} gloms on average; "
      f"{rmx3.n_empty} empty")

# %%
# --- Figure 3 ---------------------------------------------------------------
# The paper's outline: row 1 is the four odorant x odorant matrices (A the chemical MACCS
# matrix, then the three response matrices); row 2 is the recruitment matrix (E, left) and
# the four example spectra in a 2 x 2 block (F, right).
fig = plt.figure(figsize=(14.0, 10.6), dpi=200)
gs = fig.add_gridspec(2, 4, height_ratios=[1.00, 1.20], hspace=0.22, wspace=0.30,
                      left=0.07, right=0.91, top=0.965, bottom=0.045)

top_panels = (("A", mats[3], "MACCS cosine"),
              ("B", mats[0], "'most sensitive'"),
              ("C", mats[1], "'low' (1x)"),
              ("D", mats[2], "'high' (25x)"))
axes_top4 = []
for j, (L, mat, title) in enumerate(top_panels):
    ax = fig.add_subplot(gs[0, j])
    # The x-axis chemical-class strip, with its family names, is drawn on all
    # four panels (they share the same class order).
    im = P.correlation_matrix(ax, mat, title=title, boundaries=False,
                              strip_labels=True,
                              box_boundaries=FIG3_BOX_BOUNDARIES,
                              box_color=FIG3_BOX_COLOR)
    hs.panel_label(ax, L, x=-0.16)
    axes_top4.append(ax)
fig.canvas.draw()                          # aspect="equal" fixes the box size only now
pos_d = axes_top4[-1].get_position()
cax = fig.add_axes([pos_d.x1 + 0.012, pos_d.y0, 0.010, pos_d.height])
P.similarity_colorbar(fig, im, cax)

# E: the recruitment matrix, square, on the left half of row 2.
axE = fig.add_subplot(gs[1, 0:2])
axE.set_box_aspect(1)
imE = P.recruitment.matrix(axE, rmx3, "observed",
                           vmax=P.recruitment.shared_vmax(rmx3), ylabel=True,
                           title="binary recruitment matrix (observed)")
hs.panel_label(axE, "E", x=-0.16, y=1.05)
fig.canvas.draw()
pos_e = axE.get_position()
caxE = fig.add_axes([pos_e.x1 + 0.012, pos_e.y0, 0.010, pos_e.height])
cbE = fig.colorbar(imE, cax=caxE, extend="max")
cbE.set_label("recruitment probability", fontsize=hs.FS_LABEL)
cbE.ax.tick_params(labelsize=hs.FS_TICK)

# F: one 2 x 2 block of example spectra, filled row by row (173, 211 / 97, 155).
gs_F = gs[1, 2:4].subgridspec(2, 2, hspace=0.85, wspace=0.42)
LIM = P.shared_spectrum_limit(spectra) if SPECTRUM_SCALE == "shared" else None
for k, s in enumerate(spectra):
    ax = fig.add_subplot(gs_F[k // 2, k % 2])
    P.response_spectrum(ax, s, limit=LIM,
                        ylabel="Resp. amplitude" if k % 2 == 0 else "")
    ax.set_title(f"glomerulus {s.glom} · cluster {s.cluster}",
                 fontsize=hs.FS_LABEL)
    if k == 0:
        hs.panel_label(ax, "F", x=-0.30)

save(fig, "fig3_correlational_structure")
# plt.show()
# export_figure("fig3_correlational_structure")   # -> notebooks/outputs/fig3_correlational_structure.{png,pdf,svg}
# export_figure("fig3_correlational_structure", formats=("pdf",), outdir="~/Desktop")

# %% [markdown]
# ### Supplement to Figure 3 — the recruitment matrix at five cutoffs
#
# Panel E depends on one choice: how strongly an odorant must drive a glomerulus before that glomerulus counts as
# assigned to it. This supplement redraws the matrix at five values of that cutoff, observed above and held-out
# predicted below, and the cell after it prints the observed-versus-predicted correlation at each cutoff so you can
# see how much it moves.
#
# The paper's figure PDF has no such figure, so it is an extra here, and it saves under its own name.

# %%
SWEEP_FRACS_FIG3 = (0.20, 0.35, 0.50, 0.65, 0.80)

sweep3 = rc.matrix_sweep(m, fracs=SWEEP_FRACS_FIG3, level="odorant",
                         row_rule=ROW_RULE_FIG3, order=CLASS_ORDER_FIG3,
                         prediction_threshold=PRED_THRESHOLD_FIG3,
                         recruited=RECRUITED_FIG3, lo_frac=LO_FRAC_FIG3,
                         hi_frac=HI_FRAC_FIG3)
VMAX_S3 = P.recruitment.shared_vmax(*sweep3)

fig, ax = plt.subplots(2, len(sweep3), figsize=(3.0 * len(sweep3), 6.8), dpi=150)
for j, s in enumerate(sweep3):
    for i, which in enumerate(("observed", "predicted")):
        im = P.recruitment.matrix(ax[i, j], s, which, vmax=VMAX_S3,
                                  tick_labels=None, xlabel=False,
                                  ylabel=False, title="", boundaries=False)
    P.recruitment.matrix_sweep_label(ax[0, j], s)
ax[0, 0].set_ylabel("observed")
ax[1, 0].set_ylabel("predicted (held-out)")
fig.colorbar(im, ax=ax, fraction=0.014, pad=0.02, extend="max",
             label="recruitment probability")
fig.suptitle("Supplement to Figure 3 — the recruitment matrix does not depend "
             "on where the assignment rule is cut", fontsize=13, weight="bold")
save(fig, "fig3_supp_threshold_robustness")
# plt.show()
# export_figure("fig3_supp_threshold_robustness")   # -> notebooks/outputs/fig3_supp_threshold_robustness.{png,pdf,svg}
# export_figure("fig3_supp_threshold_robustness", formats=("pdf",), outdir="~/Desktop")

# %%
for s in sweep3:
    print(f"frac={s.frac:.2f}  empty rows={s.n_empty:2d}  "
          f"mean row n={s.n_rows.mean():5.1f}  r_offdiag={s.r_offdiag:.3f}")

# %% [markdown]
# ## 4. Figure 4 — the seven glomerular groups, and how they carry over when high-only glomeruli are folded in
#
# Panels A–C show the same 376 glomeruli three ways; panel D asks whether the grouping survives adding more.
#
# - **A** — every pair of glomeruli, scored on how similar their odorant preferences are, ordered by group.
#   Within-group correlation runs 0.12–0.51 against a between-group mean of −0.03. (Grey diagonal = each
#   glomerulus versus itself.) A green box outlines each of the seven diagonal within-group squares.
# - **B** — the response matrix behind that: odorants down, glomeruli across. Brightness is each glomerulus's
#   response **relative to its own strongest odorant** (0–1, linear), not an absolute ΔF/F — so the panel answers
#   "which odorants drive this glomerulus, relative to its own best one," not "how big is this response." Glomerulus
#   columns are sorted by group (the vertical lines mark where one group ends and the next begins); odorant rows
#   run in the same chemical-class order as Figure 3, with each of the ten classes named on the left, so a group's
#   odorants are scattered rather than forming a block. A response below the noise floor (1 ΔF/F) stays grey even
#   after normalising, so a glomerulus whose best response never clears the floor reads as empty, not as a falsely
#   "bright" 1.0.
# - **C** — the same data one mouse at a time, for groups 1 and 3, in the same row order and with the same
#   per-glomerulus normalisation as B, so a glomerulus looks the same in both. Thin white lines separate the
#   chemical classes.
# - **D** — what happens to the groups when the glomeruli that only respond at *high* concentration are folded in.
#   The 118 such glomeruli (silent at low, responding at high) are added to the low-concentration data, using their
#   high response, and the clustering is run again with the same seven groups. Each ribbon follows glomeruli from
#   their original group (left, 1–7; grey = the 118 newcomers, "recruited at high") to the group they land in (right,
#   numbered 1–7 by which original group each new one most resembles). Glomeruli that already responded at low
#   mostly keep their group — from 100% of group 1 and 98% of group 2 down to 67% of group 7 — and the newcomers
#   spread across all seven.
#
# **Which two groups are in panel C.** Group 1, which responds to a mix of acids and aldehydes, and group 3, which
# responds almost only to furans. The paper contrasts a mixed group with a narrowly selective one, and these two are
# the clearest such pair that is present in all four mice (group 5, the amine group, has no glomeruli in mice 3 and
# 4, which would leave two panels empty).
#
# The paper's panel for chemical class by group is now Figure 5B, and the per-odorant tuning bars are Figure 5A.

# %%
# Every constant here is prefixed and owned by this section, so that no two
# sections can share a bare name like `VMIN`/`VMAX` and silently rebind it.
FIG4_VMIN, FIG4_VMAX = 1, 140      # dF/F colour range of the response heatmaps
FIG4_CMAP = "magma"
FIG4_EMPTY = "gainsboro"           # below vmin -> grey
FIG4_BOUNDARY, FIG4_BOUNDARY_ALPHA = "blue", 0.5

FIG4_CONC_A = "low"          # 4A is the low-concentration matrix
FIG4C_CLUSTERS = (1, 3)      # mixed acids/aldehydes cluster, furan cluster
# Reuses §3's own CLASS_ORDER_FIG3 directly (not an independent copy), so "the
# same broad chemical order as Fig 3" is guaranteed rather than kept in sync
# by hand -- same coupling §5 and §7 have to §3's class order.
FIG4_CHEM_ORDER = CLASS_ORDER_FIG3
# Toggle for B/C's odorant order *within* each cluster block (columns stay
# glomerulus-cluster-sorted either way -- this only ever touches rows).
# True: chemical class (FIG4_CHEM_ORDER). False: native odorant order, i.e.
# `bicluster_order`'s own chem_order=None default -- ties in the stable sort
# keep their original order, so this is not "no sort", it's "no secondary
# sort key".
FIG4B_SORT_BY_CLASS = True

# Panel B's own row order -- independent of `bic` below, which panel C still
# reads. "cluster": the manuscript default, rows bi-clustered like the
# columns (diagonal block highlight makes sense). "odorant": rows follow
# Fig 3's own odorant order (FIG4_CHEM_ORDER, WITHIN_CLASS_ORDER) instead, ignoring
# cluster membership entirely -- a cluster's rows are then scattered, not
# contiguous, so the diagonal highlight is turned off for this mode (wired
# in below via `diagonal_blocks`). Vertical glomerulus-cluster lines draw
# either way -- only the rows move.
FIG4B_ROW_ORDER = "odorant"

# Panels B and C both: peak-normalise each glomerulus (column) to [0, 1] by its
# own max response,
# on the raw dF/F values, and draw on a LINEAR scale instead of log -- the
# two go together (see `biclustered`'s `col_normalize` docstring). `FIG4_VMIN`
# /`FIG4_VMAX` (tuned for raw dF/F on a log scale) do not apply in this mode;
# `FIG4B_NORM_VMIN`/`FIG4B_NORM_VMAX` are panel B's own limits instead.
FIG4B_COL_NORMALIZE = True
FIG4B_NORM_VMIN, FIG4B_NORM_VMAX = 0.0, 1.0
# Restores the "nonzero cutoff" (the same FIG4_VMIN=1 floor the non-normalised
# log view uses) even though column-normalising otherwise erases it: a cell
# whose RAW dF/F was below this never draws on the ramp, however large its
# normalised value would have been -- so a glomerulus whose peak never
# clears the noise floor greys out entirely instead of showing a misleadingly
# "bright" 1.0 for its own weak best odorant. `None` disables this and lets
# every normalised value draw, including noise-floor-weak ones.
FIG4B_RAW_FLOOR = FIG4_VMIN

gsim = similarity.glomerulus_correlation(d, FIG4_CONC_A)
# The block boundaries (and the diagonal alignment with the glomerulus
# columns) are the same regardless of FIG4B_SORT_BY_CLASS; only the order
# inside a block moves. Shared by C (and, when FIG4B_ROW_ORDER="cluster", by
# B too); B reads its own `bic_b` below whenever FIG4B_ROW_ORDER="odorant",
# so C's row order never moves when B's does.
bic = P.clustering.bicluster_order(
    d, chem_order=(FIG4_CHEM_ORDER if FIG4B_SORT_BY_CLASS else None),
    within_order=WITHIN_CLASS_ORDER)
bic_b = P.clustering.bicluster_order(
    d, chem_order=FIG4_CHEM_ORDER, within_order=WITHIN_CLASS_ORDER,
    row_order=FIG4B_ROW_ORDER)

bm = gsim.block_means()
within = np.diag(bm.values)
between = bm.values[~np.eye(len(bm), dtype=bool)]
print(f"panel A: within-cluster correlation {within.min():.2f}-{within.max():.2f}, "
      f"between-cluster mean {between.mean():+.3f}")
for k in FIG4C_CLUSTERS:
    n = [int((d.subject_mask(s) & d.cluster_mask(k)).sum()) for s in d.subjects]
    print(f"panel C: cluster {k} has {sum(n)} glomeruli, {n} per subject")

# %%
# Scale the whole figure by one factor rather than rebuild its gridspec
# math -- every gridspec call below is fraction-based, so it shrinks
# proportionally for free. Solved so the rendered width is MAIN_TARGET_W.
FIG4_SCALE = MAIN_TARGET_W / 14.0
FIG4_W, FIG4_H = 16.0 * FIG4_SCALE, 11.0 * FIG4_SCALE

fig = plt.figure(figsize=(FIG4_W, FIG4_H), dpi=FIG_DPI)

# --- top row: A (square correlation matrix) and B (wide response heatmap) ---
gs_top = fig.add_gridspec(1, 2, width_ratios=[1, 3.1], wspace=0.30,
                          left=0.06, right=0.93, top=0.95, bottom=0.60)
ax_a = fig.add_subplot(gs_top[0, 0])
ax_b = fig.add_subplot(gs_top[0, 1])

im_a = P.clustering.glomerulus_correlation_matrix(
    ax_a, gsim, title={"low": "low", "high": "high"}[FIG4_CONC_A],
    boundaries=False, cluster_boxes=True, cluster_box_color="green")
mesh_b = P.clustering.biclustered(d, ax_b, conc=FIG4_CONC_A, order=bic_b,
                                  diagonal_blocks=(FIG4B_ROW_ORDER == "cluster"),
                                  class_names=True,
                                  class_name_fontsize=CLASS_STRIP_FONTSIZE,
                                  col_normalize=FIG4B_COL_NORMALIZE,
                                  raw_floor=(FIG4B_RAW_FLOOR
                                            if FIG4B_COL_NORMALIZE else None),
                                  vmin=(FIG4B_NORM_VMIN if FIG4B_COL_NORMALIZE
                                       else FIG4_VMIN),
                                  vmax=(FIG4B_NORM_VMAX if FIG4B_COL_NORMALIZE
                                       else FIG4_VMAX),
                                  cmap=FIG4_CMAP, empty_color=FIG4_EMPTY,
                                  boundary_color=FIG4_BOUNDARY,
                                  boundary_alpha=FIG4_BOUNDARY_ALPHA)
# House type sizes, so the labels between A and B and under A stay clear of each other.
ax_a.title.set_fontsize(hs.FS_LABEL)
ax_a.xaxis.label.set_fontsize(hs.FS_LABEL)
ax_a.yaxis.label.set_fontsize(hs.FS_LABEL)
ax_b.set_ylabel("")          # the ten class names already label B's rows
ax_b.xaxis.label.set_fontsize(hs.FS_LABEL)

pos_a = ax_a.get_position()
cax_a = fig.add_axes([pos_a.x0, pos_a.y0 - 0.055, pos_a.width * 0.6, 0.014])
cb_a = fig.colorbar(im_a, cax=cax_a, orientation="horizontal",
                    ticks=[0.0, 0.5, 1.0])
# The label goes beside the colourbar, not under it: below it would run into C's titles.
cb_a.ax.text(1.06, 0.5, r"correlation $\rho$", transform=cb_a.ax.transAxes,
             ha="left", va="center", fontsize=hs.FS_LABEL)
cb_a.ax.tick_params(labelsize=hs.FS_TICK)
cb_a.outline.set_visible(False)

pos_b = ax_b.get_position()
cax_b = fig.add_axes([pos_b.x1 + 0.012, pos_b.y0, 0.010, pos_b.height * 0.7])
if FIG4B_COL_NORMALIZE:
    # Linear 0-1: each glomerulus's response relative to its own peak, not
    # an absolute ΔF/F -- pairs with `col_normalize` in the draw call above.
    cb_b = fig.colorbar(
        plt.cm.ScalarMappable(
            norm=plt.Normalize(vmin=FIG4B_NORM_VMIN, vmax=FIG4B_NORM_VMAX),
            cmap=hs.sequential_cmap(under=FIG4_EMPTY, bad=FIG4_EMPTY)),
        cax=cax_b)
    cb_b.set_label("resp., rel. to glom. peak", fontsize=hs.FS_LABEL)
    cb_b.set_ticks([0.0, 0.5, 1.0])
    cb_b.set_ticklabels(["0", "0.5", "1"])
else:
    cb_b = fig.colorbar(
        plt.cm.ScalarMappable(norm=LogNorm(vmin=FIG4_VMIN, vmax=FIG4_VMAX),
                              cmap=hs.sequential_cmap(under=FIG4_EMPTY,
                                                      bad=FIG4_EMPTY)),
        cax=cax_b)
    cb_b.set_label("ΔF/F", fontsize=hs.FS_LABEL)
    cb_b.set_ticks([1, 10, 100])
    cb_b.set_ticklabels(["1", "10", "100"])
cb_b.ax.minorticks_off()
cb_b.ax.tick_params(labelsize=hs.FS_TICK)

# --- bottom row: C (two clusters x four subjects) and D (class x cluster) ---
gs_bot = fig.add_gridspec(
    1, 3, width_ratios=[1, 1, 1.15], wspace=0.30,
    left=0.06, right=0.93, top=0.46, bottom=0.06)

# C's rows are B's rows: Fig 3's exact 59-odorant order (CLASS_ORDER_FIG3
# blocks, WITHIN_CLASS_ORDER within them), not the cluster-sorted rows `bic`
# carries. Only the rows move: columns stay one cluster's glomeruli per
# subject. C also takes B's peak-normalisation and raw floor, so a glomerulus
# reads identically in both and B's colorbar is C's too.
perm_c, _, _, class_bounds_c = similarity.class_sorted(
    d, CLASS_ORDER_FIG3, WITHIN_CLASS_ORDER)
assert FIG4B_ROW_ORDER == "odorant" and np.array_equal(perm_c, bic_b.rows), \
    "C's rows must be B's rows (FIG4B_ROW_ORDER='odorant', same class/within order)"
bic_c = dataclasses.replace(bic,
                            rows=perm_c,
                            row_boundaries=np.array([], dtype=int),
                            row_edges=np.array([0, len(perm_c)]))
FIG4C_LINE_COLOR, FIG4C_LINE_LW, FIG4C_LINE_ALPHA = "white", 0.6, 0.8

im_c = None
for g, k in enumerate(FIG4C_CLUSTERS):
    inner = GridSpecFromSubplotSpec(1, len(d.subjects), subplot_spec=gs_bot[0, g],
                                    wspace=0.18)
    axes_c = [fig.add_subplot(inner[0, j]) for j in range(len(d.subjects))]
    for j, subj in enumerate(d.subjects):
        im_c = P.clustering.per_subject_heatmap(
            d, axes_c[j], subj, cluster=k, conc=FIG4_CONC_A, order=bic_c,
            col_normalize=FIG4B_COL_NORMALIZE,
            raw_floor=FIG4B_RAW_FLOOR if FIG4B_COL_NORMALIZE else None,
            vmin=FIG4B_NORM_VMIN if FIG4B_COL_NORMALIZE else FIG4_VMIN,
            vmax=FIG4B_NORM_VMAX if FIG4B_COL_NORMALIZE else FIG4_VMAX,
            cmap=FIG4_CMAP, empty_color=FIG4_EMPTY, cluster_labels=False,
            title_fontsize=hs.FS_LABEL, ylabel_fontsize=hs.FS_LABEL,
            ylabel=("odorants (by chemical class)" if j == 0 else None))
        for b in class_bounds_c:               # line between chemical classes
            axes_c[j].axhline(b - 0.5, color=FIG4C_LINE_COLOR,
                              lw=FIG4C_LINE_LW, alpha=FIG4C_LINE_ALPHA)
    mid = (axes_c[0].get_position().x0 + axes_c[-1].get_position().x1) / 2
    fig.text(mid, axes_c[0].get_position().y1 + 0.045, f"cluster {k}",
             ha="center", va="bottom", fontsize=hs.FS_LABEL,
             color=hs.cluster_color(k))
    fig.text(mid, axes_c[0].get_position().y0 - 0.032, "glomeruli",
             ha="center", va="top", fontsize=hs.FS_TICK)

# D — where each glomerulus goes when the recruited non-responders are folded in.
# Glomeruli silent at low but responding at high have their high column folded into the
# low matrix and NMF is run again. That fit is cached by `dev/fold_nonresponders.py` (seeded,
# same settings as the original clustering) and only READ here -- never refit at draw time.
# Left: the original clusters 1-7 (grey = the recruited glomeruli, which were not in the
# original clustering). Right: the new clusters, numbered 1-7 by best overlap with the originals.
fig4_flows = folding.flows(d)
ax_d = fig.add_subplot(gs_bot[0, 2])
P.sankey.flow_diagram(ax_d, fig4_flows, left_title="Low-Conc.\nOnly",
                      right_title="Low-Conc. +\nHigh-Conc. recruited",
                      unmatched_label="recruited\nat high")
print(f"D: {fig4_flows.n_matched} original + {fig4_flows.n_added} recruited glomeruli; "
      f"share of each original cluster that stays in its matched new cluster: "
      f"{', '.join(f'{c + 1}: {s:.0%}' for c, s in enumerate(fig4_flows.stayed))}")

hs.panel_label(ax_a, "A", x=-0.22, y=1.12)
hs.panel_label(ax_b, "B", x=-0.035, y=1.06)
pos_d = ax_d.get_position()
fig.text(0.02, 0.48, "C", fontsize=hs.FS_LABEL + 4, weight="bold", va="top")
fig.text(pos_d.x0 - 0.055, 0.48, "D", fontsize=hs.FS_LABEL + 4,
         weight="bold", va="top")

save(fig, "fig4_tuning_by_cluster")
# plt.show()
# export_figure("fig4_tuning_by_cluster")   # -> notebooks/outputs/fig4_tuning_by_cluster.{png,pdf,svg}

# %% [markdown]
# ## 5. Figure 5 — which chemical classes each glomerular group responds to
#
# Panels A and B of the paper's Figure 5.
#
# **A — tuning by cluster.** One column per glomerular group, one row per odorant, grouped into the ten chemical
# classes. Each bar is the fraction of that group's glomeruli that respond to the odorant: low concentration is the
# blue bar to the left of zero and high concentration is the orange bar to the right. A is a plotly figure saved on its
# own (`fig5_tuning_by_cluster_A`); B is a separate matplotlib figure. The number at the bottom of each
# column is how many glomeruli the group has. Rows follow Figure 3's odorant order, including the reordered pyrazine
# block. A glomerulus "responds" when its response is above zero.
#
# **B — chemical classes.** For each of the seven glomerular groups, how many of its glomeruli respond to
# each chemical class, scaled so the group's top class is 1 — at low concentration (left) and at high concentration
# (right). Each group has its own favourites at low concentration. At high concentration the bright cells spread for
# several groups, while others — the amine group, for one — stay selective.
#
# - **Rows** are chemical classes, in Figure 3's order; **columns** are the seven groups. Class names and group
#   numbers are in their own colours.
# - **A cell** counts the glomeruli in that group with a nonzero response to at least one odorant of that class
#   (a glomerulus that responds to five odorants of a class still counts once), then divides each column by its own
#   largest cell. So columns are *relative* — a group's top class is always 1 — and do not sum to 1.
# - **Low and high are each scaled within themselves**, so a bright cell means "this class reaches many of this
#   group's glomeruli," not "more than at the other concentration." One colourbar serves both.

# %%
# Own constants — this section inherits no bare name from anywhere.
FIG5_CHEM_ORDER = CLASS_ORDER_FIG3     # same broad chemical order as Figs 3, 4 and 7 ("matt")
FIG5_WITHIN = WITHIN_CLASS_ORDER       # odorant order inside a class ("matt": the pyrazine block reordered)
FIG5_CMAP = "magma"
FIG5_EMPTY = "gainsboro"
FIG5_BOUNDARY, FIG5_BOUNDARY_ALPHA = "blue", 0.5
FIG5_A_H, FIG5_B_H = 5.5, 2.9          # inches; the width is MAIN_TARGET_W
FIG5_A_XLIM = 1.0                      # bars run to +-this fraction of a cluster's glomeruli
FIG5_A_LEFT, FIG5_A_RIGHT = 0.95, 1.60 # inches reserved for class names (left) and odorant names (right)
FIG5_A_GAP = 0.9                       # extra rows of white space between chemical classes
FIG5_A_NAME_FS = 4.6

# A's numbers: per odorant and cluster, the fraction of the cluster's glomeruli with a response above zero,
# at low and at high. Rows follow Figure 3's odorant order (class blocks, then FIG5_WITHIN inside each).
_perm5, _classes5, _blocks5, _ = similarity.class_sorted(d, FIG5_CHEM_ORDER, FIG5_WITHIN)
fig5a_low = cl.fraction_responding(d, "low", _perm5)
fig5a_high = cl.fraction_responding(d, "high", _perm5)
fig5a_names = np.asarray(d.odorants["Name"].values)[_perm5]
fig5a_table = pd.concat(
    {conc: pd.DataFrame(mat, index=fig5a_names, columns=[f"cluster_{k}" for k in range(1, 8)])
     for conc, mat in (("low", fig5a_low), ("high", fig5a_high))},
    names=["conc", "odorant"])
fig5a_table.insert(0, "chem_class", np.tile(np.asarray(_classes5), 2))
assert max(fig5a_low.max(), fig5a_high.max()) <= FIG5_A_XLIM, "a bar would be clipped; raise FIG5_A_XLIM"
if SAVE:
    fig5a_csv = OUT / "fig5_class_by_cluster_A_fraction_responsive.csv"
    fig5a_table.to_csv(fig5a_csv)
    print(f"saved: {fig5a_csv.relative_to(ROOT)}")

# B's numbers, one row per chemical class and one column per cluster, computed once and reused
# by the drawing calls and the CSV export, so the saved numbers are exactly what the panel draws.
# Each concentration is normalised within itself (every cluster column / its own top class).
_fig5_low, _fig5_letters = P.clustering.class_cluster_glom_composition(d, "low", FIG5_CHEM_ORDER)
_fig5_high, _ = P.clustering.class_cluster_glom_composition(d, "high", FIG5_CHEM_ORDER)
fig5_table = pd.concat(
    {conc: pd.DataFrame(matrix, index=_fig5_letters,
                        columns=[f"cluster_{k}" for k in range(1, 8)])
     for conc, matrix in (("low", _fig5_low), ("high", _fig5_high))},
    names=["conc", "chem_class"])
if SAVE:
    fig5_csv = OUT / "fig5_class_by_cluster_B_glom_composition.csv"
    fig5_table.to_csv(fig5_csv)
    print(f"saved: {fig5_csv.relative_to(ROOT)}")

# %%
# --- Figure 5A: tuning by cluster, mirrored bars (plotly) --------------------------
# Plotly draws its own canvas, so A is exported by itself -- PNG, SVG and PDF through kaleido -- and B below
# is a separate matplotlib figure. The two are put side by side on the page when the figure is assembled.
fig5a = P.tuning_plotly.mirrored_tuning_figure(
    fig5a_low, fig5a_high, fig5a_names, _classes5,
    class_names=cl.CHEM_CLASS_NAMES, class_colors=cl.CHEM_CLASS_COLORS,
    cluster_colors=[hs.cluster_color(n) for n in range(1, 8)], n_glom=cl.SIZES_BY_JOEL,
    width_in=MAIN_TARGET_W, height_in=FIG5_A_H, xlim=FIG5_A_XLIM,
    low_color=hs.C_LOW, high_color=hs.C_HIGH, class_gap=FIG5_A_GAP,
    left_in=FIG5_A_LEFT, right_in=FIG5_A_RIGHT, name_font_pt=FIG5_A_NAME_FS)
if SAVE:
    for _ext in ("png", "svg", "pdf"):
        fig5a.write_image(str(OUT / f"fig5_tuning_by_cluster_A.{_ext}"), scale=(FIG_DPI / 100 if _ext == "png" else 1))
    print("saved: fig5_tuning_by_cluster_A.{png,svg,pdf}")
fig5a.show() if not SAVE else None

# %%
# --- Figure 5B: low concentration beside high concentration ----------------------
fig = plt.figure(figsize=(MAIN_TARGET_W, FIG5_B_H), dpi=FIG_DPI)
gs5 = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.045], wspace=0.55,
                       left=0.11, right=0.92, top=0.84, bottom=0.17)
summary_kw_5 = dict(class_order=FIG5_CHEM_ORDER, cmap=FIG5_CMAP, empty_color=FIG5_EMPTY,
                    boundary_color=FIG5_BOUNDARY, boundary_alpha=FIG5_BOUNDARY_ALPHA,
                    xlabel=None, label_fontsize=hs.FS_LABEL, tick_fontsize=hs.FS_TICK)
ax_low = fig.add_subplot(gs5[0, 0])
ax_high = fig.add_subplot(gs5[0, 1])
mesh_low = P.clustering.cluster_summary(d, ax_low, conc="low", **summary_kw_5)
mesh_high = P.clustering.cluster_summary(d, ax_high, conc="high", **summary_kw_5)
for ax_s, title in ((ax_low, "low concentration"), (ax_high, "high concentration")):
    ax_s.set_title(title, fontsize=hs.FS_LABEL)
    for lbl in ax_s.get_xticklabels():                      # cluster numbers, in cluster colours
        lbl.set_color(hs.cluster_color(int(lbl.get_text())))
    for lbl, letter in zip(ax_s.get_yticklabels(), _fig5_letters):   # class names, in class colours
        lbl.set_color(cl.class_color(letter))

cax5 = fig.add_subplot(gs5[0, 2])
cb5 = fig.colorbar(mesh_high, cax=cax5)
cb5.set_label("responding glomeruli,\nrelative to cluster's top class", fontsize=hs.FS_LABEL)
cb5.ax.tick_params(labelsize=hs.FS_TICK)
hs.panel_label(ax_low, "B", x=-0.42, y=1.10)

save(fig, "fig5_class_by_cluster")
# plt.show()
# export_figure("fig5_class_by_cluster")   # -> notebooks/outputs/fig5_class_by_cluster.{png,pdf,svg}

# %% [markdown]
# ## 6. Figure 6 — the seven groups on the dorsal bulb
#
# Panel A of the paper's Figure 6. Each map is the same view of the dorsal bulb, with all four mice overlaid (the
# field of view is consistent between them). One group is highlighted per map: its glomeruli are the solid red
# points, and every other glomerulus is a small hollow grey circle behind them. The headings are the group numbers,
# in each group's own colour. Groups differ in how spatially compact they are.
#
# The red points show **membership only** — which glomeruli belong to the group — and carry no response value, so
# there is no colour scale and no concentration to choose. The groups have 45, 41, 40, 75, 30, 78 and 67 glomeruli
# (printed under the figure).
#
# **Only panel A is built here.** Panel B in the paper — the bulb images, maximum-projected by odorant class — comes
# from the imaging data and is drawn elsewhere.

# %%
# Membership only: every cluster member is one solid colour, and no response
# value is drawn -- so there is no concentration to choose and no colourbar.
FIG6_MEMBER_COLOR = "red"
FIG6_NONRESPONDER = "hollow"        # or "context"/"requested", or a dict of scatter kwargs
FIG6_CLUSTER_EDGE = "black"         # see below

# This section has no heatmap panel to key off (see §1), so it gets its own
# scale rather than reading FIG_SCALE directly, solved for MAIN_TARGET_W.
# `_fig6_k` is how much smaller that makes the page than it was under the
# old bare-FIG_SCALE sizing -- marker size (area, so k**2) and linewidth
# (linear, so k) are rescaled by it too, or a page this much smaller would
# leave 7 giant, overlapping dots per row instead of 7 legible clusters.
FIG6_SCALE = MAIN_TARGET_W / 16.0
_fig6_k = FIG6_SCALE / FIG_SCALE
FIG6_CLUSTER_SIZE = 62 * _fig6_k ** 2   # the cluster's own points, drawn on top
FIG6_CLUSTER_LW = 0.6 * _fig6_k
FIG6_CONTEXT_SIZE = 16 * _fig6_k ** 2   # "hollow"'s own base size (NONRESPONDER_MARKERS), scaled
FIG6_CONTEXT_LW = 0.4 * _fig6_k

fig = plt.figure(figsize=(16.0 * FIG6_SCALE, 3.6 * FIG6_SCALE), dpi=FIG_DPI)
gs5 = fig.add_gridspec(1, cl.N_CLUSTERS, wspace=0.06,
                       left=0.05, right=0.97, top=0.82, bottom=0.06)
xlim5, ylim5 = P.clustering.bulb_limits(d)

for k in range(1, cl.N_CLUSTERS + 1):
    ax = fig.add_subplot(gs5[0, k - 1])
    P.clustering.topography(
        d, ax, k, color=FIG6_MEMBER_COLOR,
        size=FIG6_CLUSTER_SIZE, edgecolor=FIG6_CLUSTER_EDGE,
        linewidth=FIG6_CLUSTER_LW, nonresponder_marker=FIG6_NONRESPONDER,
        context_size=FIG6_CONTEXT_SIZE, context_linewidth=FIG6_CONTEXT_LW,
        title=str(k), title_fontsize=hs.FS_LABEL + 2)
    ax.title.set_color(hs.cluster_color(k))
    ax.set_xlim(*xlim5)
    ax.set_ylim(*ylim5)
    if k == 1:
        hs.panel_label(ax, "A", x=-0.30, y=1.12)

save(fig, "fig6_spatial_by_cluster")
# plt.show()
# export_figure("fig6_spatial_by_cluster")   # -> notebooks/outputs/fig6_spatial_by_cluster.{png,pdf,svg}

# %%
for k in range(1, cl.N_CLUSTERS + 1):
    mask = d.cluster_mask(k)
    xy = d.xy[mask]
    spread = float(np.sqrt(((xy - xy.mean(0)) ** 2).sum(1).mean()))
    print(f"cluster {k}: n = {int(mask.sum()):3d}, "
          f"rms spread {spread:6.1f} px about its centroid")

# %% [markdown]
# ## 7. Figure 7 — a linear model of the low → high change
#
# One linear model, fitted once across all odorants, predicting high-concentration responses from the
# low-concentration ones. Every prediction shown is **held out**: the model is fitted on some odorants and asked
# about others.
#
# - **B** — two individual glomeruli, 211 (group 2) and 97 (group 6). Measured high-concentration response points
#   up; the model's prediction for the same glomerulus points down. The light line across the lower half marks the
#   smallest nonzero response ever measured, so a predicted bar under the line is smaller than anything seen in the
#   data. The correlation between measured and predicted is in each title (0.78 and 0.64).
# - **D** — the weight matrix W itself: what a glomerulus at low concentration contributes to each glomerulus at
#   high, ordered by group. **Every weight is drawn.** The colour scale clips each sign at its own 99th percentile,
#   so a few extreme weights do not wash out the rest, but nothing is left out of the picture.
#
# **A caveat about B.** The bars are in **log1p ΔF/F** — the space the model was actually fitted and scored on, not
# the raw units Figure 3F uses. In raw ΔF/F, one acid's held-out prediction for glomerulus 211 reaches **193.6**
# against a measured **66.2**, which flattens every other bar once the axis autoscales to fit it; `FIG7_SPACE = "raw"`
# restores those units, at that cost. The correlations depend on the space too — the same predictions score higher on
# the log scale the model was scored in — so the scale is stated on the axis.
#
# Predictions near zero still draw a visible bar: above a threshold of 0.4 the model has nonzero predictions for 25
# and 51 odorants (glomeruli 211 and 97) against 11 and 29 measured. A cross-validated linear model essentially never
# predicts exactly zero, so this is the model's behaviour, not a drawing artifact.
#
# The paper's Figure 7 also has the X → W → Y schematic (A) and three cartoons of the weight matrix (C), which are
# drawn elsewhere.

# %%
FIG7_SPACE = "log"              # "log" (the model's own space) | "sqrt" | "raw" (dF/F)
FIG7_GLOMERULI = (211, 97)      # B's two example glomeruli (dataset glomerulus ids)
FIG7_LIMIT = None               # None = autoscale each panel
FIG7_SUPPORT_THRESHOLD = 0.4    # value above which a bar counts as a response
FIG7_CANONICAL = "baseline"
# D shows EVERY weight of W: there is no top-K cut. Only the colour scale is
# clipped (each sign at its own 99th percentile, see below), never the data.

# Same odorant order as Fig 3's F, so B reads against it.
fig7_specs = pls.eg_predicted_spectra(m, d, FIG7_GLOMERULI, space=FIG7_SPACE,
                                      order=CLASS_ORDER_FIG3,
                                      within=WITHIN_CLASS_ORDER)
fig7_w = pls.load_weights(FIG7_CANONICAL)
# The smallest nonzero observed high response in the whole dataset, in
# FIG7_SPACE: the line drawn across C's predicted half.
fig7_obs_floor = pls.observed_high_floor(m, space=FIG7_SPACE)

for s in fig7_specs:
    n_obs, n_pred, n_both = s.support(FIG7_SUPPORT_THRESHOLD)
    print(f"7B: glomerulus {s.glom} (cluster {s.cluster}) - "
          f"r = {s.r:.3f} in {s.space} space; above "
          f"{FIG7_SUPPORT_THRESHOLD}: {n_obs} measured, {n_pred} predicted, "
          f"{n_both} both")
print(f"7B: smallest nonzero observed high response ({FIG7_SPACE}): "
      f"{fig7_obs_floor:.4f}")
print(f"7D: W is {fig7_w.real.shape[0]} x {fig7_w.real.shape[1]}, "
      f"all weights shown")

# 7B's own numbers, one row per odorant per glomerulus -- built from the same
# `fig7_specs` the panels draw, so the saved numbers are exactly what's shown.
fig7_c_table = pd.concat(
    {s.glom: pd.DataFrame({
        "odorant": s.names, "class": s.classes,
        "observed": s.observed, "predicted": s.predicted,
     }).assign(cluster=s.cluster, space=s.space, r=s.r)
     for s in fig7_specs},
    names=["glom", "row"])
if SAVE:
    fig7_c_csv = OUT / "fig7_pls_B_prediction_spectra.csv"
    fig7_c_table.to_csv(fig7_c_csv)
    print(f"saved: {fig7_c_csv.relative_to(ROOT)}")

# %%
# One row, in inches: B (two example glomeruli, side by side) | D (W, square, + colourbar).
# Panels A (the X -> W -> Y schematic) and C (the three weight-matrix cartoons) of the paper's
# Figure 7 are drawn separately; the letters here are the paper's. The width is MAIN_TARGET_W
# exactly; the right margin is whatever is left, so changing any other width cannot silently
# change the page width.
FIG7_MARGIN_L, FIG7_MARGIN_TOP, FIG7_MARGIN_BOT = 0.58, 0.55, 0.60
FIG7_B_W, FIG7_B_H, FIG7_B_GAP = 1.40, 1.05, 0.50   # each spectrum, and the gap between them
FIG7_COL_GAP = 0.45                                  # between B and D
FIG7_D_SIDE = 1.75                                   # D is square
FIG7_CBAR_PAD, FIG7_CBAR_W = 0.30, 0.12
FIG7_MARGIN_R = (MAIN_TARGET_W - FIG7_MARGIN_L - 2 * FIG7_B_W - FIG7_B_GAP
                 - FIG7_COL_GAP - FIG7_D_SIDE - FIG7_CBAR_PAD - FIG7_CBAR_W)
assert FIG7_MARGIN_R > 0.30, f"Fig 7 row does not fit: right margin {FIG7_MARGIN_R:.2f} in"
FIG7_W = MAIN_TARGET_W
FIG7_H = FIG7_MARGIN_TOP + FIG7_D_SIDE + FIG7_MARGIN_BOT


def x6(inches):
    return inches / FIG7_W


def y6(inches):
    return inches / FIG7_H


fig = plt.figure(figsize=(FIG7_W, FIG7_H), dpi=FIG_DPI)
b_bottom = FIG7_MARGIN_BOT + (FIG7_D_SIDE - FIG7_B_H) / 2     # B centred on D's height

# --- B: two example glomeruli, experimental high up / predicted high down ----
for j, s in enumerate(fig7_specs):
    x0 = FIG7_MARGIN_L + j * (FIG7_B_W + FIG7_B_GAP)
    ax = fig.add_axes([x6(x0), y6(b_bottom), x6(FIG7_B_W), y6(FIG7_B_H)])
    P.similarity.prediction_spectrum(
        ax, s, limit=FIG7_LIMIT, strip_labels=True, annotate_r=False,
        observed_floor=fig7_obs_floor, label_observed="measured",
        label_predicted="predicted", **({} if j == 0 else {"ylabel": ""}))
    # r goes in the title, as in Fig 3's F, not inside the axes
    ax.set_title(f"glomerulus {s.glom}  (cluster {s.cluster})    r = {s.r:.2f}",
                 fontsize=HEATMAP_TITLE_FONTSIZE)
    if j == 0:
        hs.panel_label(ax, "B", x=-0.30, y=1.18)

# --- D: the full W, every weight ---------------------------------------------
d_left = FIG7_MARGIN_L + 2 * FIG7_B_W + FIG7_B_GAP + FIG7_COL_GAP
W6 = fig7_w.real
W6_reord, _, bnd6, centers6 = pls.reorder_by_cluster(W6, fig7_w.cluster_number)
# Each sign scaled to its own 99th percentile -- negatives are smaller, and a
# symmetric range washes them out. TwoSlopeNorm keeps 0 on gainsboro. This only
# sets the colour range; every weight is still drawn.
pos6, neg6 = W6_reord[W6_reord > 0], -W6_reord[W6_reord < 0]
# norm6 = TwoSlopeNorm(vmin=-(np.percentile(neg6, 99) if neg6.size else 1.0),
#                      vcenter=0.0,
#                      vmax=(np.percentile(pos6, 99) if pos6.size else 1.0))
norm6 = TwoSlopeNorm(vmin=-0.03, vcenter=0.0, vmax=0.04)
ax_d = fig.add_axes([x6(d_left), y6(FIG7_MARGIN_BOT),
                     x6(FIG7_D_SIDE), y6(FIG7_D_SIDE)])
P.pls.weight_blocks(ax_d, W6_reord, bnd6, centers6, norm=norm6,
                    title="W (low → high)")
ax_d.title.set_fontsize(HEATMAP_TITLE_FONTSIZE)
ax_d.xaxis.label.set_fontsize(HEATMAP_TITLE_FONTSIZE)
ax_d.yaxis.label.set_fontsize(HEATMAP_TITLE_FONTSIZE)
ax_d.tick_params(labelsize=hs.FS_TICK)
hs.panel_label(ax_d, "D", x=-0.30, y=1.10)
cax_d = fig.add_axes([x6(d_left + FIG7_D_SIDE + FIG7_CBAR_PAD),
                      y6(FIG7_MARGIN_BOT + FIG7_D_SIDE * 0.22),
                      x6(FIG7_CBAR_W), y6(FIG7_D_SIDE * 0.56)])
cb_d = fig.colorbar(plt.cm.ScalarMappable(cmap=hs.DIV, norm=norm6), cax=cax_d,
                    extend="both")
cb_d.set_label("W weight", fontsize=hs.FS_LABEL)
cb_d.ax.tick_params(labelsize=hs.FS_TICK)

save(fig, "fig7_pls")
# plt.show()
# export_figure("fig7_pls")   # -> notebooks/outputs/fig7_pls.{png,pdf,svg}

# %% [markdown]
# ## 8. Figure 8 — directional transfer of tuning between groups
#
# **The test.** Predict a group's high-concentration responses from its own low-concentration average — that is the
# baseline. Then add a second group's low-concentration average and see whether the prediction improves. Repeat for
# every ordered pair, in both directions. Prediction is scored on held-out odorants, so an improvement is not simply
# the result of adding a parameter.
#
# - **B** — the permutation *p* for each pair, against 10,000 shuffles. A **green box** marks the pairs that survive
#   Benjamini–Hochberg correction at q < 0.05, and only those boxed squares carry a number.
# - **C** — the size of each improvement, as a percentage of the error the target's own average left behind. Most
#   cells are near zero; again only the boxed pairs carry a number.
#
# Rows are the **source** group whose low-concentration average is added ("added non-source cluster"); columns are the
# **target** group being predicted at high concentration. The matrix is not symmetric and is not meant to be.
#
# **A control that is no longer drawn.** The same test can group glomeruli by the chemical class each one prefers
# instead of by the NMF clusters, which does not use the clustering at all. That analysis is still run, and the
# cells after the figure print the cluster edges next to their chemical-class counterparts (seven of the eight
# significant cluster results have one) — but the paper's Figure 8 shows only the cluster version, so only that is
# drawn.
#
# **The explanatory header in the paper (panel A) is drawn elsewhere.**

# %%
FIG8_METRIC = "pct_err"      # "pct_err" (default) | "pct_var" | "dq2"
FIG8_ALPHA = 0.05            # BH-FDR level; the black boxes
FIG8_CLIP = 0.1              # p above this is the palest cell — colour only
FIG8_NULL = "broken_pairing"  # the published null; see the markdown above
FIG8_MIN_CLASS = 15          # glomeruli a chemical class needs to be tested
FIG8_BOX_COLOR = "#2b9743"   # the box on significant cells (green)
FIG8_P_CMAP = "OrRd"         # A keeps the original sequential ramp (hs.SEQ is magma now)
FIG8_XLABEL = "cluster (predicted at high)"
FIG8_YLABEL = "added non-source cluster\nprediction improvement"

fig8_clusters = transfer.pairwise(m, "clusters", alpha_fdr=FIG8_ALPHA,
                                  null=FIG8_NULL)
fig8_classes = transfer.pairwise(m, "chem_class", alpha_fdr=FIG8_ALPHA,
                                 null=FIG8_NULL)

for t in (fig8_clusters, fig8_classes):
    s = t.summary(FIG8_ALPHA)
    print(f"{s['grouping']:>11}: {s['n_significant']} of {s['n_pairs']} "
          f"ordered pairs at q < {FIG8_ALPHA} "
          f"({s['n_one_way']} one-way), N = {s['n_perm']:,}")
print()
print(P.transfer.edge_table(fig8_clusters, alpha=FIG8_ALPHA,
                            metric=FIG8_METRIC))

# %%
# One row: B and C (the paper's letters; its A is the procedure schematic, drawn elsewhere). The chemical-class control (what
# were C and D) is no longer drawn; it is still computed, printed and checked
# in the cells below.
# Scale the whole figure by one factor, solved so the rendered width is
# MAIN_TARGET_W.
FIG8_SCALE = MAIN_TARGET_W / 13.0
FIG8_W, FIG8_H = 13.0 * FIG8_SCALE, 6.0 * FIG8_SCALE
fig = plt.figure(figsize=(FIG8_W, FIG8_H), dpi=FIG_DPI)

gs7 = fig.add_gridspec(1, 2, left=0.09, right=0.90, top=0.82, bottom=0.17,
                       wspace=0.42)

# Only the significant squares carry a number, and they are boxed in green.
ax_p = fig.add_subplot(gs7[0, 0])
P.transfer.pvalue_matrix(ax_p, fig8_clusters, alpha=FIG8_ALPHA, clip=FIG8_CLIP,
                         xlabel=FIG8_XLABEL, ylabel=FIG8_YLABEL, cmap=FIG8_P_CMAP,
                         box_color=FIG8_BOX_COLOR, text_significant_only=True)
hs.panel_label(ax_p, "B", x=-0.26, y=1.12)

ax_m = fig.add_subplot(gs7[0, 1])
P.transfer.improvement_matrix(ax_m, fig8_clusters, metric=FIG8_METRIC,
                              alpha=FIG8_ALPHA, xlabel=FIG8_XLABEL, ylabel="",
                              box_color=FIG8_BOX_COLOR,
                              text_significant_only=True)
hs.panel_label(ax_m, "C", x=-0.16, y=1.12)
P.transfer.significance_legend(ax_m, fig8_clusters, alpha=FIG8_ALPHA,
                               loc=(1.02, -0.16), box_color=FIG8_BOX_COLOR)

fig.text(0.5, 0.98, "Motif transfer between glomerular groups "
                    f"(broken-pairing null, N = {fig8_clusters.n_perm:,})",
         ha="center", va="top", fontsize=hs.FS_TITLE + 1)

save(fig, "fig8_motif_transfer")
# plt.show()
# export_figure("fig8_motif_transfer")   # -> notebooks/outputs/fig8_motif_transfer.{png,pdf,svg}

# %%
# The control's actual claim, checked rather than asserted: each significant
# cluster edge is looked up against the class edge between the two clusters'
# dominant preferred classes. Seven of eight should match. If this ever prints
# fewer, the two analyses have drifted apart — which is exactly what the
# control exists to detect.
fig8_match = transfer.matching_edges(fig8_clusters, fig8_classes, m,
                                     alpha=FIG8_ALPHA)
print(f"{'cluster edge':>16}  {'ΔQ²':>8}  {'via':>26}  chemical-class counterpart")
for row in fig8_match:
    e, c = row["cluster_edge"], row["counterpart"]
    via = f"{row['via'][0]} → {row['via'][1]}"
    tail = (f"{c['dq2']:+.4f}  (q = {c['q']:.4f})" if c else "— none —")
    print(f"{e['source'] + ' → ' + e['target']:>16}  {e['dq2']:+8.4f}  "
          f"{via:>26}  {tail}")
n_match = sum(r["counterpart"] is not None for r in fig8_match)
print(f"\n{n_match} of {len(fig8_match)} cluster edges have a class counterpart")

# %%
# Per-group intra baselines — the denominator of both percentage readings, and
# the reason `pct_var` is not the default: cluster 5 is predicted from its own
# mean at Q² = 0.046 and cluster 2 at 0.438, so the same absolute increment
# would read ten times larger on one than on the other.
for t in (fig8_clusters, fig8_classes):
    print(f"{t.grouping.id:>11} intra Q²: " +
          "  ".join(f"{lab}={q:.3f}"
                    for lab, q in zip(t.labels, t.q2_intra)))
print(f"\nchemical classes tested: {list(fig8_classes.labels)}")
print(f"glomeruli excluded (class below n = {FIG8_MIN_CLASS}): "
      f"{fig8_classes.grouping.dropped} of {m.n_gloms}")

# %% [markdown]
# ## 9. Figure S1 — how the number of groups was chosen
#
# Seven groups at low concentration. Three panels:
#
# - **A** — cumulative variance explained against the number of groups, for PCA, PCA on shuffled data, and NMF. The
#   gap between the real and shuffled curves is the structure not attributable to chance.
# - **B** — the slope of A: how much each additional group adds.
# - **C** — adjusted Rand index across 50 NMF runs from different random starts, measuring how stable the grouping is
#   to the seed.
#
# The number chosen (7) is marked by the dotted vertical line and printed under the figure. The paper's version also shows consensus matrices at the chosen
# number and one past it, and the matching analysis at high concentration; neither is drawn here.
#
# **This one is slow and is cached.** The underlying sweep takes minutes, so it is stored on disk and read back. If
# the cache is missing the notebook stops with an error rather than recomputing, because a recomputed number might
# differ from the paper's.

# %%
# Own constants — this section inherits no bare name from anywhere.
SUPP_DIM_K_VALUES = list(range(1, 21))   # the sweep on disk; changing this invalidates it
SUPP_DIM_N_SEEDS, SUPP_DIM_N_SHUFFLES = 50, 20
SUPP_DIM_ROWS = [
    # (conc, cache dir, consensus k's, chosen k, panel letters)
    ("1x",  FIGURES / "nmf_dim_choice_1x" / "data" / "cache",  [7, 9], 7, "ABC"),
    # The 25x block (panels D-F) is no longer drawn. Its sweep cache is still on
    # disk under figures/nmf_dim_choice_25x/ and is simply not read here.
]
# One block now, so the 12.5in-tall two-block hand layout is cut to 5.5in; the
# whole thing is scaled so the width is exactly SUPP_TARGET_W, proportions held.
SUPP_DIM_SCALE = SUPP_TARGET_W / 9.0
SUPP_DIM_W, SUPP_DIM_H, SUPP_DIM_DPI = 9.0 * SUPP_DIM_SCALE, 5.5 * SUPP_DIM_SCALE, 300

fig = plt.figure(figsize=(SUPP_DIM_W, SUPP_DIM_H), dpi=SUPP_DIM_DPI)
# One block (the 1x concentration). Inside a block the layout is the retired
# per-concentration figure's: the
# two line panels on top, then a bottom strip holding the ARI curve and the
# stacked consensus matrices, whose height matches the original 7 x 1.75-inch
# rectangle.
gs_rows = fig.add_gridspec(len(SUPP_DIM_ROWS), 1, hspace=0.32)

for row, (conc, cache_dir, consensus_ks, chosen, letters) in enumerate(SUPP_DIM_ROWS):
    sweep, consensus = dimsweep.load_cached(cache_dir, conc=conc, ks=consensus_ks)

    gs_block = gs_rows[row].subgridspec(2, 1, height_ratios=[2, 1], hspace=0.40)
    gs_top = gs_block[0].subgridspec(1, 2, wspace=0.40)
    gs_bot = gs_block[1].subgridspec(1, 6, width_ratios=[0.65] * 5 + [1.0],
                                     wspace=0.20)

    ax_var = fig.add_subplot(gs_top[0, 0])
    P.dimchoice.variance_explained(sweep, ax_var, n_nmf_seeds=SUPP_DIM_N_SEEDS,
                                   n_shuffles=SUPP_DIM_N_SHUFFLES)
    ax_slope = fig.add_subplot(gs_top[0, 1])
    P.dimchoice.variance_slope(sweep, ax_slope)
    # The ARI curve stops before the consensus column so it does not run under
    # the matrices.
    ax_ari = fig.add_subplot(gs_bot[0, 0:5])
    P.dimchoice.stability(sweep, ax_ari, n_nmf_seeds=SUPP_DIM_N_SEEDS)

    # The consensus-matrix insets that used to sit beside C are omitted (a reviewer
    # asked what they show, and the paper's Figure S1 can do without them).

    # The chosen number of groups, as a dotted vertical line across all three panels
    # (the paper marks it the same way; the value is also printed below the figure).
    for ax_k in (ax_var, ax_slope, ax_ari):
        ax_k.axvline(chosen, color="0.35", ls=":", lw=1.0, zorder=0)

    for ax, letter, (lx, ly) in zip((ax_var, ax_slope, ax_ari), letters,
                                    ((-0.30, 1.10), (-0.30, 1.10), (-0.13, 1.24))):
        hs.panel_label(ax, letter, x=lx, y=ly)

    print(f"{conc}: k = {sweep.k_values.tolist()} | chosen k = {chosen} | "
          f"consensus matrices not drawn | panels {'/'.join(letters)}")

save(fig, "supp_nmf_dim_choice", dpi=SUPP_DIM_DPI)
# plt.show()
# export_figure("supp_nmf_dim_choice")   # -> notebooks/outputs/supp_nmf_dim_choice.{png,pdf,svg}

# %% [markdown]
# ## 10. Figure S2 — the same structure within each mouse
#
# Figure 4 pools glomeruli from four mice. This shows the response matrix computed within each animal separately, at
# low concentration: one heatmap per mouse (S1–S4), glomeruli across and sorted into the seven groups (the numbers
# along the top, in each group's colour; the light vertical lines separate them), and odorants down in the same
# chemical-class order as Figures 3 and 4, with the classes named on the left. Brightness is ΔF/F on a log scale from
# 1 to 140, as the colourbar on the right shows.
#
# The same block structure appears in every mouse. Group 5 is missing from mice 3 and 4, which have no glomeruli in
# that group.
#
# **If you are comparing against an older version of this figure:** the source page labels the groups A–G and this one
# numbers them 1–7. The letters do not map onto the numbers in order.

# %%
# Own constants — this section inherits no bare name from anywhere.
# Two separate figures: A (per-subject heatmaps) and B (per-subject maps).
SUPP_SUBJ_DPI = 200
SUPP_SUBJ_LABEL_FS = 7
SUPP_SUBJ_SCALE = SUPP_TARGET_W / 12.5
_supp_subj_k = SUPP_SUBJ_SCALE / (HEATMAP_IN / 2.352)   # vs. this section's old scale
SUPP_SUBJ_TITLE_PAD = 22 * _supp_subj_k   # lifts 'S1'..'S4' clear of the 1-7 cluster row
# A: raw dF/F on a log scale, odorant rows in Fig 3 / Fig 4's chemical-class order.
SUPP_SUBJ_HEAT_VMIN, SUPP_SUBJ_HEAT_VMAX = 1.0, 140.0
SUPP_SUBJ_CMAP = "magma"
SUPP_SUBJ_ODOR_ORDER = CLASS_ORDER_FIG3        # "matt", the same blocks as Fig 3 and Fig 4C
SUPP_SUBJ_ODOR_WITHIN = WITHIN_CLASS_ORDER     # "matt", within each block

# Shared odorant row sort: chemical-class blocks in the Fig 3 order, so these
# panels line up with the Fig 3 matrices. The columns stay sorted by cluster.
# The row boundaries are cleared -- they were odorant-cluster lines, which no
# longer correspond to the rows.
_perm, *_ = similarity.class_sorted(d, SUPP_SUBJ_ODOR_ORDER, SUPP_SUBJ_ODOR_WITHIN)
supp_subj_order = dataclasses.replace(
    P.clustering.bicluster_order(d),
    rows=_perm,
    row_boundaries=np.array([], dtype=int),
    row_edges=np.array([0, len(_perm)]),
)

# %%
# --- Figure S2A: per-subject heatmaps ---------------------------------------
# Canvas in the original 12.5in-wide units (3.6in tall: the old A block plus its
# title room), scaled so the width is exactly SUPP_TARGET_W.
SUPP_SUBJ_W, SUPP_SUBJ_H = 12.5 * SUPP_SUBJ_SCALE, 3.6 * SUPP_SUBJ_SCALE
fig = plt.figure(figsize=(SUPP_SUBJ_W, SUPP_SUBJ_H), dpi=SUPP_SUBJ_DPI)

gs_heat = fig.add_gridspec(1, 4, left=0.13, right=0.90, top=0.767, bottom=0.14,
                           wspace=0.14)
ax_heat = [fig.add_subplot(gs_heat[0, i]) for i in range(4)]
for ax, s in zip(ax_heat, d.subjects):
    P.clustering.per_subject_heatmap(d, ax, s, conc="low",
                                     order=supp_subj_order,
                                     cmap=SUPP_SUBJ_CMAP,
                                     vmin=SUPP_SUBJ_HEAT_VMIN,
                                     vmax=SUPP_SUBJ_HEAT_VMAX,
                                     boundary_color="white")
    # The subject title and the coloured cluster numbers both live above the
    # axes and were landing on the same line. The numbers are written by the
    # panel just above the data; the title is lifted clear of them here.
    ax.set_title(ax.get_title(), fontsize=7, pad=SUPP_SUBJ_TITLE_PAD)
fig.text(0.515, 0.03, "glomerulus #", ha="center", va="bottom", fontsize=hs.FS_LABEL)

# Chemical-class strip, with the class names, left of the first heatmap, in the plotted row order.
_row_cls = d.chem_class[supp_subj_order.rows]
_cls_name = {c: cl.class_name(c) for c in dict.fromkeys(_row_cls)}
hs.class_strip(ax_heat[0], [_cls_name[c] for c in _row_cls],
               colors={_cls_name[c]: cl.class_color(c) for c in _cls_name},
               orientation="y", side="left", labels=True, fontsize=CLASS_STRIP_FONTSIZE)

# The colourbar is shortened to 70% of the panel height and centred on it, right of the last heatmap.
p0 = ax_heat[-1].get_position()
_cb_h = p0.height * 0.70
cax = fig.add_axes([p0.x1 + 0.012, p0.y0 + (p0.height - _cb_h) / 2, 0.012, _cb_h])
cb = fig.colorbar(plt.cm.ScalarMappable(
    norm=LogNorm(SUPP_SUBJ_HEAT_VMIN, SUPP_SUBJ_HEAT_VMAX),
    cmap=SUPP_SUBJ_CMAP), cax=cax)
cb.set_label("\u0394F (log)", fontsize=hs.FS_LABEL)
cb.set_ticks([1, 10, 100])
cb.set_ticklabels(["1", "10", "100"])
cb.ax.tick_params(labelsize=8)
cb.ax.minorticks_off()
cb.outline.set_visible(False)

save(fig, "supp_per_subject_heatmaps", dpi=SUPP_SUBJ_DPI)
# plt.show()
# export_figure("supp_per_subject_heatmaps")   # -> notebooks/outputs/supp_per_subject_heatmaps.{png,pdf,svg}

# %% [markdown]
# ## 11. Figure S3 — where each mouse's groups sit on the bulb
#
# For each of the four mice (rows) and each of the seven glomerular groups (columns), where that mouse's members of
# the group lie on the dorsal bulb. Members are filled circles in the group's own colour with a thin black outline;
# every glomerulus from every mouse is drawn behind as a small hollow grey circle, so all 28 maps share one frame and
# can be compared by position. It is drawn exactly like Figure 6, which pools the four mice into one map per group.
#
# Where a mouse has no glomeruli in a group — group 5 in mice 3 and 4 — the map shows a short grey dash instead.

# %%
# Own constants -- this section inherits no bare name from anywhere (S3).
# --- Figure S3: where each mouse's members of each cluster sit -------------
# Drawn exactly like Fig 6: the same page scale (so each map is Fig 6's size),
# the same marker sizes, black outline on the members and the same pale hollow
# context. Restated here, not read from §6, so each section runs on its own.
# Only the fill differs: it changes with the cluster, as Fig 6's red does not.
SUPP_TOPO_DPI = 200
SUPP_SUBJ_TOPO_SCALE = MAIN_TARGET_W / 16.0
_supp_topo_k = SUPP_SUBJ_TOPO_SCALE / FIG_SCALE
SUPP_SUBJ_CLUSTER_EDGE = "black"
SUPP_SUBJ_CLUSTER_SIZE = 62 * _supp_topo_k ** 2
SUPP_SUBJ_CLUSTER_LW = 0.6 * _supp_topo_k
SUPP_SUBJ_CONTEXT_COLOR = P.clustering.CONTEXT_COLOR   # "darkgray", as Fig 5's hollow dots
SUPP_SUBJ_CONTEXT_ALPHA = 0.5
SUPP_SUBJ_CONTEXT_SIZE = 16 * _supp_topo_k ** 2
SUPP_SUBJ_CONTEXT_LW = 0.4 * _supp_topo_k
SUPP_SUBJ_TOPO_H = 5.2      # 4 mice, each map about as tall as Fig 5's, plus the title row

supp_subj_limits = P.clustering.bulb_limits(d)      # shared spatial frame

fig = plt.figure(figsize=(MAIN_TARGET_W, SUPP_SUBJ_TOPO_H), dpi=SUPP_TOPO_DPI)
gs_topo = fig.add_gridspec(4, cl.N_CLUSTERS, left=0.07, right=0.97, top=0.92,
                           bottom=0.02, wspace=0.06, hspace=0.06)
for si, s in enumerate(d.subjects):
    for j, number in enumerate(range(1, cl.N_CLUSTERS + 1)):
        P.clustering.per_subject_topography(
            d, fig.add_subplot(gs_topo[si, j]), s, number,
            limits=supp_subj_limits,
            size=SUPP_SUBJ_CLUSTER_SIZE, edgecolor=SUPP_SUBJ_CLUSTER_EDGE,
            linewidth=SUPP_SUBJ_CLUSTER_LW,
            context_color=SUPP_SUBJ_CONTEXT_COLOR,
            context_alpha=SUPP_SUBJ_CONTEXT_ALPHA,
            context_size=SUPP_SUBJ_CONTEXT_SIZE,
            context_linewidth=SUPP_SUBJ_CONTEXT_LW,
            empty_fontsize=8,
            title=(si == 0),                     # cluster number on row 1 only
            ylabel=(f"S{s}" if j == 0 else None))

save(fig, "supp_per_subject_topography", dpi=SUPP_TOPO_DPI)
# plt.show()
# export_figure("supp_per_subject_topography")   # -> notebooks/outputs/supp_per_subject_topography.{png,pdf,svg}

# %% [markdown]
# ## 12. Figure S4 — predicted vs measured, over every responding pair
#
# The horizontal axis is what was measured and the vertical axis is what the model predicted, both as log(1 + ΔF/F);
# perfect prediction would put every point on the diagonal. Every individual glomerulus–odorant pair that responded
# at high concentration is one point — 3,405 of them — drawn as a density plot, darker where more pairs fall. *Q²* is
# printed on the panel (0.27 in the paper).
#
# The dotted horizontal line is the smallest nonzero response ever measured. Predictions below it are smaller than
# any real response, which is where the model's small nonzero guesses for non-responding pairs pile up.
#
# The paper's figure also had two averaged versions of this plot (by group and odorant, r = 0.80, and by group and
# chemical class, r = 0.89); only this unaveraged one is drawn now. The fit is loosest here, since nothing is
# averaged away.

# %%
# Own constants — this section inherits no bare name from anywhere.
# Only panel C of the old three-panel figure remains: the parity of every
# responding glomerulus x odorant pair, observed vs held-out predicted.
SUPP_PAR_C_CMAP = "RdPu"             # hexbin colormap: pink to purple, as in the paper
SUPP_PAR_W, SUPP_PAR_H, SUPP_PAR_DPI = SUPP_TARGET_W / 2, 2.4, 200

# The smallest nonzero observed high response in the whole dataset (log1p
# space, which is this panel's space): the dotted horizontal line. A predicted
# value under it is lower than any response the experiment ever measured.
supp_par_floor = pls.observed_high_floor(m, space="log")

fig, ax = plt.subplots(1, 1, figsize=(SUPP_PAR_W, SUPP_PAR_H),
                       dpi=SUPP_PAR_DPI, constrained_layout=True)
# Bare call = the pinned framing at its defaults, which is the draft's panel.
supp_par_hb = P.pls.responding_parity(
    ax, m, cbar=True, title="all responding pairs", title_fontsize=11.5,
    cmap=SUPP_PAR_C_CMAP, floor=supp_par_floor)

print(f"responding pairs   Q2 = {pls.baseline_q2(m):.3f}  "
      f"(n = {pls.responding_pairs(m)[0].size:,}); "
      "the draft page labels this same value R2")
print(f"dotted line at the smallest nonzero observed response: "
      f"{supp_par_floor:.4f} (log1p dF/F)")

save(fig, "supp_parity", dpi=SUPP_PAR_DPI)
# plt.show()
# export_figure("supp_parity")   # -> notebooks/outputs/supp_parity.{png,pdf,svg}

# %% [markdown]
# ---
#
# **Note.** Every figure is now in `notebooks/outputs/` as PNG, PDF and SVG, and
# a second, SVG-only copy is in `joel_figures/` at the repo root.
#
# The PDFs and SVGs have transparent backgrounds and live, editable text — open
# one in Illustrator or Inkscape and the labels are still text objects you can
# restyle or move.
#
# To put copies somewhere else, uncomment an `export_figure(...)` line — there
# is one under each figure, and a general one near the top of the notebook.
#
# **If you want to change something:** a panel's contents live in
# `highconc/panels/`, a figure's layout lives in the cell that builds it, and
# fonts and colours live in `highconc/style.py`.
