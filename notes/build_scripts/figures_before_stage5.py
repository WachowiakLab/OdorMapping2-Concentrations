# %% [markdown] tags=["human"]
# # High-concentration odor mapping — the figures
#
# This notebook draws every figure in the paper, except for a small handful
# that only Matt has the data for. Run it top to bottom and it rebuilds all of
# them into `notebooks/outputs/` as PNG, PDF and SVG.
#
# ### The basic story the figures tell
#
# - Raising concentration changes which odorants evoke similar glomerular
#   response patterns, not only how strongly they respond (**Fig 3**).
# - Glomeruli group by shared odorant preference, and tuning broadens with
#   concentration (**Fig 4**); the groups map onto the bulb (**Fig 5**).
# - A linear model fitted across all odorants predicts held-out
#   high-concentration responses (**Fig 6**).
# - The broadening is selective and asymmetric: particular groups take on the
#   tuning of particular others, and not the reverse (**Fig 7**).
#
# ### What's in here
#
# Section numbers are the paper's figure numbers, so §3 draws Figure 3.
#
# | § | Figure | In one line |
# |---|---|---|
# | 1 | setup | Loads the data and fits the model. Run this first. |
# | 2 | **Fig 2** | Odorant tuning across Burton's own glomeruli, both concentrations |
# | 3 | **Fig 3** | Odorant–odorant response similarity at both concentrations |
# | 4 | **Fig 4** | The seven glomerular groups, and tuning by chemical class |
# | 5 | **Fig 5** | Where those groups sit on the bulb |
# | 6 | **Fig 6** | A linear model predicting the low → high change |
# | 7 | **Fig 7** | Directional transfer of tuning between groups |
# | 8 | **Fig S1** | How the number of groups was chosen |
# | 9 | **Fig S2** | The same structure within each mouse |
# | 10 | **Fig S3** | Direction of change through tuning space |
# | 11 | **Fig S4** | Predicted vs measured, at three levels of averaging |
#
# §2 used to be empty — a model cartoon that lived there is drawn separately
# in Illustrator now — but Figure 2 (added 2026-09) has since filled it, so
# §*n* still means Figure *n* with no gap left in the numbering.
#
# A few panels in the paper are Matt's — the bulb photographs, the per-cluster
# spectra in Fig 4E, and Fig 7's explanatory header — and are not drawn here.
#
# ### How to change something
#
# - **A panel's contents** (what's plotted, its colours): edit the function in
#   `highconc/panels/`. Every figure using it follows.
# - **A figure's layout** (grid, colourbars, panel letters): edit the cell that
#   builds it, right here in the notebook.
# - **Fonts, colours, export settings**: `highconc/style.py`, once, for
#   everything. No cell sets its own.
#
# ### Two things that will bite you
#
# 1. **Run it in order.** §1 loads the data and fits the model that later
#    sections share. Running a later section on a fresh kernel refits it, which
#    takes about a minute.
# 2. **Some results are cached on disk and are not recomputed** — the model fit
#    and the dimensionality sweeps in §8. If a cache is missing the notebook
#    stops with an error instead of quietly recomputing, because a recomputed
#    number might not match the one in the paper.
#
# Glomerular groups are numbered **1–7** throughout. Older figures and drafts
# label them A–G; the letters do not map onto the numbers in order, so don't
# assume cluster C is cluster 3. `highconc.clusters` holds the translation.

# %% [markdown] tags=["llm"]
# # High-concentration odor-mapping — every figure, one notebook
#
# This notebook draws every figure in the manuscript. It is the only place
# composition lives: the analysis is in the `highconc` package, each panel is a
# function in `highconc.panels.*` that takes an `ax` and draws into it, and the
# grids, colorbars and panel letters below are the part you are meant to edit.
#
# To change a panel's *content*, edit `highconc/panels/<module>.py`.
# To change a figure's *layout*, edit the cell that builds it.
# To change the type scale, colours, or export contract, edit `highconc/style.py` —
# no cell here sets rcParams.
#
# **Run it top to bottom.** Later sections reuse the model fitted in section 1;
# running a PLS section on its own re-fits it (about a minute) and then caches.
#
# ## Contents
#
# This notebook draws every figure in the manuscript except the handful only
# Matt has the data for — the bulb photographs, Fig 4E, Fig 5B–D and Fig 7's
# explanatory header. The numbering below follows the manuscript: section 2
# is Figure 2, sections 3–7 are Figures 3–7, and 8–11 are Supplemental
# Figures S1–S4. `dev/CONCORDANCE.yaml` is the manifest that fixes this list,
# and `dev/concordance.py --check` enforces that the notebook's output set is
# exactly these stems.
#
# | § | Figure | Stem | Panels |
# |---|---|---|---|
# | 1 | Setup — style, paths, shared data and model fit | — | — |
# | 2 | **Figure 2** — odorant x glomerulus response, Burton's own dataset, both concentrations | `fig2_response_matrix` | `panels.similarity` |
# | 3 | **Figure 3** — correlational structure vs. concentration, and its threshold supplement | `fig3_correlational_structure`, `fig3_supp_threshold_robustness` | `panels.similarity`, `panels.recruitment` |
# | 4 | **Figure 4** — change in tuning with concentration, by NMF cluster | `fig4_tuning_by_cluster` | `panels.clustering` |
# | 5 | **Figure 5** — spatial mapping of glomeruli by cluster | `fig5_spatial_by_cluster` | `panels.clustering` |
# | 6 | **Figure 6** — PLS-R: predicted spectra, recruitment matrix, retuning hypotheses, W | `fig6_pls` | `panels.schematic`, `panels.similarity`, `panels.recruitment`, `panels.pls` |
# | 7 | **Figure 7** — motif transfer between glomerular groups, with the chemical-class control | `fig7_motif_transfer` | `panels.transfer` |
# | 8 | **Figure S1** — NMF dimensionality choice, 1x and 25x as one figure | `supp_nmf_dim_choice` | `panels.dimchoice` |
# | 9 | **Figure S2** — per-subject reproducibility, heatmaps + topography | `supp_per_subject` | `panels.clustering` |
# | 10 | **Figure S3** — concentration flow through tuning space, flow + coherence | `supp_tuning_flow` | `panels.flow` |
# | 11 | **Figure S4** — parity at three granularities: cluster×odorant, cluster×class, all responding pairs | `supp_parity` | `panels.pls` |
#
# **§2 used to be empty** — the numbering was left with a gap so that §*n*
# drew Figure *n*, and the model cartoon that used to sit at Figure 2 is drawn
# separately in Illustrator — but Figure 2 (added 2026-09) has since filled
# it, so there is no longer a gap. `panels.schematic.cartoon_map` and
# `.vector_strip` are still exported by the package and are unused here;
# `weight_glyph` is live, since `retuning_glyph` builds on it for Figure 6E.
#
# Figure 4E, 5B–D and Figure 7's explanatory header are drawn elsewhere and
# are not built here. `CONCORDANCE.yaml` records each as `skipped` with its
# reason. Figure 3A used to be on this list too — it needed a "most sensitive"
# condition this dataset didn't have — until `data/burton_et_al/` supplied it
# (2026-09-02); `CONCORDANCE.yaml` now records that panel `done`.
#
# Each section opens with a markdown cell describing what the figure shows and
# how it was made. That text is the figure legend; it sits next to the code so
# the two cannot drift apart.
#
# ## Two things not to break
#
# 1. **The NMF dimensionality sweeps are cached** (`figures/nmf_dim_choice_*/data/cache/*.npz`)
#    and take minutes to recompute. `dimsweep.load_cached` raises on a miss
#    rather than silently recomputing — leave it that way.
# 2. **Nothing is refitted at draw time.** `highconc.pls.load_weights`,
#    `.component_selection` and `.lv_loadings` read the canonical run
#    pls-explorer wrote to disk, and raise if it is missing rather than
#    quietly producing different numbers.
#
# Cluster numbering is Joel's **1–7** everywhere. The legacy A–G letters survive
# only as a mapping in `highconc.clusters`; a letter does not tell you the number.

# %% [markdown] tags=["human"]
# ## 1. Setup
#
# Loads the imaging data and fits the model that most of the later sections use.
# Run this before anything else. It takes about a minute the first time and is
# then cached, so re-running costs nothing.
#
# The cell after the imports also explains **how to export figures** — the short
# version is that a normal run already writes PNG, PDF and SVG of everything.

# %% [markdown] tags=["llm"]
# ## 1. Setup
#
# One style call, one data load, one model fit — shared by every section below.

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
#: `HEATMAP_IN` so that Fig 2 A-C, Fig 3 A-D/H, Fig 4 A-B, Fig 6 D, Fig 7 A,
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

# %% [markdown] tags=["human"]
# ## 2. Figure 2 — odorant tuning across Burton's own glomeruli, three conditions
#
# Three response matrices over the same 59 odorants, this time against
# Burton et al.'s own glomeruli rather than this dataset's 376 — panel A at
# high concentration, B at low, C at the "most sensitive" condition (a
# different, smaller subject/preparation, `omp114L`/`omp8x`, 76 glomeruli
# instead of the other two panels' 113). Rows are grouped by chemical class
# the same way Figure 3's panels are, so the two figures read the same way at
# a glance.
#
# - **A, B, C** — each cell is one odorant-glomerulus pair; brighter means a
#   bigger response. All three panels share one colour scale. Odorants are
#   grouped down the rows by chemical class, with the coloured strip and
#   class names on the left — the same convention as Figure 3's chemical-class
#   strip, just turned 90° onto the row axis since odorants are rows here, not
#   columns. Columns (glomeruli) carry no identity or ordering information
#   supplied with the data and are left unlabeled. Panel C is narrower than A
#   and B only because its source glomerulus population is smaller (76 vs.
#   113) — not a scaling artifact.
#
# High concentration is visibly denser than low — printed under the figure.

# %% [markdown] tags=["llm"]
# ## 2. Figure 2 — odorant x glomerulus response, Burton's own dataset
#
# **What it shows.** The same 59 odorants against Burton's own glomeruli — not
# this repo's 376 — at three conditions. Panel A is
# `data/burton_et_al/sortedmatrix.csv` ("high"), panel B is
# `sortedmatrixlo.csv` ("low"), panel C is
# `sortedmatrix_omp114L_omp8x.csv` ("most sensitive" — a different
# subject/preparation than A/B, hence its 76 rather than 113 glomeruli).
# All three were supplied 2026-09 with no header or odorant labels, so rows
# are matched by position against `odorant_order.txt`, the same convention
# `most_sensitive_correlation` (Figure 3A) uses.
#
# **How to read it.** Rows are the 59 odorants, permuted into the same
# class-block order Figure 3 uses (`class_sorted`, `order="matt"`,
# `within="burton"`), with the chemical-class strip carried on the row axis
# instead of the column axis — `highconc.style.class_strip`, the same helper
# Figure 3 uses, just with `orientation="y"`. Columns are Burton's own
# glomeruli in whatever order each source file ships them; there is no
# identity or ordering information for them, so none is imposed, and they are
# left unlabeled. Panel C's narrower column count (76 vs. A/B's 113) reflects
# that it comes from a different, smaller population, not a shared axis that
# was cropped.
#
# **Parameters.** `FIG2_ORDER`/`FIG2_WITHIN` are this section's own
# row-ordering parameters, independent of Figure 3's `CLASS_ORDER_FIG3`/
# `WITHIN_CLASS_ORDER` even though the legal values and defaults are the
# same — every section owns its constants, and Figure 2 doesn't depend on
# Figure 3 having been run first.
#
# **Section number.** This is the first figure this notebook draws at §2 —
# previously the gap was deliberate (a model cartoon that used to live there
# is now drawn separately in Illustrator), but §*n* still means Figure *n*,
# and there is now a real Figure 2 to draw.

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
plt.show()

# %% [markdown] tags=["human"]
# ## 3. Figure 3 — odorant–odorant response similarity at both concentrations
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# Panels A, B and C ask, for every pair of the 59 odorants, how similar the
# two patterns of glomerular activity are — A at the "most sensitive"
# condition, B at low concentration, C at high. Average similarity between
# different odorants is 0.03 at low and 0.10 at high, and the low/high
# matrices agree with each other at a rank correlation of 0.55.
#
# Panel D is a chemical-structure similarity for the same 59 odorants, which
# does not change with concentration. Chemical similarity tracks response
# similarity more closely at high concentration than at low (rank correlation
# 0.21 → 0.31).
#
# **Reading order is D, A, B, C, left to right** — the fixed chemical ruler
# shown first, then the three concentration-dependent conditions it's judged
# against. A panel's letter stays tied to its content regardless of where
# it's drawn (`FIG3_PANEL_ORDER`), so D is still the chemical-similarity panel
# even though it's leftmost.
#
# **Reading the panels.**
#
# - **A, B, C** — every cell is one odorant pair; brighter means their
#   glomerular response patterns are more alike. All three share one colour
#   scale, so differences between them are not a scaling artifact. Odorants
#   are grouped along both edges by chemical class, with the aromatics split
#   into three finer groups (Ar. Ketone, Ar. Ketone/Ester, Phenol/Methoxy —
#   see below). The grey diagonal is each odorant compared with itself, which
#   is 1 by definition; four odorants (D15, E1, P1, P4) have no data at the
#   "most sensitive" condition, so panel A reads them at the colour floor
#   rather than as a gap. Only **A** carries the chemical-class names under
#   its strip — D, B and C repeat the same coloured strip with no text, so
#   the row of four doesn't say "acids, aldehydes, ..." four times over.
# - **D** — the same pairs, scored on chemical structure (MACCS keys, cosine)
#   instead of on response. Carries the same aromatics breakdown as A–C.
# - **E, F, G** — three individual glomeruli, one per panel. Each bar is one
#   odorant: the low-concentration response points up, the high-concentration
#   response points down. Bar colours are chemical class — the plain, ten-class
#   version; the aromatics split is a panel A–D thing only. The downward bars
#   are magnitudes, not negative numbers.
# - **H** — recruitment. Each row is an odorant that drove a glomerulus at low
#   concentration; each column an odorant that newly drives it at high. The
#   cell is the fraction of that row's glomeruli recruited. Grey is undefined,
#   not zero.
#
# **Panel A's data.** `data/burton_et_al/most_sensitive.csv` and
# `odorant_order.txt`, supplied 2026-09-02 — the "most sensitive" condition
# this dataset didn't otherwise have. The odorant order in that file is now
# also the order every panel here uses (A–H), in place of the odorant table's
# own order.
#
# **The three example glomeruli** are numbers 173, 211 and 97, and Figure 6
# shows the same three. Nothing on record said which glomeruli the original
# panels used, so they were identified by measuring the bars in the draft PDF
# and matching them against all 376; they match to within the PDF's own
# rounding.
#
# **Panel letters follow the figure PDF, not the manuscript text**, which
# disagree — the text calls the recruitment matrix 3E, the PDF calls it 3H.

# %% [markdown] tags=["llm"]
# ## 3. Figure 3 — correlational structure of odor representations vs. concentration
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# **What it shows.** That raising concentration does not just make responses
# bigger, it makes the *relationships between odorants* different. Panels B and
# C are the same 59 × 59 odorant × odorant correlation matrix at the two
# concentrations, and the high-concentration matrix is visibly denser: the mean
# off-diagonal correlation goes from 0.03 to 0.10, and the two matrices agree
# with each other only at Spearman 0.55. Panel A is the same matrix at a third,
# "most sensitive" condition, supplied 2026-09-02 in
# `data/burton_et_al/most_sensitive.csv` — previously out of scope here because
# this dataset didn't have that condition. Panel D is the chemical similarity of
# the same 59 odorants, which of course does not change — so it is the fixed
# ruler the other three are measured against. Chemistry predicts response
# similarity **better at high concentration than at low** (Spearman 0.21 →
# 0.31).
#
# **Panel lettering follows the figure PDF, not the manuscript text.** The two
# disagree: the text calls the recruitment matrix "Fig. 3E", the figure PDF
# puts it at 3H with E/F/G as the response spectra. The figure PDF is followed
# here, and the discrepancy is noted so it can be reconciled in the text.
#
# **How to read each panel.**
#
# - **A, B, C** — odorant × odorant Pearson correlation between response maps
#   (each odorant is a 376-glomerulus vector), at the "most sensitive"
#   condition, at 1x and at 25x. All three are pinned to the same 0–1 scale,
#   and they must stay pinned: A and B are sparse and C is not, so autoscaling
#   them separately would erase the comparison the figure exists to make. The
#   grey diagonal is *not a measurement* — self-similarity is 1 by
#   construction — while a pale cell is a real, near-zero correlation. Panel A
#   ships with two extra, unnamed rows/columns beyond the 59 in
#   `odorant_order.txt`; `most_sensitive_correlation` drops them. Four odorants
#   (D15, E1, P1, P4) have no data at this condition — filled at 0.0 rather
#   than left NaN, so they read at the colour floor instead of as a hole, while
#   the true diagonal is still re-masked to NaN like every other panel.
# - **D** — MACCS cosine similarity from the odorants' SMILES. Cosine rather
#   than Tanimoto so that it is the same normalised-inner-product quantity as
#   the correlations beside it. The 59 × 59 matrix is cached to `data/chem/`,
#   so drawing this needs no RDKit; only rebuilding the cache does.
# - **Draw order is D, A, B, C** (`FIG3_PANEL_ORDER`), not alphabetical: the
#   concentration-invariant "ruler" panel leads, then the three conditions
#   measured against it. Letters index content, not position — `mats_by_letter`
#   keys the row's data by its established letter, and `FIG3_PANEL_ORDER` only
#   picks which column of the gridspec each one lands in.
# - **Strip labels are on A only** (`strip_labels=(L == "A")` in the draw
#   loop). Every panel still gets the coloured x-axis class strip -- all four
#   share one class order, so the strip is structurally identical across the
#   row -- but repeating the ten class names under D, B and C too would be
#   pure duplication.
# - **A–D's aromatics breakdown.** The single "R" (aromatics) class is split
#   into three display subclasses for the axis strip and boundaries only —
#   `highconc.clusters.AROMATIC_SUBCLASS_BY_CODE`, applied via
#   `similarity.split_block`: Ar. Ketone (R27), Ar. Ketone/Ester (R11, R12,
#   R21, R31, R3, R29, R40), Phenol/Methoxy (R10, R23), coloured `#a32e49` /
#   `#00bbb4` / `#d59537`. It touches nothing but these four panels — Fig 4/6
#   and panels E–H still see plain "R". Toggle with `FIG3_SPLIT_AROMATICS`.
# - **A–D's diagonal boxes**, drawn the same way Figure 4A's `cluster_boxes`
#   are (`correlation_matrix`'s `box_boundaries`, an unfilled square between
#   each pair of consecutive edges): five groups, not the full ten/twelve
#   chemical-class blocks — acids through alcohols, aromatics+furans, amines,
#   pyrazines, thiazoles. `FIG3_BOX_BOUNDARIES` is the cut points; a different
#   `CLASS_ORDER_FIG3` or `FIG3_SPLIT_AROMATICS` setting changes what block
#   sizes those land on, so it isn't derived automatically.
# - **E, F, G** — three individual glomeruli, low plotted up and high plotted
#   down about a shared axis, bars coloured by chemical class (the plain
#   ten-class version, not the A–D aromatics split). The lower half is
#   a magnitude, not a decrease: the tick labels are made positive again. Each
#   panel autoscales, as the draft figure does — see `SPECTRUM_SCALE`.
# - **H** — the recruitment matrix, observed. Rows are the low-concentration
#   ("source") odorant, columns the high-concentration ("target") odorant that
#   recruits a glomerulus tuned to the source, and each row is rescaled to its
#   own [0, 1] range, so a cell reads as that source odorant's *relative*
#   pull toward the target, not an absolute probability. Grey is undefined
#   (a row with no signal to rescale), not zero. §6 draws this same matrix
#   beside its held-out prediction (Figure 6E); here only the observed half
#   is Figure 3's.
#
# **The odorant order.** `WITHIN_CLASS_ORDER = "burton"` — the full 59-odorant
# order in `data/burton_et_al/odorant_order.txt`, covering every class block
# rather than only the pyrazine reshuffle `"matt"` makes. Its block sequence is
# `CLASS_ORDER_MATT`'s, so it pairs with `CLASS_ORDER_FIG3 = "matt"` unchanged.
# §6 inherits this directly; Fig 3H/Fig 6D/E's recruitment matrices honour it
# too via `recruitment_matrix`'s new `within`/`codes` arguments. Figure 4, 5,
# S1 and S2 are deliberately untouched: Fig 4/S2's heatmaps are ordered by NMF
# bicluster structure, not chemical class — reordering them by class would
# erase the block-diagonal structure those figures exist to show — and Fig 5/S1
# have no odorant axis at all.
#
# **Which three glomeruli, and how they were identified.** The draft figure's
# E/F/G panels are three individual glomeruli, and no record states which. They
# were recovered by measurement: the panels were reduced to bar geometry from
# the figure PDF's vector content, and all 376 glomeruli were scored against
# them. Glomeruli **173, 211 and 97** reproduce every bar under a single
# per-panel scale factor, to within 0.008 pt. As an independent check, the
# scale fitted to panel E's bar *heights* predicts that panel's axis *tick
# spacing* — 11.98 pt for a 20-unit tick, against 11.980 measured — and nothing
# in the fit used the ticks. They live in `similarity.EG_GLOMERULI` because
# **Figure 6A–C must draw the same three**.
#
# **Parameters.**
#
# - `CLASS_ORDER_FIG3` — the chemical-class block order. `"matt"`
#   (C D E K L R F N P S) is the draft figure's order, verified against the
#   figure PDF rather than inferred. `"matt_esters_first"` is the variant with
#   esters immediately after acids, annotated on that page; `"table"` is the
#   odorant table's own order.
# - `WITHIN_CLASS_ORDER` — order *inside* a block. `"table"` is the odorant
#   table's own order; `"matt"` additionally reshuffles just the pyrazine
#   block to match the draft figure; `"burton"` is the full 59-odorant order
#   from `data/burton_et_al/odorant_order.txt`, covering every block. Default
#   is now `"burton"`.
# - `FIG3_SPLIT_AROMATICS` — whether panels A-D split the "R" (aromatics)
#   class into Ar. Ketone / Ar. Ketone/Ester / Phenol/Methoxy for the axis
#   strip and boundaries. `True` by default; set `False` to see plain "R"
#   like every other figure.
# - `SPECTRUM_SCALE` — `"per-panel"` (the draft's; each spectrum autoscales) or
#   `"shared"` (one y-limit across E/F/G). Per-panel reads shape; shared reads
#   magnitude. Note that per-panel means E, F and G are at roughly 50, 65 and
#   110 units of the same quantity drawn the same height.
# - `ROW_FRAC_FIG3`, `LO_FRAC_FIG3` — panel H's tally thresholds, and the same
#   values §6 restates for Figure 6E. Each glomerulus is first peak-normalised
#   to [0, 1] separately per concentration; a glomerulus counts as
#   "recruited" by an odorant when it sits below `LO_FRAC_FIG3` (0.20) of its
#   own peak at low concentration and at or above `ROW_FRAC_FIG3` (0.50) at
#   high. The same 0.50 also decides which odorants are "effective" for a
#   glomerulus at low concentration — the rows a glomerulus's recruitment
#   gets tallied into. Ported directly from the lab's original MATLAB; see
#   `highconc/recruitment.py` for the exact algorithm.

# %%
CLASS_ORDER_FIG3 = "matt"      # "table" | "matt" | "matt_esters_first"
WITHIN_CLASS_ORDER = "matt"    # "table" | "matt" (pyrazine only) | "burton" (full order)
FIG3_SPLIT_AROMATICS = True    # split "R" into 3 display subclasses on panels A-D
SPECTRUM_SCALE = "per-panel"   # "per-panel" | "shared"
EG_GLOMERULI = similarity.EG_GLOMERULI     # panel F's first three example glomeruli
# Fourth example glomerulus (panel F), by dataset glomerulus id like EG_GLOMERULI.
# Kept out of EG_GLOMERULI, which other code imports.
EXTRA_GLOM = 155

# Odorant CODEs (the leading token of `Odorant`), positionally aligned with
# `d`'s native odorant order -- `recruitment_matrix`'s `within` needs these,
# and §6 reuses this same array for its own recruitment call for the same
# reason it reuses CLASS_ORDER_FIG3/WITHIN_CLASS_ORDER.
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
# `threshold_effodor` in the lab's original algorithm -- one shared fraction,
# used both to call a glomerulus "tuned" to a low odorant and, on the high
# side, to call it "recruited". `LO_FRAC_FIG3` is the separate, lower bar for
# "not really active at low" that recruitment requires.
ROW_FRAC_FIG3 = 0.20
LO_FRAC_FIG3 = 0.20

rmx3 = rc.recruitment_matrix(m, level="odorant", lo_threshold=LO_FRAC_FIG3,
                             frac=ROW_FRAC_FIG3, order=CLASS_ORDER_FIG3,
                             within=WITHIN_CLASS_ORDER, codes=ODORANT_CODES)
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
                              strip_labels=True)
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
plt.show()
# export_figure("fig3_correlational_structure")   # -> notebooks/outputs/fig3_correlational_structure.{png,pdf,svg}
# export_figure("fig3_correlational_structure", formats=("pdf",), outdir="~/Desktop")

# %% [markdown] tags=["human"]
# ### Supplement to Figure 3 — the recruitment matrix at five cutoffs
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# Panel H depends on one choice: how strongly an odorant must drive a
# glomerulus before that glomerulus counts as assigned to it.
#
# This supplement redraws the matrix at five values of that cutoff. The block
# structure and the observed-vs-predicted correlation change little across
# them.
#
# In the paper's numbering this is its own figure rather than a panel of
# Figure 3, which is why it saves under its own name.

# %% [markdown] tags=["llm"]
# ### Supplement to Figure 3 — recruitment-matrix robustness to threshold
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# Panel H depends on one choice: what fraction of its own maximum an odorant
# must drive a glomerulus to before that glomerulus counts as "assigned" to it.
# This supplement draws the matrix at five values of that fraction. The block
# structure and the observed-vs-predicted correlation barely move, so the
# headline matrix is not an artifact of where the rule is cut.
#
# It is its own figure in the manuscript's numbering, not a panel of Figure 3,
# which is
# why it gets its own stem.

# %%
SWEEP_FRACS_FIG3 = (0.20, 0.35, 0.50, 0.65, 0.80)

sweep3 = rc.matrix_sweep(m, fracs=SWEEP_FRACS_FIG3, level="odorant",
                         lo_threshold=LO_FRAC_FIG3, order=CLASS_ORDER_FIG3)
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
plt.show()
# export_figure("fig3_supp_threshold_robustness")   # -> notebooks/outputs/fig3_supp_threshold_robustness.{png,pdf,svg}
# export_figure("fig3_supp_threshold_robustness", formats=("pdf",), outdir="~/Desktop")

# %%
for s in sweep3:
    print(f"frac={s.frac:.2f}  empty rows={s.n_empty:2d}  "
          f"mean row n={s.n_rows.mean():5.1f}  r_offdiag={s.r_offdiag:.3f}")

# %% [markdown] tags=["human"]
# ## 4. Figure 4 — the seven glomerular groups, and tuning by chemical class
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# Every panel is the same 376 glomeruli, shown four ways.
#
# - **A** — every pair of glomeruli, scored on how similar their odorant
#   preferences are, ordered by group. Within-group correlation runs 0.12–0.51
#   against a between-group mean of −0.03. (Grey diagonal = each glomerulus
#   versus itself.) A green box outlines each of the seven diagonal
#   within-group squares.
# - **B** — the response matrix behind that: odorants down, glomeruli across.
#   Brightness is now each glomerulus's response **relative to its own
#   strongest odorant** (0-1, linear), not an absolute ΔF/F on a log scale —
#   so this panel answers "which odorants drive this glomerulus, relative to
#   its own best one," not "how big is this response in absolute terms."
#   Glomerulus columns are still sorted by group (the vertical lines mark
#   where one group's block ends and the next begins); odorant rows run in
#   the plain order Figure 3/6 use instead of being grouped with the columns,
#   so a group's block of columns no longer lines up with a matching block of
#   rows — there's no shaded diagonal square any more, just the vertical
#   group lines. Each of the ten chemical-class row blocks is labelled with
#   its name, the same ten names panel D uses. Most glomerulus–odorant pairs
#   do not respond, so most of the panel is empty; a response below the same
#   floor used elsewhere (1 ΔF/F) still falls through to the grey background,
#   even after normalising — so a glomerulus whose best response never really
#   clears the noise floor reads as empty, not as a falsely "bright" 1.0.
# - **C** — the same data one mouse at a time, for two of the groups, sharing
#   the odorant sort so the panels line up row for row.
# - **D** — rows are chemical classes, columns the seven groups. Each cell
#   counts the glomeruli in that group with a nonzero response to at least
#   one odorant of that class, then scales the column so its top chemical
#   class reads as 1 and the rest are a fraction of it. Low concentration on
#   top, high underneath, each with its own colour scale.
#
# **Which two groups are in panel C.** Group 1, which responds to a mix of
# acids and aldehydes, and group 3, which responds almost only to furans. The
# paper's caption describes this panel as contrasting a mixed group with a
# narrowly selective one, and these two are the clearest such pair that is
# present in all four mice.
#
# **Note for the caption:** it currently names amines as the selective group
# (group 5). Group 5 has no glomeruli at all in mice 3 and 4, so that version
# of the panel had two empty columns; group 3 was substituted. The caption
# needs updating to match, or the panel reverting.
#
# **Panel E is Matt's**, not drawn here.

# %% [markdown] tags=["llm"]
# ## 4. Figure 4 — change in tuning with concentration by NMF clusters
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# **What it shows.** That the seven NMF clusters are real structure and not a
# sorting artifact, and that raising concentration broadens what each of them
# responds to. Every panel is the same 376 glomeruli looked at a different way:
# as a correlation matrix (A), as the raw response matrix (B), one subject at a
# time (C), and summarised by chemical class (D).
#
# **How to read each panel.**
#
# - **A** — glomerulus × glomerulus Pearson correlation of tuning curves at
#   1x. Each glomerulus is a 59-odorant vector; the matrix is permuted into
#   cluster order, so the seven bright squares on the diagonal *are* the
#   clusters. It is Figure 3B's picture one axis over, and it is on the same
#   0–1 scale for that reason. Within-cluster correlation runs 0.12–0.51
#   against a between-cluster mean of −0.03, which is the panel's whole claim.
#   The grey diagonal is not a measurement. Each of the seven diagonal blocks
#   gets an unfilled green box (`cluster_boxes=True` in
#   `glomerulus_correlation_matrix`) so the block-diagonal structure the claim
#   depends on is visible without the full boundary grid (`boundaries=False`
#   here, separate from `cluster_boxes` and independently toggleable).
# - **B** — the response heatmap: odorants down, glomeruli across.
#   `FIG4B_COL_NORMALIZE=True` (the current default) peak-normalises each
#   glomerulus (column) to [0, 1] by its own maximum response, on the raw
#   dF/F values, and draws on a **linear** scale (`FIG4B_NORM_VMIN`/
#   `FIG4B_NORM_VMAX`, `0`-`1`) instead of log — deliberately paired, since a
#   [0, 1] quantity has no absolute floor for a log scale's dynamic range to
#   do anything with. Reading it: brightness is "how close is this response
#   to this glomerulus's own best odorant," not an absolute ΔF/F, so two
#   glomeruli that respond very differently in magnitude can still both show
#   a bright cell for their respective best odorant. Set
#   `FIG4B_COL_NORMALIZE=False` to go back to raw ΔF/F on a log scale
#   (`FIG4_VMIN`/`FIG4_VMAX`), matching panel C's own scale as it did before.
#   Either way the background is drawn as an explicit grey rectangle
#   (`empty_color`), not the axes facecolor, so it renders the same way
#   regardless of backend, and either way cells below a colour floor fall
#   through to it -- but what that floor is made of differs. On the log/raw
#   scale it's `vmin` acting directly on the displayed dF/F. On the
#   linear/normalised scale, `vmin=0.0` matches the data's own floor exactly
#   (nothing is ever "below" 0 after normalising), so normalising alone would
#   erase the floor entirely -- a glomerulus whose peak response never
#   clears the noise floor would otherwise show a misleadingly "bright" 1.0
#   for its own weak best odorant. `FIG4B_RAW_FLOOR` (the same `FIG4_VMIN=1`
#   used by the log view) restores it: any cell whose *raw*, pre-normalised
#   dF/F was below that floor is greyed out regardless of its normalised
#   value, so a column that never clears the floor greys out completely
#   rather than lighting up relative to its own noise.
#   **Rows and columns are no longer sorted the same way.** Columns
#   (glomeruli) are always cluster-sorted — the vertical lines
#   (`col_boundaries`) mark where one glomerulus cluster's block ends and the
#   next begins, unconditionally. Rows (odorants) follow `FIG4B_ROW_ORDER`
#   (default `"odorant"`): Figure 3/6's own odorant order
#   (`FIG4_CHEM_ORDER`/`bicluster_order`'s `within_order="burton"`), with no
#   reference to cluster membership at all — so a cluster's odorant rows are
#   scattered through the panel rather than forming a block, and the lighter
#   "paired block" rectangles that highlight the diagonal in the bi-clustered
#   view (`diagonal_blocks=True`) are switched off (`diagonal_blocks=
#   (FIG4B_ROW_ORDER == "cluster")`), since they would otherwise shade the
#   wrong cells entirely. Thin white lines still mark each chemical class's
#   own row boundary, and now read as one clean block per class rather than a
#   scattered set of sub-block divisions, since rows are chemical-class-first
#   throughout. Set `FIG4B_ROW_ORDER = "cluster"` to restore the earlier
#   bi-clustered view (rows sorted by cluster first, `FIG4B_SORT_BY_CLASS`
#   only then choosing whether chemical class breaks ties *within* a
#   cluster's block, and the diagonal highlight back on). Panel C is
#   unaffected either way — it reads its own, always cluster-ordered `bic`,
#   never `FIG4B_ROW_ORDER`'s `bic_b`.
#   `class_names=True` writes each contiguous chemical-class row-run's name
#   (`CHEM_CLASS_NAMES`, the same ten panel D labels its own rows with) in
#   the middle of its run, in place of a numeric y-tick — under
#   `FIG4B_ROW_ORDER="odorant"` that's exactly ten labels, one per class;
#   under `"cluster"` a class split across cluster blocks gets one label per
#   block, so a name can repeat.
# - **C** — the same data broken out by mouse, for two clusters: one glomerular
#   cluster per group, four subjects side by side, sharing the odorant sort so
#   the panels line up row for row. The block repeats within each animal, so it
#   is not a product of pooling four fields of view.
# - **D** — chemical class × cluster, from `class_cluster_glom_composition`.
#   Each cell counts the **distinct glomeruli** in that cluster with a
#   nonzero response to at least one odorant of that class (a glomerulus
#   responding to several odorants of a class still counts once), then
#   scales each column by its own maximum: the chemical class reaching the
#   most responding glomeruli in a cluster is 1, the rest read as a fraction
#   of it. Columns do not sum to 1 — this is a relative magnitude within the
#   cluster, not a share of a whole. Rows run in `FIG4_CHEM_ORDER` (§3's
#   broad chemical-class order), not alphabetically. Low above, high below,
#   **each on its own colour scale** (`cluster_summary`'s default `vmax=None`
#   autoscales per call, and the two calls below are no longer pinned to a
#   shared ceiling).
# - **E** is deliberately absent — the per-cluster tuning spectra are Matt's
#   panel, like Figure 5B–D.
#
# **Which two clusters in C, and why those.** The panel draws **clusters 1 and
# 3**, and this departs from the draft caption, which describes the pair by
# content: one "characterized by mixed selectivity for acids and aldehydes",
# the other showing "near complete selectivity for amines". Read off panel D's
# numbers, those are cluster 1 and cluster 5.
#
# **Cluster 5 has no glomeruli in subjects 3 and 4** (19, 11, 0, 0), so drawing
# it leaves two of the eight sub-panels empty. It is the only cluster missing
# from any subject; the other six are present in all four.
#
# Cluster 3 replaces it on the caption's own criterion — *mixed* selectivity
# against *near-complete* selectivity for one chemical class. Of the clusters
# present in all four subjects, **cluster 3 has the sharpest single-class
# selectivity**: its top class takes 0.75 of its total activation (furans
# 0.188, next 0.045), against 0.69 for cluster 4, 0.67 for cluster 7, and 0.51
# for cluster 1, the mixed half of the pair. The panel therefore still shows
# the contrast the caption describes, with furans in place of amines.
#
# **These four numbers (0.75/0.188/0.045/0.69/0.67/0.51) were derived under an
# earlier metric implementation (an activation-rate calculation) that this
# panel no longer has.** They have not been recomputed against
# `class_cluster_glom_composition`'s responding-glomerulus basis — the two
# don't agree on what a class's "share" means, and this repo doesn't
# fabricate numbers it hasn't actually computed. Rerun this cell and reread
# these values off the current panel before trusting them for the cluster 3
# vs. 5 argument above.
#
# **The caption needs a matching change**, from "near complete selectivity for
# amines" to furans. Set `FIG4C_CLUSTERS = (1, 5)` to restore the original pair
# and accept the two empty sub-panels.
#
# The caption also identifies the two by letter, and those letters follow an
# older lettering that does not match the 4D panel; the content description is
# what was followed here.
#
# **Parameters.** `FIG4_CONC_A` is which concentration panel A correlates —
# `"low"` is the draft's; `"high"` draws the same matrix after broadening,
# which is a real comparison but not this figure's. `FIG4C_CLUSTERS` picks
# C's two clusters. Panel D has no metric switch any more — it is always
# `class_cluster_glom_composition` — and each of its two colour scales is
# autoscaled independently rather than sharing one ceiling. `FIG4_CHEM_ORDER`
# is this section's one deliberate coupling to §3 (the same one §6 already
# has): it's `CLASS_ORDER_FIG3` itself, not an independent copy, and sets
# both B's within-block odorant order and D's row order, so "the same broad
# chemical order as Figure 3/6" holds by construction rather than by hand.
# `FIG4B_SORT_BY_CLASS` toggles whether C (which always reads `bic`, the
# bi-clustered order) sorts each cluster block's odorants by chemical class
# within the block (`True`, default) or leaves them in native odorant order
# (`False`); it also governs B whenever `FIG4B_ROW_ORDER="cluster"`. Either
# way, columns stay glomerulus-cluster-sorted and the block-diagonal
# alignment is unaffected — this only ever changes the order *within* a
# block. D's row order is controlled separately, by `FIG4_CHEM_ORDER`
# directly, and doesn't follow this toggle.
#
# `FIG4B_ROW_ORDER` (default `"odorant"`) is panel B's own row order,
# independent of `bic`/`FIG4B_SORT_BY_CLASS`: `"odorant"` sorts B's rows by
# `FIG4_CHEM_ORDER` alone, ignoring cluster membership, via a second
# `BiclusterOrder` (`bic_b`) built with `bicluster_order`'s new
# `row_order="odorant"`; `"cluster"` restores the original bi-clustered rows
# (and the diagonal block highlight, via `diagonal_blocks=
# (FIG4B_ROW_ORDER == "cluster")`). Panel C always reads `bic`, never
# `bic_b`, so this toggle cannot move C's rows.

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
# same broad chemical order as Fig 3/6" is guaranteed rather than kept in sync
# by hand -- same coupling §6 already has to §3's class order.
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
# Fig 3/6's own odorant order (FIG4_CHEM_ORDER/"burton") instead, ignoring
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
gs_top = fig.add_gridspec(1, 2, width_ratios=[1, 3.1], wspace=0.16,
                          left=0.06, right=0.93, top=0.95, bottom=0.60)
ax_a = fig.add_subplot(gs_top[0, 0])
ax_b = fig.add_subplot(gs_top[0, 1])

im_a = P.clustering.glomerulus_correlation_matrix(
    ax_a, gsim, title={"low": "low", "high": "high"}[FIG4_CONC_A],
    boundaries=False, cluster_boxes=True, cluster_box_color="green")
mesh_b = P.clustering.biclustered(d, ax_b, conc=FIG4_CONC_A, order=bic_b,
                                  diagonal_blocks=(FIG4B_ROW_ORDER == "cluster"),
                                  class_names=True,
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

pos_a = ax_a.get_position()
cax_a = fig.add_axes([pos_a.x0, pos_a.y0 - 0.055, pos_a.width * 0.6, 0.014])
cb_a = fig.colorbar(im_a, cax=cax_a, orientation="horizontal",
                    ticks=[0.0, 0.5, 1.0])
cb_a.set_label(r"correlation $\rho$", fontsize=hs.FS_LABEL)
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

# C's rows are B's rows: Fig 3/6's exact 59-odorant order (CLASS_ORDER_FIG3
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
plt.show()
# export_figure("fig4_tuning_by_cluster")   # -> notebooks/outputs/fig4_tuning_by_cluster.{png,pdf,svg}

# %% [markdown] tags=["human"]
# ## 5. Figure 5 — which chemical classes each glomerular group responds to
#
# > **[Caption not yet written for the merged figure — see notes/merge_decisions.md]**
#
# Panel B of the paper's Figure 5. For each of the seven glomerular groups, how many of its
# glomeruli respond to each chemical class, scaled so the group's top class is 1 — at low
# concentration (left) and at high concentration (right). Panel A (the per-odorant tuning
# bars) is drawn separately and is not built here.

# %% [markdown] tags=["llm"]
# ## 5. Figure 5 — chemical-class composition of each NMF cluster, low vs high
#
# > **[Caption not yet written for the merged figure — see notes/merge_decisions.md]**
#
# Paper Fig 5B. `class_cluster_glom_composition(d, conc, order)`: for each (class, cluster) cell the
# number of distinct glomeruli of that cluster with a nonzero response to at least one odorant of
# that class, divided by the cluster's own maximum over classes, separately for each concentration
# (columns do not sum to 1). Class order is `CLASS_ORDER_FIG3`. Panel A is drawn elsewhere.

# %%
# Own constants — this section inherits no bare name from anywhere.
FIG5_CHEM_ORDER = CLASS_ORDER_FIG3     # same broad chemical order as Figs 3, 4 and 7 ("matt")
FIG5_CMAP = "magma"
FIG5_EMPTY = "gainsboro"
FIG5_BOUNDARY, FIG5_BOUNDARY_ALPHA = "blue", 0.5
FIG5_H = 2.9                           # inches; the width is MAIN_TARGET_W

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
# --- Figure 5B: low concentration beside high concentration ----------------------
fig = plt.figure(figsize=(MAIN_TARGET_W, FIG5_H), dpi=FIG_DPI)
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
plt.show()
# export_figure("fig5_class_by_cluster")   # -> notebooks/outputs/fig5_class_by_cluster.{png,pdf,svg}

# %% [markdown] tags=["human"]
# ## 6. Figure 6 — the seven groups on the dorsal bulb
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# Each map is the same view of the dorsal bulb, with all four mice overlaid
# (the field of view is consistent between them). One group is highlighted per
# map, coloured by mean response at low concentration; every other glomerulus
# is the pale grey background. Groups differ in how spatially compact they are.
#
# **About the colour scale.** It tops out below the brightest glomeruli, so
# about 5% of them saturate. Scaling to the single brightest glomerulus instead
# puts the median glomerulus at 8% of the ramp and the maps read as near-empty.
# It is a **low-concentration** scale: at high concentration about half the
# glomeruli would exceed it.
#
# **Only panel A is ours.** B, C and D in the paper are Matt's bulb
# photographs.

# %% [markdown] tags=["llm"]
# ## 6. Figure 6 — spatial mapping of glomeruli by cluster
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# **What it shows.** Where each of the seven clusters sits on the dorsal bulb.
# Every map is the same field of view — all four subjects superimposed in raw
# pixel coordinates, which is meaningful because the field is consistent
# between them — with one cluster's glomeruli coloured by their mean 1x
# response and every other glomerulus drawn in grey behind them. Some clusters
# are spatially compact and some are not, and that is the panel's content.
#
# **Only panel A is built here.** B, C and D are epifluorescence bulb images
# (omp114), and the "most-sensitive" condition in B is a third concentration
# this dataset does not have.
#
# **The non-responder treatment.** Non-responsive glomeruli are drawn larger
# than the cluster's own points are small, grey, without a border, and pale and
# transparent (`0.85` at `alpha=0.5`). That is
# `nonresponder_marker="requested"`; the preset lives in
# `panels.clustering.NONRESPONDER_MARKERS` beside the alternative `"context"`
# treatment, so the two are one keyword apart. `alpha=0.35` was tried and
# rejected — the bulb's outline dissolves and the map stops reading as a bulb.
# Set `FIG6_NONRESPONDER` to a dict of scatter keywords to overrule any of it.
# The cluster's own points keep a thin white outline (`FIG6_CLUSTER_EDGE`):
# the bottom of the sequential ramp is nearly white, so a weakly responding
# cluster member would otherwise read as a hole rather than as a point.
#
# **The colour ceiling is 5.0, and it is a LOW-concentration ceiling.** At 15.0
# it sat *above* the low-concentration data maximum of 13.88, so nothing
# reached full scale and the median glomerulus rendered at 8% of the ramp. At
# 5.0 the median sits at 24% of scale and 18 of 376 glomeruli (4.8%) saturate.
# **Do not reuse this ceiling for `conc="high"`** — 51% of glomeruli exceed it
# there; `P.clustering.TOPO_VMAX_HIGH` clips the same ~5% at that
# concentration.
#
# **Parameters.** `FIG6_CONC` — the draft page is labelled "low-clustered, 1x
# concentrations", so `"low"`. `FIG6_TOPO_VMIN`/`FIG6_TOPO_VMAX` must stay
# shared across all seven maps, or the comparison between them stops meaning
# anything.

# %%
# Membership only: every cluster member is one solid colour, and no response
# value is drawn -- so there is no concentration to choose and no colourbar.
FIG6_MEMBER_COLOR = "red"
FIG6_NONRESPONDER = "hollow"        # or "context"/"requested", or a dict of scatter kwargs
FIG6_CLUSTER_EDGE = "black"         # see below

# This section has no heatmap panel to key off (see §1), so it gets its own
# scale rather than reading FIG_SCALE directly, solved for MAIN_TARGET_W.
# `_fig5_k` is how much smaller that makes the page than it was under the
# old bare-FIG_SCALE sizing -- marker size (area, so k**2) and linewidth
# (linear, so k) are rescaled by it too, or a page this much smaller would
# leave 7 giant, overlapping dots per row instead of 7 legible clusters.
FIG6_SCALE = MAIN_TARGET_W / 16.0
_fig5_k = FIG6_SCALE / FIG_SCALE
FIG6_CLUSTER_SIZE = 62 * _fig5_k ** 2   # the cluster's own points, drawn on top
FIG6_CLUSTER_LW = 0.6 * _fig5_k
FIG6_CONTEXT_SIZE = 16 * _fig5_k ** 2   # "hollow"'s own base size (NONRESPONDER_MARKERS), scaled
FIG6_CONTEXT_LW = 0.4 * _fig5_k

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
        title=f"Cluster {k}")
    ax.set_xlim(*xlim5)
    ax.set_ylim(*ylim5)
    if k == 1:
        hs.panel_label(ax, "A", x=-0.30, y=1.12)

save(fig, "fig6_spatial_by_cluster")
plt.show()
# export_figure("fig6_spatial_by_cluster")   # -> notebooks/outputs/fig6_spatial_by_cluster.{png,pdf,svg}

# %%
for k in range(1, cl.N_CLUSTERS + 1):
    mask = d.cluster_mask(k)
    xy = d.xy[mask]
    spread = float(np.sqrt(((xy - xy.mean(0)) ** 2).sum(1).mean()))
    print(f"cluster {k}: n = {int(mask.sum()):3d}, "
          f"rms spread {spread:6.1f} px about its centroid")

# %% [markdown] tags=["human"]
# ## 7. Figure 7 — a linear model of the low → high change
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# One linear model, fitted once across all odorants, predicting
# high-concentration responses. Every prediction shown is **held out**: the
# model is fitted on some odorants and asked about others.
#
# - **A, B, C** — three individual glomeruli, the same three as Figure 3E–G.
#   Measured high-concentration response points up; the model's prediction for
#   the same glomerulus points down. Each panel carries its own x-axis
#   chemical-class names — unlike Figure 3E–G, which only labels the middle
#   panel, all three here are labelled.
# - **D** — the recruitment matrix from Figure 3H, measured on the left and
#   predicted on the right (r = 0.54).
# - **E** — three cartoons, not data: the three forms the weight matrix could
#   take. Every glomerulus scaling only its own response (diagonal only),
#   off-diagonal weights present but unstructured, or off-diagonal weights in
#   blocks. The paper's text refers to these as i, ii and iii.
# - **F** — the weight matrix itself, ordered by group.
#
# **Two caveats about A–C.**
#
# The bars are in **log1p ΔF/F** — the space the model was actually fitted
# and scored on, not the raw units Figure 3E–G uses. That's a deliberate
# choice, not just a default: in raw ΔF/F, one acid's held-out prediction on
# panel B (glomerulus 211) reaches **193.6** against a measured **66.2**,
# which flattens every other bar toward the baseline once the axis
# autoscales to fit it. `FIG7_SPACE = "raw"` restores those units — the same
# ones 3E–G use, so the two figures compare directly — at the cost of that
# flattening; it's one keyword away.
#
# Correlations depend on that choice: the same predictions score r = 0.53 /
# 0.78 / 0.64 on the log scale shown here and r = 0.30 / 0.53 / 0.55 in raw
# ΔF/F — same model, better-looking numbers in the space it was scored in.
# The scale is stated on the axis for that reason.
#
# Predictions near zero still draw a visible bar regardless of scale: above a
# threshold of 0.4 the model has non-zero predictions for 35 / 25 / 51
# odorants against 10 / 11 / 29 measured. The exact counts print underneath
# the figure on every run — a cross-validated linear model essentially never
# predicts exactly zero, so this is the model's behaviour, not a drawing
# artifact.
#
# **The schematic that used to be panel A was dropped**, and the rest moved up
# a letter. E and F were kept because the paper's text refers to them directly
# ("Fig. 6F, i–iii" and "Fig. 6G").

# %% [markdown] tags=["llm"]
# ## 7. Figure 7 — PLS-R
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# **What it shows.** The linear-model half of the manuscript on one page: what
# the model predicts for three individual glomeruli (A–C), what it predicts for
# recruitment across all 59 odorants (D), the hypotheses W could have confirmed
# (E), and which one it actually did (F).
#
# **Panel lettering.** The `Y = X · W` schematic that headed earlier versions
# of this figure is drawn separately and is not built here, so the panels run
# **A–F**. The draft text refers to the cartoons and the weight matrix as
# "Fig. 6F, i–iii" and "Fig. 6G"; those are **E** and **F** here, and the
# content of each is unchanged.
#
# - **A, B, C — the same three glomeruli as Fig 3E/F/G**, so that the two
#   figures can be read against each other: the same response spectra with
#   predicted high compared to experimental high. The glomeruli are
#   `similarity.EG_GLOMERULI` — imported, never re-derived — and
#   the panel is `panels.similarity.mirror_spectrum`, the same function that
#   draws 3E/F/G, with the low series replaced by the held-out prediction.
#   Experimental high is **above** the line and predicted high below.
#   `strip_labels=True` on all three (`mirror_spectrum`'s x-axis chemical-class
#   strip) — unlike 3E/F/G, which label only the middle panel to avoid
#   repeating the ten class names three times in one row, each of these three
#   carries its own labels.
# - **D — the recruitment matrix**, observed beside its held-out prediction.
#   This is the emailed `recruitment_matrix_binary.pdf`, computed here from
#   this section's own parameters rather than inherited from §3, which draws
#   the observed half as Figure 3H.
# - **E — three cartoons**, not data: independent glomerular gain control
#   (diagonal only), unstructured retuning (scattered off-diagonal terms), and
#   structured retuning (coherent blocks). All three share one diagonal, so
#   they differ *only* in what is off it. They live in `highconc.schematic`
#   beside the other cartoons, never in a module that touches real data.
# - **F — the weight matrix**, reordered by cluster and labelled 1–7, the
#   numbering used throughout. Earlier renderings carry the legacy A–G labels.
#
# **Why the prediction is denser than the measurement.** A cross-validated
# regression essentially never predicts exactly zero, so the lower half of each
# spectrum has small bars where the glomerulus measured nothing.
# `PredictionSpectrum.support` counts both sides above any threshold, and the
# counts are printed below the figure.
#
# **The display space is `"log"`.** `log1p` is the space the model was
# actually fitted and scored in, so the bars and the correlation printed on
# each panel are the ones the model was optimised against — no transform
# between what was minimised and what is drawn. Its cost is comparability:
# 3E/F/G are in raw ΔF/F, so reading 6A/B/C directly against them needs
# `FIG7_SPACE = "raw"` instead (one keyword away — see below).
#
# **Why `"raw"` isn't the default.** `expm1` undoes `log1p` back into ΔF/F,
# the same units 3E/F/G use, but un-shrinks the prediction doing it: a single
# large held-out value can dominate a panel's autoscaled axis. On glomerulus
# 211 (panel B) one acid's prediction reaches **193.6** against a measured
# **66.2** in raw units, which flattens every other bar in that panel toward
# the baseline. That's a property of the display space, not of the model —
# the same prediction reads cleanly on the log axis shown here.
#
# Predictions near zero still draw a visible bar regardless of scale: above a
# threshold of 0.4 the model has non-zero predictions for **35 / 25 / 51**
# odorants against **10 / 11 / 29** measured. A cross-validated linear rule
# spreads a little predicted response everywhere the population is active, so
# this is the model's behaviour rather than a drawing artifact;
# `PredictionSpectrum.support` prints both counts below the figure.
#
# The two alternatives, both one keyword away and drawn by the same code:
#
# - `"raw"` draws ΔF/F, matching 3E/F/G's units exactly, at the cost of
#   glomerulus 211's overshoot dominating its panel (above).
# - `"sqrt"` compresses that same overshoot without expanding the floor the
#   way `log1p` does, since `sqrt` has a finite derivative at zero where
#   `log`'s is unbounded — a middle ground between the two.
#
# The correlation genuinely differs between them — **0.53 / 0.78 / 0.64 log**,
# 0.48 / 0.70 / 0.60 sqrt, 0.30 / 0.53 / 0.55 raw — and that spread is made by
# the display, not by the prediction: log scores highest because it's the
# space the model was scored in during fitting. The panel's y-axis label
# follows `FIG7_SPACE`, so it cannot claim ΔF/F while drawing something else.
#
# `FIG7_LIMIT` is `None`, so each spectrum autoscales, as 3E/F/G do.
# Set it to `P.similarity.shared_spectrum_limit(...)` to put the three on one
# axis.
#
# `FIG7_ROW_FRAC` / `FIG7_LO_FRAC` / `FIG7_ORDER` are the recruitment-matrix
# parameters, restating §3's `ROW_FRAC_FIG3` / `LO_FRAC_FIG3` for Figure 6E
# (odorant-within-block order comes from §3's
# `WITHIN_CLASS_ORDER`/`ODORANT_CODES`, same coupling as the spectra), and
# `FIG7_TOP_K` is the W display cut. The predicted half of the matrix applies
# the identical tally rule to `m.Yhat` (the held-out CV prediction) in place
# of `m.high` — there is no separate prediction threshold to set any more.
# **This section reads no bare colour constant from anywhere** — every scale
# below is a `FIG7_*` it owns.

# %%
FIG7_SPACE = "log"              # "log" (the model's own space) | "sqrt" | "raw" (dF/F)
FIG7_GLOMERULI = (211, 97)      # C's two example glomeruli (dataset glomerulus ids)
FIG7_LIMIT = None               # None = autoscale each panel
FIG7_SUPPORT_THRESHOLD = 0.4    # value above which a bar counts as a response
FIG7_CANONICAL = "baseline"
# D shows EVERY weight of W: there is no top-K cut. Only the colour scale is
# clipped (each sign at its own 99th percentile, see below), never the data.

# Same odorant order as Fig 3's E-H, so C reads against them.
fig6_specs = pls.eg_predicted_spectra(m, d, FIG7_GLOMERULI, space=FIG7_SPACE,
                                      order=CLASS_ORDER_FIG3,
                                      within=WITHIN_CLASS_ORDER)
fig6_w = pls.load_weights(FIG7_CANONICAL)
# The smallest nonzero observed high response in the whole dataset, in
# FIG7_SPACE: the line drawn across C's predicted half.
fig6_obs_floor = pls.observed_high_floor(m, space=FIG7_SPACE)

for s in fig6_specs:
    n_obs, n_pred, n_both = s.support(FIG7_SUPPORT_THRESHOLD)
    print(f"6C: glomerulus {s.glom} (cluster {s.cluster}) - "
          f"r = {s.r:.3f} in {s.space} space; above "
          f"{FIG7_SUPPORT_THRESHOLD}: {n_obs} measured, {n_pred} predicted, "
          f"{n_both} both")
print(f"6C: smallest nonzero observed high response ({FIG7_SPACE}): "
      f"{fig6_obs_floor:.4f}")
print(f"6D: W is {fig6_w.real.shape[0]} x {fig6_w.real.shape[1]}, "
      f"all weights shown")

# 6C's own numbers, one row per odorant per glomerulus -- built from the same
# `fig6_specs` the panels draw, so the saved numbers are exactly what's shown.
fig6_c_table = pd.concat(
    {s.glom: pd.DataFrame({
        "odorant": s.names, "class": s.classes,
        "observed": s.observed, "predicted": s.predicted,
     }).assign(cluster=s.cluster, space=s.space, r=s.r)
     for s in fig6_specs},
    names=["glom", "row"])
if SAVE:
    fig6_c_csv = OUT / "fig7_pls_C_prediction_spectra.csv"
    fig6_c_table.to_csv(fig6_c_csv)
    print(f"saved: {fig6_c_csv.relative_to(ROOT)}")

# %%
# One row, in inches: C (two example glomeruli, stacked) | D (W, square, + colourbar).
# The width is MAIN_TARGET_W exactly; the right margin is whatever is left, so
# changing any other width cannot silently change the page width.
FIG7_MARGIN_L, FIG7_MARGIN_TOP, FIG7_MARGIN_BOT = 0.65, 0.55, 0.55
FIG7_C_W, FIG7_C_H, FIG7_C_GAP = 2.25, 0.90, 0.60   # each spectrum, and the gap between
FIG7_COL_GAP = 0.75                                  # between C and D
FIG7_D_SIDE = 2 * FIG7_C_H + FIG7_C_GAP              # D is square, as tall as C's stack
FIG7_CBAR_PAD, FIG7_CBAR_W = 0.40, 0.14
FIG7_MARGIN_R = (MAIN_TARGET_W - FIG7_MARGIN_L - FIG7_C_W - FIG7_COL_GAP
                 - FIG7_D_SIDE - FIG7_CBAR_PAD - FIG7_CBAR_W)
assert FIG7_MARGIN_R > 0.30, f"Fig 6 row does not fit: right margin {FIG7_MARGIN_R:.2f} in"
FIG7_W = MAIN_TARGET_W
FIG7_H = FIG7_MARGIN_TOP + FIG7_D_SIDE + FIG7_MARGIN_BOT


def x6(inches):
    return inches / FIG7_W


def y6(inches):
    return inches / FIG7_H


fig = plt.figure(figsize=(FIG7_W, FIG7_H), dpi=FIG_DPI)
top_in = FIG7_H - FIG7_MARGIN_TOP          # top edge of the row, from the bottom

# --- C: two example glomeruli, experimental high up / predicted high down ----
for j, s in enumerate(fig6_specs):
    y0 = top_in - FIG7_C_H - j * (FIG7_C_H + FIG7_C_GAP)
    ax = fig.add_axes([x6(FIG7_MARGIN_L), y6(y0), x6(FIG7_C_W), y6(FIG7_C_H)])
    P.similarity.prediction_spectrum(
        ax, s, limit=FIG7_LIMIT, strip_labels=True, annotate_r=False,
        observed_floor=fig6_obs_floor)
    # r goes in the title, as in Fig 3's E-H, not inside the axes
    ax.set_title(f"glomerulus {s.glom}  (cluster {s.cluster})    r = {s.r:.2f}",
                 fontsize=HEATMAP_TITLE_FONTSIZE)
    if j == 0:
        hs.panel_label(ax, "C", x=-0.18, y=1.18)

# --- D: the full W, every weight ---------------------------------------------
d_left = FIG7_MARGIN_L + FIG7_C_W + FIG7_COL_GAP
W6 = fig6_w.real
W6_reord, _, bnd6, centers6 = pls.reorder_by_cluster(W6, fig6_w.cluster_number)
# Each sign scaled to its own 99th percentile -- negatives are smaller, and a
# symmetric range washes them out. TwoSlopeNorm keeps 0 on gainsboro. This only
# sets the colour range; every weight is still drawn.
pos6, neg6 = W6_reord[W6_reord > 0], -W6_reord[W6_reord < 0]
norm6 = TwoSlopeNorm(vmin=-(np.percentile(neg6, 99) if neg6.size else 1.0),
                     vcenter=0.0,
                     vmax=(np.percentile(pos6, 99) if pos6.size else 1.0))
ax_d = fig.add_axes([x6(d_left), y6(FIG7_MARGIN_BOT),
                     x6(FIG7_D_SIDE), y6(FIG7_D_SIDE)])
P.pls.weight_blocks(ax_d, W6_reord, bnd6, centers6, norm=norm6,
                    title="W (low → high concentration mapping)")
hs.panel_label(ax_d, "D", x=-0.11, y=1.05)
cax_d = fig.add_axes([x6(d_left + FIG7_D_SIDE + FIG7_CBAR_PAD),
                      y6(FIG7_MARGIN_BOT + FIG7_D_SIDE * 0.22),
                      x6(FIG7_CBAR_W), y6(FIG7_D_SIDE * 0.56)])
cb_d = fig.colorbar(plt.cm.ScalarMappable(cmap=hs.DIV, norm=norm6), cax=cax_d,
                    extend="both")
cb_d.set_label("W weight", fontsize=hs.FS_LABEL)
cb_d.ax.tick_params(labelsize=hs.FS_TICK)

save(fig, "fig7_pls")
plt.show()
# export_figure("fig7_pls")   # -> notebooks/outputs/fig7_pls.{png,pdf,svg}

# %% [markdown] tags=["human"]
# ## 8. Figure 8 — directional transfer of tuning between groups
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# **The test.** Predict a group's high-concentration responses from its own
# low-concentration average — that is the baseline. Then add a second group's
# low-concentration average and see whether the prediction improves. Repeat for
# every ordered pair, in both directions. Prediction is scored on held-out
# odorants, so an improvement is not simply the result of adding a parameter.
#
# - **Left panels** — the permutation *p* for each pair, against 10,000
#   shuffles. A **black box** marks pairs surviving Benjamini–Hochberg
#   correction at q < 0.05.
# - **Right panels** — the size of each improvement, as a percentage of the
#   error the target's own average left behind. Most cells are near zero.
#
# **Two groupings.** The top row groups glomeruli by the seven NMF clusters;
# the bottom row groups them by the chemical class each glomerulus prefers,
# which does not use the clustering. Seven of the eight significant cluster
# results have a counterpart in the class analysis.
#
# **One thing not to misread.** In the bottom-right panel the furan row is
# strongly coloured in the negative direction. Furan-preferring glomeruli have
# a negative intra-group baseline (Q² = −0.08), so the percentage runs against
# a negative denominator. None of those cells is boxed.
#
# **The explanatory header in the paper's version is Matt's**, not drawn here.

# %% [markdown] tags=["llm"]
# ## 8. Figure 8 — motif transfer between glomerular groups
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# **What it shows.** The direct test that replaces reading the off-diagonal
# weights of `W`. For every *ordered* pair of glomerular groups `F → T`: does
# adding group **F**'s mean low-concentration response to the prediction of
# group **T**'s glomeruli at high improve held-out performance, beyond what
# **T**'s own low-concentration mean already achieves?
#
#     intra    = predict T's high from  [T's low mean]
#     intra+F  = predict T's high from  [T's low mean, F's low mean]
#     ΔQ²      = Q²(intra+F) − Q²(intra)        # out-of-fold, pooled over T
#
# A positive ΔQ² says a tuning motif carried by F at low concentration
# reappears *inside* T at high. That is the manuscript's claim of structured
# retuning, stated as a prediction rather than as an inspection of weights.
#
# **Why not just read W.** Because `W`'s off-diagonal weights cannot identify
# which groups inherit whose tuning: the diagonal and off-diagonal predictors
# are 94% redundant (`cov(X)` is not diagonal), so the same held-out
# performance can be assembled many ways. The pairwise test avoids that by
# asking whether the *prediction* improves, which is a property of performance
# rather than of a particular parameterisation.
#
# - **Left — permutation p per pair**, clipped at 0.1. Colour separates the
#   transfers that happened from the ones that did not; the printed number is
#   the unclipped p, and `<1e-4` means the 10,000 draws cannot resolve it
#   further. The **black box** is the claim: BH q < 0.05 across all ordered
#   off-diagonal pairs (42 for clusters, 72 for classes).
# - **Right — the size of the transfer**, as a percentage of the error the
#   target's own mean left behind, on the diverging ramp with zero at grey.
#   Most cells are grey.
#
# **Two groupings, the same test.** The upper row groups glomeruli by NMF
# cluster (renumbered to Joel's 1–7); the lower row is the **chemical-class
# control** — each glomerulus assigned to the class evoking its largest mean
# low-concentration response, keeping the nine classes with at least 15
# glomeruli (372 of 376; the ester-preferring group, n = 4, is excluded). The
# control matters because it never saw the NMF factorisation: if the transfers
# survive it, they are a property of the population and not of the clustering.
# They do — **seven of the eight cluster edges have a direct class counterpart**,
# and the pairing is printed below the figure rather than asserted.
#
# **One group is predicted worse than its own mean.** Furan-preferring
# glomeruli have an intra baseline of Q² = −0.084: predicting them from their
# own low-concentration group mean is *worse* than predicting the mean of the
# data. Their row in the right-hand panel is correspondingly deep teal — adding
# furans as a source costs the target held-out variance almost everywhere. That
# is a real property of the group, not a plotting artifact, and it is a second
# reason `FIG8_METRIC` defaults to `"pct_err"`: `"pct_var"` divides by that
# baseline and would flip the sign of every furan cell while inflating its
# magnitude. Switch metrics knowing this.
#
# **Read the axes, not the shape.** Rows are the source, columns are the
# target. The matrix is *not* symmetric and is not meant to be: 11 of the 13
# class-level transfers are one-way. The exception is worth reporting rather
# than smoothing over — **amines ↔ thiazoles is reciprocal** at class
# resolution, whereas NMF cluster 5 is a strict sink. The amine → thiazole
# increment is the smallest significant transfer in either analysis (+0.013).
#
# **The null, and that it is a choice.** Significance comes from 10,000
# *broken-pairing* permutations: the odorant order of Y is shuffled against X,
# which preserves the marginal distribution of high responses and all
# correlation structure within X and within Y, while destroying the low-to-high
# correspondence. This is the null the manuscript publishes. A size-matched
# control run on 2026-08-20 leaves only three of the eight cluster edges
# standing and removes the cluster-5 sink example — a real caveat, recorded in
# the methods audit, and deliberately not what this figure shows
# (`FIG8_NULL`).
#
# **Nothing here recomputes a permutation.** Both nulls are cached in
# `data/pairwise/` and ship with the clone; `transfer.pairwise` reads them and
# recomputes only the intra-group baselines, which are seconds. Rebuilding a
# cache is an explicit `transfer.compute_null(...)` call — 11 minutes for the
# clusters, 16 for the classes serially — and `pairwise` raises rather than
# silently starting one.
#
# **The three-step explanatory header** on Matt's version of this page is
# Illustrator work and is out of scope here (2026-08-27).

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

fig7_clusters = transfer.pairwise(m, "clusters", alpha_fdr=FIG8_ALPHA,
                                  null=FIG8_NULL)
fig7_classes = transfer.pairwise(m, "chem_class", alpha_fdr=FIG8_ALPHA,
                                 null=FIG8_NULL)

for t in (fig7_clusters, fig7_classes):
    s = t.summary(FIG8_ALPHA)
    print(f"{s['grouping']:>11}: {s['n_significant']} of {s['n_pairs']} "
          f"ordered pairs at q < {FIG8_ALPHA} "
          f"({s['n_one_way']} one-way), N = {s['n_perm']:,}")
print()
print(P.transfer.edge_table(fig7_clusters, alpha=FIG8_ALPHA,
                            metric=FIG8_METRIC))

# %%
# One row: A and B, the cluster analysis. The chemical-class control (what
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
P.transfer.pvalue_matrix(ax_p, fig7_clusters, alpha=FIG8_ALPHA, clip=FIG8_CLIP,
                         xlabel=FIG8_XLABEL, ylabel=FIG8_YLABEL, cmap=FIG8_P_CMAP,
                         box_color=FIG8_BOX_COLOR, text_significant_only=True)
hs.panel_label(ax_p, "A", x=-0.26, y=1.12)

ax_m = fig.add_subplot(gs7[0, 1])
P.transfer.improvement_matrix(ax_m, fig7_clusters, metric=FIG8_METRIC,
                              alpha=FIG8_ALPHA, xlabel=FIG8_XLABEL, ylabel="",
                              box_color=FIG8_BOX_COLOR,
                              text_significant_only=True)
hs.panel_label(ax_m, "B", x=-0.16, y=1.12)
P.transfer.significance_legend(ax_m, fig7_clusters, alpha=FIG8_ALPHA,
                               loc=(1.02, -0.16), box_color=FIG8_BOX_COLOR)

fig.text(0.5, 0.98, "Motif transfer between glomerular groups "
                    f"(broken-pairing null, N = {fig7_clusters.n_perm:,})",
         ha="center", va="top", fontsize=hs.FS_TITLE + 1)

save(fig, "fig8_motif_transfer")
plt.show()
# export_figure("fig8_motif_transfer")   # -> notebooks/outputs/fig8_motif_transfer.{png,pdf,svg}

# %%
# The control's actual claim, checked rather than asserted: each significant
# cluster edge is looked up against the class edge between the two clusters'
# dominant preferred classes. Seven of eight should match. If this ever prints
# fewer, the two analyses have drifted apart — which is exactly what the
# control exists to detect.
fig7_match = transfer.matching_edges(fig7_clusters, fig7_classes, m,
                                     alpha=FIG8_ALPHA)
print(f"{'cluster edge':>16}  {'ΔQ²':>8}  {'via':>26}  chemical-class counterpart")
for row in fig7_match:
    e, c = row["cluster_edge"], row["counterpart"]
    via = f"{row['via'][0]} → {row['via'][1]}"
    tail = (f"{c['dq2']:+.4f}  (q = {c['q']:.4f})" if c else "— none —")
    print(f"{e['source'] + ' → ' + e['target']:>16}  {e['dq2']:+8.4f}  "
          f"{via:>26}  {tail}")
n_match = sum(r["counterpart"] is not None for r in fig7_match)
print(f"\n{n_match} of {len(fig7_match)} cluster edges have a class counterpart")

# %%
# Per-group intra baselines — the denominator of both percentage readings, and
# the reason `pct_var` is not the default: cluster 5 is predicted from its own
# mean at Q² = 0.046 and cluster 2 at 0.438, so the same absolute increment
# would read ten times larger on one than on the other.
for t in (fig7_clusters, fig7_classes):
    print(f"{t.grouping.id:>11} intra Q²: " +
          "  ".join(f"{lab}={q:.3f}"
                    for lab, q in zip(t.labels, t.q2_intra)))
print(f"\nchemical classes tested: {list(fig7_classes.labels)}")
print(f"glomeruli excluded (class below n = {FIG8_MIN_CLASS}): "
      f"{fig7_classes.grouping.dropped} of {m.n_gloms}")

# %% [markdown] tags=["human"]
# ## 9. Figure S1 — how the number of groups was chosen
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# Seven at low concentration (top row, **A–C**) and five at high (bottom row,
# **D–F**).
#
# - **A / D** — cumulative variance explained against the number of groups, for
#   PCA, PCA on shuffled data, and NMF. The gap between the real and shuffled
#   curves is the structure not attributable to chance.
# - **B / E** — the slope of A: how much each additional group adds.
# - **C / F** — adjusted Rand index across 50 NMF runs from different random
#   starts, measuring how stable the grouping is to the seed. The insets are
#   consensus matrices at the chosen number and at one past it.
#
# **This one is slow and is cached.** The underlying sweep takes minutes, so it
# is stored on disk and read back. If the cache is missing the notebook stops
# with an error rather than recomputing, because a recomputed number might
# differ from the paper's.

# %% [markdown] tags=["llm"]
# ## 9. Figure S1 — NMF dimensionality choice (companion to Figure 4)
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# **What it shows.** Why the clustering uses *k* = 7 at 1x and *k* = 5 at 25x.
# The manuscript page carries both concentrations as **one** figure, so this
# section
# draws them as one: the 1x row is **A–C**, the 25x row is **D–F**.
#
# - **A / D** cumulative variance explained vs. *k*, for PCA, PCA on shuffled
#   data, and NMF. The gap between real and shuffled is the dimensionality
#   worth keeping.
# - **B / E** the slope of A — where the curve stops paying for itself.
# - **C / F** adjusted Rand index across 50 NMF initialisations vs. *k*: how
#   stable the clustering is to the seed, with consensus matrices (insets) at
#   two *k* values — the chosen one, and one past it where the block coherence
#   starts to go.
#
# The two concentrations are drawn as one figure, to match the manuscript.
# Both halves come from `panels.dimchoice`; this section is composition only.
#
# **This is the expensive one, and it is cached.** The sweep is NMF at *k* =
# 1…20 x 50 seeds, plus PCA, plus shuffled PCA over 20 permutations — minutes
# of compute. `dimsweep.load_cached` reads the `.npz` files under
# `figures/nmf_dim_choice_*/data/cache/` **in place** and raises if one is
# missing rather than recomputing, so that a silent recompute cannot put the
# notebook out of step with the manuscript.

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

    inner = GridSpecFromSubplotSpec(len(consensus_ks), 1,
                                    subplot_spec=gs_bot[0, 5], hspace=0.30)
    for i, k in enumerate(consensus_ks):
        P.dimchoice.consensus(consensus[k], fig.add_subplot(inner[i, 0]))

    for ax, letter in zip((ax_var, ax_slope, ax_ari), letters):
        hs.panel_label(ax, letter)
    ax_var.set_title(f"{conc} concentration — chosen k = {chosen}",
                     fontsize=hs.FS_PANEL, loc="left")

    print(f"{conc}: k = {sweep.k_values.tolist()} | chosen k = {chosen} | "
          f"consensus shown at {consensus_ks} | panels {'/'.join(letters)}")

save(fig, "supp_nmf_dim_choice", dpi=SUPP_DIM_DPI)
plt.show()
# export_figure("supp_nmf_dim_choice")   # -> notebooks/outputs/supp_nmf_dim_choice.{png,pdf,svg}

# %% [markdown] tags=["human"]
# ## 10. Figure S2 — the same structure within each mouse
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# Figure 4 pools glomeruli from four mice. This shows the same structure
# computed within each animal separately.
#
# - **A** — each mouse's own glomeruli, sorted the same way.
# - **B** — where each group sits on the bulb, one mouse per row and one group
#   per column, over a grey backdrop of every glomerulus from every animal so
#   all 28 maps share a frame.
#
# Group 5's maps are empty for two of the mice, which have no glomeruli in that
# group.
#
# **If you are comparing against an older version of this figure:** the source
# page labels the groups A–G and this one numbers them 1–7. The letters do not
# map onto the numbers in order.

# %% [markdown] tags=["llm"]
# ## 10. Figure S2 — per-subject reproducibility (companion to Figure 3)
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# **What it shows.** That the bi-cluster block structure is not an artifact of
# pooling glomeruli across the four subjects. Everything here is at **low**
# concentration.
#
# - **A** — each subject's own glomeruli, sorted by cluster, against the shared
#   odorant-cluster row order. The same diagonal block structure appears
#   independently within every subject.
# - **B** — topography grid, subjects x clusters. Each small map is one
#   subject's glomeruli of one cluster, in that cluster's colour, over a grey
#   backdrop of every glomerulus from every subject, so all 28 maps share a
#   frame.
#
# **Two panels, matching the manuscript's Figure S2 page**, which carries the
# heatmaps and the topography grid. An earlier version also drew per-subject
# block averages; `panels.clustering.per_subject_blocks` remains in the package
# for anyone who wants it.
#
# **The 1-7 numbering is applied in the panel layer.** `per_subject_heatmap`
# and `per_subject_topography` both write the cluster number in the cluster's
# colour, so this section inherits it. The source page labels the clusters A-G,
# and a letter does not tell you the number — `highconc.clusters` holds the
# mapping.
#
# Subjects 3 and 4 have **no cluster-5 glomeruli**, so their cluster-5 maps are
# empty. (That cluster is "D" under the legacy letters, which is why the source
# page's D column is the empty one.)
#
# **Note on type size.** `SUPP_SUBJ_LABEL_FS`/`SUPP_SUBJ_TITLE_PAD` scale with
# the canvas's own width (`SUPP_SUBJ_SCALE`) rather than reading the house
# `hs.FS_PANEL`, which would read smaller than the subject titles it heads;
# the region letters are set from these, not `FS_PANEL`, for that reason.

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

gs_heat = fig.add_gridspec(1, 4, left=0.11, right=0.97, top=0.767, bottom=0.097,
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
ax_heat[0].set_ylabel("odorants\n(by chemical class)", fontsize=5)
ax_heat[0].yaxis.set_label_coords(-0.16, 0.5)   # clear of the class strip

# Chemical-class strip on the leftmost heatmap only, in the plotted row order.
_row_cls = d.chem_class[supp_subj_order.rows]
_cls_name = {c: cl.class_name(c) for c in dict.fromkeys(_row_cls)}
hs.class_strip(ax_heat[0], [_cls_name[c] for c in _row_cls],
               colors={_cls_name[c]: cl.class_color(c) for c in _cls_name},
               orientation="y", labels=False)

# The colourbar is shortened to 70% of the panel height and centred on it. It
# sits left of the y-label, which sits left of the class strip.
p0 = ax_heat[0].get_position()
_cb_h = p0.height * 0.70
cax = fig.add_axes([0.030, p0.y0 + (p0.height - _cb_h) / 2, 0.012, _cb_h])
cb = fig.colorbar(plt.cm.ScalarMappable(
    norm=LogNorm(SUPP_SUBJ_HEAT_VMIN, SUPP_SUBJ_HEAT_VMAX),
    cmap=SUPP_SUBJ_CMAP), cax=cax)
cax.yaxis.set_ticks_position("left")
cax.yaxis.set_label_position("left")
cb.set_label("\u0394F/F", fontsize=10)
cb.set_ticks([1, 10, 100])
cb.set_ticklabels(["1", "10", "100"])
cb.ax.tick_params(labelsize=8)
cb.ax.minorticks_off()
cb.outline.set_visible(False)

save(fig, "supp_per_subject_heatmaps", dpi=SUPP_SUBJ_DPI)
plt.show()
# export_figure("supp_per_subject_heatmaps")   # -> notebooks/outputs/supp_per_subject_heatmaps.{png,pdf,svg}

# %%
# --- Figure S2B: where each mouse's members of each cluster sit -------------
# Drawn exactly like Fig 5: the same page scale (so each map is Fig 5's size),
# the same marker sizes, black outline on the members and the same pale hollow
# context. Restated here, not read from §5, so each section runs on its own.
# Only the fill differs: it changes with the cluster, as Fig 5's red does not.
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

fig = plt.figure(figsize=(MAIN_TARGET_W, SUPP_SUBJ_TOPO_H), dpi=SUPP_SUBJ_DPI)
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

save(fig, "supp_per_subject_topography", dpi=SUPP_SUBJ_DPI)
plt.show()
# export_figure("supp_per_subject_topography")   # -> notebooks/outputs/supp_per_subject_topography.{png,pdf,svg}

# %% [markdown] tags=["human"]
# ## 12. Figure S4 — predicted vs measured, at three levels of averaging
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# In each panel the horizontal axis is what was measured and the vertical axis
# is what the model predicted; perfect prediction would put every point on the
# diagonal.
#
# - **A** — averaged by group and odorant. *r* = 0.80.
# - **B** — averaged by group and chemical class. *r* = 0.89.
# - **C** — every individual glomerulus–odorant pair, unaveraged, as a density
#   plot. *Q²* = 0.27 over 3,405 responding pairs.
#
# The three are the same predictions at three levels of averaging, and the fit
# is looser the less averaging is applied.
#
# A fourth view, averaging by glomerulus, is not shown: its *r* is near 1 by
# construction, since mean-centring guarantees it.
#
# **Where this figure comes from.** These three panels sat loose on a page of
# the draft under "Not sure about where these figure panels should go ???".
# They are referenced in the text, so they were collected into this supplement.
#
# The draft's version of panel C annotates the statistic as R²; it is computed
# out-of-fold, so it is written here as Q². The value is the same to three
# decimals.

# %% [markdown] tags=["llm"]
# ## 12. Figure S4 — parity of the linear model at three granularities
#
# > **[Caption not yet updated for the merged figure — see notes/merge_decisions.md]**
#
#
# **What it shows.** How well the PLS-R model's held-out predictions track the
# observed high-concentration responses, at three levels of aggregation:
#
# - **A** — cluster (n = 7) x odorant (n = 59) means, one colour per cluster.
#   *r* = 0.80.
# - **B** — cluster (n = 7) x chemical-class (n = 10) block means: does the
#   model reproduce the chemical *rules*, not just the totals? *r* = 0.89.
# - **C** — every responding glomerulus x odorant pair (observed > 0), as a
#   hexbin. *Q²* = 0.27 over n = 3,405 pairs. Unaggregated, and the weakest of
#   the three; the three panels are the same predictions at three levels of
#   averaging.
#
# **Where this figure comes from.** These three panels appear on page 7 of the
# draft figure PDF without an assigned figure number. They are referenced in
# the manuscript, so they are collected here as a supplement.
#
# **Three panels and no more.** The rules matrices, the residual split, and the
# per-glomerulus and per-odorant mean parities are not on the draft page and
# are not drawn. Per-glomerulus-mean parity in particular is left out: its
# *r* ≈ 1.0 follows from mean-centering rather than from the model.
#
# **Panel C uses the pinned framing.** `responding_parity` offers two:
# `frame="pinned"` fixes the axes to 1–5 and raises the colour ceiling above
# the densest bin so the cloud's shape reads, and `frame="square"` scales a
# full parity frame to the data. The draft panel is the pinned one — 1–5 axes,
# a 1–10 "# pairs" bar — so this section calls the panel at its defaults.
#
# **One naming difference, flagged rather than silently resolved.** The draft
# panel annotates "R² = 0.27". The statistic is `pls.baseline_q2`, computed
# out-of-fold, so it is written here as **Q²**. The value is identical to three
# decimals — 0.272 — so this is a label difference, not a numerical one.
#
# All three numbers reproduce the draft page: *r* = 0.803, *r* = 0.887, and
# Q² = 0.272 over 3,405 responding pairs.

# %%
# Own constants — this section inherits no bare name from anywhere.
# Only panel C of the old three-panel figure remains: the parity of every
# responding glomerulus x odorant pair, observed vs held-out predicted.
SUPP_PAR_C_CMAP = "Purples"          # hexbin colormap
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
plt.show()
# export_figure("supp_parity")   # -> notebooks/outputs/supp_parity.{png,pdf,svg}

# %% [markdown] tags=["human"]
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

# %% [markdown] tags=["llm"]
# ---
#
# **Note.** Every figure above is in `notebooks/outputs/` as PNG, PDF and SVG,
# plus a second, SVG-only copy in `joel_figures/` at the repo root.
# The PDFs and SVGs carry editable text and a transparent background, so either
# opens in Illustrator with the labels still real text objects.
#
# To change what a panel draws, edit its function in `highconc/panels/`; every
# figure that uses it follows. To change a layout, edit the cell that builds it.
