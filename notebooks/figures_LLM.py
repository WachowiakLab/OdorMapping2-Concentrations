# %% [markdown]
# <!-- GENERATED FILE — do not edit. Source: notebooks/figures.py.
#      Rebuild with `python dev/build_notebooks.py`. -->
#
# *This is the **llm** edition of the figure notebook — the exact one — provenance, parameters and values, stated in full.*
# *Same code, same figures, different words: the other edition is `figures_human`.*
# *Edits go in `notebooks/figures.py`, never here — this file is overwritten.*

# %% [markdown]
# # High-concentration odor-mapping — every figure, one notebook
#
# This notebook draws every figure in the manuscript that is built from this data. It is the only place
# composition lives: the analysis is in the `highconc` package, each panel is a function in
# `highconc.panels.*` that takes an `ax` and draws into it, and the grids, colorbars and panel letters below
# are the part you are meant to edit.
#
# To change a panel's *content*, edit `highconc/panels/<module>.py`.
# To change a figure's *layout*, edit the cell that builds it.
# To change the type scale, colours, or export contract, edit `highconc/style.py` — no cell here sets rcParams.
#
# **Run it top to bottom.** Later sections reuse the model fitted in section 1; running a PLS section on its
# own re-fits it (about a minute) and then caches.
#
# ## Contents
#
# The numbering follows the manuscript PDF (`manuscript/allfigs_hiconcms_penultimate_102126.pdf`, Figures 1–8
# and S1–S4) and uses the paper's panel letters: §*n* draws Figure *n* for 2 ≤ *n* ≤ 8 and §9–§12 draw Figures
# S1–S4. A figure here draws only the panels the paper builds from this data; the rest are drawn elsewhere and
# are recorded as `skipped` in `dev/CONCORDANCE.yaml`: Figure 1 entirely, Fig 2B–E, Fig 6B, Fig 7A
# (the X → W → Y schematic) and 7C (the three weight-matrix cartoons), and Fig 8A (the procedure schematic).
# `dev/concordance.py --check` enforces that the notebook's output set is exactly the manifest's stems.
#
# | § | Figure | Stem | Panels |
# |---|---|---|---|
# | 1 | Setup — style, paths, shared data and model fit | — | — |
# | 2 | **Figure 2** — odorant x glomerulus response, Burton's own dataset, three concentrations | `fig2_response_matrix` | `panels.similarity` |
# | 3 | **Figure 3** — correlational structure vs. concentration, recruitment, example spectra; and its threshold supplement | `fig3_correlational_structure`, `fig3_supp_threshold_robustness` | `panels.similarity`, `panels.recruitment` |
# | 4 | **Figure 4** — NMF clusters, per-subject detail, and the folded-in Sankey | `fig4_tuning_by_cluster` | `panels.clustering`, `panels.sankey` (analysis: `folding`) |
# | 5 | **Figure 5** — chemical-class composition of each cluster, low vs high | `fig5_class_by_cluster` | `panels.clustering` |
# | 6 | **Figure 6** — spatial mapping of glomeruli by cluster | `fig6_spatial_by_cluster` | `panels.clustering` |
# | 7 | **Figure 7** — PLS-R: predicted spectra and the W matrix | `fig7_pls` | `panels.similarity`, `panels.pls` |
# | 8 | **Figure 8** — motif transfer between glomerular groups | `fig8_motif_transfer` | `panels.transfer` |
# | 9 | **Figure S1** — NMF dimensionality choice (low concentration) | `supp_nmf_dim_choice` | `panels.dimchoice` |
# | 10 | **Figure S2** — per-subject response heatmaps | `supp_per_subject_heatmaps` | `panels.clustering` |
# | 11 | **Figure S3** — per-subject topography of the clusters | `supp_per_subject_topography` | `panels.clustering` |
# | 12 | **Figure S4** — parity over every responding pair | `supp_parity` | `panels.pls` |
#
# **Not built any more.** The old concentration-flow supplement (`supp_tuning_flow`) is not in the paper PDF; it
# is `retired` in the manifest and its section is gone, which leaves `highconc.flow` and `panels.flow` unused.
# `panels.schematic` is unused too (the retuning cartoons that sat at Fig 6E are now the paper's Fig 7C, drawn
# elsewhere). `fig3_supp_threshold_robustness` is kept although the PDF has no such figure.
#
# Each section opens with a markdown cell describing what the figure shows and how it was made. That text is
# the figure legend; it sits next to the code so the two cannot drift apart.
#
# ## Things not to break
#
# 1. **The NMF dimensionality sweep is cached** (`figures/nmf_dim_choice_1x/data/cache/*.npz`; the 25x cache is
#    no longer read) and takes minutes to recompute. `dimsweep.load_cached` raises on a miss rather than silently
#    recomputing — leave it that way.
# 2. **Nothing is refitted at draw time.** `highconc.pls.load_weights`, `.component_selection` and `.lv_loadings`
#    read the canonical run pls-explorer wrote to disk, and raise if it is missing rather than quietly
#    producing different numbers.
# 3. **The folded-in reclustering is cached** (`data/nonresponders_folded/glom_clusters.csv`, with `meta.json`).
#    `highconc.folding.load_cache` raises if it is missing. Rebuild it only deliberately, with
#    `python dev/fold_nonresponders.py`: it is a seeded NMF fit and uses the package's own seed and iteration cap.
#    Calling the clustering with `params={"n_components": 7}` instead would drop both and change on every run.
#
# Cluster numbering is Joel's **1–7** everywhere. The legacy A–G letters survive only as a mapping in
# `highconc.clusters`; a letter does not tell you the number.

# %% [markdown]
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
# ## 2. Figure 2 — odorant x glomerulus response, Burton's own dataset
#
# **What it shows.** Paper Fig 2A: the same 59 odorants against Burton's own glomeruli — not this repo's 376 —
# at three conditions, stacked top to bottom in the paper's order: "most sensitive", "low", "high". The
# "most sensitive" matrix is from a different
# subject/preparation (`sortedmatrix_omp114L_omp8x.csv`, 76 glomeruli); "low" is `sortedmatrixlo.csv` and "high"
# is `sortedmatrix.csv` (113 glomeruli each). All three were supplied 2026-09 with no header or odorant labels,
# so rows are matched by position against `odorant_order.txt`, the same convention `most_sensitive_correlation`
# (Figure 3B) uses. The share of nonzero cells printed under the figure is 3.8% / 4.9% / 12.1% at most-sensitive /
# low / high.
#
# **How to read it.** Rows are the 59 odorants, permuted into the same class-block order Figure 3 uses
# (`class_sorted`, `order="matt"`), with the chemical-class strip and class names on the *right* of each matrix —
# `highconc.style.class_strip` with `orientation="y"`. All three share one colour scale. Columns are Burton's own
# glomeruli in whatever order each source file ships them; there is no identity or ordering information for them,
# so none is imposed, and they are left unlabelled. The top matrix's narrower column count (76 vs 113) reflects a
# smaller source population, not a cropped shared axis.
#
# **Panel letters.** The paper has one panel A holding all three matrices, so only the top matrix is lettered; the
# other panels of the paper's Figure 2 (B–E: sorted spectra, histogram, cumulative lifetime sparseness, S_L low vs
# high) are not drawn here.
#
# **Parameters.** `FIG2_CONCS = ("most_sensitive", "low", "high")` sets the top-to-bottom order.
# `FIG2_ORDER = "matt"` and `FIG2_WITHIN = "burton"` are this section's own row-ordering parameters. Note that
# Figure 3 now uses `WITHIN_CLASS_ORDER = "matt"` (only the pyrazine block reshuffled), whereas this figure
# keeps the full 59-odorant `"burton"` order inside each class; the class blocks are identical, only the
# within-class order differs. Every section owns its constants, so Figure 2 does not depend on Figure 3 having
# been run.

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
# ## 3. Figure 3 — correlational structure of odor representations vs. concentration
#
# **What it shows.** That raising concentration does not just make responses bigger, it makes the
# *relationships between odorants* different. Panels C and D are the same 59 × 59 odorant × odorant correlation
# matrix at the two concentrations, and the high-concentration matrix is visibly denser: the mean off-diagonal
# correlation goes from 0.03 to 0.10, and the two matrices agree with each other only at Spearman 0.55. Panel B
# is the same matrix at a third, "most sensitive" condition, supplied 2026-09-02 in
# `data/burton_et_al/most_sensitive.csv`. Panel A is the chemical similarity of the same 59 odorants, which of
# course does not change — so it is the fixed ruler the other three are measured against. Chemistry predicts
# response similarity **better at high concentration than at low** (Spearman 0.21 → 0.31).
#
# **Panel lettering follows the manuscript figure PDF** (Figure 3, A–F): A MACCS, B most-sensitive, C low, D high,
# E recruitment, F the example spectra. The notebook's earlier lettering (A–D shuffled, E–H spectra, I
# recruitment) is gone; the draw order is now simply left to right, with no switch for it.
#
# **How to read each panel.**
#
# - **A** — MACCS cosine similarity from the odorants' SMILES. Cosine rather than Tanimoto so that it is the same
#   normalised-inner-product quantity as the correlations beside it. The 59 × 59 matrix is cached to `data/chem/`,
#   so drawing this needs no RDKit; only rebuilding the cache does.
# - **B, C, D** — odorant × odorant Pearson correlation between response maps (each odorant is a 376-glomerulus
#   vector), at the "most sensitive" condition, at 1x and at 25x. All four of A–D are pinned to the same 0–1 scale,
#   and they must stay pinned: B and C are sparse and D is not, so autoscaling them separately would erase the
#   comparison the figure exists to make. The grey diagonal is *not a measurement* — self-similarity is 1 by
#   construction — while a pale cell is a real, near-zero correlation. Panel B ships with two extra, unnamed
#   rows/columns beyond the 59 in `odorant_order.txt`; `most_sensitive_correlation` drops them. Four odorants
#   (D15, E1, P1, P4) have no data at this condition — filled at 0.0 rather than left NaN, so they read at the
#   colour floor instead of as a hole, while the true diagonal is still re-masked to NaN like every other panel.
# - **Strip labels are on all four of A–D** (`strip_labels=True`): each panel carries the coloured class strip with
#   the class names, since all four share one class order.
# - **A–D's aromatics breakdown.** The single "R" (aromatics) class is split into three display subclasses for the
#   axis strip and boundaries only — `highconc.clusters.AROMATIC_SUBCLASS_BY_CODE`, applied via
#   `similarity.split_block`: Ar. Ketone (R27), Ar. Ketone/Ester (R11, R12, R21, R31, R3, R29, R40), Phenol/Methoxy
#   (R10, R23), coloured `#a32e49` / `#00bbb4` / `#d59537`. It touches nothing but these four panels — the other
#   figures and panel F still see plain "R". Toggle with `FIG3_SPLIT_AROMATICS`. Each of A–D also carries five green diagonal
#   boxes (`box_boundaries=FIG3_BOX_BOUNDARIES`, `box_color=FIG3_BOX_COLOR`), as in the paper: acids through
#   alcohols, aromatics + furans, amines, pyrazines, thiazoles — five groups, not the full ten/twelve class
#   blocks. The cut points `(23, 35, 42, 52)` are specific to `CLASS_ORDER_FIG3 = "matt"` with the aromatics
#   split, so a different order or setting changes what block sizes they land on and they are not derived.
# - **E** — the recruitment matrix, observed. Rows are the low-concentration ("source") odorant, columns the
#   high-concentration ("target") odorant that recruits a glomerulus tuned to the source, and each row is rescaled
#   to its own [0, 1] range, so a cell reads as that source odorant's *relative* pull toward the target, not an
#   absolute probability. Grey is undefined (a row with no signal to rescale), not zero. It is computed here, not
#   inherited from §7, so the two sections cannot drift and either runs on its own. Square, on the left of row 2,
#   with its own colourbar ("recruitment probability").
# - **F** — four individual glomeruli in a 2 × 2 block, filled row by row (173, 211 / 97, 155): low plotted up and
#   high plotted down about a shared axis, bars coloured by chemical class (the plain ten-class version, not the
#   A–D aromatics split). The lower half is a magnitude, not a decrease: the tick labels are made positive again.
#   Each panel autoscales, as the draft figure does — see `SPECTRUM_SCALE`. `EG_GLOMERULI` supplies the first three
#   and `EXTRA_GLOM = 155` the fourth, by dataset glomerulus number, kept out of `EG_GLOMERULI` so other code that
#   imports it is unaffected.
#
# **The odorant order.** `WITHIN_CLASS_ORDER = "matt"` — only the pyrazine block is reshuffled to match the draft
# figure; the other classes keep the odorant table's own within-class order. Its block sequence is
# `CLASS_ORDER_MATT`'s, so it pairs with `CLASS_ORDER_FIG3 = "matt"`. Figures 4 (panels B and C), 5, S2 and the
# recruitment matrices all read these two constants, so they share one odorant order. `"burton"` — the full
# 59-odorant order in `data/burton_et_al/odorant_order.txt` — is the alternative, and is what Figure 2 uses for
# its own rows.
#
# **Which glomeruli, and how the first three were identified.** The draft figure's example spectra are individual
# glomeruli, and no record states which. They were recovered by measurement: the panels were reduced to bar
# geometry from the figure PDF's vector content, and all 376 glomeruli were scored against them. Glomeruli
# **173, 211 and 97** reproduce every bar under a single per-panel scale factor, to within 0.008 pt. As an
# independent check, the scale fitted to panel E's bar *heights* predicts that panel's axis *tick spacing* —
# 11.98 pt for a 20-unit tick, against 11.980 measured — and nothing in the fit used the ticks. (In that draft the
# panel was lettered E; it is F in the current manuscript PDF.) 211 and 97 are the two examples Figure 7B shows
# again. 155 was chosen later by the authors.
#
# **Parameters.**
#
# - `CLASS_ORDER_FIG3` — the chemical-class block order. `"matt"` (C D E K L R F N P S) is the draft figure's
#   order, verified against the figure PDF rather than inferred. `"matt_esters_first"` is the variant with esters
#   immediately after acids, annotated on that page; `"table"` is the odorant table's own order.
# - `WITHIN_CLASS_ORDER` — order *inside* a block: `"table"`, `"matt"` (pyrazine block only; the value used) or
#   `"burton"` (every block). Currently `"matt"`.
# - `FIG3_SPLIT_AROMATICS` — whether panels A–D split the "R" class into Ar. Ketone / Ar. Ketone/Ester /
#   Phenol/Methoxy. `True` by default.
# - `SPECTRUM_SCALE` — `"per-panel"` (the draft's; each spectrum autoscales) or `"shared"` (one y-limit across the
#   four). Per-panel reads shape; shared reads magnitude.
# - `ROW_RULE_FIG3 = "halfmax"`, `ROW_FRAC_FIG3 = 0.20` — panel E's rows. A row is a low-concentration odorant, and a
#   glomerulus belongs to the row of every odorant that drives it at or above `ROW_FRAC_FIG3` (a fifth) of its own
#   low-concentration peak.
# - `RECRUITED_FIG3 = "relative"`, `LO_FRAC_FIG3 = 0.20`, `HI_FRAC_FIG3 = 0.20` — what counts as recruited. Against its
#   own peak at each concentration, a glomerulus is recruited by an odorant when it is below `LO_FRAC_FIG3` of its low
#   peak and above `HI_FRAC_FIG3` of its high peak. The held-out prediction (`expm1(Yhat)`) goes through the same
#   test, so `PRED_THRESHOLD_FIG3` is not used. `"exact"` is the earlier rule: silent at low (exactly zero) and
#   responding at high, with the prediction rate-matched to the observed total. Self-pairs are excluded either way,
#   so the diagonal is blank. The supplement below sweeps the row fraction (0.20 to 0.80) with the same recruited rule.

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
# ### Supplement to Figure 3 — recruitment-matrix robustness to threshold
#
# Panel E depends on one choice: what fraction of its own maximum an odorant must drive a glomerulus to before that
# glomerulus counts as "assigned" to it. This supplement draws the matrix at five values of that fraction
# (`SWEEP_FRACS_FIG3 = 0.20, 0.35, 0.50, 0.65, 0.80`), observed on the top row and held-out predicted below, with the
# relative recruited rule (`recruitment.matrix_sweep`). The following cell
# prints, per cutoff, the number of empty rows, the mean row size and the off-diagonal observed-vs-predicted
# correlation, which is the quantity to read for robustness.
#
# It is its own figure in the repo and is not in the manuscript figure PDF (Figures 1–8, S1–S4), which is why it keeps
# its own stem, `fig3_supp_threshold_robustness`, and a `done` entry in the manifest on its own.

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
# ## 4. Figure 4 — change in tuning with concentration by NMF clusters
#
# **What it shows.** That the seven NMF clusters are real structure and not a sorting artifact, and that they
# survive adding the glomeruli that only respond at high concentration. Panels A–C are the same 376 glomeruli looked
# at three ways: as a correlation matrix (A), as the response matrix (B), and one subject at a time (C). Panel D is
# the flow of glomeruli from the original clusters to the clusters found after the high-only glomeruli are folded in.
#
# **How to read each panel.**
#
# - **A** — glomerulus × glomerulus Pearson correlation of tuning curves at 1x. Each glomerulus is a 59-odorant
#   vector; the matrix is permuted into cluster order, so the seven bright squares on the diagonal *are* the clusters.
#   It is Figure 3's picture one axis over, on the same 0–1 scale. Within-cluster correlation runs 0.12–0.51 against a
#   between-cluster mean of −0.03, which is the panel's whole claim. The grey diagonal is not a measurement. Each of the
#   seven diagonal blocks gets an unfilled green box (`cluster_boxes=True` in `glomerulus_correlation_matrix`) so the
#   block-diagonal structure is visible without the full boundary grid (`boundaries=False`).
# - **B** — the response heatmap: odorants down, glomeruli across. `FIG4B_COL_NORMALIZE=True` peak-normalises each
#   glomerulus (column) to [0, 1] by its own maximum response, on the raw dF/F values, and draws on a **linear**
#   scale (`FIG4B_NORM_VMIN`/`FIG4B_NORM_VMAX`, 0–1) instead of log — deliberately paired, since a [0, 1] quantity has
#   no absolute floor for a log scale's dynamic range to use. Reading it: brightness is "how close is this response to
#   this glomerulus's own best odorant," not an absolute ΔF/F. `FIG4B_RAW_FLOOR` (= `FIG4_VMIN`, 1 ΔF/F) restores the
#   noise floor that normalising would erase: any cell whose *raw*, pre-normalised dF/F was below it is greyed out
#   regardless of its normalised value, so a column that never clears the floor greys out completely rather than
#   lighting up relative to its own noise. Set `FIG4B_COL_NORMALIZE=False` for raw ΔF/F on a log scale. The
#   background is an explicit grey rectangle (`empty_color`), so it renders the same on every backend.
#   **Rows and columns are sorted differently.** Columns (glomeruli) are always cluster-sorted; the vertical lines mark
#   where each cluster's block ends. Rows (odorants) follow `FIG4B_ROW_ORDER = "odorant"`: Figure 3's own odorant order
#   (`FIG4_CHEM_ORDER = CLASS_ORDER_FIG3`, within-class `WITHIN_CLASS_ORDER`, both `"matt"`), with no reference to
#   cluster membership, so the diagonal highlight (`diagonal_blocks`) is off. `class_names=True` writes the ten class
#   names on the left, one per class. `FIG4B_ROW_ORDER = "cluster"` restores the bi-clustered rows and the highlight.
# - **C** — the same data by mouse, for clusters 1 and 3, four subjects side by side. Its rows are **B's rows**
#   (asserted equal in the cell) and it takes B's per-glomerulus peak-normalisation and raw floor
#   (`per_subject_heatmap(col_normalize=..., raw_floor=...)`), so a glomerulus reads identically in both panels and
#   B's colourbar is C's. White lines (`FIG4C_LINE_*`) separate the chemical classes. The cluster structure repeats
#   within each animal, so it is not a product of pooling four fields of view.
# - **D** — the Sankey. The 118 glomeruli that are silent at low (zero for every odorant) but respond to at least one
#   odorant at high have their HIGH column copied over their LOW one, and NMF is rerun (k = 7) on the 494 glomeruli
#   with any variation (73 are flat and left out). 376 of them were in the original clustering and carry that cluster;
#   the 118 are the grey "recruited at high" source. Ribbons are coloured by source cluster. The new clusters are
#   renamed 1–7 by best one-to-one overlap with the originals (Hungarian assignment on the 376 shared glomeruli) — a
#   naming choice that does not change any count. Share of each original cluster that stays in its matched new cluster:
#   1: 100%, 2: 98%, 3: 90%, 4: 96%, 5: 83%, 6: 73%, 7: 67% (printed under the figure). The reclustering is
#   **cached** by `dev/fold_nonresponders.py` in `data/nonresponders_folded/glom_clusters.csv` and only read here
#   (`highconc.folding`); it uses the package's seeded NMF defaults (seed 706231266, 1000 iterations), the same as the
#   original clustering. The drawing is `panels.sankey.flow_diagram`.
#
# **Which two clusters in C, and why those.** The panel draws **clusters 1 and 3**. The draft caption described the
# pair by content: one "characterized by mixed selectivity for acids and aldehydes", the other showing "near complete
# selectivity for amines", which by panel-D-style numbers would be clusters 1 and 5. **Cluster 5 has no glomeruli in
# subjects 3 and 4** (19, 11, 0, 0 per subject), so drawing it leaves two of the eight sub-panels empty; cluster 3
# replaces it on the caption's own criterion, mixed against near-complete selectivity for a single class (furans in
# place of amines). The selectivity figures that backed this choice (cluster 3's top class taking 0.75 of its
# activation, against 0.69, 0.67 and 0.51 for clusters 4, 7 and 1) were derived under an earlier activation-rate metric
# and have **not** been recomputed against the current `class_cluster_glom_composition`; reread them off Figure 5B
# before relying on them. Set `FIG4C_CLUSTERS = (1, 5)` to restore the original pair and accept the empty sub-panels.
#
# **Parameters.** `FIG4_CONC_A` is which concentration panel A correlates (`"low"` is the draft's). `FIG4C_CLUSTERS`
# picks C's two clusters. `FIG4_CHEM_ORDER` is this section's one deliberate coupling to §3: it is `CLASS_ORDER_FIG3`
# itself, not an independent copy. `FIG4B_SORT_BY_CLASS` only matters when `FIG4B_ROW_ORDER = "cluster"`. The
# class-by-cluster panel that used to be D now lives in §5 (Figure 5B), next to the per-cluster tuning bars (Figure 5A).

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
# ## 5. Figure 5 — chemical-class composition of each NMF cluster, low vs high
#
# **What it shows.** Paper Fig 5A and 5B.
#
# **Panel A.** `clusters.fraction_responding(d, conc, perm)`: for every odorant and cluster, the fraction of the cluster's
# glomeruli with a response above zero, at low and at high concentration. `perm` is
# `similarity.class_sorted(d, FIG5_CHEM_ORDER, FIG5_WITHIN)` — Figure 3's odorant order (`"matt"` classes, and `"matt"`
# inside a class, which reorders the pyrazine block). Drawn with **plotly** by `panels.tuning_plotly.mirrored_tuning_figure`
# — the one panel that is not matplotlib — and exported through kaleido as its own figure,
# `fig5_tuning_by_cluster_A.{png,svg,pdf}`. One axis per cluster: low is minus the fraction (left, blue), high is plus
# (right, orange), the bar axis runs to ±`FIG5_A_XLIM`, a rule separates the chemical classes, and the cluster's glomerulus count
# is printed at the lower left. The numbers are written to `fig5_class_by_cluster_A_fraction_responsive.csv`. The
# layout is modelled on `glom_explorer.grid_mirrored_horizontal(..., estimator="fraction_nonzero")`, but uses this
# repo's cluster numbering, colours and odorant order. Static export needs `kaleido` (see `pyproject.toml`).
#
# **Panel B.** `class_cluster_glom_composition(d, conc, order)`: for
# each (chemical class, cluster) cell, the number of **distinct glomeruli** of that cluster with a nonzero response to at
# least one odorant of that class (a glomerulus responding to several odorants of a class counts once; a glomerulus
# responding to more than one class is counted in each), divided by the cluster's own maximum over classes — separately
# at low and at high concentration. Columns do not sum to 1: this is a relative magnitude within the cluster, not a
# share of a whole. This used to be panel D of Figure 4 (stacked low over high); the paper moved it here and puts the
# two concentrations side by side.
#
# **How to read it.** Rows are the ten chemical classes in `FIG5_CHEM_ORDER = CLASS_ORDER_FIG3` (§3's broad order);
# columns are clusters 1–7 in the canonical numbering. The class names and cluster numbers are coloured with the house
# class and cluster colours. The two concentrations are each normalised within themselves and share one colourbar
# (0–1, "responding glomeruli, relative to cluster's top class"); because each column's maximum is 1 at both
# concentrations the two scales coincide by construction.
#
# **Numbers.** The table behind the panel is written to `fig5_class_by_cluster_B_glom_composition.csv` (one row per
# concentration × class, one column per cluster), computed once and reused for drawing and export so the saved numbers
# are exactly what the panel draws.
#
# **Parameters.** `FIG5_CHEM_ORDER`, `FIG5_WITHIN` (A's row order), `FIG5_A_H`/`FIG5_B_H`, `FIG5_A_XLIM`, `FIG5_A_LEFT`/
# `FIG5_A_RIGHT` (space for the class and odorant names), `FIG5_A_GAP` (white space between classes), `FIG5_CMAP` (`"magma"`), `FIG5_EMPTY`, `FIG5_BOUNDARY`/`FIG5_BOUNDARY_ALPHA` (the
# faint cluster-boundary rules), `FIG5_H` (height; the width is `MAIN_TARGET_W`).

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
# ## 6. Figure 6 — spatial mapping of glomeruli by cluster
#
# **What it shows.** Paper Fig 6A: where each of the seven clusters sits on the dorsal bulb. Every map is the same
# field of view — all four subjects superimposed in raw pixel coordinates, which is meaningful because the field is
# consistent between them. One cluster's glomeruli are drawn solid red and every other glomerulus is drawn behind them
# as a hollow grey circle. Some clusters are spatially compact and some are not, and that is the panel's content. The
# rms spread of each cluster about its centroid (382–512 px) is printed under the figure with the cluster sizes
# (45, 41, 40, 75, 30, 78, 67).
#
# **Membership, not response.** The members are one solid colour (`FIG6_MEMBER_COLOR = "red"`) with no ΔF/F encoded:
# `topography(..., color=...)` ignores `conc`, `vmin`, `vmax` and `cmap`, so there is no colourbar and no concentration
# parameter in this section (it used to colour members by mean low-concentration response on a capped ramp). Headings
# are the cluster numbers 1–7, coloured with `hs.cluster_color`.
#
# **Marker treatment.** `FIG6_NONRESPONDER = "hollow"` — small, hollow, pale-grey context circles
# (`NONRESPONDER_MARKERS["hollow"]`); members carry a thin black outline (`FIG6_CLUSTER_EDGE = "black"`). The figure is
# scaled to `MAIN_TARGET_W` (7.0 in) with `FIG6_SCALE`, and marker sizes and line widths are rescaled by the same
# factor (`_fig6_k`; area scales as k², widths as k), or a page this much smaller would be seven overlapping blobs.
# `FIG6_CLUSTER_SIZE`, `FIG6_CLUSTER_LW`, `FIG6_CONTEXT_SIZE` and `FIG6_CONTEXT_LW` are the scaled values. These
# numbers are restated by Figure S3, so that figure draws exactly like this one.
#
# **Only panel A is built here.** Paper panel B is the bulb images (omp114), maximum-projected by odorant class
# (acids, amines, and a class-II group) for the "most sensitive", "low" and "high" conditions, and is drawn elsewhere.
# The paper also draws a bulb outline and anterior/lateral arrows, which this panel does not.

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
# ## 7. Figure 7 — PLS-R
#
# **What it shows.** The linear-model half of the manuscript: what the model predicts for two individual glomeruli (B),
# and the weight matrix it learned (D). The paper's panel A (the `Y = X · W` schematic) and panel C (three cartoons of
# independent gain, unstructured retuning and structured retuning) are drawn separately and are not built here, so the
# lettering runs **B, D**. The section used to hold three example glomeruli, the observed/predicted recruitment matrix,
# and the cartoons; none of those are in this figure any more (the observed recruitment matrix is Figure 3E).
#
# - **B — glomeruli 211 and 97**, side by side: experimental high **above** the line and held-out predicted high
#   below, the same picture `panels.similarity.mirror_spectrum` draws for Figure 3F with the low series replaced by the
#   prediction. 211 and 97 are two of Figure 3F's examples. The glomeruli come from `FIG7_GLOMERULI`, passed to
#   `pls.eg_predicted_spectra`, in the odorant order of §3 (`CLASS_ORDER_FIG3`, `WITHIN_CLASS_ORDER`). The panel's own
#   words are "measured" and "predicted". The correlation sits in each title (r = 0.78 and 0.64 on the log scale).
#   A light line across the predicted half is `pls.observed_high_floor(m, space)`: the smallest **nonzero** observed
#   high response over the whole dataset in the drawn space (1.5386 in log1p ΔF/F). The plain minimum is always 0
#   because the display clips at zero, so the smallest *nonzero* value is the lowest response the experiment measured;
#   a predicted bar beneath the line is smaller than any real response. It is drawn by
#   `prediction_spectrum(observed_floor=...)`.
# - **D — the weight matrix W**, reordered by cluster and labelled 1–7, with the block lines. **No top-K cut**: every
#   one of the 376 × 376 weights is drawn (`pls.top_weight_gloms` is no longer called). Only the *colour range* is
#   clipped — each sign at its own 99th percentile through a `TwoSlopeNorm` centred on 0 — since negatives are smaller
#   and a symmetric range washes them out. Diagonal weights scale a glomerulus by its own response; off-diagonal
#   weights are the shared co-tuning.
#
# **Layout.** One row in inches, `MAIN_TARGET_W` (7.0 in) wide: B's two panels (`FIG7_B_W`, `FIG7_B_H`, `FIG7_B_GAP`),
# then D (`FIG7_D_SIDE`, square) and its colourbar. The right margin is derived, and an assertion fails if the row ever
# does not fit, so no width change can silently change the page width.
#
# **The display space is `"log"`.** `log1p` is the space the model was fitted and scored in, so the bars and the printed
# correlation are the quantities it was optimised against. Its cost is comparability with Figure 3F, which is in raw
# ΔF/F. `"raw"` undoes `log1p` with `expm1` and un-shrinks the prediction: for glomerulus 211 one acid's prediction
# reaches **193.6** against a measured **66.2**, which flattens every other bar in the panel — a property of the
# display space, not of the model. `"sqrt"` is a middle ground. The correlations differ by space for the same
# predictions (log is highest because it is the space the model was scored in), and the y-axis label follows
# `FIG7_SPACE`, so it cannot claim ΔF/F while drawing something else.
#
# Predictions near zero still draw a visible bar: above `FIG7_SUPPORT_THRESHOLD = 0.4` the model has nonzero
# predictions for **25 / 51** odorants (glomeruli 211 / 97) against **11 / 29** measured. A cross-validated linear rule
# spreads a little predicted response everywhere the population is active; `PredictionSpectrum.support` prints both
# counts below the figure.
#
# **Numbers.** `fig7_pls_B_prediction_spectra.csv` holds each odorant's observed and predicted value for the two
# glomeruli, built from the same spectra the panels draw.
#
# **Parameters.** `FIG7_SPACE`, `FIG7_GLOMERULI`, `FIG7_LIMIT` (`None` autoscales each spectrum; a shared limit puts them
# on one axis), `FIG7_SUPPORT_THRESHOLD`, `FIG7_CANONICAL` (which saved PLS run supplies W: `"baseline"`). This section
# reads no bare constant from another section. W comes from `pls.load_weights`, read from the saved canonical run;
# nothing is refitted here.

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
# ## 8. Figure 8 — motif transfer between glomerular groups
#
# **What it shows.** The direct test that replaces reading the off-diagonal weights of `W`. For every *ordered* pair
# of glomerular groups `F → T`: does adding group **F**'s mean low-concentration response to the prediction of group
# **T**'s glomeruli at high improve held-out performance, beyond what **T**'s own low-concentration mean already
# achieves?
#
#     intra    = predict T's high from  [T's low mean]
#     intra+F  = predict T's high from  [T's low mean, F's low mean]
#     ΔQ²      = Q²(intra+F) − Q²(intra)        # out-of-fold, pooled over T
#
# A positive ΔQ² says a tuning motif carried by F at low concentration reappears *inside* T at high. That is the
# manuscript's claim of structured retuning, stated as a prediction rather than as an inspection of weights.
#
# **Why not just read W.** Because `W`'s off-diagonal weights cannot identify which groups inherit whose tuning: the
# diagonal and off-diagonal predictors are 94% redundant (`cov(X)` is not diagonal), so the same held-out performance
# can be assembled many ways. The pairwise test avoids that by asking whether the *prediction* improves, which is a
# property of performance rather than of a particular parameterisation.
#
# **Panels (the paper's letters; its A is the procedure schematic, drawn elsewhere).**
#
# - **B — permutation p per pair**, clipped at `FIG8_CLIP = 0.1`, on the original sequential ramp (OrRd, passed as
#   `cmap`; `hs.SEQ` is magma elsewhere). Colour separates the transfers that happened from the ones that did not. The
#   **box** is the claim: BH q < `FIG8_ALPHA` (0.05) across all ordered off-diagonal pairs (42 for the clusters), drawn
#   in `FIG8_BOX_COLOR` (`#2b9743`, green). **Only the boxed squares carry a number** (`text_significant_only=True`); it
#   is the unclipped p, and `<1e-4` means the 10,000 draws cannot resolve it further. Axes: x "cluster (predicted at
#   high)", y "added non-source cluster / prediction improvement"; rows are the source, columns the target.
# - **C — the size of the transfer**, as a percentage of the error the target's own mean left behind
#   (`FIG8_METRIC = "pct_err"`), on the diverging ramp with zero at grey, same boxes and same numbers-only-where-boxed
#   rule. Most cells are grey. It has no y label of its own (B's applies). A legend in the same green explains the box.
#
# **The chemical-class control is computed but not drawn.** The same test groups each glomerulus by the class evoking
# its largest mean low-concentration response, keeping the nine classes with at least 15 glomeruli (`FIG8_MIN_CLASS`;
# 372 of 376, the ester-preferring group, n = 4, is excluded). It never saw the NMF factorisation, so if the transfers
# survive it they are a property of the population and not of the clustering. The cells after the figure print the
# cluster edges beside their class counterparts — **seven of the eight cluster edges have one** — and the per-group
# intra baselines, but the paper's Figure 8 has no class panels, so none is drawn (the earlier two-row layout with
# its C and D is gone). 72 ordered pairs are tested for classes.
#
# **One group is predicted worse than its own mean.** Furan-preferring glomeruli have an intra baseline of
# Q² = −0.084: predicting them from their own low-concentration group mean is *worse* than predicting the mean of the
# data, so any percentage of that baseline runs against a negative denominator. It is a real property of the group,
# and a second reason `FIG8_METRIC` defaults to `"pct_err"`: `"pct_var"` divides by that baseline and would flip the
# sign of every furan cell while inflating its magnitude. (It was most visible in the class row; see the printed
# baselines.)
#
# **Read the axes, not the shape.** The matrix is *not* symmetric and is not meant to be: 11 of the 13 class-level
# transfers are one-way. The exception is worth reporting rather than smoothing over — **amines ↔ thiazoles is
# reciprocal** at class resolution, whereas NMF cluster 5 is a strict sink. The amine → thiazole increment is the
# smallest significant transfer in either analysis (+0.013).
#
# **The null, and that it is a choice.** Significance comes from 10,000 *broken-pairing* permutations: the odorant
# order of Y is shuffled against X, which preserves the marginal distribution of high responses and all correlation
# structure within X and within Y, while destroying the low-to-high correspondence. This is the null the manuscript
# publishes. A size-matched control run on 2026-08-20 leaves only three of the eight cluster edges standing and removes
# the cluster-5 sink example — a real caveat, recorded in the methods audit, and deliberately not what this figure
# shows (`FIG8_NULL`).
#
# **Nothing here recomputes a permutation.** Both nulls are cached in `data/pairwise/` and ship with the clone;
# `transfer.pairwise` reads them and recomputes only the intra-group baselines, which are seconds. Rebuilding a cache
# is an explicit `transfer.compute_null(...)` call — 11 minutes for the clusters, 16 for the classes serially — and
# `pairwise` raises rather than silently starting one.
#
# **Layout and colours.** One row (B, C) scaled to `MAIN_TARGET_W` (7.0 in); the box colour, p-panel ramp and labels
# are constants of the cell (`FIG8_BOX_COLOR`, `FIG8_P_CMAP = "OrRd"`, `FIG8_XLABEL`, `FIG8_YLABEL`), not hard-coded in
# `panels.transfer`, which takes `box_color`, `cmap` and `text_significant_only`.

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
# ## 9. Figure S1 — NMF dimensionality choice
#
# **What it shows.** Why the clustering uses *k* = 7 at low (1x) concentration. One block, three panels:
#
# - **A** cumulative variance explained vs. *k*, for PCA, PCA on shuffled data, and NMF. The gap between real and
#   shuffled is the dimensionality worth keeping.
# - **B** the slope of A — where the curve stops paying for itself.
# - **C** adjusted Rand index across 50 NMF initialisations (`SUPP_DIM_N_SEEDS`) vs. *k*: how stable the clustering is
#   to the seed.
#
# **What is no longer drawn.** The consensus-matrix insets that sat beside C (at the chosen *k* and one past it) are
# omitted — a reviewer asked what they show — and so is the whole high-concentration block (the old D–F, *k* = 5):
# `SUPP_DIM_ROWS` now holds the 1x row only, and the 25x cache under `figures/nmf_dim_choice_25x/` is no longer read.
# The chosen *k* (7) is marked by a dotted vertical line across A, B and C, as in the paper, and is also printed under
# the figure.
#
# **This is the expensive one, and it is cached.** The sweep is NMF at *k* = 1–20 × 50 seeds, plus PCA, plus shuffled
# PCA over 20 permutations (`SUPP_DIM_N_SHUFFLES`) — minutes of compute. `dimsweep.load_cached` reads the `.npz` files
# under `figures/nmf_dim_choice_1x/data/cache/` **in place** and raises if one is missing rather than recomputing, so
# that a silent recompute cannot put the notebook out of step with the manuscript. `SUPP_DIM_K_VALUES` must match the
# sweep on disk.
#
# **Layout.** Scaled from the original hand layout so the width is exactly `SUPP_TARGET_W` (5.5 in), height cut to one
# block (`SUPP_DIM_SCALE`). Both halves come from `panels.dimchoice`; this section is composition only.

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
# ## 10. Figure S2 — per-subject reproducibility
#
# **What it shows.** That the bi-cluster block structure is not an artifact of pooling glomeruli across the four
# subjects. Everything is at **low** concentration. One panel per subject (S1–S4), each that subject's own glomeruli
# sorted by cluster; the same diagonal block structure appears independently within every subject.
#
# **How it is drawn.** `panels.clustering.per_subject_heatmap` for each subject with the **same odorant row order for
# all four**: chemical-class blocks in `SUPP_SUBJ_ODOR_ORDER = CLASS_ORDER_FIG3` and `SUPP_SUBJ_ODOR_WITHIN =
# WITHIN_CLASS_ORDER` (both `"matt"`, as Figures 3 and 4), applied by replacing the rows of a `BiclusterOrder` and
# clearing its odorant-cluster row boundaries (they described odorant clusters, which no longer match the rows).
# Columns stay cluster-sorted; white vertical lines mark where one cluster ends and the next begins, and each block's
# number is written above it in the cluster's colour. The chemical-class strip, with the class names, is on the left of
# the first heatmap (`hs.class_strip`, `side="left"`). Colour is raw ΔF/F on a **log** scale (`SUPP_SUBJ_HEAT_VMIN` 1,
# `SUPP_SUBJ_HEAT_VMAX` 140) in `SUPP_SUBJ_CMAP = "magma"`; cells below the floor fall to the grey background. The
# colourbar is shortened to 70% of the panel height, on the right.
#
# **Subjects 3 and 4 have no cluster-5 glomeruli**, so their cluster-5 blocks are empty. (That cluster is "D" under
# the legacy letters, which is why the source page's D column is the empty one.)
#
# **The 1–7 numbering is applied in the panel layer.** `per_subject_heatmap` writes the cluster number in the cluster's
# colour, so this section inherits it. The source page labels the clusters A–G; a letter does not tell you the number —
# `highconc.clusters` holds the mapping.
#
# **Layout and type size.** Scaled so the width is exactly `SUPP_TARGET_W` (5.5 in) with the canvas height cut to one
# row of panels; the subject titles and `SUPP_SUBJ_TITLE_PAD` scale with `SUPP_SUBJ_SCALE` rather than reading the house
# `hs.FS_PANEL`. The old two-panel version (heatmaps plus the topography grid) is split: the grid is Figure S3.

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
# ## 11. Figure S3 — per-subject topography of the seven clusters
#
# **What it shows.** Whether each cluster occupies the same territory in every mouse. `per_subject_topography(d, ax,
# subject, cluster)` over subjects × clusters (4 × 7): one subject's members of one cluster in the cluster's own colour,
# over a grey backdrop of every glomerulus from every subject, so all 28 maps share a frame (`bulb_limits`).
#
# **Drawn exactly like Figure 6.** The same page scale (`SUPP_SUBJ_TOPO_SCALE = MAIN_TARGET_W / 16.0`, so each map is
# Figure 6's size), marker size, black outline width and pale hollow context (`CONTEXT_COLOR` "darkgray" at
# `alpha = 0.5`), restated here as this section's own constants rather than read from §6. Only the fill differs: it is
# the cluster's own colour here, where Figure 6's members are all red. For that match the figure is `MAIN_TARGET_W`
# (7.0 in) wide, not the supplement's `SUPP_TARGET_W`; its height (`SUPP_SUBJ_TOPO_H = 5.2` in) holds four rows of maps.
#
# **Empty cells are a result, not a failure.** Subjects 3 and 4 have no cluster-5 glomeruli, so those two maps show the
# `—` placeholder (`empty_marker`). Cluster numbers head the columns, in the cluster colours, on the first row only; the
# subject labels S1–S4 are the row labels.
#
# The paper's version also draws the bulb outline and anterior/lateral arrows, which this panel does not.

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
# ## 12. Figure S4 — parity of the linear model
#
# **What it shows.** How well the PLS-R model's held-out predictions track the observed high-concentration responses at
# the finest grain: every responding glomerulus × odorant pair (observed > 0), as a hexbin. *Q²* = 0.27 over
# n = 3,405 pairs in the paper; the count is a property of the data and the cell prints it with the value.
#
# **One panel.** The earlier figure had three at increasing levels of averaging — cluster × odorant means (r = 0.80),
# cluster × chemical-class block means (r = 0.89), and this one. The first two are no longer drawn (`pls.cluster_by_odorant`
# and `pls.class_by_cluster` are not called). Per-glomerulus-mean parity was never shown: its *r* ≈ 1.0 follows from
# mean-centring rather than from the model.
#
# **The pinned framing.** `responding_parity` offers two framings: `frame="pinned"` fixes the axes to 1–5 and raises the
# colour ceiling above the densest bin so the cloud's shape reads, and `frame="square"` scales a full parity frame to
# the data. This panel is the pinned one, called at its defaults. Colour is `SUPP_PAR_C_CMAP = "RdPu"` (pink to purple,
# as in the paper), passed to the panel as `cmap`.
#
# **The dotted line** is `pls.observed_high_floor(m, space="log")` drawn at that predicted value
# (`responding_parity(floor=...)`): the smallest *nonzero* observed high response in the whole dataset, 1.5386 in
# log1p ΔF/F — the same value Figure 7B draws on its spectra. The plain minimum is 0 because the display clips there, so
# the smallest nonzero value is the lowest response the experiment measured.
#
# **One naming difference, flagged rather than silently resolved.** The draft panel annotates "R² = 0.27". The statistic
# is `pls.baseline_q2`, computed out-of-fold, so it is written here as **Q²**. The value is identical to three decimals.
#
# **Layout.** `SUPP_TARGET_W / 2` wide and 2.4 in tall (`SUPP_PAR_W`, `SUPP_PAR_H`).

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
# **Note.** Every figure above is in `notebooks/outputs/` as PNG, PDF and SVG,
# plus a second, SVG-only copy in `joel_figures/` at the repo root.
# The PDFs and SVGs carry editable text and a transparent background, so either
# opens in Illustrator with the labels still real text objects.
#
# To change what a panel draws, edit its function in `highconc/panels/`; every
# figure that uses it follows. To change a layout, edit the cell that builds it.
