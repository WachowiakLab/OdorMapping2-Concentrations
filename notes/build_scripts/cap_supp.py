import sys
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from cap_lib import load, save, md


def put(cells, tag, prefix, text):
    hits = [i for i, c in enumerate(cells) if c.startswith(f'# %% [markdown] tags=["{tag}"]')
            and c.split("\n")[1].startswith(prefix)]
    assert len(hits) == 1, (tag, prefix, hits)
    cells[hits[0]] = md(tag, text)


S9_HUMAN = r"""
## 9. Figure S1 — how the number of groups was chosen

Seven groups at low concentration. Three panels:

- **A** — cumulative variance explained against the number of groups, for PCA, PCA on shuffled data, and NMF. The
  gap between the real and shuffled curves is the structure not attributable to chance.
- **B** — the slope of A: how much each additional group adds.
- **C** — adjusted Rand index across 50 NMF runs from different random starts, measuring how stable the grouping is
  to the seed.

The number chosen (7) is printed under the figure. The paper's version also shows consensus matrices at the chosen
number and one past it, and the matching analysis at high concentration; neither is drawn here.

**This one is slow and is cached.** The underlying sweep takes minutes, so it is stored on disk and read back. If
the cache is missing the notebook stops with an error rather than recomputing, because a recomputed number might
differ from the paper's.
"""

S9_LLM = r"""
## 9. Figure S1 — NMF dimensionality choice

**What it shows.** Why the clustering uses *k* = 7 at low (1x) concentration. One block, three panels:

- **A** cumulative variance explained vs. *k*, for PCA, PCA on shuffled data, and NMF. The gap between real and
  shuffled is the dimensionality worth keeping.
- **B** the slope of A — where the curve stops paying for itself.
- **C** adjusted Rand index across 50 NMF initialisations (`SUPP_DIM_N_SEEDS`) vs. *k*: how stable the clustering is
  to the seed.

**What is no longer drawn.** The consensus-matrix insets that sat beside C (at the chosen *k* and one past it) are
omitted — a reviewer asked what they show — and so is the whole high-concentration block (the old D–F, *k* = 5):
`SUPP_DIM_ROWS` now holds the 1x row only, and the 25x cache under `figures/nmf_dim_choice_25x/` is no longer read.
The chosen *k* is printed under the figure; it is not marked on the curves (the paper draws a dotted vertical line at
7, which this panel does not).

**This is the expensive one, and it is cached.** The sweep is NMF at *k* = 1–20 × 50 seeds, plus PCA, plus shuffled
PCA over 20 permutations (`SUPP_DIM_N_SHUFFLES`) — minutes of compute. `dimsweep.load_cached` reads the `.npz` files
under `figures/nmf_dim_choice_1x/data/cache/` **in place** and raises if one is missing rather than recomputing, so
that a silent recompute cannot put the notebook out of step with the manuscript. `SUPP_DIM_K_VALUES` must match the
sweep on disk.

**Layout.** Scaled from the original hand layout so the width is exactly `SUPP_TARGET_W` (5.5 in), height cut to one
block (`SUPP_DIM_SCALE`). Both halves come from `panels.dimchoice`; this section is composition only.
"""

S10_HUMAN = r"""
## 10. Figure S2 — the same structure within each mouse

Figure 4 pools glomeruli from four mice. This shows the response matrix computed within each animal separately, at
low concentration: one heatmap per mouse (S1–S4), glomeruli across and sorted into the seven groups (the numbers
along the top, in each group's colour; the light vertical lines separate them), and odorants down in the same
chemical-class order as Figures 3 and 4, with the classes named on the left. Brightness is ΔF/F on a log scale from
1 to 140, as the colourbar on the right shows.

The same block structure appears in every mouse. Group 5 is missing from mice 3 and 4, which have no glomeruli in
that group.

**If you are comparing against an older version of this figure:** the source page labels the groups A–G and this one
numbers them 1–7. The letters do not map onto the numbers in order.
"""

S10_LLM = r"""
## 10. Figure S2 — per-subject reproducibility

**What it shows.** That the bi-cluster block structure is not an artifact of pooling glomeruli across the four
subjects. Everything is at **low** concentration. One panel per subject (S1–S4), each that subject's own glomeruli
sorted by cluster; the same diagonal block structure appears independently within every subject.

**How it is drawn.** `panels.clustering.per_subject_heatmap` for each subject with the **same odorant row order for
all four**: chemical-class blocks in `SUPP_SUBJ_ODOR_ORDER = CLASS_ORDER_FIG3` and `SUPP_SUBJ_ODOR_WITHIN =
WITHIN_CLASS_ORDER` (both `"matt"`, as Figures 3 and 4), applied by replacing the rows of a `BiclusterOrder` and
clearing its odorant-cluster row boundaries (they described odorant clusters, which no longer match the rows).
Columns stay cluster-sorted; white vertical lines mark where one cluster ends and the next begins, and each block's
number is written above it in the cluster's colour. The chemical-class strip, with the class names, is on the left of
the first heatmap (`hs.class_strip`, `side="left"`). Colour is raw ΔF/F on a **log** scale (`SUPP_SUBJ_HEAT_VMIN` 1,
`SUPP_SUBJ_HEAT_VMAX` 140) in `SUPP_SUBJ_CMAP = "magma"`; cells below the floor fall to the grey background. The
colourbar is shortened to 70% of the panel height, on the right.

**Subjects 3 and 4 have no cluster-5 glomeruli**, so their cluster-5 blocks are empty. (That cluster is "D" under
the legacy letters, which is why the source page's D column is the empty one.)

**The 1–7 numbering is applied in the panel layer.** `per_subject_heatmap` writes the cluster number in the cluster's
colour, so this section inherits it. The source page labels the clusters A–G; a letter does not tell you the number —
`highconc.clusters` holds the mapping.

**Layout and type size.** Scaled so the width is exactly `SUPP_TARGET_W` (5.5 in) with the canvas height cut to one
row of panels; the subject titles and `SUPP_SUBJ_TITLE_PAD` scale with `SUPP_SUBJ_SCALE` rather than reading the house
`hs.FS_PANEL`. The old two-panel version (heatmaps plus the topography grid) is split: the grid is Figure S3.
"""

S11_HUMAN = r"""
## 11. Figure S3 — where each mouse's groups sit on the bulb

For each of the four mice (rows) and each of the seven glomerular groups (columns), where that mouse's members of
the group lie on the dorsal bulb. Members are filled circles in the group's own colour with a thin black outline;
every glomerulus from every mouse is drawn behind as a small hollow grey circle, so all 28 maps share one frame and
can be compared by position. It is drawn exactly like Figure 6, which pools the four mice into one map per group.

Where a mouse has no glomeruli in a group — group 5 in mice 3 and 4 — the map shows a short grey dash instead.
"""

S11_LLM = r"""
## 11. Figure S3 — per-subject topography of the seven clusters

**What it shows.** Whether each cluster occupies the same territory in every mouse. `per_subject_topography(d, ax,
subject, cluster)` over subjects × clusters (4 × 7): one subject's members of one cluster in the cluster's own colour,
over a grey backdrop of every glomerulus from every subject, so all 28 maps share a frame (`bulb_limits`).

**Drawn exactly like Figure 6.** The same page scale (`SUPP_SUBJ_TOPO_SCALE = MAIN_TARGET_W / 16.0`, so each map is
Figure 6's size), marker size, black outline width and pale hollow context (`CONTEXT_COLOR` "darkgray" at
`alpha = 0.5`), restated here as this section's own constants rather than read from §6. Only the fill differs: it is
the cluster's own colour here, where Figure 6's members are all red. For that match the figure is `MAIN_TARGET_W`
(7.0 in) wide, not the supplement's `SUPP_TARGET_W`; its height (`SUPP_SUBJ_TOPO_H = 5.2` in) holds four rows of maps.

**Empty cells are a result, not a failure.** Subjects 3 and 4 have no cluster-5 glomeruli, so those two maps show the
`—` placeholder (`empty_marker`). Cluster numbers head the columns, in the cluster colours, on the first row only; the
subject labels S1–S4 are the row labels.

The paper's version also draws the bulb outline and anterior/lateral arrows, which this panel does not.
"""

S12_HUMAN = r"""
## 12. Figure S4 — predicted vs measured, over every responding pair

The horizontal axis is what was measured and the vertical axis is what the model predicted, both as log(1 + ΔF/F);
perfect prediction would put every point on the diagonal. Every individual glomerulus–odorant pair that responded
at high concentration is one point — 3,405 of them — drawn as a density plot, darker where more pairs fall. *Q²* is
printed on the panel (0.27 in the paper).

The dotted horizontal line is the smallest nonzero response ever measured. Predictions below it are smaller than
any real response, which is where the model's small nonzero guesses for non-responding pairs pile up.

The paper's figure also had two averaged versions of this plot (by group and odorant, r = 0.80, and by group and
chemical class, r = 0.89); only this unaveraged one is drawn now. The fit is loosest here, since nothing is
averaged away.
"""

S12_LLM = r"""
## 12. Figure S4 — parity of the linear model

**What it shows.** How well the PLS-R model's held-out predictions track the observed high-concentration responses at
the finest grain: every responding glomerulus × odorant pair (observed > 0), as a hexbin. *Q²* = 0.27 over
n = 3,405 pairs in the paper; the count is a property of the data and the cell prints it with the value.

**One panel.** The earlier figure had three at increasing levels of averaging — cluster × odorant means (r = 0.80),
cluster × chemical-class block means (r = 0.89), and this one. The first two are no longer drawn (`pls.cluster_by_odorant`
and `pls.class_by_cluster` are not called). Per-glomerulus-mean parity was never shown: its *r* ≈ 1.0 follows from
mean-centring rather than from the model.

**The pinned framing.** `responding_parity` offers two framings: `frame="pinned"` fixes the axes to 1–5 and raises the
colour ceiling above the densest bin so the cloud's shape reads, and `frame="square"` scales a full parity frame to
the data. This panel is the pinned one, called at its defaults. Colour is `SUPP_PAR_C_CMAP = "RdPu"` (pink to purple,
as in the paper), passed to the panel as `cmap`.

**The dotted line** is `pls.observed_high_floor(m, space="log")` drawn at that predicted value
(`responding_parity(floor=...)`): the smallest *nonzero* observed high response in the whole dataset, 1.5386 in
log1p ΔF/F — the same value Figure 7B draws on its spectra. The plain minimum is 0 because the display clips there, so
the smallest nonzero value is the lowest response the experiment measured.

**One naming difference, flagged rather than silently resolved.** The draft panel annotates "R² = 0.27". The statistic
is `pls.baseline_q2`, computed out-of-fold, so it is written here as **Q²**. The value is identical to three decimals.

**Layout.** `SUPP_TARGET_W / 2` wide and 2.4 in tall (`SUPP_PAR_W`, `SUPP_PAR_H`).
"""

cells = load()
put(cells, "human", "# ## 9. Figure S1", "# " + S9_HUMAN.strip("\n"))
put(cells, "llm", "# ## 9. Figure S1", "# " + S9_LLM.strip("\n"))
put(cells, "human", "# ## 10. Figure S2", "# " + S10_HUMAN.strip("\n"))
put(cells, "llm", "# ## 10. Figure S2", "# " + S10_LLM.strip("\n"))
put(cells, "human", "# ## 11. Figure S3", "# " + S11_HUMAN.strip("\n"))
put(cells, "llm", "# ## 11. Figure S3", "# " + S11_LLM.strip("\n"))
put(cells, "human", "# ## 12. Figure S4", "# " + S12_HUMAN.strip("\n"))
put(cells, "llm", "# ## 12. Figure S4", "# " + S12_LLM.strip("\n"))
save(cells)
print("S1-S4 captions replaced")
