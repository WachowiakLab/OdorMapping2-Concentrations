import sys
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from cap_lib import load, save, md


def put(cells, tag, prefix, text):
    hits = [i for i, c in enumerate(cells) if c.startswith(f'# %% [markdown] tags=["{tag}"]')
            and c.split("\n")[1].startswith(prefix)]
    assert len(hits) == 1, (tag, prefix, hits)
    cells[hits[0]] = md(tag, text)


S4_HUMAN = r"""
## 4. Figure 4 — the seven glomerular groups, and how they carry over when high-only glomeruli are folded in

Panels A–C show the same 376 glomeruli three ways; panel D asks whether the grouping survives adding more.

- **A** — every pair of glomeruli, scored on how similar their odorant preferences are, ordered by group.
  Within-group correlation runs 0.12–0.51 against a between-group mean of −0.03. (Grey diagonal = each
  glomerulus versus itself.) A green box outlines each of the seven diagonal within-group squares.
- **B** — the response matrix behind that: odorants down, glomeruli across. Brightness is each glomerulus's
  response **relative to its own strongest odorant** (0–1, linear), not an absolute ΔF/F — so the panel answers
  "which odorants drive this glomerulus, relative to its own best one," not "how big is this response." Glomerulus
  columns are sorted by group (the vertical lines mark where one group ends and the next begins); odorant rows
  run in the same chemical-class order as Figure 3, with each of the ten classes named on the left, so a group's
  odorants are scattered rather than forming a block. A response below the noise floor (1 ΔF/F) stays grey even
  after normalising, so a glomerulus whose best response never clears the floor reads as empty, not as a falsely
  "bright" 1.0.
- **C** — the same data one mouse at a time, for groups 1 and 3, in the same row order and with the same
  per-glomerulus normalisation as B, so a glomerulus looks the same in both. Thin white lines separate the
  chemical classes.
- **D** — what happens to the groups when the glomeruli that only respond at *high* concentration are folded in.
  The 118 such glomeruli (silent at low, responding at high) are added to the low-concentration data, using their
  high response, and the clustering is run again with the same seven groups. Each ribbon follows glomeruli from
  their original group (left, 1–7; grey = the 118 newcomers, "recruited at high") to the group they land in (right,
  numbered 1–7 by which original group each new one most resembles). Glomeruli that already responded at low
  mostly keep their group — from 100% of group 1 and 98% of group 2 down to 67% of group 7 — and the newcomers
  spread across all seven.

**Which two groups are in panel C.** Group 1, which responds to a mix of acids and aldehydes, and group 3, which
responds almost only to furans. The paper contrasts a mixed group with a narrowly selective one, and these two are
the clearest such pair that is present in all four mice (group 5, the amine group, has no glomeruli in mice 3 and
4, which would leave two panels empty).

The paper's panel for chemical class by group is now Figure 5B, and the per-odorant tuning bars in the paper's
Figure 5A are drawn elsewhere.
"""

S4_LLM = r"""
## 4. Figure 4 — change in tuning with concentration by NMF clusters

**What it shows.** That the seven NMF clusters are real structure and not a sorting artifact, and that they
survive adding the glomeruli that only respond at high concentration. Panels A–C are the same 376 glomeruli looked
at three ways: as a correlation matrix (A), as the response matrix (B), and one subject at a time (C). Panel D is
the flow of glomeruli from the original clusters to the clusters found after the high-only glomeruli are folded in.

**How to read each panel.**

- **A** — glomerulus × glomerulus Pearson correlation of tuning curves at 1x. Each glomerulus is a 59-odorant
  vector; the matrix is permuted into cluster order, so the seven bright squares on the diagonal *are* the clusters.
  It is Figure 3's picture one axis over, on the same 0–1 scale. Within-cluster correlation runs 0.12–0.51 against a
  between-cluster mean of −0.03, which is the panel's whole claim. The grey diagonal is not a measurement. Each of the
  seven diagonal blocks gets an unfilled green box (`cluster_boxes=True` in `glomerulus_correlation_matrix`) so the
  block-diagonal structure is visible without the full boundary grid (`boundaries=False`).
- **B** — the response heatmap: odorants down, glomeruli across. `FIG4B_COL_NORMALIZE=True` peak-normalises each
  glomerulus (column) to [0, 1] by its own maximum response, on the raw dF/F values, and draws on a **linear**
  scale (`FIG4B_NORM_VMIN`/`FIG4B_NORM_VMAX`, 0–1) instead of log — deliberately paired, since a [0, 1] quantity has
  no absolute floor for a log scale's dynamic range to use. Reading it: brightness is "how close is this response to
  this glomerulus's own best odorant," not an absolute ΔF/F. `FIG4B_RAW_FLOOR` (= `FIG4_VMIN`, 1 ΔF/F) restores the
  noise floor that normalising would erase: any cell whose *raw*, pre-normalised dF/F was below it is greyed out
  regardless of its normalised value, so a column that never clears the floor greys out completely rather than
  lighting up relative to its own noise. Set `FIG4B_COL_NORMALIZE=False` for raw ΔF/F on a log scale. The
  background is an explicit grey rectangle (`empty_color`), so it renders the same on every backend.
  **Rows and columns are sorted differently.** Columns (glomeruli) are always cluster-sorted; the vertical lines mark
  where each cluster's block ends. Rows (odorants) follow `FIG4B_ROW_ORDER = "odorant"`: Figure 3's own odorant order
  (`FIG4_CHEM_ORDER = CLASS_ORDER_FIG3`, within-class `WITHIN_CLASS_ORDER`, both `"matt"`), with no reference to
  cluster membership, so the diagonal highlight (`diagonal_blocks`) is off. `class_names=True` writes the ten class
  names on the left, one per class. `FIG4B_ROW_ORDER = "cluster"` restores the bi-clustered rows and the highlight.
- **C** — the same data by mouse, for clusters 1 and 3, four subjects side by side. Its rows are **B's rows**
  (asserted equal in the cell) and it takes B's per-glomerulus peak-normalisation and raw floor
  (`per_subject_heatmap(col_normalize=..., raw_floor=...)`), so a glomerulus reads identically in both panels and
  B's colourbar is C's. White lines (`FIG4C_LINE_*`) separate the chemical classes. The cluster structure repeats
  within each animal, so it is not a product of pooling four fields of view.
- **D** — the Sankey. The 118 glomeruli that are silent at low (zero for every odorant) but respond to at least one
  odorant at high have their HIGH column copied over their LOW one, and NMF is rerun (k = 7) on the 494 glomeruli
  with any variation (73 are flat and left out). 376 of them were in the original clustering and carry that cluster;
  the 118 are the grey "recruited at high" source. Ribbons are coloured by source cluster. The new clusters are
  renamed 1–7 by best one-to-one overlap with the originals (Hungarian assignment on the 376 shared glomeruli) — a
  naming choice that does not change any count. Share of each original cluster that stays in its matched new cluster:
  1: 100%, 2: 98%, 3: 90%, 4: 96%, 5: 83%, 6: 73%, 7: 67% (printed under the figure). The reclustering is
  **cached** by `dev/fold_nonresponders.py` in `data/nonresponders_folded/glom_clusters.csv` and only read here
  (`highconc.folding`); it uses the package's seeded NMF defaults (seed 706231266, 1000 iterations), the same as the
  original clustering. The drawing is `panels.sankey.flow_diagram`.

**Which two clusters in C, and why those.** The panel draws **clusters 1 and 3**. The draft caption described the
pair by content: one "characterized by mixed selectivity for acids and aldehydes", the other showing "near complete
selectivity for amines", which by panel-D-style numbers would be clusters 1 and 5. **Cluster 5 has no glomeruli in
subjects 3 and 4** (19, 11, 0, 0 per subject), so drawing it leaves two of the eight sub-panels empty; cluster 3
replaces it on the caption's own criterion, mixed against near-complete selectivity for a single class (furans in
place of amines). The selectivity figures that backed this choice (cluster 3's top class taking 0.75 of its
activation, against 0.69, 0.67 and 0.51 for clusters 4, 7 and 1) were derived under an earlier activation-rate metric
and have **not** been recomputed against the current `class_cluster_glom_composition`; reread them off Figure 5B
before relying on them. Set `FIG4C_CLUSTERS = (1, 5)` to restore the original pair and accept the empty sub-panels.

**Parameters.** `FIG4_CONC_A` is which concentration panel A correlates (`"low"` is the draft's). `FIG4C_CLUSTERS`
picks C's two clusters. `FIG4_CHEM_ORDER` is this section's one deliberate coupling to §3: it is `CLASS_ORDER_FIG3`
itself, not an independent copy. `FIG4B_SORT_BY_CLASS` only matters when `FIG4B_ROW_ORDER = "cluster"`. The
class-by-cluster panel that used to be D now lives in §5 (Figure 5B). The paper's per-cluster tuning bars (Figure 5A)
are drawn elsewhere.
"""

S5_HUMAN = r"""
## 5. Figure 5 — which chemical classes each glomerular group responds to

Panel B of the paper's Figure 5. For each of the seven glomerular groups, how many of its glomeruli respond to
each chemical class, scaled so the group's top class is 1 — at low concentration (left) and at high concentration
(right). Each group has its own favourites at low concentration, and at high concentration most of them respond to
more classes: the bright cells spread.

- **Rows** are chemical classes, in Figure 3's order; **columns** are the seven groups. Class names and group
  numbers are in their own colours.
- **A cell** counts the glomeruli in that group with a nonzero response to at least one odorant of that class
  (a glomerulus that responds to five odorants of a class still counts once), then divides each column by its own
  largest cell. So columns are *relative* — a group's top class is always 1 — and do not sum to 1.
- **Low and high are each scaled within themselves**, so a bright cell means "this class reaches many of this
  group's glomeruli," not "more than at the other concentration." One colourbar serves both.

Panel A, the per-odorant tuning bars for each group, is drawn elsewhere.
"""

S5_LLM = r"""
## 5. Figure 5 — chemical-class composition of each NMF cluster, low vs high

**What it shows.** Paper Fig 5B (panel A is drawn elsewhere). `class_cluster_glom_composition(d, conc, order)`: for
each (chemical class, cluster) cell, the number of **distinct glomeruli** of that cluster with a nonzero response to at
least one odorant of that class (a glomerulus responding to several odorants of a class counts once; a glomerulus
responding to more than one class is counted in each), divided by the cluster's own maximum over classes — separately
at low and at high concentration. Columns do not sum to 1: this is a relative magnitude within the cluster, not a
share of a whole. This used to be panel D of Figure 4 (stacked low over high); the paper moved it here and puts the
two concentrations side by side.

**How to read it.** Rows are the ten chemical classes in `FIG5_CHEM_ORDER = CLASS_ORDER_FIG3` (§3's broad order);
columns are clusters 1–7 in the canonical numbering. The class names and cluster numbers are coloured with the house
class and cluster colours. The two concentrations are each normalised within themselves and share one colourbar
(0–1, "responding glomeruli, relative to cluster's top class"); because each column's maximum is 1 at both
concentrations the two scales coincide by construction.

**Numbers.** The table behind the panel is written to `fig5_class_by_cluster_B_glom_composition.csv` (one row per
concentration × class, one column per cluster), computed once and reused for drawing and export so the saved numbers
are exactly what the panel draws.

**Parameters.** `FIG5_CHEM_ORDER`, `FIG5_CMAP` (`"magma"`), `FIG5_EMPTY`, `FIG5_BOUNDARY`/`FIG5_BOUNDARY_ALPHA` (the
faint cluster-boundary rules), `FIG5_H` (height; the width is `MAIN_TARGET_W`).
"""

cells = load()
put(cells, "human", "# ## 4. Figure 4", "# " + S4_HUMAN.strip("\n"))
put(cells, "llm", "# ## 4. Figure 4", "# " + S4_LLM.strip("\n"))
put(cells, "human", "# ## 5. Figure 5", "# " + S5_HUMAN.strip("\n"))
put(cells, "llm", "# ## 5. Figure 5", "# " + S5_LLM.strip("\n"))
save(cells)
print("S4 and S5 captions replaced")
