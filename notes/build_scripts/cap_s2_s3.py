import sys
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from cap_lib import load, save, md


def put(cells, tag, prefix, text):
    hits = [i for i, c in enumerate(cells) if c.startswith(f'# %% [markdown] tags=["{tag}"]')
            and c.split("\n")[1].startswith(prefix)]
    assert len(hits) == 1, (tag, prefix, hits)
    cells[hits[0]] = md(tag, text)


S2_HUMAN = r"""
## 2. Figure 2 — odorant responses across Burton's own glomeruli, at three concentrations

Panel A of the paper's Figure 2: three response matrices over the same 59 odorants, stacked top to bottom
at the "most sensitive", "low" and "high" concentrations. These are Burton et al.'s own glomeruli rather
than this dataset's 376 — the "most sensitive" matrix has 76 glomeruli (a different, smaller preparation,
`omp114L`/`omp8x`) and the other two have 113.

- **Each cell** is one odorant–glomerulus pair; brighter means a bigger response. All three matrices share
  one colour scale. Odorants run down the rows grouped by chemical class, with the coloured strip and the
  class names on the right of each matrix — the same convention as Figure 3's strip, turned onto the row
  axis since odorants are rows here. Columns (glomeruli) carry no identity or ordering information and are
  left unlabelled.
- **The matrices get denser as concentration rises**: 3.8% of cells are nonzero at "most sensitive", 4.9%
  at low and 12.1% at high (the counts are printed under the figure).
- **The top matrix is narrower** (76 columns against 113) only because its source preparation has fewer
  glomeruli, not because anything was cropped.

The rest of the paper's Figure 2 — the sorted response spectra, the histogram of glomeruli per odorant, and
the lifetime-sparseness panels (B–E) — is drawn elsewhere.
"""

S2_LLM = r"""
## 2. Figure 2 — odorant x glomerulus response, Burton's own dataset

**What it shows.** Paper Fig 2A: the same 59 odorants against Burton's own glomeruli — not this repo's 376 —
at three conditions, stacked top to bottom in the paper's order: "most sensitive", "low", "high". The
"most sensitive" matrix is `data/burton_et_al/most_sensitive.csv`-family data from a different
subject/preparation (`sortedmatrix_omp114L_omp8x.csv`, 76 glomeruli); "low" is `sortedmatrixlo.csv` and "high"
is `sortedmatrix.csv` (113 glomeruli each). All three were supplied 2026-09 with no header or odorant labels,
so rows are matched by position against `odorant_order.txt`, the same convention `most_sensitive_correlation`
(Figure 3B) uses. The share of nonzero cells printed under the figure is 3.8% / 4.9% / 12.1% at most-sensitive /
low / high.

**How to read it.** Rows are the 59 odorants, permuted into the same class-block order Figure 3 uses
(`class_sorted`, `order="matt"`), with the chemical-class strip and class names on the *right* of each matrix —
`highconc.style.class_strip` with `orientation="y"`. All three share one colour scale. Columns are Burton's own
glomeruli in whatever order each source file ships them; there is no identity or ordering information for them,
so none is imposed, and they are left unlabelled. The top matrix's narrower column count (76 vs 113) reflects a
smaller source population, not a cropped shared axis.

**Panel letters.** The paper has one panel A holding all three matrices, so only the top matrix is lettered; the
other panels of the paper's Figure 2 (B–E: sorted spectra, histogram, cumulative lifetime sparseness, S_L low vs
high) are not drawn here.

**Parameters.** `FIG2_CONCS = ("most_sensitive", "low", "high")` sets the top-to-bottom order.
`FIG2_ORDER = "matt"` and `FIG2_WITHIN = "burton"` are this section's own row-ordering parameters. Note that
Figure 3 now uses `WITHIN_CLASS_ORDER = "matt"` (only the pyrazine block reshuffled), whereas this figure
keeps the full 59-odorant `"burton"` order inside each class; the class blocks are identical, only the
within-class order differs. Every section owns its constants, so Figure 2 does not depend on Figure 3 having
been run.
"""

S3_HUMAN = r"""
## 3. Figure 3 — odorant–odorant similarity, the recruitment matrix, and four example glomeruli

Panel A is a chemical-structure similarity for the 59 odorants (MACCS fingerprints, cosine), which does not
change with concentration. Panels B, C and D ask, for every pair of odorants, how similar the two patterns of
glomerular activity are — at the "most sensitive" condition (B), at low concentration (C) and at high (D).
Average similarity between different odorants is 0.03 at low and 0.10 at high, and the low and high matrices
agree with each other at a rank correlation of 0.55. Chemical similarity tracks response similarity more
closely at high concentration than at low (rank correlation 0.21 → 0.31).

**Reading order is A, B, C, D, left to right** — the fixed chemical ruler first, then the three
concentration-dependent conditions it is judged against. Panels E and F are the second row.

**Reading the panels.**

- **A–D** — every cell is one odorant pair; brighter means more alike. All four share one colour scale, so
  differences between them are not a scaling artifact. Odorants are grouped along both edges by chemical class
  (the coloured strip carries the class names under every panel), with the aromatics split into three finer
  groups (Ar. Ketone, Ar. Ketone/Ester, Phenol/Methoxy). The grey diagonal is each odorant compared with
  itself, which is 1 by definition; four odorants (D15, E1, P1, P4) have no data at the "most sensitive"
  condition, so panel B reads them at the colour floor rather than as a gap.
- **E** — recruitment. Each row is an odorant that drove a glomerulus at low concentration; each column an
  odorant that newly drives it at high. The cell is the fraction of that row's glomeruli recruited. Grey is
  undefined, not zero.
- **F** — four individual glomeruli (173 and 211 on top, 97 and 155 below). Each bar is one odorant: the
  low-concentration response points up, the high-concentration response points down. Bar colours are chemical
  class — the plain ten-class version; the aromatics split is a panel A–D thing only. The downward bars are
  magnitudes, not negative numbers.

**Panel B's data.** `data/burton_et_al/most_sensitive.csv` and `odorant_order.txt`, supplied 2026-09-02 — the
"most sensitive" condition this dataset didn't otherwise have.

**The example glomeruli** are numbers 173, 211, 97 and 155. Nothing on record said which glomeruli the original
panels used for the first three, so they were identified by measuring the bars in the draft PDF and matching them
against all 376; they match to within the PDF's own rounding. 211 and 97 appear again as the predicted-versus-
measured examples in Figure 7B. Glomerulus 155 was added later, by its number in the dataset.

**Panel letters follow the manuscript figure PDF** (`allfigs_hiconcms_penultimate_102126.pdf`).
"""

S3_LLM = r"""
## 3. Figure 3 — correlational structure of odor representations vs. concentration

**What it shows.** That raising concentration does not just make responses bigger, it makes the
*relationships between odorants* different. Panels C and D are the same 59 × 59 odorant × odorant correlation
matrix at the two concentrations, and the high-concentration matrix is visibly denser: the mean off-diagonal
correlation goes from 0.03 to 0.10, and the two matrices agree with each other only at Spearman 0.55. Panel B
is the same matrix at a third, "most sensitive" condition, supplied 2026-09-02 in
`data/burton_et_al/most_sensitive.csv`. Panel A is the chemical similarity of the same 59 odorants, which of
course does not change — so it is the fixed ruler the other three are measured against. Chemistry predicts
response similarity **better at high concentration than at low** (Spearman 0.21 → 0.31).

**Panel lettering follows the manuscript figure PDF** (Figure 3, A–F): A MACCS, B most-sensitive, C low, D high,
E recruitment, F the example spectra. The notebook's earlier lettering (A–D shuffled, E–H spectra, I
recruitment) is gone; the draw order is now simply left to right, with no switch for it.

**How to read each panel.**

- **A** — MACCS cosine similarity from the odorants' SMILES. Cosine rather than Tanimoto so that it is the same
  normalised-inner-product quantity as the correlations beside it. The 59 × 59 matrix is cached to `data/chem/`,
  so drawing this needs no RDKit; only rebuilding the cache does.
- **B, C, D** — odorant × odorant Pearson correlation between response maps (each odorant is a 376-glomerulus
  vector), at the "most sensitive" condition, at 1x and at 25x. All four of A–D are pinned to the same 0–1 scale,
  and they must stay pinned: B and C are sparse and D is not, so autoscaling them separately would erase the
  comparison the figure exists to make. The grey diagonal is *not a measurement* — self-similarity is 1 by
  construction — while a pale cell is a real, near-zero correlation. Panel B ships with two extra, unnamed
  rows/columns beyond the 59 in `odorant_order.txt`; `most_sensitive_correlation` drops them. Four odorants
  (D15, E1, P1, P4) have no data at this condition — filled at 0.0 rather than left NaN, so they read at the
  colour floor instead of as a hole, while the true diagonal is still re-masked to NaN like every other panel.
- **Strip labels are on all four of A–D** (`strip_labels=True`): each panel carries the coloured class strip with
  the class names, since all four share one class order.
- **A–D's aromatics breakdown.** The single "R" (aromatics) class is split into three display subclasses for the
  axis strip and boundaries only — `highconc.clusters.AROMATIC_SUBCLASS_BY_CODE`, applied via
  `similarity.split_block`: Ar. Ketone (R27), Ar. Ketone/Ester (R11, R12, R21, R31, R3, R29, R40), Phenol/Methoxy
  (R10, R23), coloured `#a32e49` / `#00bbb4` / `#d59537`. It touches nothing but these four panels — the other
  figures and panel F still see plain "R". Toggle with `FIG3_SPLIT_AROMATICS`. The diagonal block boxes that
  `correlation_matrix` can draw (`box_boundaries`) are not drawn in this layout.
- **E** — the recruitment matrix, observed. Rows are the low-concentration ("source") odorant, columns the
  high-concentration ("target") odorant that recruits a glomerulus tuned to the source, and each row is rescaled
  to its own [0, 1] range, so a cell reads as that source odorant's *relative* pull toward the target, not an
  absolute probability. Grey is undefined (a row with no signal to rescale), not zero. It is computed here, not
  inherited from §7, so the two sections cannot drift and either runs on its own. Square, on the left of row 2,
  with its own colourbar ("recruitment probability").
- **F** — four individual glomeruli in a 2 × 2 block, filled row by row (173, 211 / 97, 155): low plotted up and
  high plotted down about a shared axis, bars coloured by chemical class (the plain ten-class version, not the
  A–D aromatics split). The lower half is a magnitude, not a decrease: the tick labels are made positive again.
  Each panel autoscales, as the draft figure does — see `SPECTRUM_SCALE`. `EG_GLOMERULI` supplies the first three
  and `EXTRA_GLOM = 155` the fourth, by dataset glomerulus number, kept out of `EG_GLOMERULI` so other code that
  imports it is unaffected.

**The odorant order.** `WITHIN_CLASS_ORDER = "matt"` — only the pyrazine block is reshuffled to match the draft
figure; the other classes keep the odorant table's own within-class order. Its block sequence is
`CLASS_ORDER_MATT`'s, so it pairs with `CLASS_ORDER_FIG3 = "matt"`. Figures 4 (panels B and C), 5, S2 and the
recruitment matrices all read these two constants, so they share one odorant order. `"burton"` — the full
59-odorant order in `data/burton_et_al/odorant_order.txt` — is the alternative, and is what Figure 2 uses for
its own rows.

**Which glomeruli, and how the first three were identified.** The draft figure's example spectra are individual
glomeruli, and no record states which. They were recovered by measurement: the panels were reduced to bar
geometry from the figure PDF's vector content, and all 376 glomeruli were scored against them. Glomeruli
**173, 211 and 97** reproduce every bar under a single per-panel scale factor, to within 0.008 pt. As an
independent check, the scale fitted to panel E's bar *heights* predicts that panel's axis *tick spacing* —
11.98 pt for a 20-unit tick, against 11.980 measured — and nothing in the fit used the ticks. (In that draft the
panel was lettered E; it is F in the current manuscript PDF.) 211 and 97 are the two examples Figure 7B shows
again. 155 was chosen later by the authors.

**Parameters.**

- `CLASS_ORDER_FIG3` — the chemical-class block order. `"matt"` (C D E K L R F N P S) is the draft figure's
  order, verified against the figure PDF rather than inferred. `"matt_esters_first"` is the variant with esters
  immediately after acids, annotated on that page; `"table"` is the odorant table's own order.
- `WITHIN_CLASS_ORDER` — order *inside* a block: `"table"`, `"matt"` (pyrazine block only; the value used) or
  `"burton"` (every block). Currently `"matt"`.
- `FIG3_SPLIT_AROMATICS` — whether panels A–D split the "R" class into Ar. Ketone / Ar. Ketone/Ester /
  Phenol/Methoxy. `True` by default.
- `SPECTRUM_SCALE` — `"per-panel"` (the draft's; each spectrum autoscales) or `"shared"` (one y-limit across the
  four). Per-panel reads shape; shared reads magnitude.
- `ROW_FRAC_FIG3 = 0.20`, `LO_FRAC_FIG3 = 0.20` — panel E's recruitment thresholds, passed to
  `recruitment_matrix(lo_threshold=..., frac=...)`. In the lab's original algorithm (`threshold_effodor`) there is
  one shared fraction, used both to call a glomerulus "tuned" to a low odorant and, on the high side, to call it
  "recruited"; `LO_FRAC_FIG3` is the separate bar for "not really active at low" that recruitment requires. Each
  glomerulus is peak-normalised to [0, 1] separately per concentration. The supplement below sweeps `ROW_FRAC`
  and holds `LO_FRAC` fixed. This replaces the earlier half-max / rate-matched rule.
"""

SUP_HUMAN = r"""
### Supplement to Figure 3 — the recruitment matrix at five cutoffs

Panel E depends on one choice: how strongly an odorant must drive a glomerulus before that glomerulus counts as
assigned to it. This supplement redraws the matrix at five values of that cutoff, observed above and held-out
predicted below, and the cell after it prints the observed-versus-predicted correlation at each cutoff so you can
see how much it moves.

The paper's figure PDF has no such figure, so it is an extra here, and it saves under its own name.
"""

SUP_LLM = r"""
### Supplement to Figure 3 — recruitment-matrix robustness to threshold

Panel E depends on one choice: what fraction of its own maximum an odorant must drive a glomerulus to before that
glomerulus counts as "assigned" to it. This supplement draws the matrix at five values of that fraction
(`SWEEP_FRACS_FIG3 = 0.20, 0.35, 0.50, 0.65, 0.80`), observed on the top row and held-out predicted below, with the
low-side bar `LO_FRAC_FIG3` held at 0.20, using `recruitment.matrix_sweep(lo_threshold=...)`. The following cell
prints, per cutoff, the number of empty rows, the mean row size and the off-diagonal observed-vs-predicted
correlation, which is the quantity to read for robustness.

It is its own figure in the repo and is not in the manuscript figure PDF (Figures 1–8, S1–S4), which is why it keeps
its own stem, `fig3_supp_threshold_robustness`, and a `done` entry in the manifest on its own.
"""

cells = load()
put(cells, "human", "# ## 2. Figure 2", "# " + S2_HUMAN.strip("\n"))
put(cells, "llm", "# ## 2. Figure 2", "# " + S2_LLM.strip("\n"))
put(cells, "human", "# ## 3. Figure 3", "# " + S3_HUMAN.strip("\n"))
put(cells, "llm", "# ## 3. Figure 3", "# " + S3_LLM.strip("\n"))
put(cells, "human", "# ### Supplement to Figure 3", "# " + SUP_HUMAN.strip("\n"))
put(cells, "llm", "# ### Supplement to Figure 3", "# " + SUP_LLM.strip("\n"))
save(cells)
print("S2, S3 and supplement captions replaced")
