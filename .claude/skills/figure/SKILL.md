---
name: figure
description: Regenerate or adjust one figure in this repo — "redraw Figure 4", "make Fig 6 draw membership only", "why is cluster 5 blank in two mice?", "change which clusters Fig 4C shows", "add a fifth glomerulus to Fig 3F". Knows the section-to-figure map, every FIG*_ parameter and its legal values, and the edit → rebuild → verify loop. Use for any request to change what a figure shows.
---

# Working on one figure

Every figure in the manuscript is drawn by one numbered section of
`notebooks/figures.py`. **The section number is the figure number** (§2–§8 are
Figures 2–8; §9–§12 are the supplements S1–S4), and the **panel letters are the
paper's** — taken from the manuscript PDF in the workspace's `manuscript/`. A
section draws only the panels of its figure that come from this data, so a letter
can start at B or D. There is no per-figure script.

## Find the section

| § | Figure (panels drawn) | Stem | Panel module |
|---|---|---|---|
| 1 | Setup — style, data, the PLS fit | — | — |
| 2 | **Fig 2** (A) — Burton's response matrices | `fig2_response_matrix` | `panels.similarity` |
| 3 | **Fig 3** (A–F) + a threshold supplement | `fig3_correlational_structure`, `fig3_supp_threshold_robustness` | `panels.similarity`, `panels.recruitment` |
| 4 | **Fig 4** (A–D) — clusters, per-subject, folded-in Sankey | `fig4_tuning_by_cluster` | `panels.clustering`, `panels.sankey` |
| 5 | **Fig 5** (A, B) — A: tuning bars by cluster (plotly); B: class × cluster, low beside high | `fig5_tuning_by_cluster_A`, `fig5_class_by_cluster` | `panels.tuning_plotly`, `panels.clustering` |
| 6 | **Fig 6** (A) — clusters on the dorsal bulb | `fig6_spatial_by_cluster` | `panels.clustering` |
| 7 | **Fig 7** (B, D) — PLS-R spectra and W | `fig7_pls` | `panels.pls`, `panels.similarity` |
| 8 | **Fig 8** (B, C) — motif transfer | `fig8_motif_transfer` | `panels.transfer` |
| 9 | **Fig S1** — NMF dimensionality (low concentration) | `supp_nmf_dim_choice` | `panels.dimchoice` |
| 10 | **Fig S2** — per-subject response heatmaps | `supp_per_subject_heatmaps` | `panels.clustering` |
| 11 | **Fig S3** — per-subject topography | `supp_per_subject_topography` | `panels.clustering` |
| 12 | **Fig S4** — parity over every responding pair | `supp_parity` | `panels.pls` |

Drawn elsewhere, not here: Figure 1, Fig 2B–E, Fig 6B, Fig 7A and 7C, Fig 8A.
The old concentration-flow supplement is retired (not in the paper draft).

## Decide which layer to edit

Three layers, and picking the wrong one is the usual mistake:

- **A constant at the top of the section** — the section's `FIG*_` / `SUPP_*_` block.
  Most requests land here. Every section owns its constants and shares no bare name
  with another, so changing one cannot leak into a different figure. (The couplings
  that *are* deliberate are named below.)
- **Composition, in the notebook cell** — grids, colorbars, panel letters, figure
  size. Layout requests land here.
- **A panel function in `highconc/panels/`** — what a panel *draws*. Editing one
  changes **every figure that uses it**; check before you do. Panels take an `ax`,
  never create a figure, never save, never load data.

Analysis lives in `highconc/*.py` and returns numbers. If a request needs a
different *number*, it belongs there — but see the invariants below first.

## The parameters, by section

**§2 — Figure 2**
- `FIG2_CONCS` = `("most_sensitive", "low", "high")` — top-to-bottom order of the
  three matrices, as in the paper's panel A.
- `FIG2_ORDER` = `"matt"`, `FIG2_WITHIN` = `"burton"` — this section's own row order
  (the full 59-odorant order). Note Fig 3 uses `"matt"` for the within-class order.

**§3 — Figure 3**
- `CLASS_ORDER_FIG3` — `"table"` | `"matt"` | `"matt_esters_first"`. Default
  `"matt"`. **§4, §5, §7 and §10 read it** — changing it changes all of them.
- `WITHIN_CLASS_ORDER` — `"table"` | `"matt"` (pyrazine block only) | `"burton"`
  (full order). Currently `"matt"`; also read by §4, §7 and §10.
- `FIG3_SPLIT_AROMATICS` — split the "R" class into three display subclasses on A–D.
- `FIG3_BOX_BOUNDARIES` (`(23, 35, 42, 52)`), `FIG3_BOX_COLOR` (`"green"`) — the five green diagonal
  boxes on A–D. The cut points are specific to `"matt"` order with the aromatics split; they are not derived.
- `SPECTRUM_SCALE` — `"per-panel"` | `"shared"` (panel F's four spectra).
- `EG_GLOMERULI` (from `similarity.EG_GLOMERULI`) + `EXTRA_GLOM = 155` — F's four
  glomeruli, by dataset number, in the order 173, 211 / 97, 155. Fig 7B reuses 211
  and 97.
- `ROW_RULE_FIG3` (`"halfmax"`), `ROW_FRAC_FIG3` (`0.50`), `PRED_THRESHOLD_FIG3`
  (`"rate_matched"`) — panel E's recruitment rule, the **original** one (it reproduces the paper's
  panel: three empty rows, blank diagonal). `0.20` is the comparison value; the supplement sweeps
  0.20–0.80. The lab's rewritten 0.20/0.20 rule is kept in `notes/alternatives/` and is not used.
- `SWEEP_FRACS_FIG3` — the five cutoffs in the supplement.

**§4 — Figure 4**
- `FIG4_VMIN`, `FIG4_VMAX` = `1, 140`; `FIG4_CMAP = "magma"`; `FIG4_EMPTY`.
- `FIG4_CONC_A` — `"low"` | `"high"`. Which concentration panel A shows.
- `FIG4C_CLUSTERS` — currently `(1, 3)`. Which two clusters panel C draws.
- `FIG4_CHEM_ORDER` — **is** `CLASS_ORDER_FIG3` (a deliberate coupling).
- `FIG4B_ROW_ORDER` (`"odorant"` | `"cluster"`), `FIG4B_COL_NORMALIZE` (`True`),
  `FIG4B_NORM_VMIN/VMAX`, `FIG4B_RAW_FLOOR` — B's row order and per-glomerulus
  peak-normalisation. **C reuses B's rows and normalisation** (it asserts the rows are
  equal), so change them together.
- Panel D (the Sankey) has no constant of its own: it reads the cache in
  `data/nonresponders_folded/` through `highconc.folding`.

**§5 — Figure 5**
- `FIG5_CHEM_ORDER` (= `CLASS_ORDER_FIG3`), `FIG5_WITHIN` (= `WITHIN_CLASS_ORDER`) — A's row order; `FIG5_A_XLIM`, `FIG5_A_GAP`, `FIG5_A_H`/`FIG5_B_H`; `FIG5_CMAP`.

**§6 — Figure 6**
- `FIG6_MEMBER_COLOR = "red"` — members are one solid colour; **no ΔF/F is encoded**,
  so there is no colourbar and no concentration.
- `FIG6_NONRESPONDER = "hollow"`, `FIG6_CLUSTER_EDGE = "black"`; the sizes
  `FIG6_CLUSTER_SIZE/LW`, `FIG6_CONTEXT_SIZE/LW` are scaled with the page. §11 restates
  the same numbers so S3 draws like this figure.

**§7 — Figure 7**
- `FIG7_SPACE` — `"log"` | `"sqrt"` | `"raw"`. `"log"` is the space the model is
  fitted and scored in. **Do not flip it without reading the section's caption.**
- `FIG7_GLOMERULI = (211, 97)` — B's two glomeruli (dataset numbers).
- `FIG7_LIMIT` (`None`), `FIG7_SUPPORT_THRESHOLD` (`0.4`), `FIG7_CANONICAL` (`"baseline"`).
- Layout is in inches: `FIG7_B_W/H/GAP`, `FIG7_D_SIDE`, … The right margin is derived
  and asserted, so the page stays `MAIN_TARGET_W` wide.
- D draws **every** weight of W; only the colour range is clipped (99th percentile
  per sign). The grey line across B is `pls.observed_high_floor`, not a literal.

**§8 — Figure 8**
- `FIG8_METRIC` — `"pct_err"` | `"pct_var"` | `"dq2"`; `FIG8_ALPHA` = `0.05`;
  `FIG8_CLIP` = `0.1`; `FIG8_NULL = "broken_pairing"` (the published null);
  `FIG8_MIN_CLASS = 15`.
- `FIG8_BOX_COLOR`, `FIG8_P_CMAP = "OrRd"`, `FIG8_XLABEL`, `FIG8_YLABEL` — the green
  box, the p panel's ramp, the axis labels. Only the boxed squares carry numbers.
- The chemical-class control is still computed and printed, not drawn.

**§9–§12 — supplements**
- S1: `SUPP_DIM_K_VALUES` = `range(1, 21)` — **the sweep that is on disk**; changing it
  invalidates the cache. `SUPP_DIM_ROWS` holds the 1x row only.
- S2: `SUPP_SUBJ_HEAT_VMIN/VMAX` = `1.0, 140.0` (log), `SUPP_SUBJ_CMAP`,
  `SUPP_SUBJ_ODOR_ORDER/WITHIN` (= §3's orders).
- S3: `SUPP_SUBJ_TOPO_*` — the Fig 6 marker numbers, restated; `SUPP_SUBJ_TOPO_H`.
- S4: `SUPP_PAR_C_CMAP = "RdPu"`, `SUPP_PAR_W/H`.

## The loop

1. Edit `notebooks/figures.py` — **never** `figures_human.*` or `figures_LLM.*`, which
   are generated and will be overwritten.
2. Rebuild the editions and verify:
   ```bash
   python dev/build_notebooks.py
   ```
3. Re-render and check. Use `/verify` for the full gate, or for a quick look run §1
   and then just that figure's section in a kernel.

A figure's caption lives in the source twice, tagged `human` and `llm`. **If a change
alters what the figure shows, update both captions** — they sit next to the code so
the two cannot drift.

## Invariants — say no rather than break these

- **Cluster numbering is 1–7.** The legacy A–G letters survive only as a mapping in
  `highconc.clusters` and **do not map in order**.
- **Nothing is refitted at draw time.** `dimsweep.load_cached`,
  `highconc.pls.load_weights` / `.component_selection` / `.lv_loadings` and
  `highconc.folding.load_cache` read from disk and **raise** on a miss. If you hit
  that, restore the cache; do not make it recompute. The folded-in labels are rebuilt
  deliberately with `python dev/fold_nonresponders.py` — and that script must keep the
  package's seeded NMF defaults (see `CLAUDE.md`, trap 5).
- **Figs 3A–D and 4A share a 0–1 correlation scale.** They are only comparable because
  of it.
- **Fig 7B reuses two of Fig 3F's glomeruli (211 and 97).** Do not change one without
  the other unless you mean to.
- **One odorant order.** Figs 3, 4B/C, 5, 7B and S2 all read `CLASS_ORDER_FIG3` and
  `WITHIN_CLASS_ORDER`; the figures are read against each other.
- **The panels the paper draws elsewhere are `skipped` in `dev/CONCORDANCE.yaml`.**
  Do not "fill them in".
- **The figure set is fixed by `dev/CONCORDANCE.yaml`.** Adding a figure means adding
  it to the manifest with a reason, or `concordance.py --check` fails it as scope
  creep. Removing one fails as unfinished work. Old-numbered renders left in
  `notebooks/outputs/` fail as regressions until deleted.
- **Never edit a render by hand.** If a figure needs touch-up, the touch-up belongs
  upstream in the code.

If a request is incompatible with one of these, say so and explain which — do not
quietly produce a figure that no longer matches the manuscript.

## A different skill you may also see

`style/for-collaborators/claude-skill/highconc-figure/` is a **separate, older skill
with a different job**: turning an arbitrary data file — a CSV, an `.npy` — into a
house-style panel, for someone working *outside* this repo. It is installed into
`~/.claude/skills/` by `style/for-collaborators/install-claude-skill.sh`, and it treats
`style/STYLE_GUIDE.md` as its authority.

Use **this** skill for anything inside `highconc-figures`: the figures here are already
built, and the work is editing a section of `notebooks/figures.py`, not plotting a file
from scratch.
