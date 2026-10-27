# Merge decisions

Baseline: `github/original_codebase`. Versions: `lab` = highconc-figures-lab,
`main` = highconc-figures-lab-main, `copy` = highconc-figures-lab-main - Copy.
`joel_figures` is dropped (user decision): never copied, never exported.

`final/` is seeded with the baseline. Figure sections are assembled from the
choices below once the per-figure pass is done; generated files
(`figures_human.*`, `figures_LLM.*`) are rebuilt with `dev/build_notebooks.py`,
never merged.

## §1 Setup  — decided

Source: `main` setup cell (it is a strict superset of `copy`; `lab` changed nothing).
Changes on top of it:

- Remove the `JOEL_FIGURES` folder, its `_write(...)` call, and its docstring/print mentions.
- Remove unused `import os`.
- Drop `mpl.use("Agg")` (keeps inline display as in the baseline; use `MPLBACKEND=Agg` for headless runs).
- Keep `import pandas as pd`, `MAIN_TARGET_W = 7.0`, `SUPP_TARGET_W = 5.5`,
  `CLASS_STRIP_FONTSIZE`, `HEATMAP_TITLE_FONTSIZE` and their explanatory comment
  (exact-width figure scaling is the standard; assumed, not yet explicitly confirmed).

## §2 Figure 2 (odorant x glomerulus response, Burton's own data) — decided

Only in main/copy. Source: **main's `figures_human.py`** code cells (live, newest).
main's `figures.py` has the whole section commented out; main's `figures_LLM.py` has no §2 code;
copy's `figures.py` is live but older.

- Put the live code into final `figures.py`, uncommented.
- Drop the dead `# cax = fig.add_axes([0.84, 0.15, ...])` line and the two commented `export_figure` examples.
- Prose (both editions) identical in main and copy: keep.
- Needs main's `highconc/panels/similarity.py` (`strip_fontsize` option) and `similarity.burton_response_matrix`
  (in main and copy, not in baseline), plus `data/burton_et_al/*.csv` (identical in main and copy)
  and the `fig2_response_matrix` entry in `dev/CONCORDANCE.yaml`.

## §3 Figure 3 (odorant-odorant similarity) — decided, one item open

Assembled code: `final/notes/sections/03_fig3_code.py` (built from copy's `figures.py` by
`build03.py`; syntax-checked, NOT yet run). Choices:

- Layout: **copy's** (14 x 11.6 in, 3x4 grid, `FIG3_PANEL_ORDER = "DABC"`, panel A = "most sensitive").
  NOT main's inch-based redesign, so Fig 3 does **not** hit `MAIN_TARGET_W` and has no green
  diagonal boxes (`FIG3_BOX_*`) — consequence of this choice, not yet raised with the user.
- Chemical-family strip labels on **all four** matrices A-D (`strip_labels=True`).
- `CLASS_ORDER_FIG3 = "matt"`, `WITHIN_CLASS_ORDER = "matt"` (copy and main had `"burton"`).
- 4th example glomerulus: **dataset id 155** (`EXTRA_GLOM = 155`; user's choice after lab's
  "subject 1, #8" lookup raised IndexError on real data). Spectra are E-H,
  recruitment matrix becomes **I**. Not in `EG_GLOMERULI`, so Fig 6 is unaffected.
- NOT ported from lab: `SPECTRUM_YLIM(S)` options and fixed spectrum size/fonts (not requested).
- Recruitment rule (decided by user): **main's** — `ROW_FRAC_FIG3 = 0.20`, `LO_FRAC_FIG3 = 0.20`,
  `recruitment_matrix(..., lo_threshold=LO_FRAC_FIG3, ...)` and the same in the S-sweep
  (`matrix_sweep(..., lo_threshold=...)`). `row_rule` / `prediction_threshold` are gone from §3.
  **Consequence:** this requires main's rewritten `highconc/recruitment.py` (`lo_threshold`
  kwarg); copy's and the baseline's do not take it. Check Fig 6 (and anything else calling
  `recruitment_matrix` / `matrix_sweep`) when merging that module.
- Do NOT take lab's `highconc/similarity.py`: its 3-line change splits `_clusters` into
  `_cluste` / `rs` (syntax error).
- Needs copy/main's `similarity.most_sensitive_correlation`, `similarity.split_block`,
  `cl.AROMATIC_SUBCLASS_BY_CODE`, `recruitment_matrix(..., within=, codes=)`.

## §4 Figure 4 (glomerular groups, tuning by chemical class) — decided

Assembled code: `final/notes/sections/04_fig4_code.py` (from **main's `figures.py`**, which is newer
than main's `figures_human.py`; built by `build04.py`; syntax-checked, NOT yet run).
Panel code written to `final/highconc/panels/clustering.py` (= main's file + patch, see below).

- **A, C:** main (green cluster boxes; C rows in chemical-class order with white class lines).
  lab's C change is subsumed by main's; lab's own 14-line `panels/clustering.py` edit
  (half-finished `norm`/`normalize` args) is dropped.
- **B:** rows ordered by class `"matt"` and within-class `"matt"` (`within_order=WITHIN_CLASS_ORDER`),
  column-normalised (peak per glomerulus, linear 0-1, `raw_floor = FIG4_VMIN`), class names on rows,
  diagonal highlight off. Colormap `magma` (copy and main).
- **C:** same row order as B (asserted equal to B's rows in the cell) and same normalisation.
  Needed a **new patch**: `per_subject_heatmap` gained `col_normalize` / `raw_floor` (default off, so
  S2's use is unchanged), same logic as `biclustered`. B's colorbar serves C.
- **D:** main's `class_cluster_glom_composition` metric, low and high each normalised within their own
  concentration (each cluster column / its top class) with their own colorbar; writes
  `fig4_tuning_by_cluster_D_glom_composition.csv`. Copy's `class_cluster_composition` and
  `activation_fraction` are dropped (not called anywhere in main's notebook).
- **Width:** main's exact `MAIN_TARGET_W` (7.0 in) — user did not answer; assumed, same as §1 (c).
- **To do at assembly:** add `import dataclasses` to §1 (main's notebook uses `dataclasses.replace`
  in §4 but never imports it).

## §5 Figure 5 (clusters on the dorsal bulb) — decided

No conflict: baseline, lab and copy notebook code identical; only main changed it.
Assembled code: `final/notes/sections/05_fig5_code.py` (from main's `figures.py`; `build05.py`;
syntax-checked, NOT yet run).

- Keep main's formatting: hollow non-responders, black outline, exact `MAIN_TARGET_W` width with
  marker sizes/line widths rescaled.
- **User change:** every cluster member is one solid **red** (`FIG5_MEMBER_COLOR = "red"`), showing
  membership only — no dF/F is encoded.
  - Panel patch in `final/highconc/panels/clustering.py`: `topography(..., color=None)`; when set, members
    are one solid colour and `conc/vmin/vmax/cmap` are ignored. Default behaviour unchanged.
  - Consequences, made by me: colourbar and its "mean dF/F" label removed; `FIG5_CONC` and the row label
    "low" removed (membership has no concentration); gridspec `right` 0.90 -> 0.97 since the colourbar
    margin is gone. Caption/legend text for Fig 5 still says "mean dF/F" if it does — check at prose pass.
- Related, for the analysis step: lab's `spatial_statistics.svg`, `smoothed_cluster_spatial_overlap.svg`
  and main's `dev/spatial_autocorrelation.py` look like Fig 5 spatial analyses.

## §6 Figure 6 (PLS-R) — decided: NEW DESIGN, not any version's

User decision: one row, **panels C and D only**.
Assembled: `final/notes/sections/06_fig6_code.py` (new code, based on main's `figures.py`;
`build06.py`; syntax-checked, NOT run, layout NOT seen — see environment note below).

- **C:** two example glomeruli, ids **211 and 97** (dataset glomerulus ids), experimental high up /
  predicted high down, `FIG6_SPACE = "log"`, class labels on both, stacked on the left.
  Horizontal line = smallest **nonzero** observed high response over the whole dataset in log space
  (computed by new `pls.observed_high_floor`; replaces main's hard-coded `-1.538585657976575`). Plain
  min is always 0 because `_to_space` clips at 0. User did not say per-panel vs dataset-wide: I used
  dataset-wide (one shared line). Drawn at the mirrored (negative) position, like main's line.
- **D:** the W matrix with **every weight** (no `FIG6_TOP_K` / `top_weight_gloms` cut). Colour range
  still each sign at its 99th percentile (a display scale, not a data cut) — not asked about.
- **Dropped:** old A, B (Fig 3 glomeruli), old D (recruitment matrix), E (retuning cartoons),
  `FIG6_ROW_FRAC/LO_FRAC`, the recruitment CSV exports, and the "6A/B/C vs 3E/F/G" parity cell
  (211/97 are not Fig 3's glomeruli). Only `fig6_pls_C_prediction_spectra.csv` is exported.
  `schematic` (retuning cartoons) is now unused by Fig 6 — check whether anything else uses it.
- Width: `MAIN_TARGET_W` exactly (7.0 in); layout in inches, height 3.5 in; right margin is derived.
- Code written into `final/`: `highconc/pls.py` (= main's + `observed_high_floor`) and
  `highconc/panels/similarity.py` (= main's + `observed_floor` option on `prediction_spectrum`).
- **Panels C/D letters** are C and D as the user said (not A/B). Fig 6 in the manuscript may need
  renumbering later (Step 4).
- **Manifest:** `dev/CONCORDANCE.yaml` must list only the new Fig 6 panels (C, D), with the dropped
  ones recorded as skipped/removed, or `concordance.py --check` will fail.

## §7 Figure 7 (motif transfer) — decided by user, one item open

Assembled: `final/notes/sections/07_fig7_code.py` (from main's `figures.py`; `build07.py`;
syntax-checked, NOT run). Panel code: `final/highconc/panels/transfer.py` (baseline file + patch).

- **Only A and B** (the cluster row). C and D (chemical-class control) are no longer drawn; the
  class analysis, `matching_edges` control and intra-Q² prints are still computed and printed
  (not asked to remove — can be dropped, along with `FIG7_MIN_CLASS`).
- **Significant squares are boxed in green** (`FIG7_BOX_COLOR = "#2b9743"`, main's green) — I read
  "green for significance" as the existing box outline, not a fill. lab used plain `"green"`.
- **Only significant squares carry a value** (new `text_significant_only`).
- Labels: x = `cluster (predicted at high)` on A and B; y = `added non-source cluster\nprediction
  improvement` on A only (B had no y label before; I left it). Fixed the typo "predictin".
- Width `MAIN_TARGET_W` (7.0 in), one row, base height 6.0 (so ~3.2 in tall); layout NOT seen.
- Panel patch adds `box_color` + `text_significant_only` to `cell_grid` (passed through `**kw`) and
  `box_color` to `significance_legend`. This is lab's `box_color` idea, not main's hard-coded green.
- **Correction to an earlier message of mine:** analysis `highconc/transfer.py` is identical in all
  versions and the baseline; only `panels/transfer.py` differs (lab 18 lines, main 6, copy 0).
- lab's notebook also had muted colour ramps (`FIG7_MUTE = 0.35`, p panel and B both) and
  every matrix at a fixed 1.56 in with 2-line axis labels. **User: keep the original colormap**, so none
  of lab's Fig 7 colour/size changes are taken (fixed-size matrices not addressed; not taken).
  A and B use `hs.SEQ` (p) and `hs.DIV` (B) as the style module defines them — which `SEQ` that is
  (baseline OrRd vs main/copy magma) is decided in the `style.py` merge, not here.

## §8 Figure S1 (NMF dimensionality choice) — decided

No conflict: baseline, lab and copy identical; only main changed it (exact `SUPP_TARGET_W` = 5.5 in width).
Assembled: `final/notes/sections/08_figS1_code.py` (from main's `figures.py`; `build08.py`;
syntax-checked, NOT run). `panels/dimchoice.py` and `dimsweep.py` are identical in every version.

- **User: only A, B, C** = the 1x block (variance explained, variance slope, ARI stability). The 25x
  block (D-F) is dropped; its cache stays on disk and is not read.
- The unlettered consensus matrices beside C stay (part of the same strip) — my reading, not stated.
- Height 12.5 -> 5.5 in-equivalent (x `SUPP_DIM_SCALE`); the `ax_var` title still says
  "1x concentration - chosen k = 7". Layout NOT seen.
- **OPEN (global):** `hs.SEQ` is OrRd in baseline/lab and magma in copy/main. Fig 4 asks for magma; Fig 7
  A uses `hs.SEQ` and the user said "keep the original colormap" — unclear if that means OrRd. Resolve
  in the `style.py` merge; pinning Fig 7 A to OrRd explicitly is possible.

## §9 Figure S2 (per-subject structure) — decided; SPLIT INTO TWO FIGURES

Assembled: `final/notes/sections/09_figS2_code.py` (`build09.py`; syntax-checked, NOT run, layout NOT seen).
Outputs two figures, no panel letters (each stands alone): `supp_per_subject_heatmaps` (was A) and
`supp_per_subject_topography` (was B). Needs `dev/CONCORDANCE.yaml` to list both and drop `supp_per_subject`.

- Sources: baseline = copy; lab changed A; main changed B + scaling to 5.5 in.
- **A (heatmaps):** lab's design — rows in chemical-class order, linear raw dF/F 0-100 in magma, white
  boundaries, class strip on the first heatmap, y label "odorants (by chemical class)". User: ordering
  **'matt' / 'matt'** like Fig 4 (`SUPP_SUBJ_ODOR_ORDER = CLASS_ORDER_FIG3`, `..._WITHIN = WITHIN_CLASS_ORDER`),
  not lab's hard-coded `("matt", "table")`. **Scale: user chose LOG (dF/F 1-140, ticks 1/10/100), not lab's
  linear 0-100.** Colormap stays magma (lab's; matches Fig 4) — user did not comment on the colormap.
  Width `SUPP_TARGET_W` (5.5 in), fonts as main (7 pt title, 5 pt y label).
- **B (topography):** "exact formatting as Figure 5, colours changing by cluster" — I read it as: Fig 5's marker
  size, black outline width, hollow pale-grey context (darkgray at alpha 0.5), same page scale; members keep
  their per-cluster fill (not Fig 5's red). Consequence: B is **7.0 in wide** (`MAIN_TARGET_W`), not the
  supplement's 5.5 in, so the maps are Fig 5's size; height 5.2 in is a guess.
  Also: main's "marker size / 2" styling is replaced by Fig 5's numbers.
- Panel patch in `final/highconc/panels/clustering.py`: `per_subject_topography(context_alpha=None)` only.
  (An earlier `per_subject_heatmap(norm=...)` patch for lab's linear scale was removed once log was chosen.)
- **To do at assembly (§1):** add `import dataclasses` (used by §4 and §9).

## §10 Figure S3 (direction of change through tuning space) — SKIPPED

User: not in the current draft. Nothing assembled. Open consequences, to handle at assembly:
- Decide what `final/notebooks/figures.py` does with §10: remove the section (keeping §11 numbered as is,
  since "section number = figure number") or leave the baseline cell. Not decided; user's call.
- `dev/CONCORDANCE.yaml` must record S3's outputs as skipped/removed or `concordance.py --check` fails.
- `highconc/flow.py` and `highconc/panels/flow.py` would be unused. main/copy changed `flow.py` (17 lines);
  that edit is therefore moot unless S3 returns. Keep the modules; do not delete.

## §11 Figure S4 (parity of the linear model) — decided on my reading of an earlier user message

Assembled: `final/notes/sections/11_figS4_code.py` (`build11.py`; syntax-checked, NOT run, layout NOT seen).
Panel patch: `final/highconc/panels/pls.py` (baseline file + `responding_parity(cmap=None, floor=None, floor_color)`).

- **Only panel C** (parity of all responding pairs, hexbin); A and B (cluster x odorant, cluster x class)
  are dropped, with their `cluster_by_odorant` / `class_by_cluster` calls. Single panel, no letter, stem
  `supp_parity` kept, `SUPP_TARGET_W / 2` x 2.4 in (a guess, layout not seen).
- **Purples** hexbin colormap, via `SUPP_PAR_C_CMAP = "Purples"` in the cell. main hard-coded `"Purples"` inside
  the panel; I did NOT take that (colours should not be hard-coded in panels) — made it a `cmap` argument.
- **Dotted horizontal line** at `pls.observed_high_floor(m, "log")` (smallest nonzero observed high response,
  the same value as Fig 6C's line) — `floor=` option on the panel. Predictions below it are below anything
  measured. I read "dotted horizontal line for under the observed minimum" this way.
- **Provenance of this decision:** the user's earlier message "First, only part C. Purples color map and a dotted
  horizontal line for under the observed minimum" arrived while S2 was on screen; I first took it as Fig 6.
  Finding main's hard-coded Purples in S4's panel suggests it was about S4. If it was about Fig 6, redo.
- Copy and main (`SUPP_PAR_CLASS_ORDER = "matt"`) only affected A/B, now moot. main's 5.5 in scaling not used
  (single panel). lab's change: C colormap magma (`cmap` param) and `ax.set_aspect("equal", "box")` in the
  panel — the aspect change is NOT taken (not requested).
- **lab's extra ~270 lines in §11** (Ripley's K, Moran/Geary autocorrelation, k-sweep; new deps `esda`,
  `libpysal`; placeholder `gini`; saves nothing) are NOT part of S4 and NOT assembled here. Deferred to the
  `analysis` step together with main's `dev/spatial_autocorrelation.py` (a separate implementation of the
  same statistics). `esda`/`libpysal` are not installed in this Python.

## Shared modules and tooling — decided and APPLIED to `final/` (script `shared.py`)

Taken from **main** unchanged (copy == main; lab unchanged or corrupt): `highconc/similarity.py`,
`clusters.py`, `flow.py`, `recruitment.py`, `panels/recruitment.py`, `panels/__init__.py` (drops the
`activation_fraction` export, matching main's `panels/clustering.py`), `dev/build_notebooks.py`
(explicit UTF-8 — the likely cause of lab's mojibake). Also `data/burton_et_al/*` (4 csv + odorant_order.txt;
NOT main's stray `test.ipynb`).

- `style.py` = main's (**`SEQ = magma` globally**, `class_strip(side=...)`) + lab's `class_strip(offset=...)`.
- Fig 4B's class strip now passes `offset=0` (seaborn heatmap: row i spans i..i+1). Main's default would have
  put it half a row off — my inference, confirm when rendered.
- **Fig 7 A pinned to OrRd** (`FIG7_P_CMAP`, new `cmap` option on `pvalue_matrix`): I read "keep the original
  colormap" as the baseline's. B uses `hs.DIV` (unchanged in every version).
- `pyproject.toml`: private deps `git+ssh://git@github.com/...` -> `git+https://github.com/...` (main wrote
  `git+https://git@github.com/...`; I dropped the odd `git@` user). README wording SSH -> HTTPS. `lab`'s `spatial`
  extra (esda/libpysal/pointpats) NOT added — analysis step.
- `.claude/settings.json`: baseline's exact commands + PowerShell mirrors. Dropped main's wide `:*` prefixes
  and `dev/where.py` (that file exists in no version).
- `.gitignore`: left as the baseline (no `joel_figures/` line).
- `data/pls_results/baseline/component_selection.json` and all other data/cache files: identical in every
  version (lab's differs only by CRLF). Nothing to merge.
- `environment.yml`: main's only change was a missing trailing newline; not taken.
- **Check run:** every `.py` in `final/` parses, and `highconc` + all its modules and panels import with
  `glom_explorer`/`pls_explorer`/`plotly` stubbed. This is an import check, NOT a figure render.

## STITCHED: `final/notebooks/figures.py` written (2026-10-26)

`stitch.py` built it: main's `figures.py` as the skeleton (its Fig 2 section, intro and section layout), code
cells replaced cell-for-cell from `notes/sections/*.py`, §1 edited (import `dataclasses`; removed
`mpl.use("Agg")`, `import os` and all `JOEL_FIGURES` code), Fig 2 taken live from main's `figures_human.py`
(dead `cax` line and export examples dropped), **§10 (S3) removed** (my recommendation; the user had not yet
answered "remove or keep the baseline cell" — easy to restore from the baseline `figures.py`).

- 53 cells. Compiles as one script; no name is read without being bound (AST check). Not executed.
- `figures_human.py` and `figures_LLM.py` regenerated with the repo's own `parse_cells`/`render`
  (`build_notebooks.py --check` says both are current).
- **STALE:** `notebooks/figures_human.ipynb` and `figures_LLM.ipynb` are still the BASELINE's (with the old
  figures' outputs). Regenerating them needs `jupytext`, which is not installed here; I did not install it.
  Run `python -m pip install jupytext` then `python dev/build_notebooks.py`.
- **Captions are stale.** Every changed section's markdown starts with a visible line
  "[Caption not yet updated for the merged figure ...]" (sections 3, 4, 5, 6, 7, 8, 9, 11). The intro/table of
  contents (cells 0-1) also still describe the old figure set; §3's supplement headings too.
- Build scripts that produced all of this are archived in `notes/build_scripts/` (they read scratchpad
  split files and the original folders, so they document the process rather than re-run standalone).

## Still to do (not yet applied)
1. Stitch `notes/sections/*` into `final/notebooks/figures.py` (+ §1: `import dataclasses`, drop `JOEL_FIGURES`,
   `mpl.use("Agg")`, `import os`). §10 (S3): remove or keep — user to say.
2. `dev/CONCORDANCE.yaml`: baseline ids are `fig3 fig3_supp fig4 fig5 fig6 fig7 supp_dim_choice supp_per_subject
   supp_tuning_flow supp_parity`; needs fig2, the two S2 stems, panel lists for 6/7/S1/S4, S3 skipped.
3. Captions (human + LLM) for every changed figure; `.claude/skills/*/SKILL.md` and `style/` docs mention
   removed constants and OrRd.
4. Rebuild the two notebook editions; render; `concordance.py --check`; `build_notebooks.py --check`
   (needs an environment where the private packages install).
5. `analysis\` folder, then manuscript renumbering.

## Environment note

`glom_explorer` and `pls_explorer` (private repos) are **not installed** in this Python, so nothing
can be rendered or run here yet. All assembled code is syntax-checked only. Nothing about layout,
overlap or visual results has been seen. First full run must be done where those packages install.

## Findings to remember

- `lab`'s `figures_human.py` is not valid UTF-8 (mixed cp1252 / mojibake). Do not merge from it; write UTF-8 in `final/`.
- Every folder's `figures_human.py` / `figures_LLM.py` code differs from its own `figures.py`
  (stale or hand-edited), so check them per figure before discarding.

## Pending

- Fig 2 (only in main/copy), Fig 3 … S4, shared modules, dev tooling, pyproject/environment.

## Fixes from the first real runs (user-reported tracebacks)
1. `ImportError`-like `NameError: rs` — the user's Python had `highconc` installed editable from the Box lab
   folder, whose `similarity.py` has lab's corrupted line. Not a `final/` bug. Run with
   `PYTHONPATH=...\final`, or `pip install -e ...\final --no-deps`.
2. Fig 3: `IndexError` on lab's "subject 1, #8" lookup -> `EXTRA_GLOM = 155` set by dataset id.
3. Fig 7 A: `cmap="OrRd"` was a string but `pvalue_matrix` calls it -> panel now accepts a name or a Colormap
   (my bug). Tested both panels with a stand-in `Transfer` object (not real data).

## MILESTONE: the full notebook runs (user report, 2026-10-26)
After the three fixes above, `python final\notebooks\figures.py` ran to the end on the user's machine
(with `PYTHONPATH` pointing at `final`). Renders were not inspected by me. A snapshot was taken at this point
(`notes/snapshots/final_<stamp>.zip`, newer of the two). Still open: captions, `CONCORDANCE.yaml`, `.ipynb`
editions (need jupytext), analysis step (Sankey / nonresponders-folded), manuscript renumbering.

## Fig 4E — Sankey of the non-responders-folded reclustering (added 2026-10-26)

User: fold the recruited non-responders in (no KL/stability diagnostics); Sankey only against the ORIGINAL
(non-responders-dropped) clustering, not "retained"; add it to Figure 4 as a matplotlib panel; cluster
numbers on the sides. Built (three layers, as the repo asks):
- `dev/fold_nonresponders.py` — runs the folding + NMF ONCE, writes `data/nonresponders_folded/glom_clusters.csv`
  (494 rows: subject, gloms, new_label, recruited) and `meta.json`. Rule: silent at low (all zero) AND responding
  at high -> the HIGH column replaces the LOW column (118 glomeruli); 73 flat glomeruli left out; 494 clustered.
  Needs `glom_explorer` with `load_processed_data` / `assign_cluster_indices`; run here from the user's
  `analysis/GlomExplorer` copy (read-only). Reproducible: two runs gave identical cluster sizes.
- **Seeding finding:** the notebook's `cluster_odorants` passes `params={'n_components': 7}`, which REPLACES
  the package defaults and drops the seed and the 1000-iteration cap, so the notebook's folded labels change
  on every run. The cache uses the defaults (seed 706231266, max_iter 1000 — the same as the original clustering).
  Labels will therefore differ from any earlier unseeded run.
- `highconc/folding.py` (analysis): reads the cache (raises if missing; never refits), matches on
  (subject, gloms), counts flows; new clusters renamed 1-7 by best one-to-one overlap (Hungarian) — a naming
  choice, counts unaffected. Result: 376 original + 118 added; share staying in the matched cluster
  1: 100%, 2: 98%, 3: 90%, 4: 96%, 5: 83%, 6: 73%, 7: 67%.
- `highconc/panels/sankey.py` (`flow_diagram`): matplotlib rectangles + Bezier ribbons, cluster colours, grey for
  "not in original", numbers on both sides. Registered in `panels/__init__.py`.
- Fig 4: new third row E; A-D remapped by `_y4`/`_k4` so their size is unchanged. VERIFIED by rendering with
  real data: A-D identical with and without E; E renders cleanly.
- **Found while rendering (first real look at Fig 4):** A-D have text collisions at the 7 in width — B's left
  class names overlap A and B's y label; A's colorbar label meets C's cluster titles; D's two long colorbar
  labels overlap each other and its tick labels; the "D" letter touches "acids". Pre-existing in main's
  design (not caused by E). NOT fixed — user to decide.
- Not done: caption for E; `CONCORDANCE.yaml` entry for panel E; docs (CLAUDE.md should mention the cache
  and `dev/fold_nonresponders.py` like it does for the NMF sweep caches).
- Note: the user's machine locks `notebooks/outputs/fig4_tuning_by_cluster.pdf` when it is open in a viewer;
  re-rendering then fails with PermissionError at the PDF step (PNG is written first).

## MANUSCRIPT PDF READ (2026-10-26): `manuscript/allfigs_hiconcms_penultimate_102126.pdf`, 12 pages, Illustrator, 10/21
Paper = Figs 1-8 + S1-S4 (the older legends .docx, 09/14, has Figs 1-7 + S1-S5; the PDF is newer and governs).
Method: pymupdf installed in an isolated venv under %TEMP%\pv (user's Python untouched) to render pages;
page images are in the scratchpad only.

| Paper | Panels in the paper | Repo (current) | Gap |
|---|---|---|---|
| Fig 1 | A PCA, B concentrations, C response maps | none (drawn elsewhere) | - |
| Fig 2 | A 3 stacked matrices (most-sens, low, high); B-E sparseness | §2 (3 matrices A/B/C, order high,low,most) | one panel A, reorder; B-E not in repo |
| Fig 3 | A MACCS, B most-sens, C low, D high, E recruitment, F 4 spectra (2x2: 173,211,97,155) | §3 (letters A=most-sens,B=low,C=high,D=MACCS,E-H spectra,I recruitment) | relabel; layout: row1 four matrices, row2 recruitment + 2x2 spectra |
| Fig 4 | A corr matrix, B biclustered, C clusters 1&3 per subject, D SANKEY | §4 A,B,C + my E (Sankey, third row) + D (class x cluster) | Sankey becomes D in bottom-right; old D leaves; no third row |
| Fig 5 | A tuning-by-cluster bars (Matt's), B class x cluster low | high side by side | §4 D (stacked low/high) | new figure from old 4D; horizontal; one colorbar |
| Fig 6 | A bulb maps by cluster (dark-red members), B max projection (elsewhere) | §5 | letter A; arrows ant./lat. |
| Fig 7 | A schematic, B 2 glomeruli (211, 97) side by side, C 3 cartoons, D W matrix | §6 (C stacked, D) | letters C->B, side by side |
| Fig 8 | A procedure schematic, B p-values, C % error reduction | §7 (A, B) | letters A->B, B->C |
| S1 | A, B, C (+ consensus insets, reviewer: "can omit") | §8 | ok; insets? |
| S2 | per-subject heatmaps (log, class strips, cluster numbers) | §9 S2A | ok |
| S3 | per-subject bulb maps (subject x cluster) | §9 S2B | renumber |
| S4 | parity hexbin Q2=0.27 n=3,405, dashed floor lines | §11 | ok; colormap pink |
Repo S3 (concentration flow) is not in the paper (consistent with the skip). No renumbering or relayout has been applied yet.

## RENUMBERED / RELAID OUT TO THE MANUSCRIPT (2026-10-26) — applied to `final/notebooks/figures.py`
User: use the paper's panel letters even where a neighbouring panel is not drawn here; drop S1's consensus
insets; S4 pink-to-purple. Section numbering (§1 setup, §2-8 = Figs 2-8, §9-12 = S1-S4) was my proposal, accepted
by silence. Done in stages (scripts in `notes/build_scripts/`: renumber, stage2..7); each checked by compiling
and, where possible, rendering with the real data.

| § | Paper | Stem | Letters drawn |
|---|---|---|---|
| 2 | Fig 2 | fig2_response_matrix | A (3 matrices stacked: most-sensitive, low, high) |
| 3 | Fig 3 | fig3_correlational_structure | A MACCS, B most-sens, C low, D high, E recruitment, F 2x2 spectra (173, 211, 97, 155); layout = paper's |
| 4 | Fig 4 | fig4_tuning_by_cluster | A, B, C, D = Sankey (bottom right; class x cluster summary removed; no third row) |
| 5 | Fig 5 | fig5_class_by_cluster (NEW section) | B only: class x cluster, low beside high |
| 6 | Fig 6 | fig6_spatial_by_cluster (was fig5_...) | A (cluster numbers in cluster colours as headings) |
| 7 | Fig 7 | fig7_pls (was fig6_pls) | B (211 and 97 side by side, "measured"/"predicted"), D (W, all weights) |
| 8 | Fig 8 | fig8_motif_transfer (was fig7_...) | B (p-values), C (effect size) |
| 9 | S1 | supp_nmf_dim_choice | A, B, C; consensus insets dropped; title removed |
| 10 | S2 | supp_per_subject_heatmaps | heatmaps; class names left, colourbar right |
| 11 | S3 | supp_per_subject_topography (S2B split out) | topography grid |
| 12 | S4 | supp_parity | hexbin in RdPu; dotted line at the smallest nonzero observed response |
Constants renamed FIG5_->FIG6_, FIG6_->FIG7_, FIG7_->FIG8_ (and lowercase fig6_/fig7_ variables); new FIG5_ for §5.
`dev/CONCORDANCE.yaml` rewritten to this set (backup in scratchpad); S3 (concentration flow) is `retired`.
Intro contents tables rewritten (both editions); the intro/caption PROSE is still stale (marked in each section).

Verified by rendering with real data: Figs 2, 4, 5, 6, S1, S2, S3, S4 (S4 reports n = 3,405, matching the paper;
its floor line = 1.5386, equal to the user's earlier hard-coded value).
NOT verified with real data (need the private `pls_explorer`): Fig 3E (recruitment values), Fig 7 (spectra, W) and Fig 8 —
Fig 3E and Fig 7 layouts were rendered with a TEST STAND-IN model (plain sklearn PLS / random W; harness
`notes/build_scripts/render_cells.py`, env FAKE_M), so only their layout is checked. Fig 8 not rendered at all.

Known gaps vs the paper (not done, not asked): bulb outline + ant./lat. arrows on Figs 6 and S3; Fig 4A/B text
collisions (class names over A, colourbar label vs C titles); S1's dotted vertical line at k = 7; S4's extra dashed
lines; paper Fig 2B-E, 5A, 6B, 7A/C, 8A and Fig 1 are drawn elsewhere. `fig3_supp` (threshold sweep) is not in the
paper PDF but is kept. Old-named renders in `notebooks/outputs/` (fig5_spatial_by_cluster, fig6_pls, fig7_motif_transfer,
supp_per_subject*) make `concordance.py --check` report REGRESSION until deleted / re-rendered clean.

## CAPTIONS REWRITTEN (2026-10-26)
All markdown cells of `figures.py` (intro + every section, human and LLM editions, plus the Fig 3 supplement) were
rewritten for the merged figures; the "caption not yet updated" markers are gone (0 left). Scripts:
`notes/build_scripts/cap_*.py`. Facts in them come from real runs where possible: Fig 2 nonzero shares 3.8/4.9/12.1%
(76/113/113 glomeruli); Fig 3 mean off-diagonal r 0.0287 -> 0.1003, Spearman 0.553 over 1,711 pairs; Fig 4D 376 original
+ 118 recruited, 494 clustered, 73 flat left out, share staying 100/98/90/96/83/73/67%; Fig 6 cluster sizes
45/41/40/75/30/78/67; S4 n = 3,405 and floor 1.5386. Carried over from the earlier captions (NOT recomputed here, since
they need the private model): Fig 3 chemistry-vs-response Spearman 0.21 -> 0.31; Fig 7B r = 0.78 / 0.64 and the support counts
25/51 vs 11/29; the 193.6-vs-66.2 raw-space overshoot for glomerulus 211; Fig 8's "seven of eight", 94% redundancy, null
discussion; S4's Q2 = 0.27 (stand-in model gave 0.276). Honest gaps stated in captions: S1's chosen k is printed, not
drawn; Figs 6 and S3 lack the paper's bulb outline/arrows; Fig 3 has no green diagonal boxes (the Copy layout).
Cross-reference comments in code cells were updated (§ numbers, Fig numbers, `_fig6_k`). Checked: every constant a
caption names exists in the notebook or the package.
NOT yet updated: `final/CLAUDE.md`, `README.md`, `.claude/skills/*/SKILL.md`, `style/` docs still describe the old
numbering and constants.

## DOCS PASS + FIG 4 TEXT FIXES (2026-10-27)
**Bug found and fixed:** the caption pass had doubled the leading `#` on every section heading (`# # ## 2. ...`), so
headings would have rendered as a literal "## 2." inside an H1 and section detection failed. Fixed (24 heading lines
+ the human intro title); verified: generated human edition now has `# High-concentration...` then `## 1.`..`## 12.`.
**Fig 4 text collisions fixed** (rendered with real data and looked at): A/B gap widened (wspace 0.30), B's class names at
`CLASS_STRIP_FONTSIZE` and its redundant y label removed, A's colorbar label moved beside the bar, A/C/B label fonts to the house
size, C's titles/labels to `FS_LABEL`. Layout otherwise unchanged.
**Docs rewritten to the paper's numbering / current constants:** `CLAUDE.md` (numbering statement, layout tree with
`folding.py`, `dev/fold_nonresponders.py`, new trap 5 on the seeded cache, the wrong-folder editable-install trap, HTTPS,
"not built here" list), `README.md` (figure table, cache note, examples), `.claude/skills/{figure,restyle,verify}/SKILL.md`
(figure skill rewritten from the real constants; every named constant checked against the notebook).
**NOT rewritten (open decision):** `style/STYLE_GUIDE.md`, `style/CHECKLIST.md`, `style/palette/*`, `style/for-collaborators/*`
still teach "sequential = OrRd on white". The paper and `highconc.style.SEQ` are magma. A dated status note was added to
`style/README.md` and `style/STYLE_GUIDE.md` saying so; `style/generate_assets.py` has its own OrRd and does not read the code.
Deciding whether the guide teaches magma (and regenerating palette files, collaborator prompt/skill, `style_guide.html`,
`style_reference.ipynb`) is the user's call.
Still not done: Fig 3 green diagonal boxes (paper has them; Copy layout does not); bulb outline/arrows (Figs 6, S3); S1's dotted line at k=7;
Fig 8 never rendered; the `.ipynb` editions (need jupytext); deleting old-named renders in `notebooks/outputs/`.

## MAGMA IS THE HOUSE SEQUENTIAL STANDARD; FIG 3 GREEN BOXES (2026-10-27)
User: "do the magma as the default one"; "try the green box for 3"; don't bother with bulb outlines or old-name renders.
- **Fig 3:** five green diagonal boxes on A-D (`FIG3_BOX_BOUNDARIES = (23, 35, 42, 52)`, `FIG3_BOX_COLOR = "green"`, from main's
  cut points; they are specific to class order "matt" + the aromatics split, not derived). Rendered with the stand-in model and looked at:
  matches the paper's boxes. Captions (human + LLM) and the figure skill updated.
- **style/ now teaches magma** (replaces my earlier "open decision" status note, which was removed): `STYLE_GUIDE.md` (scope,
  rule 3, §2.1, §2.5, §6 map bullet, §10), `CHECKLIST.md`, `style/README.md`, `style_guide.html` (text + the CSS gradient, with true
  magma stops), `generate_assets.py` (`SEQ = magma`, output `seq_magma_256.csv`), the collaborator prompt (`build_prompt.py` -> regenerated
  `FIGURE_PROMPT.md`), the collaborator Claude skill, and `style_reference.ipynb` (source edited and RE-EXECUTED, so its embedded
  images are magma). Palette files regenerated (`highconc_palette.csv`, `.ase`, `swatches.png/.pdf`, `seq_magma_256.csv`).
- **Documented exception kept:** Fig 8B's p-value matrix is OrRd (`FIG8_P_CMAP`); raw bulb activity maps keep the acquisition colormap.
  The diverging map is unchanged (its name still contains "orrd").
- **Left in place, not deleted (asked first):** `style/palette/seq_orrd_256.csv`, now superseded by `seq_magma_256.csv`.
- Not done by choice: bulb outlines/arrows; deleting old-name renders. Still open: S1 dotted line at k=7, Fig 8 never rendered,
  the two .ipynb editions need jupytext.

## JUPYTEXT, S1 LINE, FIG 8 RENDER (2026-10-27)
- **jupytext 1.19.6 installed** into the Python used for the build (`pip install jupytext`; it is already in the repo's `notebook` extra).
  `python dev/build_notebooks.py` now writes the two `.ipynb` editions: 45 cells, 12 sections, no outputs (they are unexecuted — running
  them needs the user's environment with the private packages). They previously were the BASELINE's, with the old figures' outputs.
  Note: `build_notebooks.py --check` compares only the `.py` editions, so it cannot tell a stale `.ipynb`; rebuild with the full command.
- **S1:** dotted vertical line at the chosen k (7) on A, B and C; captions updated; panel letters moved off the y-labels. Rendered, looked at.
- **Fig 8 rendered** (cells run with the TEST STAND-IN model, real cached nulls): layout matches the paper — B p-values (OrRd, green
  boxes, numbers only in boxes), C effect size, legend, the user's axis labels. The boxed pairs and values are the stand-in model's
  (they happen to resemble the paper's), so confirm with the real model. Class control printed (13 of 72 pairs at q < 0.05 in this run).
- Everything now renders here except what needs the real PLS model for its VALUES: Fig 3E, Fig 7 and Fig 8.

## RECRUITMENT MATRIX BACK TO THE ORIGINAL RULE (2026-10-27) — reverses the earlier "0.2 and 0.2" decision for Fig 3
User asked why the recruitment matrix differed. Measured on real data (observed side, same order):
merged rule (main's rewrite, 0.20/0.20) vs original half-max 0.50: r = 0.88, 3.5% of cells differ by > 0.2; empty rows 1 vs 3;
mean glomeruli per row 20.8 vs 14.7; the merged rule fills the diagonal (original excludes self-pairs). At the same nominal 0.50/0.50
the two ALGORITHMS disagree (r = 0.67), so the rewrite itself matters, not only the thresholds.
The paper PDF's Fig 3E shows two grey bars near the amines (3 empty rows), a blank diagonal and a colour scale to ~0.6 — i.e. the
ORIGINAL rule. The older legends doc also says "half-max response cutoff". The 0.2/0.2 merged rule did not reproduce it.
**Done:** `final/highconc/recruitment.py` = the Copy's (baseline algorithm + the within/codes ordering); restored module at 0.50 is
numerically identical to the untouched baseline (max |diff| 0.0). Fig 3 cells back to `row_rule="halfmax"`, `ROW_FRAC_FIG3 = 0.50`,
`prediction_threshold="rate_matched"`; supplement sweep back to the baseline call. main's rewrite kept at
`notes/alternatives/recruitment_main_rewrite.py`. Rendered (stand-in model for the predicted half): panel E now shows the paper's
features; colourbar tops out ~0.6.
**Comparison requested (original rule, 0.50 vs 0.20):** r = 0.913, mean |diff| 0.035, 2.2% of cells differ by > 0.2, concentrated in
the amine rows (3 empty rows at 0.50, 1 at 0.20). Picture: scratchpad `rec_compare_thresholds.png`; the Fig 3 supplement draws 0.20-0.80
with predictions (empty rows 1/2/3/5/7 at 0.20/0.35/0.50/0.65/0.80).
Captions and the figure skill updated. Note: Fig 3E's rule now matches the paper; the manuscript Methods should say half-max 0.50.
