import sys
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from cap_lib import load, save, replace_md

HUMAN = r"""
# High-concentration odor mapping — the figures

This notebook draws the figures of the paper that can be built from the data, using the paper's own
numbering and panel letters. Run it top to bottom and it rebuilds all of them into
`notebooks/outputs/` as PNG, PDF and SVG.

### The basic story the figures tell

- At three odorant concentrations — most-sensitive, low and high — the same 59 odorants drive very
  different numbers of glomeruli, and more of them as the concentration rises (**Fig 2**).
- Raising concentration changes which odorants evoke similar glomerular response patterns, not only
  how strongly they respond (**Fig 3**).
- Glomeruli group into seven clusters by shared odorant preference (**Fig 4**). The groups are stable
  when glomeruli that only respond at high concentration are folded in: most keep their group, and the
  newcomers spread across all seven (**Fig 4D**).
- Each group favours particular chemical classes, and its reach widens at high concentration (**Fig 5**),
  and the groups occupy different parts of the bulb (**Fig 6**).
- A linear model fitted across all odorants predicts held-out high-concentration responses from the
  low-concentration ones (**Fig 7**).
- The broadening is selective and asymmetric: particular groups take on the tuning of particular others,
  and not the reverse (**Fig 8**).

### What's in here

Section numbers are the paper's figure numbers, so §3 draws Figure 3; the supplements follow, §9 to §12
being Figures S1 to S4. A section draws only the panels of its figure that come from this data and code.
Everything else in the paper — Figure 1, Figure 2B–E, Fig 5A, Fig 6B, the schematics in Fig 7A and 7C, and
the procedure diagram in Fig 8A — is drawn elsewhere, so a figure here can start at B or D.

| § | Figure | In one line |
|---|---|---|
| 1 | setup | Loads the data and fits the model. Run this first. |
| 2 | **Fig 2** | Odorant responses across Burton's own glomeruli, at three concentrations |
| 3 | **Fig 3** | Odorant–odorant similarity, the recruitment matrix, and four example glomeruli |
| 4 | **Fig 4** | The seven glomerular groups, and how they carry over when high-only glomeruli are folded in |
| 5 | **Fig 5** | Which chemical classes each group responds to |
| 6 | **Fig 6** | Where those groups sit on the bulb |
| 7 | **Fig 7** | A linear model predicting the low → high change |
| 8 | **Fig 8** | Directional transfer of tuning between groups |
| 9 | **Fig S1** | How the number of groups was chosen |
| 10 | **Fig S2** | The same structure within each mouse |
| 11 | **Fig S3** | Where each mouse's groups sit on the bulb |
| 12 | **Fig S4** | Predicted vs measured, over every responding pair |

The paper has no Figure S3 from the old concentration-flow analysis, so that figure is no longer built.

### How to change something

- **A panel's contents** (what's plotted, its colours): edit the function in `highconc/panels/`. Every
  figure using it follows.
- **A figure's layout** (grid, colourbars, panel letters): edit the cell that builds it, right here in the
  notebook.
- **Fonts, colours, export settings**: `highconc/style.py`, once, for everything. No cell sets its own.

### Two things that will bite you

1. **Run it in order.** §1 loads the data and fits the model that later sections share. Running a later
   section on a fresh kernel refits it, which takes about a minute.
2. **Some results are cached on disk and are not recomputed** — the model fit, the dimensionality sweep in
   §9, and the re-clustering with high-only glomeruli folded in (§4D, built once by
   `dev/fold_nonresponders.py`). If a cache is missing the notebook stops with an error instead of quietly
   recomputing, because a recomputed number might not match the one in the paper.

Glomerular groups are numbered **1–7** throughout. Older figures and drafts label them A–G; the letters do
not map onto the numbers in order, so don't assume cluster C is cluster 3. `highconc.clusters` holds the
translation.
"""

LLM = r"""
High-concentration odor-mapping — every figure, one notebook

This notebook draws every figure in the manuscript that is built from this data. It is the only place
composition lives: the analysis is in the `highconc` package, each panel is a function in
`highconc.panels.*` that takes an `ax` and draws into it, and the grids, colorbars and panel letters below
are the part you are meant to edit.

To change a panel's *content*, edit `highconc/panels/<module>.py`.
To change a figure's *layout*, edit the cell that builds it.
To change the type scale, colours, or export contract, edit `highconc/style.py` — no cell here sets rcParams.

**Run it top to bottom.** Later sections reuse the model fitted in section 1; running a PLS section on its
own re-fits it (about a minute) and then caches.

## Contents

The numbering follows the manuscript PDF (`manuscript/allfigs_hiconcms_penultimate_102126.pdf`, Figures 1–8
and S1–S4) and uses the paper's panel letters: §*n* draws Figure *n* for 2 ≤ *n* ≤ 8 and §9–§12 draw Figures
S1–S4. A figure here draws only the panels the paper builds from this data; the rest are drawn elsewhere and
are recorded as `skipped` in `dev/CONCORDANCE.yaml`: Figure 1 entirely, Fig 2B–E, Fig 5A, Fig 6B, Fig 7A
(the X → W → Y schematic) and 7C (the three weight-matrix cartoons), and Fig 8A (the procedure schematic).
`dev/concordance.py --check` enforces that the notebook's output set is exactly the manifest's stems.

| § | Figure | Stem | Panels |
|---|---|---|---|
| 1 | Setup — style, paths, shared data and model fit | — | — |
| 2 | **Figure 2** — odorant x glomerulus response, Burton's own dataset, three concentrations | `fig2_response_matrix` | `panels.similarity` |
| 3 | **Figure 3** — correlational structure vs. concentration, recruitment, example spectra; and its threshold supplement | `fig3_correlational_structure`, `fig3_supp_threshold_robustness` | `panels.similarity`, `panels.recruitment` |
| 4 | **Figure 4** — NMF clusters, per-subject detail, and the folded-in Sankey | `fig4_tuning_by_cluster` | `panels.clustering`, `panels.sankey` (analysis: `folding`) |
| 5 | **Figure 5** — chemical-class composition of each cluster, low vs high | `fig5_class_by_cluster` | `panels.clustering` |
| 6 | **Figure 6** — spatial mapping of glomeruli by cluster | `fig6_spatial_by_cluster` | `panels.clustering` |
| 7 | **Figure 7** — PLS-R: predicted spectra and the W matrix | `fig7_pls` | `panels.similarity`, `panels.pls` |
| 8 | **Figure 8** — motif transfer between glomerular groups | `fig8_motif_transfer` | `panels.transfer` |
| 9 | **Figure S1** — NMF dimensionality choice (low concentration) | `supp_nmf_dim_choice` | `panels.dimchoice` |
| 10 | **Figure S2** — per-subject response heatmaps | `supp_per_subject_heatmaps` | `panels.clustering` |
| 11 | **Figure S3** — per-subject topography of the clusters | `supp_per_subject_topography` | `panels.clustering` |
| 12 | **Figure S4** — parity over every responding pair | `supp_parity` | `panels.pls` |

**Not built any more.** The old concentration-flow supplement (`supp_tuning_flow`) is not in the paper PDF; it
is `retired` in the manifest and its section is gone, which leaves `highconc.flow` and `panels.flow` unused.
`panels.schematic` is unused too (the retuning cartoons that sat at Fig 6E are now the paper's Fig 7C, drawn
elsewhere). `fig3_supp_threshold_robustness` is kept although the PDF has no such figure.

Each section opens with a markdown cell describing what the figure shows and how it was made. That text is
the figure legend; it sits next to the code so the two cannot drift apart.

## Things not to break

1. **The NMF dimensionality sweep is cached** (`figures/nmf_dim_choice_1x/data/cache/*.npz`; the 25x cache is
   no longer read) and takes minutes to recompute. `dimsweep.load_cached` raises on a miss rather than silently
   recomputing — leave it that way.
2. **Nothing is refitted at draw time.** `highconc.pls.load_weights`, `.component_selection` and `.lv_loadings`
   read the canonical run pls-explorer wrote to disk, and raise if it is missing rather than quietly
   producing different numbers.
3. **The folded-in reclustering is cached** (`data/nonresponders_folded/glom_clusters.csv`, with `meta.json`).
   `highconc.folding.load_cache` raises if it is missing. Rebuild it only deliberately, with
   `python dev/fold_nonresponders.py`: it is a seeded NMF fit and uses the package's own seed and iteration cap.
   Calling the clustering with `params={"n_components": 7}` instead would drop both and change on every run.

Cluster numbering is Joel's **1–7** everywhere. The legacy A–G letters survive only as a mapping in
`highconc.clusters`; a letter does not tell you the number.
"""

cells = load()
replace_md(cells, "human", "# # High-concentration odor mapping", "# " + HUMAN.strip("\n").split("\n")[0] + "\n" + "\n".join(HUMAN.strip("\n").split("\n")[1:])) \
    if False else None
# the first line of each cell is "# # <title>" for human and llm alike, so match on the title text
from cap_lib import md
def put(tag, title_start, text):
    hits = [i for i, c in enumerate(cells) if c.startswith(f'# %% [markdown] tags=["{tag}"]')
            and c.split("\n")[1].startswith(title_start)]
    assert len(hits) == 1, (tag, hits)
    cells[hits[0]] = md(tag, text)

put("human", "# # High-concentration odor mapping", "# " + HUMAN.strip("\n"))
put("llm", "# # High-concentration odor-mapping", "# " + LLM.strip("\n"))
save(cells)
print("intro captions replaced")
