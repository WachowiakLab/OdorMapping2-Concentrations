# OdorMapping2- concentration
Figure code for Odorant representations are structured by heterogeneous tuning and
concentration-dependence of functionally distinct olfactory sensory neuron
subpopulations.

*Jason Castro 1, *Shawn Burton2, 3, Joel Hallkaj 3, Matt Wachowiak3
1 Program in Neuroscience, and Program in Digital and Computational Studies, Bates
College, Lewiston, ME, 04240, USA.
2 Department of Biological Sciences, Lehigh University, Bethlehem, PA, 18015, USA.
3 Department of Neurobiology, University of Utah School of Medicine, Salt Lake City,
UT, 84112, USA.
* the first two authors contributed equally

Principal Code Credit: Jason Castro
From original repostory: https://github.com/WachowiakLab/highconc-figures-lab.git

Everything is drawn from one notebook, which ships in **two editions**:

| | for | |
|---|---|---|
| **`notebooks/figures_human.ipynb`** | people | what each figure shows, in plain language, and how it was made |
| **`notebooks/figures_LLM.ipynb`** | machines and close readers | full provenance: every parameter, every value, every judgement call and its reason |

**They draw the same figures from the same code.** Only the prose differs — both
are generated from `notebooks/figures.py`, which holds each figure's code once
and each caption twice, so the two editions cannot drift apart. Run either top
to bottom and it rebuilds every figure in the paper into `notebooks/outputs/` as
PNG, PDF and SVG.

Open the human edition unless you specifically want the provenance.

## Install

Needs Python ≥ 3.11 and HTTPS access to the two private workhorse repos (a GitHub account with access to CastroLab, signed in through a personal access token or the GitHub CLI credential helper).

```bash
pip install -e ".[pls,notebook]"
```

- The base install pulls in **glom-explorer**, which carries the clustering
  data *inside the package* — there is no data path to configure and nothing to
  download.
- The `pls` extra pulls in **pls-explorer**. §1 fits the baseline PLS model
  and shares it, so the notebook needs it end to end; the clustering figures
  (2, 4, 5, 6 and S1–S3) do not.
- The `chem` extra pulls in **RDKit**, needed only to *rebuild* the
  chemical-similarity cache behind Fig 3A. The 59×59 MACCS matrix is cached to
  `data/chem/` and ships with the repo, so drawing the figure never needs RDKit.
- Both dependencies are pinned to a commit in `pyproject.toml`, so a fresh
  clone rebuilds the same figures. Bump those SHAs deliberately; never float
  them.

Then:

```bash
jupyter lab notebooks/figures_human.ipynb
```

**Before running a notebook headlessly, read the kernel note in
[`CLAUDE.md`](CLAUDE.md).** The saved kernelspec is `python3`, which is
probably not the environment you just installed into, and the failure looks
like a broken repo rather than a wrong kernel.

## Quickstart

```bash
python notebooks/figures.py
```

That runs the whole notebook as a plain script — same figures, no kernel. It
takes a few minutes; the PLS fit in §1 and the permutation caches are the bulk
of it, and nothing is refit afterwards.

Two results are computed once and read back, never refit at draw time: the NMF
dimensionality sweep behind Fig S1, and the folded-in reclustering behind Fig 4D
(`python dev/fold_nonresponders.py` rebuilds that one — see `CLAUDE.md`, trap 5).

To render a single figure, run the setup cells (§1) and then that figure's
section. Every section owns its constants and inherits no bare name from
another, so a section runs correctly on its own or in order.

## Where the figures are

The notebook's section numbers follow the manuscript's figure numbers (§2–§8 are
Figures 2–8, §9–§12 are S1–S4), and each figure uses the paper's own panel letters —
so a figure can start at B or D where the paper's other panels are drawn elsewhere.

| § | Manuscript figure | Stem in `notebooks/outputs/` | What it shows |
|---|---|---|---|
| 1 | — | — | Setup: style, paths, shared data, the PLS fit |
| 2 | **Fig 2** (A) | `fig2_response_matrix` | Odorant × glomerulus response, Burton's own data, at three concentrations |
| 3 | **Fig 3** (A–F) | `fig3_correlational_structure` | MACCS chemical similarity (A), odorant × odorant correlation at most-sensitive / low / high (B–D), the recruitment matrix (E), four example spectra (F) |
| 3 | extra | `fig3_supp_threshold_robustness` | The recruitment matrix at five thresholds (not in the paper PDF) |
| 4 | **Fig 4** (A–D) | `fig4_tuning_by_cluster` | Glomerulus × glomerulus correlation, bi-clustered heatmap, two clusters by subject, and the Sankey of how clusters carry over when high-only glomeruli are folded in |
| 5 | **Fig 5** (B) | `fig5_class_by_cluster` | Chemical class × cluster, low beside high |
| 5 | **Fig 5** (A) | `fig5_tuning_by_cluster_A` | Per-odorant tuning bars by cluster — plotly, exported with kaleido |
| 6 | **Fig 6** (A) | `fig6_spatial_by_cluster` | The seven clusters on the dorsal bulb |
| 7 | **Fig 7** (B, D) | `fig7_pls` | PLS-R: measured vs. predicted spectra for two glomeruli (B), and the full W matrix (D) |
| 8 | **Fig 8** (B, C) | `fig8_motif_transfer` | Motif transfer between clusters: permutation p, and size of the effect |
| 9 | **Fig S1** | `supp_nmf_dim_choice` | NMF dimensionality at low concentration — *k* = 7 |
| 10 | **Fig S2** | `supp_per_subject_heatmaps` | Per-subject response heatmaps |
| 11 | **Fig S3** | `supp_per_subject_topography` | Per-subject topography of the clusters |
| 12 | **Fig S4** | `supp_parity` | Model parity over every responding pair |

Figure 1, Fig 2B–E, Fig 6B, the schematics in Fig 7A and 7C, and the
procedure diagram in Fig 8A are drawn elsewhere and are not built here.
[`dev/CONCORDANCE.yaml`](dev/CONCORDANCE.yaml) is the manifest that fixes this
list — which figures ship, panel by panel, with every omission recorded — and
`python dev/concordance.py --check` enforces that the notebook emits exactly these
stems and no others. (The old concentration-flow supplement is retired: it is not
in the current paper draft.)

## Exporting figures

A top-to-bottom run **already writes** `png`, `pdf` and `svg` of every figure
into `notebooks/outputs/`. You do not have to do anything to get them.

For anything else, each figure cell ends with a commented `export_figure(...)`
line, and there is a general one near the top of the notebook. Uncomment and run
after the notebook has been run:

```python
export_figure()                                    # the figure just drawn
export_figure("fig7_pls")                          # one, by name
export_figure("all")                               # everything
export_figure("all", outdir="~/Desktop/highconc")  # somewhere else
export_figure("fig4_tuning_by_cluster", formats=("pdf",), dpi=600)
```

Nothing is recomputed — it re-writes figures already drawn in this kernel
session. On a fresh kernel it tells you to run the notebook first rather than
silently doing nothing.

## Editing the notebooks

**Edit `notebooks/figures.py`. Never edit `figures_human.*` or `figures_LLM.*`
— they are generated and will be overwritten.**

Captions live in the source twice, tagged:

```python
# %% [markdown] tags=["human"]
# ## 6. Figure 6 — where the seven groups sit on the bulb
...
# %% [markdown] tags=["llm"]
# ## 6. Figure 6 — spatial mapping of glomeruli by cluster
...
```

Markdown with no tag goes into both. **Code cells always go into both, verbatim
— that is the point of the generator.** Then:

```bash
python dev/build_notebooks.py           # rebuild both editions
python dev/build_notebooks.py --check   # are they current? exit 1 if not
```

Run `--check` before sharing anything, so an edit that was never rebuilt fails
loudly instead of leaving a stale notebook in someone else's clone.

## The shareable site

`python dev/build_site.py` writes **one self-contained file**,
`notebooks/site/index.html` (~9 MB): a contents page, one page per figure with
its plain-language caption and a collapsible block carrying the full
provenance, and a downloadable PDF and SVG of every figure embedded in the
page. No assets directory, no server logic, nothing else to upload. It reads
`notebooks/outputs/`, so run the notebook first; if a figure is missing it says
which and stops rather than producing a page with a hole in it. The site is
gitignored — rebuild it rather than committing it.

**Treat that file as unpublished.** It embeds the underlying values along with
the figures. It works from disk, downloads included, so it can be sent as a
file rather than hosted.

## How a figure is built

**Analysis** lives in the `highconc` package — `highconc/clusters.py`,
`highconc/pls.py`, `highconc/similarity.py`, `highconc/recruitment.py`,
`highconc/transfer.py`, `highconc/flow.py`, `highconc/dimsweep.py`.

**Drawing** lives in `highconc/panels/*`. Every panel is a function that takes
an existing `ax` and draws into it. A panel never creates a figure, never calls
`plt.show` or `savefig`, and never loads data.

**Composition** — the grids, the colorbars, the panel letters — is in the
notebook, where it is visible and editable. That is the part you are meant to
change.

So: to change what a panel *draws*, edit its function in `highconc/panels/`, and
every figure using it follows. To change a *layout*, edit the cell that builds it.

```
highconc/            the package: analysis + panels + style
notebooks/
├── figures.py            THE SOURCE — edit this; code once, captions twice
├── figures_human.ipynb   generated: the friendly edition — open this
├── figures_LLM.ipynb     generated: the exact edition
└── outputs/              renders, gitignored
data/                vendored inputs: chemical cache, PLS results, pairwise nulls
figures/             the NMF sweep caches — load-bearing
dev/                 CONCORDANCE.yaml + concordance.py, build_notebooks.py, build_site.py
style/               the style guide and its assets
```

## Principles

1. **Panels take an `ax`.** They never create a figure, call `plt.show`, or save.
2. **Style is central.** `highconc.style.apply()` pins the export contract
   (editable SVG text, TrueType in PDF) and the type scale for every figure. No
   cell sets its own rcParams.
3. **Spec-dependent choices are parameters** with documented defaults —
   activation threshold, normalisation mode, chemical-class order, cluster
   numbering. They are set in the notebook cell, not buried in a panel.
4. **Nothing is refitted at draw time.** Cached NMF sweeps and the canonical PLS
   run are read from disk, and a missing one **raises** rather than silently
   recomputing numbers that would differ from the manuscript's.
5. **Never edit a render by hand.** Every output is reproducible from the code.
   If a figure needs touch-up, the touch-up belongs upstream.

## Four things not to break

- **The kernel.** The notebooks are saved with kernelspec `python3`, which is
  probably not the environment you installed into. Register a distinctly named
  kernel — *not* `python3`, which a user-level kernelspec will outrank — and
  pass `--ExecutePreprocessor.kernel_name=` when executing headlessly. Use
  `python -m nbconvert`, not `python -m jupyter nbconvert`: the launcher runs
  whichever `jupyter-nbconvert` is first on your PATH and then cannot see your
  kernel. Full detail in [`CLAUDE.md`](CLAUDE.md).
- **The NMF dimensionality sweeps are cached** to
  `figures/nmf_dim_choice_*/data/cache/*.npz` and take minutes to recompute.
  `dimsweep.load_cached` reads them in place and raises on a miss. Leave it that
  way — a silent recompute is how a notebook stops matching the paper.
- **Cluster numbering is 1–7** everywhere. The legacy A–G letters survive only
  as a mapping in `highconc.clusters`, and they do not map in order; a letter
  does not tell you the number.
- **The canonical PLS run's result files ship with this repo**, in
  `data/pls_results/baseline/` (3.3 MB: `component_selection.json`,
  `real/loadings.npz`, and `B.npy` for `real`, `null1`, `null2`). They are
  byte-for-byte copies of the pls-explorer `baseline` run. They live here
  because `pls_explorer.config.RESULTS_DIR` points inside a *source checkout*
  of pls-explorer, and the pinned dependency installs as a wheel with no
  `results/` beside it — so on a fresh clone that path does not exist.
  `highconc.pls` reads them through `pls.results_root()`, which resolves in
  this order: `$HIGHCONC_PLS_RESULTS` → `data/pls_results/` → a pls-explorer
  source checkout. Set `HIGHCONC_PLS_RESULTS` to render the figures against a
  different run; otherwise leave it unset and you get the manuscript's.

## Editable vector output

The PDFs and SVGs are saved transparent with `svg.fonttype="none"` and
`pdf.fonttype=42`, so either drops straight into Illustrator or Inkscape with
the labels still real text objects you can restyle and reposition. That contract
comes from `highconc.style.apply()`, once, for every figure.

## The `figures/` folder

It holds **no figure code** — every figure is drawn by
[`notebooks/figures.py`](notebooks/figures.py). What lives there is
`figures/nmf_dim_choice_*/data/cache/*.npz`: the NMF dimensionality sweeps,
which are **load-bearing** and are read rather than recomputed. See "Four
things not to break".
