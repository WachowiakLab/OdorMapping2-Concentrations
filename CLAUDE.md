# highconc-figures

Figure code for the high-concentration glomerular odor-mapping manuscript
(Burton / Wachowiak / Castro et al.).

**Every figure in the paper is drawn by one notebook.** There is no per-figure
script. If you are looking for the code behind a figure, it is a numbered
section of `notebooks/figures.py`, and the section number is the figure number.
The numbering and the panel letters are the manuscript's (see the PDF in the
workspace's `manuscript/`): §2–§8 draw Figures 2–8 and §9–§12 draw the
supplements S1–S4. A figure draws only the panels of its paper figure that come
from this data, so a letter can start at B or D.

Read this file before editing anything here. The layout is deliberate and two
of the traps below will waste an afternoon if you meet them by surprise.

## The shape of the repo

```
highconc/            the package
├── clusters.py          analysis: NMF clustering, cluster numbering, colours
├── pls.py               analysis: the PLS-R model, loaded from disk, never refit
├── similarity.py        analysis: correlation and chemical-similarity matrices
├── recruitment.py       analysis: the recruitment matrix and its thresholds
├── transfer.py          analysis: motif transfer between glomerular groups
├── flow.py              analysis: motion through tuning space (no current figure uses it)
├── dimsweep.py          analysis: NMF dimensionality sweeps, read from cache
├── folding.py           analysis: the folded-in reclustering (Fig 4D), read from cache
├── style.py             the export contract and the type scale — the only rcParams
└── panels/              drawing: one function per panel, each takes an `ax`
notebooks/
├── figures.py           THE SOURCE. Edit this. Code once, captions twice.
├── figures_human.ipynb  generated — the friendly edition
├── figures_LLM.ipynb    generated — the exact edition, full provenance
└── outputs/             renders (png/pdf/svg), gitignored, rebuilt by a run
data/                 vendored inputs: chemical cache, PLS results, pairwise nulls,
                      Burton's response matrices, the folded-in cluster labels
figures/              the NMF sweep caches — load-bearing, see trap 2
dev/                  concordance.py, build_notebooks.py, build_site.py, fold_nonresponders.py
                      + the manifest
style/                the style guide and its assets
```

**Three layers, and the boundary between them is the point:**

- **Analysis** is in `highconc/*.py`. It returns numbers. It never draws.
- **Drawing** is in `highconc/panels/*.py`. Every panel is a function that
  takes an existing `ax` and draws into it. A panel never creates a figure,
  never calls `plt.show` or `savefig`, and never loads data.
- **Composition** — grids, colorbars, panel letters, scales — is in the
  notebook, where it is visible and editable. **That is the part you are meant
  to change.**

So: to change what a panel *draws*, edit its function in `highconc/panels/` and
every figure using it follows. To change a *layout*, edit the notebook cell that
builds it.

## Running it end to end

```bash
pip install -e ".[pls,notebook]"
python notebooks/figures.py
```

That runs the whole notebook as a plain script — same figures, no kernel — and
writes PNG, PDF and SVG of every figure into `notebooks/outputs/`. It takes a
few minutes; the PLS fit in §1 is the bulk of it and nothing is refit
afterwards.

To work interactively instead:

```bash
jupyter lab notebooks/figures_human.ipynb
```

Run §1 first, then any figure's section. Every section owns its constants and
inherits no bare name from another, so a section runs correctly on its own.

### Verifying a change

```bash
python dev/build_notebooks.py --check    # are the two editions current?
python dev/concordance.py --check        # does the notebook draw exactly the agreed set?
```

`concordance.py --check` reads `notebooks/outputs/`, so run the notebook first.
It does not run the notebook itself — rendering is minutes, and a harness that
silently triggers it is the same trap trap 2 exists to avoid. Exit 0 means the
rendered set is exactly the manifest: an extra render is a figure nobody asked
for, a missing one is unfinished work, and both are otherwise silent.

## The five traps

**1. The notebook's kernelspec is `python3`, which is probably not the
environment you installed into.** A notebook saved with that name will happily
start some other interpreter that has no `highconc` in it, and you get
`ModuleNotFoundError: No module named 'highconc'` from a repo that is correctly
installed. Two rules:

- Register a **distinctly named** kernel. Do not call it `python3` — a
  user-level `python3` kernelspec outranks one registered into an environment,
  so yours would lose silently:
  ```bash
  python -m ipykernel install --user --name highconc --display-name highconc
  ```
- When executing headlessly, say which kernel, and invoke nbconvert **as a
  module**:
  ```bash
  python -m nbconvert --to notebook --execute --inplace \
      --ExecutePreprocessor.kernel_name=highconc \
      notebooks/figures_LLM.ipynb
  ```
  `python -m jupyter nbconvert` is not the same thing: the `jupyter` launcher
  runs whichever `jupyter-nbconvert` is first on your PATH, which may be a
  different Python entirely, and it then reads *its* kernel list and fails with
  `NoSuchKernel`.

A close cousin: **an editable install that points at a different copy of the
repo.** `pip install -e` records one folder; if that is an older clone, the
notebook silently imports *its* `highconc` and you debug code you did not edit
(the traceback names the other path). Check which one you are getting:
`python -c "import highconc; print(highconc.__file__)"`. Reinstall with
`pip install -e . --no-deps` from the repo you mean, or set `PYTHONPATH` to it.

**2. The NMF dimensionality sweeps are cached and will not recompute
themselves.** `figures/nmf_dim_choice_*/data/cache/*.npz` takes minutes of NMF
to rebuild. `dimsweep.load_cached` reads those files in place and **raises** on
a miss rather than recomputing. Leave it that way and do not "fix" the
exception — a silent recompute is how a notebook quietly stops matching the
paper. The same rule holds for the PLS model: `highconc.pls.load_weights`,
`.component_selection` and `.lv_loadings` read the canonical run from
`data/pls_results/baseline/` and raise if it is missing. So does the Fig 4D
reclustering — see trap 5.

**3. Edit `notebooks/figures.py`, never `figures_human.*` or `figures_LLM.*`.**
The two editions are generated from that one source and will be overwritten.
Code cells are copied into both byte for byte — that is the whole point, so the
two editions cannot draw different figures. Only prose differs, routed by a cell
tag:

```python
# %% [markdown] tags=["human"]
# ## 5. Figure 5 — where the seven groups sit on the bulb
# %% [markdown] tags=["llm"]
# ## 5. Figure 5 — spatial mapping of glomeruli by cluster
```

Untagged markdown goes into both. After editing:

```bash
python dev/build_notebooks.py            # rebuild both editions
python dev/build_notebooks.py --check    # exit 1 if they are stale
```

**4. Clusters are numbered 1–7 everywhere.** The legacy A–G letters survive
only as a mapping in `highconc.clusters`, and **they do not map in order** — a
letter does not tell you the number. If you are reading an older figure or note
that uses letters, convert through that mapping rather than assuming C is 3.

**5. The folded-in reclustering (Fig 4D) is a cached fit, and it is seeded.**
Folding the recruited non-responders into the low-concentration data and
reclustering is an NMF fit, so it is done once by `python dev/fold_nonresponders.py`
and the labels live in `data/nonresponders_folded/glom_clusters.csv`
(`meta.json` records the run). `highconc.folding.load_cache` reads them and
**raises** if they are missing — do not make it recompute. When you rebuild it,
keep the clustering call at the package's own defaults
(`assign_cluster_indices(df, method="nmf")`, seed and iteration cap included):
passing `params={"n_components": 7}` *replaces* those defaults, drops the seed,
and gives different labels on every run.

## Repo-local skills

`.claude/skills/` ships with this repo, so these are available as soon as you
clone it — no install step:

- **`/figure`** — regenerate or adjust one figure. Carries the
  section-to-figure map, every `FIG*_` parameter and its legal values, and the
  invariants a change must not break.
- **`/verify`** — re-render everything and run both checks, so a stale render
  cannot be mistaken for a current one.
- **`/restyle`** — colour, scale and marker changes routed through
  `highconc/style.py` and the panel presets rather than hand-edits.

Separately, `style/for-collaborators/` holds an older skill with a different
job — turning an arbitrary data file into a house-style panel from outside this
repo — installed by its own shell script. For work *inside* this repo, use
`/figure`.

## What `.claude/settings.json` pre-allows

The repo ships a permission allowlist so you are not approving a prompt every
thirty seconds. **It is deliberately narrow** — exact commands, plus two
tightly bounded prefixes for `nbconvert --to notebook --execute` and
`ipykernel install --user --name highconc`. Rendering the figures and running
the checks proceed without asking; everything else still prompts.

**It grants nothing outside this repo's own build.** Installing packages,
`git`, and anything writing outside the working tree are all absent on purpose.

If you widen it, keep entries specific. `Bash(python:*)` and `Bash(python -m:*)`
read like "let it render the figures" and actually mean "run arbitrary Python
without asking" — that is the whole permission surface, granted silently at
clone time. Delete the file if you would rather approve each command yourself;
nothing depends on it.

## Conventions that are settled

Do not reopen these without a reason; the code depends on all of them holding
uniformly.

1. **One notebook** for every figure, not one per figure.
2. **Panels take an `ax`.** They never create a figure, call `plt.show`, or
   save.
3. **Style comes from `highconc.style`.** `hs.apply()` pins the export contract
   — editable SVG text, TrueType in PDF — and the type scale for every figure.
   No cell sets its own rcParams.
4. **Spec-dependent choices are parameters** with documented defaults, set in
   the notebook cell rather than buried in a panel: activation threshold,
   normalisation mode, chemical-class order, cluster numbering, colour scales.
   Change the constant at the top of a section, not the panel function.
5. **Nothing is refitted at draw time.** See trap 2.
6. **Never edit a render by hand.** Every output is reproducible from the code;
   if a figure needs touch-up, the touch-up belongs upstream.

## What this repo does not build

Figure 1, Fig 2B–E, Fig 6B, the schematics in Fig 7A and 7C, and the
procedure diagram in Fig 8A are drawn elsewhere and are deliberately not built
here. `dev/CONCORDANCE.yaml` is the manifest that fixes the list — which
figures ship, panel by panel, with each omission recorded as `skipped` — and
`dev/concordance.py --check` enforces it. If you add a figure, add it to the
manifest too, or the check will call it scope creep and fail. The check reads
`notebooks/outputs/`, so a render left over from an older numbering is reported
as a regression until you delete it.

## Dependencies worth knowing about

- **glom-explorer** carries the clustering data *inside the package*. There is
  no data path to configure and nothing to download.
- **pls-explorer** comes with the `pls` extra. §1 fits the baseline model and
  shares it, so the notebook needs it end to end; the clustering figures (2, 4, 5,
  6 and S1–S3) do not.
- **RDKit** comes with the `chem` extra and is needed only to *rebuild* the
  chemical-similarity cache behind Fig 3A. The 59×59 MACCS matrix ships in
  `data/chem/`, so drawing that figure never needs RDKit.
- Both workhorse packages are **pinned to a commit** in `pyproject.toml` so a
  fresh clone rebuilds the same figures. Bump those SHAs deliberately; never
  float them. They are private GitHub repos, installed over **HTTPS**, so you need
  a GitHub login with access (a personal access token or the GitHub CLI's
  credential helper) — not an SSH key.

Note the import names: `glom_explorer.topography` and `glom_explorer.tuning`
are **submodules** of `glom-explorer`. There is no separate `glom-topography`
or `glom-tuning` package.
