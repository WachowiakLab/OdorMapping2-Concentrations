---
name: verify
description: Re-render every figure and check the repo is self-consistent — execute the notebook, then concordance.py --check and build_notebooks.py --check, optionally rebuilding the shareable site. Use after any figure edit, before sharing anything, or when asked "is this current?", "re-run everything", "check the figures still match".
---

# Verifying the figure set

Three questions, three checks. Run them in this order — the first produces what
the second reads.

First make sure Python is importing **this** repo, not an older clone (an editable
install points at one folder):

```bash
python -c "import highconc; print(highconc.__file__)"   # must be inside this repo
```

```bash
# 1. Re-render everything (minutes). Writes notebooks/outputs/*.{png,pdf,svg}
python notebooks/figures.py

# 2. Does the notebook draw exactly the agreed set?
python dev/concordance.py --check

# 3. Are the two reader editions current with the source?
python dev/build_notebooks.py --check
```

All three must exit 0. Report the actual output, not a summary of what you
expected.

## What each one catches

**`python notebooks/figures.py`** runs the whole notebook as a plain script —
no kernel, no kernelspec problem. Prefer it for verification. It takes a few
minutes; the PLS fit in §1 is the bulk, and nothing is refit afterwards. The
folded-in reclustering behind Fig 4D is read from `data/nonresponders_folded/`
and is **not** rebuilt by a run; if it is missing the run stops with an error.
On Windows, a figure's PDF that is open in a viewer is locked and the render
fails at that file (the PNG is written first).

**`concordance.py --check`** compares the stems in `notebooks/outputs/` against
`dev/CONCORDANCE.yaml`. Exit 0 means the rendered set is exactly the manifest:
12 figures, 12 stems. It reads the directory and does **not** render, so a
stale `outputs/` will pass a check it should fail — always re-render first. It
also fails the other way: renders left over from an older numbering
(`fig5_spatial_by_cluster`, `fig6_pls`, `fig7_motif_transfer`, `supp_per_subject`)
are reported as regressions, so clear `notebooks/outputs/` before the run.
What its failure modes mean:

- **MISSING** — in the manifest, not rendered. Either the run failed partway or
  a section was removed.
- **UNACCOUNTED** — rendered and in no list. Either add it to the manifest with
  a reason, or stop drawing it.
- **REGRESSION** — a retired stem is being drawn again. Usually an old file
  surviving in `notebooks/outputs/`; delete the directory and re-render before
  concluding anything else.

**`build_notebooks.py --check`** rebuilds both editions in memory and compares.
Exit 1 means `notebooks/figures.py` was edited and never rebuilt, so the
`.ipynb` files someone else opens are stale. Fix with `python
dev/build_notebooks.py` — no arguments.

## Optionally: the shareable site

```bash
python dev/build_site.py     # -> notebooks/site/index.html
```

One self-contained file: a contents page, one page per figure with its caption
and a collapsible block carrying the full provenance, plus a downloadable PDF
and SVG of every figure embedded in the page. It reads `notebooks/outputs/` and
stops with the name of any missing figure rather than producing a page with a
hole in it. It recomputes nothing, so **run the notebook first or the site is a
snapshot of an older render**. The output is gitignored — rebuild rather than
commit it, and treat the file as unpublished: it embeds the underlying values.

## Executing a notebook instead of the script

Only if you specifically need the `.ipynb` with outputs embedded. Two traps,
both of which look like a broken repo:

```bash
python -m ipykernel install --user --name highconc --display-name highconc
python -m nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=highconc \
    notebooks/figures_LLM.ipynb
```

- **Do not name the kernel `python3`.** The notebooks carry that kernelspec
  name and a user-level `python3` kernel outranks an environment one, so yours
  loses silently and you get `ModuleNotFoundError: No module named 'highconc'`
  from a correctly installed repo.
- **Use `python -m nbconvert`, not `python -m jupyter nbconvert`.** The
  `jupyter` launcher runs whichever `jupyter-nbconvert` is first on PATH, which
  may be a different Python that cannot see your kernel — it fails with
  `NoSuchKernel`.

## If a check fails

Stop and report it. Do not re-run hoping for a different answer, and do not
"fix" a raised exception from `dimsweep.load_cached` or `highconc.pls` by
making it recompute — those raise on purpose, because a silent recompute
produces numbers that differ from the manuscript's.
