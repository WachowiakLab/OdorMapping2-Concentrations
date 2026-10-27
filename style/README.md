# style/

Figure conventions for the `highconcmss` manuscript, written for collaborators
who are plotting in other software.

| File | |
|---|---|
| **`for-collaborators/`** | **Give collaborators this folder.** Everything they need to produce house-style figures from a bare CSV — a paste-anywhere AI prompt, an importable Python module, and a Claude Code skill. Start at its `README.md`. |
| **`STYLE_GUIDE.md`** | The spec. Colour, type, geometry, rules by figure type, export requirements, Illustrator gotchas. Start here. |
| **`CHECKLIST.md`** | One-page tick list to run before sending a panel over. |
| **`style_reference.ipynb`** | The guide as runnable code: every convention demonstrated on synthetic data, one section per rule. Runs on numpy + matplotlib alone — no lab packages, no real data. |
| `examples_output/` | Figures the notebook writes (PNG + PDF + SVG for each). Gitignored — regenerate by running the notebook. |

**Note on the notebook:** this repo runs `nbstripout` (see `.gitattributes`), so the committed copy of `style_reference.ipynb` carries **no outputs** by design — history stays clean. Run it once to see the figures:

```bash
cd style && python -m nbconvert --to notebook --execute --inplace style_reference.ipynb
```

Copies sent to collaborators outside the repo should be executed first, so they can read it without running anything.
| `palette/highconc_palette.csv` | Every named colour → hex / RGB255 / RGB01 / role, incl. the three concentration conditions |
| `palette/seq_magma_256.csv` | Sequential map (magma), 256 steps — loadable as a colormap in any tool |
| `palette/div_teal_gainsboro_orrd_256.csv` | House diverging map, 256 steps |
| `palette/highconc_swatches.ase` | Adobe swatch file — drag into Illustrator |
| `palette/swatches.png` / `.pdf` | Visual reference sheet |
| `generate_assets.py` | Regenerates everything in `palette/` |

The spec was derived from the existing generator scripts in `figures/*/`,
which remain the working reference implementation. Where the repo is
internally inconsistent, `STYLE_GUIDE.md` §10 says so and names the standard.

Regenerate the palette assets:

```bash
python style/generate_assets.py
```

Re-run the reference notebook end to end:

```bash
cd style && python -m nbconvert --to notebook --execute --inplace style_reference.ipynb
```
