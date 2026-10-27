---
name: highconc-figure
description: Turn a data file into a manuscript figure panel in the highconcmss house style (Burton/Wachowiak/Castro, concentration-dependent glomerular odour coding). Use when the user wants to plot a CSV/TSV/Excel/npy for the paper, restyle an existing figure to match the house style, or check whether a panel is ready to hand over. Triggers on "make a figure for the paper", "plot this in the house style", "style this panel", "is this figure ready to send", "/highconc-figure".
---

# highconc-figure

Produce a figure panel that matches the manuscript's house style and exports
cleanly into Illustrator.

## Before you plot

**Read `STYLE_GUIDE.md`** from this skill's directory (or from the repo's
`style/` directory if the user is working inside `highconc-figures`). It is the
authority. Everything here is workflow; the rules live there.

**Look at the data first.** Load it, print shape, columns, dtypes, head, and the
fraction missing per column. State in one or two sentences what you think each
column is. If the structure is genuinely ambiguous — which column is observed vs
predicted, whether values are already log-transformed, whether rows are
glomeruli or odorants — **ask, do not guess**. Guessing here produces a figure
that is confidently wrong, which is worse than a delay.

**Ask about comparability if it is not already clear:** will this panel sit
beside another one showing the same quantity? If so the two must share a colour
scale, and you need the shared limits before you plot, not after.

## Writing the figure

Use `highconc_style.py` — it is in this skill's directory and is the reference
implementation:

```python
import highconc_style as hs
hs.apply()

fig, ax = plt.subplots(figsize=(hs.WIDTH_1COL, 2.8))
hs.parity(ax, observed, predicted)      # or hs.heatmap / hs.topomap / hs.vector_matrix
hs.panel_label(ax, "A")
hs.save(fig, "panel_name", outdir="figures")
```

Helpers: `apply`, `despine`, `panel_label`, `identity_line`, `parity`,
`heatmap`, `topomap`, `vector_matrix`, `sequential_cmap`, `ecdf`, `box`,
`rank_align`, `class_strip`, `save`, `check_export`. Palette: `SEQ`, `DIV`,
`CLUSTER_COLORS`, `CONDITION_COLORS`, `condition_color`, `GAINSBORO`,
`IDENTITY`, `C_PRIMARY`, `C_OVERFIT`, `C_WIN`, `C_LOSE`, `C_SHUFFLED`.

If the module is not importable from where the user is working, copy it next to
their script rather than reimplementing it inline.

**Write a script, not throwaway cells.** The script is part of the deliverable:
a figure that cannot be regenerated cannot survive review.

## The rules you are most likely to get wrong

These are the ones that matter and that a general-purpose plotting instinct gets
backwards. The full set is in `STYLE_GUIDE.md`.

- **Anything split by concentration uses the condition colours.**
  `'most-sensitive'` amber, `'low'` (1x) blue, `'high'` (25x) vermilion — via
  `hs.condition_color(name)`, which accepts `"1x"`, `"25x"`, `"most_sensitive"`
  and similar spellings. A condition is the same colour in every figure in the
  paper, and legends run weakest → strongest. This is the paper's primary
  categorical scheme; do not let a plotting library assign its own colours to it.
- **Zero is not missing.** A zero response is a measurement and belongs at the
  floor of the ramp (black in magma); grey means *not measured*. Keep `vmin=0` on sparse
  response matrices — colouring zeros grey turns most of the panel into "no
  data". A floor above zero hides data and must be stated in the caption.
- **Diverging maps must be centred on zero.** `vmin=-v, vmax=+v`. An autoscaled
  diverging map invents a sign flip that is not in the data.
- **Panels that will be compared share one colour scale and one colorbar.**
  Autoscaling each panel separately makes the figure lie.
- **Grey must be opaque** — never transparent (a transparent raster becomes a
  PDF soft mask that Illustrator renders wrong), never white, never silently
  dropped.
- **`interpolation="nearest"` on every heatmap.** These are discrete
  measurements; smoothing invents structure.
- **Small matrices (≲20×20) are drawn as vector rectangles**, not `imshow` — a
  tiny raster can be dropped entirely on Illustrator import. Use
  `hs.vector_matrix`.
- **Nulls and shuffles are grey**, never coloured.
- **`figsize` is the final printed size** (3.5 in single-column, 7.0 in
  full-width), so 8–9 pt text is really 8–9 pt. Never plan to shrink afterwards.
- **Top/right spines off on data axes; matrices keep their frame.**

## Finishing

1. Run `hs.check_export("<stem>.svg")` and report the result. If it fails, fix
   the export settings — do not hand over a panel whose text has been outlined.
2. Show the PNG to the user.
3. Give a one-line caption: what is plotted, n, and any fixed colour-scale
   limits.
4. If the user is preparing a handover, walk `CHECKLIST.md` and report anything
   outstanding — particularly whether the underlying numbers are going with the
   figure.

## Converting one of Matt's draft panels

The older Fig 1–3/5/6 drafts used a black-background hot colormap (black → red →
yellow). **The standard is now magma**, which sits on the same dark ground. When
converting one:

- Rebuild from the data, not the image.
- The black ground used to hide both zeros and gaps. With magma, zero is the
  ramp's black floor (`vmin=0`) and a missing cell must be drawn explicitly in
  `gainsboro`. Ask if the data does not distinguish them.
- Raw photometric bulb activity maps are the documented exception and keep their
  acquisition colormap. Do not convert those; do not copy that colormap onto
  anything you are plotting from numbers.

## Restyling an existing figure

When the user has a figure made elsewhere and wants it to match: ask for the
**data**, not the image. Rebuild the panel from numbers. Recolouring an exported
image cannot fix scale limits, interpolation, or text that was already
outlined — and those are usually the actual problems.
