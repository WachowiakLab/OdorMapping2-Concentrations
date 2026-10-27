---
name: restyle
description: Colour, scale, marker and typography changes across the figures — "make the cluster maps paler", "use a different colormap for Fig 4", "the non-responder dots are too loud", "bump the axis label size". Routes changes through highconc/style.py and the panel presets instead of hand-edits, and states the scale constraints that must not break.
---

# Restyling

**Style is central.** `highconc.style.apply()` is called once in §1 and pins the
export contract and the type scale for every figure. No cell sets its own
rcParams, and none should start.

So a styling request is answered in one of three places, in this order of
preference:

1. **`highconc/style.py`** — anything that should be true of every figure: the
   type scale, the palette, `sequential_cmap`, `CLUSTER_COLORS`,
   `cluster_color`, `condition_color`, `despine`, `panel_label`. `style.SEQ` is
   **magma** (the paper's ramp); a panel that must differ takes a `cmap=` argument
   from its notebook cell — Fig 8's p-value panel pins OrRd that way (`FIG8_P_CMAP`).
2. **A named preset in the panel module** — `P.clustering.TOPO_VMIN/VMAX`,
   `TOPO_VMAX_HIGH`, `NONRESPONDER_MARKERS`. Change the preset, not the call
   site, when the change should hold everywhere the preset is used.
3. **The section's `FIG*_` constant in the notebook** — when it is genuinely
   one figure's choice.

Never hand-edit a rendered PDF or SVG. Every output is reproducible from the
code; a touch-up belongs upstream.

## The presets you will reach for

**Cluster colours.** `style.CLUSTER_COLORS` is `tab10` indexed 0–6, fixed per
component **across the whole paper**. `style.cluster_color(n)` takes the 1–7
number. Changing a cluster's colour changes it in Figures 4, 5, 6, S2 and S3 at
once — that is the point, so do not work around it by overriding one call site.

**Topography scale.** `P.clustering.TOPO_VMIN, TOPO_VMAX = 0.0, 5.0`, with
`TOPO_VMAX_HIGH = 20.0` for the high-concentration case, apply when
`topography` colours members by their mean response. **Fig 6 no longer does**:
its members are one solid colour (`FIG6_MEMBER_COLOR`), membership only, so it has
no scale. If you bring a response colouring back, **all seven maps must share
one scale**, or comparing them means nothing.

**Non-responder markers.** `P.clustering.NONRESPONDER_MARKERS` has three named
presets — `"context"` (`s=16`, faint, `alpha=0.35`), `"requested"` (`s=46`,
`c="0.85"`, `alpha=0.5`) and `"hollow"` (the one Fig 6 uses: `s=16`, no fill, a
pale-grey edge). `FIG6_NONRESPONDER` takes any name or a dict of scatter kwargs. Prefer adding a named preset over passing a one-off dict, so
the next figure can reuse it.

**Sequential ramp.** `style.sequential_cmap(under=..., bad=...)` is the house
ramp with an explicit under/bad colour. On a white ground the zero/missing
distinction has to be made explicitly — that is what `FIG4_EMPTY = "gainsboro"`
is doing under `FIG4_VMIN = 1`. If you change a colormap, check what happens to
below-range and missing cells.

## Scale constraints that must not break

These are correctness, not taste. A change that violates one produces a figure
that is still pretty and no longer means what the caption says.

- **Figs 3A–D and 4A share a 0–1 correlation scale.** `correlation_matrix`
  defaults to `vmin=0.0, vmax=1.0` for exactly this reason: the panels are only
  comparable if they share a scale, and autoscaling 3B alone breaks the
  comparison with 3C and 3D.
- **`TOPO_VMAX` is low on purpose** and shared across all seven maps. See above.
- **Fig 5B normalises each concentration within itself** (every cluster column
  is divided by its own top class, separately for low and high), so both are
  0–1 and share one colourbar. That is what makes the two halves comparable;
  do not pin them to a count scale.
- **Fig 4B and 4C are peak-normalised (0–1, linear)** with the raw floor
  `FIG4B_RAW_FLOOR` (= 1 dF/F), and C takes B's settings so a glomerulus reads
  the same in both. **Fig S2 is different on purpose**: raw dF/F on a log scale,
  `SUPP_SUBJ_HEAT_VMIN/VMAX = 1.0, 140.0`. They are different quantities — do
  not "make them match".
- **Fig 7's `FIG7_SPACE` is not a display choice.** `"log"` is the space the
  model is fitted and scored in, so the residuals drawn are the ones the model
  minimised; `"raw"` and `"sqrt"` are after-the-fact transforms with no
  modelling meaning. The section's markdown has the full argument. Do not flip
  it as part of a styling pass.
- **Figs 6 and S3 share marker numbers.** S3 restates Fig 6's marker size,
  outline and context-circle values (scaled to its own page), so the two read as
  the same drawing. Change them together.
- **Page widths are fixed**: `MAIN_TARGET_W = 7.0` in for Figs 2–8 and
  `SUPP_TARGET_W = 5.5` in for the supplements (S3 is 7.0 so its maps are Fig 6's
  size). Marker sizes and line widths are rescaled with the page; if you change a
  width, rescale them with it or the maps turn into blobs.

## Vector output contract

PDFs and SVGs are saved transparent with `svg.fonttype="none"` and
`pdf.fonttype=42`, so both open in Illustrator or Inkscape with labels as real,
editable text. That comes from `style.apply()`, once, for every figure. **Do
not set `fonttype` anywhere else, and do not rasterise text to "fix" a font
issue** — `style.check_export()` exists to confirm the contract holds.

For font problems, use `style.resolved_font()` to see what actually got picked
before changing anything.

## After a restyle

Style changes touch many figures at once, so re-render and check everything
rather than eyeballing one panel:

```bash
python notebooks/figures.py
python dev/concordance.py --check
python dev/build_notebooks.py --check
```

Or use `/verify`, which runs the same gate. The style guide in `style/` is the
reference for what the house conventions are; if a change makes the guide wrong,
update the guide in the same pass.
