# Panel handover checklist

One page. Run through it before sending a figure panel for the manuscript.
Full reasoning in `STYLE_GUIDE.md`.

## Colour
- [ ] Sequential quantity → **magma**; signed quantity → **teal↔red diverging**; categories → **cluster colours** (`palette/highconc_palette.csv`)
- [ ] Diverging map is centred on **zero** (symmetric limits)
- [ ] Colour scale is **fixed across every panel that will be compared**, and the limits are in the caption
- [ ] Split by concentration uses the **condition colours** — 'most-sensitive' amber, 'low' blue, 'high' vermilion, legend weakest → strongest
- [ ] **Zeros are the ramp's floor (black), not grey** (`vmin = 0` on sparse matrices) — grey means *not measured*
- [ ] Any floor above zero is **stated in the caption**
- [ ] Missing values are **opaque `gainsboro`**, not transparent, not white
- [ ] Nulls and shuffles are **grey**, not coloured
- [ ] No jet, no rainbow, no red-vs-green as the only cue

## Type
- [ ] **Arial** throughout (check: no other font sneaked in via a default)
- [ ] Tick labels 8 pt · axis labels 9–10 pt · panel letters 11–12 pt bold
- [ ] Sizes are correct **at final print size** — not shrunk after the fact
- [ ] Panel letter is bold capital, upper-left, no parentheses

## Axes
- [ ] Top and right spines **off**
- [ ] Grid off unless it carries information
- [ ] Parity/comparison panels are **square** with equal limits
- [ ] Heatmaps use **nearest-neighbour**, never smoothing
- [ ] One **shared** colorbar per group of same-quantity panels

## Export
- [ ] Vector **PDF or SVG** — plus a white-background PNG for preview
- [ ] Text is still **live, editable text** (open it and click a label to check)
- [ ] Fonts embedded **TrueType/Type 42**, not Type 3
- [ ] Vector background **transparent**
- [ ] Cropped tight to content
- [ ] Any raster element is **opaque** and ≥300 dpi at final size
- [ ] Small matrices drawn as **vector cells**, not a tiny image

## Provenance
- [ ] Caption text included
- [ ] Source of the data stated — paths, derivation, version
- [ ] **The numbers behind the panel included** (CSV), so it can be rebuilt for a reviewer
- [ ] Nothing was hand-edited that changes what a reader would conclude

## Send
- [ ] Panel at final size, named descriptively (`pls_recruitment`, not `fig3v2_final`)
