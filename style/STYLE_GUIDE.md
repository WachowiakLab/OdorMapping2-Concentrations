# highconcmss — figure style guide

For collaborators porting figures made in other software into the manuscript
figure set (Burton/Wachowiak/Castro et al.).

This guide is **tool-agnostic**. Nothing here requires Python, matplotlib, or
the `highconc-figures` repo. It states the target; how you hit it in MATLAB,
R, Igor, Prism, or by hand in Illustrator is up to you. Where a rule exists
for a non-obvious reason, the reason is given — several of them are scars from
Illustrator misrendering things, and skipping them costs an afternoon.

**Scope.** These conventions were derived from the analysis figures in
`figures/` (the NMF/PLS work — Fig 4 and supplements) and the manuscript figure
set. **The standard sequential map is magma** (§2.1): black at zero, pale yellow
at the ceiling, which reads on the dark ground the paper's response matrices
use. The earlier drafts of Figs 1–3, 5 and 6 used a black-background hot
colormap for ΔF; magma replaces it with a perceptually even, colourblind-safe
map that sits on the same ground. (An earlier version of this guide set OrRd on
white as the standard; the repo's `highconc.style.SEQ` and the paper are magma
now.) One documented exception: raw photometric bulb activity maps keep the
acquisition colormap — see §10.

**If you just need to make a figure, go to `for-collaborators/`** — it has a
paste-anywhere prompt for any AI assistant, an importable Python module, and a
Claude Code skill, all generated from the same source. This document is the
reasoning behind what they do.

**`style_reference.ipynb` is the worked version of this document** — every rule below,
demonstrated on synthetic data, section by section. If you'd rather read code than
prose, start there; it runs anywhere with numpy and matplotlib.

Companion files in `style/palette/`:

| File | What it's for |
|---|---|
| `highconc_palette.csv` | Every named colour → hex, RGB 0–255, RGB 0–1, and its role |
| `seq_magma_256.csv` | The sequential map (magma), 256 steps. Load directly as a MATLAB/R/ImageJ colormap |
| `div_teal_gainsboro_orrd_256.csv` | The house diverging map, 256 steps |
| `highconc_swatches.ase` | Adobe swatch file — drag into Illustrator, all 20 colours appear |
| `swatches.png` / `.pdf` | Visual reference sheet |

---

## 1. The short version

If you do only five things:

1. **Export vector (PDF or SVG) with text kept as text**, not outlines, not raster. §6.
2. **Arial**, 8–9 pt for tick labels and body text, 11–12 pt bold for panel letters. §3.
3. **Sequential data → magma. Signed data → the teal↔red diverging map. Categories → the cluster colours.** Never a rainbow, never jet. §2.
4. **Top and right spines off.** Ticks small or absent. §4.
5. **Fix the colour scale across every panel that will be compared**, and say the limits in the caption. §2.4.

---

## 2. Colour

### 2.1 Sequential — activation and magnitudes

**magma** (`plt.cm.magma`), full range. Used for ΔF/F, |values|, densities,
response matrices — anything that runs from "nothing" to "a lot" with no
meaningful zero-crossing.

Black at the floor, pale yellow at the ceiling, through purple and orange-pink.
Because zero is black, **missing data must be drawn explicitly** in `gainsboro`,
never left to the axes background (§2.5). Load `palette/seq_magma_256.csv` if
your tool doesn't ship magma (R has it via `viridisLite::magma(256)`; Fiji ships
the `mpl-magma` LUT).

One panel keeps a different ramp on purpose: Figure 8's p-value matrix is
OrRd, passed to that panel explicitly (`FIG8_P_CMAP` in the notebook). OrRd on
white, the earlier standard, is retired for everything else.

### 2.2 Diverging — signed weights and differences

The house map, `teal_gainsboro_orrd`, interpolated through five anchors:

| Position | Hex | Meaning |
|---|---|---|
| −max | `#01665e` | strong negative (mixing / loss) |
| −mid | `#5ab4ac` | |
| 0 | `gainsboro` `#dcdcdc` | **zero, and it must land on zero** |
| +mid | `#fc8d59` | |
| +max | `#b30000` | strong positive (gain) |

Used for PLS weight matrices, block structure, observed − predicted
differences. The neutral midpoint is deliberately grey rather than white, so
that zero reads as "nothing here" and stays distinguishable from the page.

**The midpoint must be zero.** Set your limits symmetrically (`vmin = −vmax`)
or use a norm that pins the centre. A diverging map centred anywhere but zero
is a bug, not a style choice.

*(Where a diverging map is needed but the quantity is non-negative — block
norms, for instance — the same palette is reused from its midpoint upward.
That's intentional, for visual family resemblance across panels.)*

### 2.3 Categorical — the concentration conditions

The paper's **primary** categorical scheme, and the one that appears in the most
panels: `'most-sensitive'`, `'low'` (1x), `'high'` (25x 'low').

| Condition | Hex | |
|---|---|---|
| `'most-sensitive'` | `#e69f00` | Okabe-Ito amber |
| `'low'` (1x) | `#0072b2` | Okabe-Ito blue — same value as the accent blue below |
| `'high'` (25x) | `#d55e00` | Okabe-Ito vermilion — same value as the accent vermilion |

Colourblind-safe, and ordered weakest → strongest, which is the order they should
appear in every legend. **A condition is the same colour in every figure in the
paper.** If a panel splits by concentration, it uses these three and no others.

### 2.4 Categorical — cluster and series identity

**Glomerulus/odorant clusters A–G** use matplotlib `tab10`, positions 0–6.
Hex values are in `highconc_palette.csv` as `cluster_A` … `cluster_G`.
Cluster ↔ colour assignment is **fixed across every figure in the paper** —
cluster D is the same red in Fig 3 as in Fig 7. If you introduce a panel with
clusters in it, use these, in this order.

**Model/condition series** use Okabe-Ito, which is colourblind-safe:

| Colour | Hex | Convention |
|---|---|---|
| Blue | `#0072b2` | primary metric (Q²) |
| Vermilion | `#d55e00` | overfit / canonical comparison |
| Green | `#009e73` | "wins" |
| Grey | `#888888` | "loses" / neutral |
| Light grey | `#bbbbbb` | shuffled / null distributions |
| Purple | `#9b5bd0` | NMF |
| Mid grey | `#777777` | PCA (reference against NMF) |

Nulls and shuffles are **always grey**, never a colour. The eye should read
"this is the boring comparison" before reading the axis label.

### 2.5 Greys, floors, and missing data

- **Zero is not missing, and they must not share a colour.** A zero response is
  a *measurement*: it belongs at the floor of the sequential ramp (black in magma). `gainsboro`
  means **not measured**. On a sparse response matrix — and most of ours are
  sparse — colouring zeros grey turns the majority of the panel into "no data"
  and makes the figure claim something false. Keep `vmin = 0` on sparse
  matrices and reserve grey for genuine gaps.
- **A floor above zero is a separate, deliberate choice.** Example glomerular
  maps clip at `vmin = 0.5, vmax = 3.5` in log(1+ΔF/F) so that responders
  saturate and shrunken predictions become visible — a full-range scale hides
  that effect. But a floor *hides data*, so say in the caption that you set one
  and what it was.
- **Identity/unity lines**: grey `#737373`, dashed, ~1.2 pt, drawn *under* the
  data (`zorder` below the points).

### 2.6 What not to use

No jet, no rainbow, no hsv. No red/green as the only distinction between two
conditions. Don't recolour a panel to "make it pop" — panels that share a
quantity share a colour scale.

---

## 3. Typography

| Element | Size | Weight |
|---|---|---|
| Tick labels | 8 pt (7 pt when crowded) | normal |
| Axis labels, in-panel annotation | 9–10 pt | normal |
| Panel titles | 9–11 pt | normal |
| Panel letters (A, B, C…) | 11–12 pt | **bold** |
| Figure-level titles | 12–13 pt | bold |

**Font: Arial**, with Helvetica and DejaVu Sans as fallbacks. Not because
Arial is beautiful — because it is installed on the Macs where these get
assembled, so Illustrator doesn't flag imported text as a missing font and
silently substitute it. If your tool defaults to something else (Helvetica
Neue, Calibri, DejaVu), change it before exporting.

Sizes above are **final, on-the-page sizes**. If a figure will be reduced to
fit a column, render it at the final size rather than shrinking it later — see
§7.4.

Panel letters sit outside the axes, upper-left, roughly at axes-fraction
`(−0.10, 1.05)`. Capital, bold, no parentheses or periods: `A`, not `(a)`.

---

## 4. Geometry and axes

- **Widths**: single-column figures ~3.5 in, full-width ~7 in at final print
  size. Working/exploratory renders in the repo are much larger (11–17 in) and
  scaled down at assembly — if you're handing over a finished panel, hand it
  over at final size.
- **DPI**: 300 minimum for any raster; 600 for final submission. Working
  previews are 200.
- **Spines**: top and right **off** on panels with data axes — scatters, lines,
  bars. Left and bottom stay. **Matrices and images keep their full frame**: there
  the frame is the extent of the data, not decoration, and dropping it leaves the
  cells floating.
- **Ticks**: label size 8. On heatmaps, tick length 0 (labels only, no marks).
- **Grids**: off, unless a grid is carrying real information.
- **Aspect**: matrices where both axes are the same kind of thing (odorant ×
  odorant, cluster × cluster) are drawn square. Matrices where they aren't
  (glomerulus × odorant) fill their cell.
- **Colorbars are shared.** When several panels show the same quantity, draw
  **one** colorbar for the group rather than one per panel, positioned just
  outside the block it serves. Per-panel colorbars imply per-panel scales.

---

## 5. By figure type

### 5.1 Scatter and parity plots

- Dense clouds: small markers (s ≈ 1–7 pt²), **no marker edge**
  (`linewidths=0`), alpha 0.2–0.6. Let overplotting show density.
- Emphasis points: large (s ≈ 40–150 pt²), **white edge** ~1 pt, drawn on top.
  The white ring is what separates a highlighted point from the cloud.
- Parity/identity line: dashed, grey `#737373`, ~1.2 pt, **underneath** the data.
- Parity panels are **square**, with equal x and y limits. A parity plot with
  unequal axes misrepresents the deviation.
- Put the statistic (r, Q², n) in the panel as text, not in the title — it
  survives cropping and reflows better in Illustrator.

### 5.2 Heatmaps and matrices

- `interpolation = nearest` — **never** bilinear/smooth. These are discrete
  measurements; smoothing invents structure.
- Cluster/block boundaries: thin lines (~0.5–0.8 pt) drawn over the cells.
  Black at 0.8 pt in the scaffolded figure; a thin blue at 50% alpha in some
  working panels. Prefer **black, 0.8 pt**.
- Row/column order is meaningful — state the ordering (seriated? by cluster?
  by chemical class?) in the caption. Within-cluster, responders are ordered
  before recruited/non-responders, contiguously.
- Missing cells: opaque `gainsboro`. **Not** transparent — see §7.2.
- Small matrices (say, 7×7): draw as **vector rectangles, one per cell**, not
  as an image. §7.3.

### 5.3 Spatial maps (glomerular topography)

- Glomeruli as filled circles on a fixed spatial frame, magma, generous marker
  size (they should read as a map, not a scatter).
- **Fixed `vmin`/`vmax` across every map in a comparison** — low vs high
  concentration, observed vs predicted, subject vs subject. If two maps in the
  same figure have different scales, the figure is making a false claim.
- Sub-floor glomeruli in `gainsboro` so they read as context.
- One shared colorbar for the whole array of maps.

### 5.4 Line plots and sweeps

- Real series in colour (§2.3); nulls/shuffles in grey.
- Line weight ~1.5–1.8 pt; markers only when the x-values are few and discrete.
- Where a sweep has a chosen operating point (e.g. the selected NMF rank),
  mark it explicitly — a vertical rule or an annotated point — rather than
  leaving the reader to infer it.
- Error/spread bands: filled, same hue as the line, low alpha, no edge.

### 5.5 Bar charts

Used sparingly in this figure set, and that's deliberate — where a bar chart
would show a handful of values, prefer a dot/strip plot with the individual
observations visible. If a bar is right:

- Zero baseline, always.
- Single hue from the palette; don't colour bars just to differentiate them
  when the axis already labels them.
- Show individual points over the bar when n is small.
- No 3D bars, no gradients, no drop shadows.

### 5.6 Cumulative distributions

Preferred over histograms for comparing conditions: no bin-width choice to
defend, and the whole distribution is visible at once.

- One line per condition, in the condition colours (§2.3), drawn in canonical
  weakest → strongest order so the legend order matches the physical order.
- Step interpolation, ~1.6 pt.
- y from 0 to 1, labelled `cumul. prob.`
- Where a specific quantile matters (the fraction at S_L = 1, say), mark it
  explicitly with a rule rather than leaving the reader to read it off.

### 5.7 Box plots

- One box per condition, filled with the condition colour, black median and
  whiskers at ~0.8–1.0 pt.
- **Log axis where the quantity spans orders of magnitude** — delivered
  concentration does, and a linear axis there is unreadable. The "zero baseline
  always" rule of §5.5 is about *bars*; it does not apply to a box plot on a log
  scale, where zero does not exist.
- Show outliers as small diamonds rather than hiding them.
- Label the axis with the unit, and if two units are in play (ppb and ppm),
  say which is which rather than relying on colour alone.

### 5.8 Rank-aligned response heatmaps

Each row sorted by response magnitude and centred on its own peak, so tuning
*breadth* is comparable across rows regardless of which odorant drove each one.

- Sort descending, place the peak at the centre column, mirror out to ±N.
- Label the axis `+N | 1 | +N` — the axis is a rank, not an odorant identity.
- Rows ordered by peak magnitude or by breadth; say which in the caption.
- Pad short rows with NaN → grey, so a row with fewer odorants is visibly
  different from a row of zeros.

### 5.9 Category strips on matrix axes

Chemical-family (or cluster) membership along a matrix edge, as a coloured
segment per contiguous run, with the family name in the same colour.

- Categories must be **grouped** in the plotting order. The strip makes that
  obvious, which is the point of drawing it.
- Use one mapping across every panel that shows the same categories — pass the
  returned mapping to sibling panels rather than regenerating it.
- Rotate labels 90° when the categories are narrow; keep them horizontal when
  they are not.

### 5.10 Schematics

- Same fonts and palette as the data panels, so a schematic doesn't look
  imported from a different paper.
- Large glyph letters (X, W, Y) at ~20 pt; captions at ~11 pt.
- Build them so the pieces stay separable in Illustrator — grouped, labelled,
  not flattened.
- 3D/WebGL scenes rasterise on export. Render those at high resolution as a
  static image and treat them as raster elements deliberately, rather than
  hoping for vector.

---

## 6. Export: the Illustrator contract

Every panel gets exported three ways:

| Format | Background | Purpose |
|---|---|---|
| PNG | white | preview, slides, emailing around |
| PDF | transparent | assembly in Illustrator |
| SVG | transparent | assembly, web |

Non-negotiable properties of the vector exports:

1. **Text stays live text.** Not outlines, not paths, not pixels. You must be
   able to click a label in Illustrator and retype it.
2. **Fonts embedded as TrueType (Type 42), not Type 3.** Type 3 text arrives
   in Illustrator as unselectable garbage.
3. **Transparent background**, so panels layer without white boxes.
4. **Cropped tight** to the content — no acres of whitespace to trim by hand.

How to get there in common tools:

- **matplotlib**: `svg.fonttype="none"`, `pdf.fonttype=42`, `ps.fonttype=42`;
  save with `bbox_inches="tight"`, `transparent=True` for pdf/svg.
- **MATLAB**: `exportgraphics(fig, 'panel.pdf', 'ContentType', 'vector',
  'BackgroundColor', 'none')`. Avoid `print -dpdf` with painters on figures
  containing transparency. Check the result — MATLAB sometimes outlines text.
- **R**: `svglite::svglite()` keeps text as text (base `svg()` is less
  reliable); `ggsave(..., device = cairo_pdf)` embeds fonts properly.
- **Prism**: export as PDF or EPS, *not* "high-resolution TIFF". Under
  export options, keep "Text as text".
- **Igor**: export as PDF with "Embed fonts" on.

If your tool simply cannot emit editable-text vector, say so and send the
highest-resolution raster you can (600 dpi at final size) plus the underlying
numbers — a panel can be rebuilt from data, but not from a JPEG.

---

## 7. Illustrator gotchas

These are all real failures that cost real time in this project.

### 7.1 Missing-font substitution
Text exported in a font the assembling machine lacks gets flagged (pink
highlight) and silently substituted, changing metrics and breaking layout.
Arial avoids this. Check every imported panel for the missing-font warning
before you trust it.

### 7.2 Transparent rasters become soft masks
A raster with an alpha channel exports into PDF as a soft-masked image, and
Illustrator renders those wrong. **Make raster elements opaque.** For empty or
sub-threshold cells, set an opaque grey (`gainsboro`) rather than leaving them
transparent.

### 7.3 Tiny rasters get dropped
A small image — a 7×7 heatmap, say — can be dropped entirely on import;
the panel arrives with its cells simply missing. For any small matrix, draw
each cell as an explicit vector rectangle. It looks identical in the PNG and
survives the round trip.

### 7.4 Never resize in Illustrator
Scaling a vector panel after import corrupts stroke weights and throws off
font sizes relative to the rest of the page. **Re-render at the target size**
instead, scaling figure dimensions, font sizes, line widths, and marker areas
together. (Marker *area* scales as the square of the linear factor.)

### 7.5 Use vector shading, not raster overlays
To wash out or highlight a region, draw a semi-transparent vector rectangle
with a visible edge — it stays selectable and editable, and it doesn't trigger
§7.2. A raster overlay does the opposite of both.

---

## 8. What gets hand-tuned, and what doesn't

The repo's stated principle is that outputs are reproducible from
(config + data), and that touch-ups belong upstream in the generating code
rather than on the exported file.

In practice a small, bounded set of things are done by hand in Illustrator:

- **Multi-panel assembly** — arranging panels into Fig 3, Fig 5, etc.
- **Panel letters** placed at assembly time.
- **A few labels deliberately omitted from the scripts** because hand
  placement is better than any automatic rule (cluster labels on the
  concentration-flow figure are the example).
- **Callouts, arrows, and brackets** that annotate rather than encode data.

Everything else — colours, scales, tick values, marker sizes, what's plotted —
belongs in the code that made the panel. **The rule of thumb: if a change
would alter what a reader concludes from the data, it must happen upstream and
be reproducible. If it only affects arrangement and annotation, hand-work is
fine.**

Practical consequence for collaborators: **send the panel and the numbers
behind it.** A hand-edited `.ai` file with no path back to data can't be
regenerated when a reviewer asks for a change, and reviewers always ask.

---

## 9. Naming and provenance

Each figure lives in its own folder with a short, descriptive name (not
`fig3`, which changes as the paper is reordered — `pls_performance`,
`nmf_per_subject`, `concentration_flow`). Alongside the images:

- `caption.md` — the caption text, version-controlled
- `source.md` — where the data came from: paths, derivation, any SHAs
- the panel-level numbers as CSV where practical, with a column dictionary

When you hand a panel over, include the equivalent: what it shows, what
produced it, and the numbers. That's the part that makes the figure survive
revision.

---

## 10. Known inconsistencies in the current repo

Stated plainly so nobody ports a bug:

- **Font is not uniform.** Nine of twelve generator scripts set Arial
  explicitly; `nmf_clusters_main`, `pls_recruitment`, and `pls_result` fall
  through to DejaVu Sans (`pls_result` says it does so "to match the other
  figures," which was true of an earlier state of the repo and isn't now).
  **Arial is the standard.** Those three should be brought into line.
- **Working DPI varies** (120–300 across scripts). 300 is the floor for
  anything final.
- **Cluster-boundary colour varies** (black in the scaffolded figure, blue at
  50% alpha in some working panels). Black, 0.8 pt is the standard.
- **The older Fig 1–3/5/6 drafts used black-background hot** (black → red →
  yellow) for ΔF heatmaps, ranked response maps and correlation matrices. They
  are on magma now, which sits on the same dark ground. Because zero is black,
  the zero/missing distinction of §2.5 has to be made explicitly — missing
  cells in opaque `gainsboro` — where the old black ground hid both.
- **Raw bulb activity maps are a deliberate exception.** The panoramic
  dorsal-surface activity overlays keep the acquisition colormap (a
  jet-like blue→green→red map) rather than magma. They are photometric images
  reproduced as acquired, not re-plotted data, and the field reads them that
  way. This is the *only* sanctioned use of a rainbow map in the paper; do not
  generalise it to any panel you are plotting from numbers.
- **The house diverging map is redefined in four separate files** rather than
  imported from one place. It's identical in all four — but if you're
  transcribing it, transcribe from `palette/` here, which is the reference.

---

*Questions, or a case this doesn't cover: ask before you build the panel
rather than after.*
