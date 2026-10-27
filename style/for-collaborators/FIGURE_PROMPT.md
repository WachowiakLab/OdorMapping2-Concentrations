# Figure prompt — highconcmss house style

**What this is.** Paste the block below into any AI assistant (ChatGPT, Claude,
Gemini, Copilot — it does not matter which), attach or paste your CSV, and say
what you want plotted. It will write Python that produces a figure in the
manuscript's house style, exported so it drops straight into Illustrator.

**You do not need this repo, a Python environment you have configured, or any
prior setup** — only numpy, pandas and matplotlib, which every scientific Python
install already has. If you *do* have the repo, `highconc/style.py` does the same
job with less pasting; see `README.md`.

**Two things to tell it up front**, because it cannot see them and will otherwise
guess wrong:

1. **What the columns mean** — units, and whether values are raw or transformed.
2. **Whether this panel will sit beside another one showing the same quantity.**
   If it will, they must share a colour scale, and the assistant needs to be told
   the shared limits.

---

Everything from here to the end is the prompt. Copy it all.

---

````text
You are helping prepare a figure panel for a neuroscience manuscript
(Burton/Wachowiak/Castro, concentration-dependent glomerular odour coding). The
figure must match an existing house style exactly. Follow these instructions over
your own defaults, including when you think a different choice would look better.

## Your workflow

1. **Look at the data first.** Print the shape, the column names, dtypes, and the
   first few rows. State what you think each column is. If the structure is
   ambiguous, ask before plotting rather than guessing.
2. **Say which figure type fits and why**, in one or two sentences, before you
   write plotting code. If the user asked for a bar chart of few values, say so
   and offer a dot/strip plot instead (see below).
3. **Write one self-contained script.** It must start with the style block below,
   verbatim. Do not restyle it, do not "modernise" it, do not substitute seaborn.
4. **Export the trio** (PNG + PDF + SVG) with the `save()` helper.
5. **Verify the export** with the check at the end, and report the result.
6. **Hand back**: the script, the figure, and a one-line caption stating what is
   plotted, the n, and any colour-scale limits you fixed.

## The style block — copy verbatim into every script

```python
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# --- Illustrator contract: text stays editable, fonts embed as TrueType ---
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

# --- type scale (points), final on-the-page sizes ---
FS_TICK, FS_LABEL, FS_TITLE, FS_PANEL = 8, 9, 10, 11
mpl.rcParams.update({
    "font.size": FS_LABEL, "axes.labelsize": FS_LABEL, "axes.titlesize": FS_TITLE,
    "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
    "legend.fontsize": FS_TICK, "legend.frameon": False, "axes.grid": False,
})

# --- palette ---
SEQ = plt.cm.magma                                    # sequential: magnitudes
DIV = LinearSegmentedColormap.from_list(              # diverging: signed values
    "teal_gainsboro_orrd",
    ["#01665e", "#5ab4ac", "gainsboro", "#fc8d59", "#b30000"])
CLUSTER_COLORS = plt.cm.tab10(np.arange(7))           # clusters A-G, fixed
GAINSBORO  = "#dcdcdc"    # NOT MEASURED (not "zero" - zero is data)
IDENTITY   = "0.45"       # unity lines (== #737373)
BOUNDARY, BOUNDARY_LW = "black", 0.8
C_PRIMARY  = "#0072b2"    # primary metric
C_OVERFIT  = "#d55e00"    # overfit / canonical comparison
C_WIN      = "#009e73"    # "wins"
C_SHUFFLED = "#bbbbbb"    # nulls and shuffles are ALWAYS grey

# --- concentration conditions: the paper's primary categorical scheme ---
# Same colour for a condition in EVERY figure. Order = weakest to strongest.
CONDITION_COLORS = {
    "most-sensitive": "#e69f00",   # Okabe-Ito amber
    "low":            "#0072b2",   # 1x   (== C_PRIMARY)
    "high":           "#d55e00",   # 25x  (== C_OVERFIT)
}
CONDITIONS = ("most-sensitive", "low", "high")

def despine(ax, keep=("left", "bottom")):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)

def panel_label(ax, letter, x=-0.10, y=1.05):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=FS_PANEL,
            fontweight="bold", va="bottom", ha="right")

def save(fig, stem, dpi=300):
    """White PNG to look at; transparent PDF + SVG with live text for Illustrator."""
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{stem}.{ext}", dpi=dpi, bbox_inches="tight",
                    facecolor="white" if ext == "png" else "none",
                    transparent=(ext != "png"))
    print(f"saved {stem}.png / .pdf / .svg")
```

## Rules

**Colour**
- Magnitudes that run from nothing to a lot → `SEQ` (magma).
- Signed quantities (differences, weights, correlations) → `DIV`, and the scale
  MUST be centred on zero: `vmin=-v, vmax=+v`. Never autoscale a diverging map.
- **Split by concentration condition → `CONDITION_COLORS`.** This is the
  paper's primary categorical scheme: 'most-sensitive' amber, 'low' (1x) blue,
  'high' (25x) vermilion. A condition is the SAME colour in every figure, and
  legends list them weakest → strongest. Column values like "1x", "25x" or
  "most_sensitive" refer to these same three conditions.
- Other unordered categories → `CLUSTER_COLORS` (clusters A–G:
  #1f77b4, #ff7f0e, #2ca02c, #d62728, …). If the data has named clusters A–G, the
  colour for each letter is fixed and must not be reassigned.
- Nulls, shuffles, controls, chance levels → grey (`C_SHUFFLED`). Never coloured.
- **Zero is not missing.** A zero response is a measurement and belongs at the
  floor of the ramp (black in magma); `GAINSBORO` grey means NOT MEASURED. On a sparse matrix
  keep `vmin=0` — colouring zeros grey turns most of the panel into "no data"
  and makes the figure claim something false.
- Genuinely missing values → opaque `GAINSBORO`. Never transparent, never white,
  never dropped silently. A floor above zero (`vmin > 0`) is a separate,
  deliberate choice that hides data — say so in the caption if you set one.
- Never jet, rainbow, hsv, or viridis-by-default. Never red-vs-green as the only
  distinction between two conditions.

**Comparability — the rule that matters most**
- Any panels that will be looked at side by side must share `vmin`/`vmax` and
  ideally one colorbar. If you autoscale each panel separately, the figure makes
  a false claim. When in doubt, fix the limits and say what they are.

**Axes**
- Top and right spines off (`despine(ax)`) on scatters, lines and bars.
- Matrices and images KEEP their full frame — there the frame is the data extent.
- No gridlines unless the grid carries information.
- Parity/observed-vs-predicted panels: square, equal x and y limits, dashed grey
  unity line at `lw=1.2, color=IDENTITY`, drawn UNDER the data (`zorder=1`).
- Put statistics (r, n, p) inside the panel as text, not in the title.
- Panel letters: bold capital, upper-left, outside the axes. `A`, not `(a)`.

**Per figure type**
- *Scatter*: dense clouds get tiny edgeless translucent markers (`s=1–6`,
  `alpha=0.2–0.6`) so overplotting reads as density; highlighted points get large
  white-edged markers on top.
- *Heatmap*: `interpolation="nearest"` ALWAYS — never smooth discrete
  measurements. Draw cluster/block boundaries as thin black lines
  (`BOUNDARY_LW`). State the row/column ordering in the caption.
- *Small matrix* (up to ~20×20): draw each cell as a `plt.Rectangle`, not
  `imshow`. A tiny raster can be dropped entirely on Illustrator import.
- *Spatial map*: fixed `vmin`/`vmax` shared across every map being compared;
  sub-threshold points in `GAINSBORO`; one shared colorbar.
- *Line/sweep*: real series coloured, nulls grey, spread as a filled band with no
  edge; mark any chosen operating point explicitly.
- *Bar*: zero baseline always; overlay the individual points when n is small; no
  3D, no gradients, no shadows. Prefer a dot/strip plot for few values.
- *Cumulative distribution* (preferred over a histogram for comparing
  conditions): one step line per condition in the condition colours, ~1.6 pt,
  y from 0 to 1 labelled "cumul. prob.".
- *Box plot*: one box per condition, filled with the condition colour, black
  median and whiskers. Use a LOG y-axis when the quantity spans orders of
  magnitude (delivered concentration does). The "zero baseline" rule is about
  bars and does not apply here.
- *Rank-aligned heatmap*: sort each row descending, centre it on its own peak,
  mirror out to ±N, and label the axis "+N | 1 | +N" — the axis is a rank, not
  an odorant identity. Pad short rows with NaN so they read as missing, not zero.
- *Category strip*: to show chemical-family membership along a matrix edge, draw
  one coloured segment per contiguous run with the family name in the same
  colour. Categories must be grouped in the plotting order. Use one mapping
  across all sibling panels.

**Sizing**
- Single-column ≈ 3.5 in wide, full-width ≈ 7.0 in. Set `figsize` to the FINAL
  printed size so the 8–9 pt text is genuinely 8–9 pt on the page.
- Never plan to shrink the figure afterwards. If it needs to be smaller,
  re-render with every dimension scaled together (marker *area* scales as the
  square of the linear factor).
- Use `layout="constrained"` on multi-panel figures.

## Export check — run this and report the output

```python
import re
svg = open("YOUR_STEM.svg").read()
n = len(re.findall(r"<text[^>]*>", svg))
print(f"<text> elements: {n}")
print("PASS - text is live and editable" if n else
      "FAIL - text was outlined or rasterised")
```

If it fails, fix the export settings rather than shipping the figure: text that
arrives in Illustrator as outlines or pixels cannot be corrected at assembly.

## Do not

- Do not invent data, fill gaps, or drop rows silently. Report what is missing.
- Do not apply a transform (log, z-score, normalisation) without saying so in
  your reply and in the caption.
- Do not tidy the figure by hand-editing the exported file.
````

---

## After you get the figure back

Send the collaborator group **three things**: the PDF or SVG, the script that
made it, and the CSV it was made from. A figure can be rebuilt from data when a
reviewer asks for a change; it cannot be rebuilt from a PDF.

Run the handover checklist in `CHECKLIST.md` if you have it.
