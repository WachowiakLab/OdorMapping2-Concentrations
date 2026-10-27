import sys
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from cap_lib import load, save, md


def put(cells, tag, prefix, text):
    hits = [i for i, c in enumerate(cells) if c.startswith(f'# %% [markdown] tags=["{tag}"]')
            and c.split("\n")[1].startswith(prefix)]
    assert len(hits) == 1, (tag, prefix, hits)
    cells[hits[0]] = md(tag, text)


S6_HUMAN = r"""
## 6. Figure 6 — the seven groups on the dorsal bulb

Panel A of the paper's Figure 6. Each map is the same view of the dorsal bulb, with all four mice overlaid (the
field of view is consistent between them). One group is highlighted per map: its glomeruli are the solid red
points, and every other glomerulus is a small hollow grey circle behind them. The headings are the group numbers,
in each group's own colour. Groups differ in how spatially compact they are.

The red points show **membership only** — which glomeruli belong to the group — and carry no response value, so
there is no colour scale and no concentration to choose. The groups have 45, 41, 40, 75, 30, 78 and 67 glomeruli
(printed under the figure).

**Only panel A is built here.** Panel B in the paper — the bulb images, maximum-projected by odorant class — comes
from the imaging data and is drawn elsewhere.
"""

S6_LLM = r"""
## 6. Figure 6 — spatial mapping of glomeruli by cluster

**What it shows.** Paper Fig 6A: where each of the seven clusters sits on the dorsal bulb. Every map is the same
field of view — all four subjects superimposed in raw pixel coordinates, which is meaningful because the field is
consistent between them. One cluster's glomeruli are drawn solid red and every other glomerulus is drawn behind them
as a hollow grey circle. Some clusters are spatially compact and some are not, and that is the panel's content. The
rms spread of each cluster about its centroid (382–512 px) is printed under the figure with the cluster sizes
(45, 41, 40, 75, 30, 78, 67).

**Membership, not response.** The members are one solid colour (`FIG6_MEMBER_COLOR = "red"`) with no ΔF/F encoded:
`topography(..., color=...)` ignores `conc`, `vmin`, `vmax` and `cmap`, so there is no colourbar and no concentration
parameter in this section (it used to colour members by mean low-concentration response on a capped ramp). Headings
are the cluster numbers 1–7, coloured with `hs.cluster_color`.

**Marker treatment.** `FIG6_NONRESPONDER = "hollow"` — small, hollow, pale-grey context circles
(`NONRESPONDER_MARKERS["hollow"]`); members carry a thin black outline (`FIG6_CLUSTER_EDGE = "black"`). The figure is
scaled to `MAIN_TARGET_W` (7.0 in) with `FIG6_SCALE`, and marker sizes and line widths are rescaled by the same
factor (`_fig6_k`; area scales as k², widths as k), or a page this much smaller would be seven overlapping blobs.
`FIG6_CLUSTER_SIZE`, `FIG6_CLUSTER_LW`, `FIG6_CONTEXT_SIZE` and `FIG6_CONTEXT_LW` are the scaled values. These
numbers are restated by Figure S3, so that figure draws exactly like this one.

**Only panel A is built here.** Paper panel B is the bulb images (omp114), maximum-projected by odorant class
(acids, amines, and a class-II group) for the "most sensitive", "low" and "high" conditions, and is drawn elsewhere.
The paper also draws a bulb outline and anterior/lateral arrows, which this panel does not.
"""

S7_HUMAN = r"""
## 7. Figure 7 — a linear model of the low → high change

One linear model, fitted once across all odorants, predicting high-concentration responses from the
low-concentration ones. Every prediction shown is **held out**: the model is fitted on some odorants and asked
about others.

- **B** — two individual glomeruli, 211 (group 2) and 97 (group 6). Measured high-concentration response points
  up; the model's prediction for the same glomerulus points down. The light line across the lower half marks the
  smallest nonzero response ever measured, so a predicted bar under the line is smaller than anything seen in the
  data. The correlation between measured and predicted is in each title (0.78 and 0.64).
- **D** — the weight matrix W itself: what a glomerulus at low concentration contributes to each glomerulus at
  high, ordered by group. **Every weight is drawn.** The colour scale clips each sign at its own 99th percentile,
  so a few extreme weights do not wash out the rest, but nothing is left out of the picture.

**A caveat about B.** The bars are in **log1p ΔF/F** — the space the model was actually fitted and scored on, not
the raw units Figure 3F uses. In raw ΔF/F, one acid's held-out prediction for glomerulus 211 reaches **193.6**
against a measured **66.2**, which flattens every other bar once the axis autoscales to fit it; `FIG7_SPACE = "raw"`
restores those units, at that cost. The correlations depend on the space too — the same predictions score higher on
the log scale the model was scored in — so the scale is stated on the axis.

Predictions near zero still draw a visible bar: above a threshold of 0.4 the model has nonzero predictions for 25
and 51 odorants (glomeruli 211 and 97) against 11 and 29 measured. A cross-validated linear model essentially never
predicts exactly zero, so this is the model's behaviour, not a drawing artifact.

The paper's Figure 7 also has the X → W → Y schematic (A) and three cartoons of the weight matrix (C), which are
drawn elsewhere.
"""

S7_LLM = r"""
## 7. Figure 7 — PLS-R

**What it shows.** The linear-model half of the manuscript: what the model predicts for two individual glomeruli (B),
and the weight matrix it learned (D). The paper's panel A (the `Y = X · W` schematic) and panel C (three cartoons of
independent gain, unstructured retuning and structured retuning) are drawn separately and are not built here, so the
lettering runs **B, D**. The section used to hold three example glomeruli, the observed/predicted recruitment matrix,
and the cartoons; none of those are in this figure any more (the observed recruitment matrix is Figure 3E).

- **B — glomeruli 211 and 97**, side by side: experimental high **above** the line and held-out predicted high
  below, the same picture `panels.similarity.mirror_spectrum` draws for Figure 3F with the low series replaced by the
  prediction. 211 and 97 are two of Figure 3F's examples. The glomeruli come from `FIG7_GLOMERULI`, passed to
  `pls.eg_predicted_spectra`, in the odorant order of §3 (`CLASS_ORDER_FIG3`, `WITHIN_CLASS_ORDER`). The panel's own
  words are "measured" and "predicted". The correlation sits in each title (r = 0.78 and 0.64 on the log scale).
  A light line across the predicted half is `pls.observed_high_floor(m, space)`: the smallest **nonzero** observed
  high response over the whole dataset in the drawn space (1.5386 in log1p ΔF/F). The plain minimum is always 0
  because the display clips at zero, so the smallest *nonzero* value is the lowest response the experiment measured;
  a predicted bar beneath the line is smaller than any real response. It is drawn by
  `prediction_spectrum(observed_floor=...)`.
- **D — the weight matrix W**, reordered by cluster and labelled 1–7, with the block lines. **No top-K cut**: every
  one of the 376 × 376 weights is drawn (`pls.top_weight_gloms` is no longer called). Only the *colour range* is
  clipped — each sign at its own 99th percentile through a `TwoSlopeNorm` centred on 0 — since negatives are smaller
  and a symmetric range washes them out. Diagonal weights scale a glomerulus by its own response; off-diagonal
  weights are the shared co-tuning.

**Layout.** One row in inches, `MAIN_TARGET_W` (7.0 in) wide: B's two panels (`FIG7_B_W`, `FIG7_B_H`, `FIG7_B_GAP`),
then D (`FIG7_D_SIDE`, square) and its colourbar. The right margin is derived, and an assertion fails if the row ever
does not fit, so no width change can silently change the page width.

**The display space is `"log"`.** `log1p` is the space the model was fitted and scored in, so the bars and the printed
correlation are the quantities it was optimised against. Its cost is comparability with Figure 3F, which is in raw
ΔF/F. `"raw"` undoes `log1p` with `expm1` and un-shrinks the prediction: for glomerulus 211 one acid's prediction
reaches **193.6** against a measured **66.2**, which flattens every other bar in the panel — a property of the
display space, not of the model. `"sqrt"` is a middle ground. The correlations differ by space for the same
predictions (log is highest because it is the space the model was scored in), and the y-axis label follows
`FIG7_SPACE`, so it cannot claim ΔF/F while drawing something else.

Predictions near zero still draw a visible bar: above `FIG7_SUPPORT_THRESHOLD = 0.4` the model has nonzero
predictions for **25 / 51** odorants (glomeruli 211 / 97) against **11 / 29** measured. A cross-validated linear rule
spreads a little predicted response everywhere the population is active; `PredictionSpectrum.support` prints both
counts below the figure.

**Numbers.** `fig7_pls_B_prediction_spectra.csv` holds each odorant's observed and predicted value for the two
glomeruli, built from the same spectra the panels draw.

**Parameters.** `FIG7_SPACE`, `FIG7_GLOMERULI`, `FIG7_LIMIT` (`None` autoscales each spectrum; a shared limit puts them
on one axis), `FIG7_SUPPORT_THRESHOLD`, `FIG7_CANONICAL` (which saved PLS run supplies W: `"baseline"`). This section
reads no bare constant from another section. W comes from `pls.load_weights`, read from the saved canonical run;
nothing is refitted here.
"""

S8_PLACEHOLDER = None
cells = load()
put(cells, "human", "# ## 6. Figure 6", "# " + S6_HUMAN.strip("\n"))
put(cells, "llm", "# ## 6. Figure 6", "# " + S6_LLM.strip("\n"))
put(cells, "human", "# ## 7. Figure 7", "# " + S7_HUMAN.strip("\n"))
put(cells, "llm", "# ## 7. Figure 7", "# " + S7_LLM.strip("\n"))
save(cells)
print("S6 and S7 captions replaced")
