"""Fig 6 (new design): C = two example glomeruli, D = full W. Writes
  final/highconc/pls.py               (main's + observed_high_floor)
  final/highconc/panels/similarity.py (main's + observed_floor on prediction_spectrum)
  final/notes/sections/06_fig6_code.py
"""
from pathlib import Path

ROOT = Path(r"C:\Users\u6068031\Documents\highconc-merge")
MAIN = ROOT / "highconc-figures-lab-main"
FINAL = ROOT / "final"


def rd(p):
    return p.read_bytes().decode("utf-8").replace("\r\n", "\n")


def edit(text, edits):
    for old, new in edits:
        assert text.count(old) == 1, (text.count(old), old[:70])
        text = text.replace(old, new)
    return text


# ---- analysis: pls.observed_high_floor --------------------------------------
t = edit(rd(MAIN / "highconc" / "pls.py"), [(
    '''    return float(_to_space(m.Yhat, space).min())
''',
    '''    return float(_to_space(m.Yhat, space).min())


def observed_high_floor(m: PLSModel, space: str = "log") -> float:
    """Smallest NONZERO observed high-concentration response over the WHOLE
    dataset (every odorant x every glomerulus), in whichever space Fig 6C draws.

    `_to_space` clips at zero, so the plain minimum is always 0 (a glomerulus
    that did not respond) and says nothing; the smallest value that is an
    actual response is the lowest response the experiment ever measured. Fig 6C
    draws it as a horizontal line so the predicted bars can be read against it.
    """
    if space not in SPECTRUM_SPACES:
        raise ValueError(f"unknown space {space!r}; expected one of "
                         f"{list(SPECTRUM_SPACES)}")
    x = _to_space(m.Ylog, space)
    return float(x[x > 0].min())
''')])
(FINAL / "highconc" / "pls.py").write_bytes(t.encode("utf-8"))

# ---- panel: prediction_spectrum(observed_floor=...) -------------------------
t = edit(rd(MAIN / "highconc" / "panels" / "similarity.py"), [
    ('''                        r_fontsize=None, predicted_floor=None,
                        floor_color="0.4", **kw):''',
     '''                        r_fontsize=None, predicted_floor=None,
                        observed_floor=None, floor_color="0.4", **kw):'''),
    ('''    rather than each other's own minimum. `None` (default) draws nothing.
    """''',
     '''    rather than each other's own minimum. `None` (default) draws nothing.

    `observed_floor`, if given, draws a solid horizontal line at that value on
    the same (bottom) half -- meant for `highconc.pls.observed_high_floor`, the
    smallest nonzero response the experiment ever measured, so a predicted bar
    can be read against the lowest real response. It is drawn at the mirrored
    (negative) position, like `predicted_floor`.
    """'''),
    ('''    if annotate_r:
        ax.text(0.985, 0.97, f"r = {spec.r:.2f}", transform=ax.transAxes,''',
     '''    if observed_floor is not None:
        ax.axhline(-observed_floor, color="#cdcdcd", lw=1.0, zorder=2)
    if annotate_r:
        ax.text(0.985, 0.97, f"r = {spec.r:.2f}", transform=ax.transAxes,'''),
])
(FINAL / "highconc" / "panels" / "similarity.py").write_bytes(t.encode("utf-8"))

# ---- notebook cells ---------------------------------------------------------
cell = '''# %%
FIG6_SPACE = "log"              # "log" (the model's own space) | "sqrt" | "raw" (dF/F)
FIG6_GLOMERULI = (211, 97)      # C's two example glomeruli (dataset glomerulus ids)
FIG6_LIMIT = None               # None = autoscale each panel
FIG6_SUPPORT_THRESHOLD = 0.4    # value above which a bar counts as a response
FIG6_CANONICAL = "baseline"
# D shows EVERY weight of W: there is no top-K cut. Only the colour scale is
# clipped (each sign at its own 99th percentile, see below), never the data.

# Same odorant order as Fig 3's E-H, so C reads against them.
fig6_specs = pls.eg_predicted_spectra(m, d, FIG6_GLOMERULI, space=FIG6_SPACE,
                                      order=CLASS_ORDER_FIG3,
                                      within=WITHIN_CLASS_ORDER)
fig6_w = pls.load_weights(FIG6_CANONICAL)
# The smallest nonzero observed high response in the whole dataset, in
# FIG6_SPACE: the line drawn across C's predicted half.
fig6_obs_floor = pls.observed_high_floor(m, space=FIG6_SPACE)

for s in fig6_specs:
    n_obs, n_pred, n_both = s.support(FIG6_SUPPORT_THRESHOLD)
    print(f"6C: glomerulus {s.glom} (cluster {s.cluster}) - "
          f"r = {s.r:.3f} in {s.space} space; above "
          f"{FIG6_SUPPORT_THRESHOLD}: {n_obs} measured, {n_pred} predicted, "
          f"{n_both} both")
print(f"6C: smallest nonzero observed high response ({FIG6_SPACE}): "
      f"{fig6_obs_floor:.4f}")
print(f"6D: W is {fig6_w.real.shape[0]} x {fig6_w.real.shape[1]}, "
      f"all weights shown")

# 6C's own numbers, one row per odorant per glomerulus -- built from the same
# `fig6_specs` the panels draw, so the saved numbers are exactly what's shown.
fig6_c_table = pd.concat(
    {s.glom: pd.DataFrame({
        "odorant": s.names, "class": s.classes,
        "observed": s.observed, "predicted": s.predicted,
     }).assign(cluster=s.cluster, space=s.space, r=s.r)
     for s in fig6_specs},
    names=["glom", "row"])
if SAVE:
    fig6_c_csv = OUT / "fig6_pls_C_prediction_spectra.csv"
    fig6_c_table.to_csv(fig6_c_csv)
    print(f"saved: {fig6_c_csv.relative_to(ROOT)}")

# %%
# One row, in inches: C (two example glomeruli, stacked) | D (W, square, + colourbar).
# The width is MAIN_TARGET_W exactly; the right margin is whatever is left, so
# changing any other width cannot silently change the page width.
FIG6_MARGIN_L, FIG6_MARGIN_TOP, FIG6_MARGIN_BOT = 0.65, 0.55, 0.55
FIG6_C_W, FIG6_C_H, FIG6_C_GAP = 2.25, 0.90, 0.60   # each spectrum, and the gap between
FIG6_COL_GAP = 0.75                                  # between C and D
FIG6_D_SIDE = 2 * FIG6_C_H + FIG6_C_GAP              # D is square, as tall as C's stack
FIG6_CBAR_PAD, FIG6_CBAR_W = 0.40, 0.14
FIG6_MARGIN_R = (MAIN_TARGET_W - FIG6_MARGIN_L - FIG6_C_W - FIG6_COL_GAP
                 - FIG6_D_SIDE - FIG6_CBAR_PAD - FIG6_CBAR_W)
assert FIG6_MARGIN_R > 0.30, f"Fig 6 row does not fit: right margin {FIG6_MARGIN_R:.2f} in"
FIG6_W = MAIN_TARGET_W
FIG6_H = FIG6_MARGIN_TOP + FIG6_D_SIDE + FIG6_MARGIN_BOT


def x6(inches):
    return inches / FIG6_W


def y6(inches):
    return inches / FIG6_H


fig = plt.figure(figsize=(FIG6_W, FIG6_H), dpi=FIG_DPI)
top_in = FIG6_H - FIG6_MARGIN_TOP          # top edge of the row, from the bottom

# --- C: two example glomeruli, experimental high up / predicted high down ----
for j, s in enumerate(fig6_specs):
    y0 = top_in - FIG6_C_H - j * (FIG6_C_H + FIG6_C_GAP)
    ax = fig.add_axes([x6(FIG6_MARGIN_L), y6(y0), x6(FIG6_C_W), y6(FIG6_C_H)])
    P.similarity.prediction_spectrum(
        ax, s, limit=FIG6_LIMIT, strip_labels=True, annotate_r=False,
        observed_floor=fig6_obs_floor)
    # r goes in the title, as in Fig 3's E-H, not inside the axes
    ax.set_title(f"glomerulus {s.glom}  (cluster {s.cluster})    r = {s.r:.2f}",
                 fontsize=HEATMAP_TITLE_FONTSIZE)
    if j == 0:
        hs.panel_label(ax, "C", x=-0.18, y=1.18)

# --- D: the full W, every weight ---------------------------------------------
d_left = FIG6_MARGIN_L + FIG6_C_W + FIG6_COL_GAP
W6 = fig6_w.real
W6_reord, _, bnd6, centers6 = pls.reorder_by_cluster(W6, fig6_w.cluster_number)
# Each sign scaled to its own 99th percentile -- negatives are smaller, and a
# symmetric range washes them out. TwoSlopeNorm keeps 0 on gainsboro. This only
# sets the colour range; every weight is still drawn.
pos6, neg6 = W6_reord[W6_reord > 0], -W6_reord[W6_reord < 0]
norm6 = TwoSlopeNorm(vmin=-(np.percentile(neg6, 99) if neg6.size else 1.0),
                     vcenter=0.0,
                     vmax=(np.percentile(pos6, 99) if pos6.size else 1.0))
ax_d = fig.add_axes([x6(d_left), y6(FIG6_MARGIN_BOT),
                     x6(FIG6_D_SIDE), y6(FIG6_D_SIDE)])
P.pls.weight_blocks(ax_d, W6_reord, bnd6, centers6, norm=norm6,
                    title="W (low \u2192 high concentration mapping)")
hs.panel_label(ax_d, "D", x=-0.11, y=1.05)
cax_d = fig.add_axes([x6(d_left + FIG6_D_SIDE + FIG6_CBAR_PAD),
                      y6(FIG6_MARGIN_BOT + FIG6_D_SIDE * 0.22),
                      x6(FIG6_CBAR_W), y6(FIG6_D_SIDE * 0.56)])
cb_d = fig.colorbar(plt.cm.ScalarMappable(cmap=hs.DIV, norm=norm6), cax=cax_d,
                    extend="both")
cb_d.set_label("W weight", fontsize=hs.FS_LABEL)
cb_d.ax.tick_params(labelsize=hs.FS_TICK)

save(fig, "fig6_pls")
plt.show()
# export_figure("fig6_pls")   # -> notebooks/outputs/fig6_pls.{png,pdf,svg}
'''
out = FINAL / "notes" / "sections" / "06_fig6_code.py"
out.write_text(cell, encoding="utf-8", newline="\n")
print("wrote pls.py, panels/similarity.py,", out, len(cell.splitlines()), "lines")
