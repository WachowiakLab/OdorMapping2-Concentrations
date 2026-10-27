# %%
CLASS_ORDER_FIG3 = "matt"      # "table" | "matt" | "matt_esters_first"
WITHIN_CLASS_ORDER = "matt"    # "table" | "matt" (pyrazine only) | "burton" (full order)
FIG3_SPLIT_AROMATICS = True    # split "R" into 3 display subclasses on panels A-D
FIG3_PANEL_ORDER = "DABC"      # left-to-right draw order for the top row -- a
                               # panel's letter stays tied to its content
                               # (D is still the MACCS matrix) regardless of
                               # where it's drawn
SPECTRUM_SCALE = "per-panel"   # "per-panel" | "shared"
EG_GLOMERULI = similarity.EG_GLOMERULI     # Fig 6B-D must reuse these
# Fourth example glomerulus (panel H), by dataset glomerulus id like EG_GLOMERULI.
# Kept out of EG_GLOMERULI, which other code imports.
EXTRA_GLOM = 155

# Odorant CODEs (the leading token of `Odorant`), positionally aligned with
# `d`'s native odorant order -- `recruitment_matrix`'s `within` needs these,
# and §6 reuses this same array for its own recruitment call for the same
# reason it reuses CLASS_ORDER_FIG3/WITHIN_CLASS_ORDER.
ODORANT_CODES = np.asarray([s.split()[0] for s in d.odorants["Odorant"].values])

mats = [similarity.most_sensitive_correlation(d, order=CLASS_ORDER_FIG3,
                                              within=WITHIN_CLASS_ORDER),
        similarity.response_correlation(d, "low", order=CLASS_ORDER_FIG3,
                                        within=WITHIN_CLASS_ORDER),
        similarity.response_correlation(d, "high", order=CLASS_ORDER_FIG3,
                                        within=WITHIN_CLASS_ORDER),
        similarity.chemical_similarity(d, order=CLASS_ORDER_FIG3,
                                       within=WITHIN_CLASS_ORDER)]
if FIG3_SPLIT_AROMATICS:
    mats = [similarity.split_block(s, d, "R", cl.AROMATIC_SUBCLASS_BY_CODE)
           for s in mats]
spectra = similarity.eg_spectra(d, [*EG_GLOMERULI, EXTRA_GLOM],
                                order=CLASS_ORDER_FIG3, within=WITHIN_CLASS_ORDER)

shift = similarity.concentration_shift(mats[1], mats[2])
print(f"mean off-diagonal correlation: {shift['mean_low']:.4f} (low) -> "
      f"{shift['mean_high']:.4f} (high), delta {shift['delta']:+.4f}")
print(f"the two matrices agree at Spearman {shift['spearman']:.3f} "
      f"over {shift['n_pairs']} odorant pairs")
for s in spectra:
    print(f"glom {s.glom:3d} (cluster {s.cluster}): "
          f"{int(s.responded.sum()):2d} odorants responding, max {s.limit:.1f}")

# %%
# Panel I's matrix is computed here rather than inherited from §6, so the two
# sections cannot drift and either can be run without the other.
# `recruitment.fit` is memoized, so the second call costs nothing.
# `threshold_effodor` in the lab's original algorithm -- one shared fraction,
# used both to call a glomerulus "tuned" to a low odorant and, on the high
# side, to call it "recruited". `LO_FRAC_FIG3` is the separate, lower bar for
# "not really active at low" that recruitment requires.
ROW_FRAC_FIG3 = 0.20
LO_FRAC_FIG3 = 0.20

rmx3 = rc.recruitment_matrix(m, level="odorant", lo_threshold=LO_FRAC_FIG3,
                             frac=ROW_FRAC_FIG3, order=CLASS_ORDER_FIG3,
                             within=WITHIN_CLASS_ORDER, codes=ODORANT_CODES)
print(f"panel I: rows built from {rmx3.n_rows.mean():.1f} gloms on average; "
      f"{rmx3.n_empty} empty")

# %%
# --- Figure 3 ---------------------------------------------------------------
fig = plt.figure(figsize=(14.0, 11.6), dpi=200)
gs = fig.add_gridspec(3, 4, height_ratios=[1.00, 0.58, 1.55],
                      hspace=0.38, wspace=0.30,
                      left=0.07, right=0.91, top=0.965, bottom=0.045)

# Letter -> (matrix, title), so FIG3_PANEL_ORDER can reorder where each is
# drawn without disturbing which letter names which content.
titles_by_letter = dict(zip("ABCD", ["'most sensitive'", "'low' (1x)",
                                     "'high' (25x)", "MACCS cosine"]))
mats_by_letter = dict(zip("ABCD", mats))
axes_top4 = []
for j, L in enumerate(FIG3_PANEL_ORDER):
    ax = fig.add_subplot(gs[0, j])
    # The x-axis chemical-class strip, with its family names, is drawn on all
    # four panels (they share the same class order).
    im = P.correlation_matrix(ax, mats_by_letter[L], title=titles_by_letter[L],
                              boundaries=False, strip_labels=True)
    hs.panel_label(ax, L, x=-0.16)
    axes_top4.append(ax)
cax = fig.add_axes([0.925, 0.775, 0.010, 0.14])
P.similarity_colorbar(fig, im, cax)

# E-H: one column per spectrum, nested inside the outer row's full span rather
# than the top row's 4 columns.
gs_efg = GridSpecFromSubplotSpec(1, len(spectra), subplot_spec=gs[1, :],
                                 wspace=0.30)
LIM = P.shared_spectrum_limit(spectra) if SPECTRUM_SCALE == "shared" else None
for j, (s, L) in enumerate(zip(spectra, "EFGH")):
    ax = fig.add_subplot(gs_efg[0, j])
    P.response_spectrum(ax, s, limit=LIM,
                        ylabel="Resp. amplitude" if j == 0 else "")
    ax.set_title(f"glomerulus {s.glom} · cluster {s.cluster}",
                 fontsize=hs.FS_LABEL)
    hs.panel_label(ax, L, x=-0.22)

# H gets exactly A-D's rendered box size (not gs[2, :]'s full-width span +
# set_box_aspect(1), which made it much bigger) -- measured off a real A-D
# axes after a draw, since aspect="equal" shrinks their box from the raw
# gridspec cell to a square, and that final size isn't known until then.
fig.canvas.draw()
pos_top = axes_top4[0].get_position()
row2_bbox = gs[2, :].get_position(fig)
h_x0 = (row2_bbox.x0 + row2_bbox.x1) / 2 - pos_top.width / 2
h_y0 = (row2_bbox.y0 + row2_bbox.y1) / 2 - pos_top.height / 2
axH = fig.add_axes([h_x0, h_y0, pos_top.width, pos_top.height])
imH = P.recruitment.matrix(axH, rmx3, "observed",
                           vmax=P.recruitment.shared_vmax(rmx3), ylabel=True,
                           title="binary recruitment matrix (observed)")
hs.panel_label(axH, "I", x=-0.16, y=1.05)
caxH = fig.add_axes([h_x0 + pos_top.width + 0.012, h_y0, 0.010, pos_top.height])
cbH = fig.colorbar(imH, cax=caxH, extend="max")
cbH.set_label("recruitment probability", fontsize=hs.FS_LABEL)
cbH.ax.tick_params(labelsize=hs.FS_TICK)

save(fig, "fig3_correlational_structure")
plt.show()
# export_figure("fig3_correlational_structure")   # -> notebooks/outputs/fig3_correlational_structure.{png,pdf,svg}
# export_figure("fig3_correlational_structure", formats=("pdf",), outdir="~/Desktop")

# %%
SWEEP_FRACS_FIG3 = (0.20, 0.35, 0.50, 0.65, 0.80)

sweep3 = rc.matrix_sweep(m, fracs=SWEEP_FRACS_FIG3, level="odorant",
                         lo_threshold=LO_FRAC_FIG3, order=CLASS_ORDER_FIG3)
VMAX_S3 = P.recruitment.shared_vmax(*sweep3)

fig, ax = plt.subplots(2, len(sweep3), figsize=(3.0 * len(sweep3), 6.8), dpi=150)
for j, s in enumerate(sweep3):
    for i, which in enumerate(("observed", "predicted")):
        im = P.recruitment.matrix(ax[i, j], s, which, vmax=VMAX_S3,
                                  tick_labels=None, xlabel=False,
                                  ylabel=False, title="", boundaries=False)
    P.recruitment.matrix_sweep_label(ax[0, j], s)
ax[0, 0].set_ylabel("observed")
ax[1, 0].set_ylabel("predicted (held-out)")
fig.colorbar(im, ax=ax, fraction=0.014, pad=0.02, extend="max",
             label="recruitment probability")
fig.suptitle("Supplement to Figure 3 — the recruitment matrix does not depend "
             "on where the assignment rule is cut", fontsize=13, weight="bold")
save(fig, "fig3_supp_threshold_robustness")
plt.show()
# export_figure("fig3_supp_threshold_robustness")   # -> notebooks/outputs/fig3_supp_threshold_robustness.{png,pdf,svg}
# export_figure("fig3_supp_threshold_robustness", formats=("pdf",), outdir="~/Desktop")

# %%
for s in sweep3:
    print(f"frac={s.frac:.2f}  empty rows={s.n_empty:2d}  "
          f"mean row n={s.n_rows.mean():5.1f}  r_offdiag={s.r_offdiag:.3f}")
