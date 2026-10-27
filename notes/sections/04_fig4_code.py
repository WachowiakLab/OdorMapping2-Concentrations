# %%
# Every constant here is prefixed and owned by this section, so that no two
# sections can share a bare name like `VMIN`/`VMAX` and silently rebind it.
FIG4_VMIN, FIG4_VMAX = 1, 140      # dF/F colour range of the response heatmaps
FIG4_CMAP = "magma"
FIG4_EMPTY = "gainsboro"           # below vmin -> grey
FIG4_BOUNDARY, FIG4_BOUNDARY_ALPHA = "blue", 0.5

FIG4_CONC_A = "low"          # 4A is the low-concentration matrix
FIG4C_CLUSTERS = (1, 3)      # mixed acids/aldehydes cluster, furan cluster
# Reuses §3's own CLASS_ORDER_FIG3 directly (not an independent copy), so "the
# same broad chemical order as Fig 3/6" is guaranteed rather than kept in sync
# by hand -- same coupling §6 already has to §3's class order.
FIG4_CHEM_ORDER = CLASS_ORDER_FIG3
# Toggle for B/C's odorant order *within* each cluster block (columns stay
# glomerulus-cluster-sorted either way -- this only ever touches rows).
# True: chemical class (FIG4_CHEM_ORDER). False: native odorant order, i.e.
# `bicluster_order`'s own chem_order=None default -- ties in the stable sort
# keep their original order, so this is not "no sort", it's "no secondary
# sort key".
FIG4B_SORT_BY_CLASS = True

# Panel B's own row order -- independent of `bic` below, which panel C still
# reads. "cluster": the manuscript default, rows bi-clustered like the
# columns (diagonal block highlight makes sense). "odorant": rows follow
# Fig 3/6's own odorant order (FIG4_CHEM_ORDER/"burton") instead, ignoring
# cluster membership entirely -- a cluster's rows are then scattered, not
# contiguous, so the diagonal highlight is turned off for this mode (wired
# in below via `diagonal_blocks`). Vertical glomerulus-cluster lines draw
# either way -- only the rows move.
FIG4B_ROW_ORDER = "odorant"

# Panels B and C both: peak-normalise each glomerulus (column) to [0, 1] by its
# own max response,
# on the raw dF/F values, and draw on a LINEAR scale instead of log -- the
# two go together (see `biclustered`'s `col_normalize` docstring). `FIG4_VMIN`
# /`FIG4_VMAX` (tuned for raw dF/F on a log scale) do not apply in this mode;
# `FIG4B_NORM_VMIN`/`FIG4B_NORM_VMAX` are panel B's own limits instead.
FIG4B_COL_NORMALIZE = True
FIG4B_NORM_VMIN, FIG4B_NORM_VMAX = 0.0, 1.0
# Restores the "nonzero cutoff" (the same FIG4_VMIN=1 floor the non-normalised
# log view uses) even though column-normalising otherwise erases it: a cell
# whose RAW dF/F was below this never draws on the ramp, however large its
# normalised value would have been -- so a glomerulus whose peak never
# clears the noise floor greys out entirely instead of showing a misleadingly
# "bright" 1.0 for its own weak best odorant. `None` disables this and lets
# every normalised value draw, including noise-floor-weak ones.
FIG4B_RAW_FLOOR = FIG4_VMIN

gsim = similarity.glomerulus_correlation(d, FIG4_CONC_A)
# The block boundaries (and the diagonal alignment with the glomerulus
# columns) are the same regardless of FIG4B_SORT_BY_CLASS; only the order
# inside a block moves. Shared by C (and, when FIG4B_ROW_ORDER="cluster", by
# B too); B reads its own `bic_b` below whenever FIG4B_ROW_ORDER="odorant",
# so C's row order never moves when B's does.
bic = P.clustering.bicluster_order(
    d, chem_order=(FIG4_CHEM_ORDER if FIG4B_SORT_BY_CLASS else None),
    within_order=WITHIN_CLASS_ORDER)
bic_b = P.clustering.bicluster_order(
    d, chem_order=FIG4_CHEM_ORDER, within_order=WITHIN_CLASS_ORDER,
    row_order=FIG4B_ROW_ORDER)

# Panel D's underlying numbers, one row per chemical class, one column per
# cluster -- computed once and reused for both the drawing calls below and
# the CSV export, so the saved numbers are exactly what the panel draws.
_fig4d_low, _fig4d_letters = P.clustering.class_cluster_glom_composition(
    d, "low", FIG4_CHEM_ORDER)
_fig4d_high, _ = P.clustering.class_cluster_glom_composition(
    d, "high", FIG4_CHEM_ORDER)

fig4d_table = pd.concat(
    {conc: pd.DataFrame(matrix, index=_fig4d_letters,
                        columns=[f"cluster_{k}" for k in range(1, 8)])
     for conc, matrix in (("low", _fig4d_low), ("high", _fig4d_high))},
    names=["conc", "chem_class"])
if SAVE:
    fig4d_csv = OUT / "fig4_tuning_by_cluster_D_glom_composition.csv"
    fig4d_table.to_csv(fig4d_csv)
    print(f"saved: {fig4d_csv.relative_to(ROOT)}")

bm = gsim.block_means()
within = np.diag(bm.values)
between = bm.values[~np.eye(len(bm), dtype=bool)]
print(f"panel A: within-cluster correlation {within.min():.2f}-{within.max():.2f}, "
      f"between-cluster mean {between.mean():+.3f}")
for k in FIG4C_CLUSTERS:
    n = [int((d.subject_mask(s) & d.cluster_mask(k)).sum()) for s in d.subjects]
    print(f"panel C: cluster {k} has {sum(n)} glomeruli, {n} per subject")

# %%
# Scale the whole figure by one factor rather than rebuild its gridspec
# math -- every gridspec call below is fraction-based, so it shrinks
# proportionally for free. Solved so the rendered width is MAIN_TARGET_W.
FIG4_SCALE = MAIN_TARGET_W / 14.0
# A-D keep their original 11.0-unit hand layout; E (the flow diagram) is a third row of
# FIG4_E_H units underneath. Every fraction-based call for A-D is remapped into the top part
# of the taller canvas -- `_y4` for positions, `_k4` for lengths along the height -- so A-D
# come out the same size as before.
FIG4_H0, FIG4_E_H = 11.0, 3.6
FIG4_W, FIG4_H = 16.0 * FIG4_SCALE, (FIG4_H0 + FIG4_E_H) * FIG4_SCALE
_k4 = FIG4_H0 / (FIG4_H0 + FIG4_E_H)           # lengths along the old height
_y0_4 = FIG4_E_H / (FIG4_H0 + FIG4_E_H)        # where the old canvas starts, from the bottom


def _y4(y):
    return y * _k4 + _y0_4


fig = plt.figure(figsize=(FIG4_W, FIG4_H), dpi=FIG_DPI)

# --- top row: A (square correlation matrix) and B (wide response heatmap) ---
gs_top = fig.add_gridspec(1, 2, width_ratios=[1, 3.1], wspace=0.16,
                          left=0.06, right=0.93, top=_y4(0.95), bottom=_y4(0.60))
ax_a = fig.add_subplot(gs_top[0, 0])
ax_b = fig.add_subplot(gs_top[0, 1])

im_a = P.clustering.glomerulus_correlation_matrix(
    ax_a, gsim, title={"low": "low", "high": "high"}[FIG4_CONC_A],
    boundaries=False, cluster_boxes=True, cluster_box_color="green")
mesh_b = P.clustering.biclustered(d, ax_b, conc=FIG4_CONC_A, order=bic_b,
                                  diagonal_blocks=(FIG4B_ROW_ORDER == "cluster"),
                                  class_names=True,
                                  col_normalize=FIG4B_COL_NORMALIZE,
                                  raw_floor=(FIG4B_RAW_FLOOR
                                            if FIG4B_COL_NORMALIZE else None),
                                  vmin=(FIG4B_NORM_VMIN if FIG4B_COL_NORMALIZE
                                       else FIG4_VMIN),
                                  vmax=(FIG4B_NORM_VMAX if FIG4B_COL_NORMALIZE
                                       else FIG4_VMAX),
                                  cmap=FIG4_CMAP, empty_color=FIG4_EMPTY,
                                  boundary_color=FIG4_BOUNDARY,
                                  boundary_alpha=FIG4_BOUNDARY_ALPHA)

pos_a = ax_a.get_position()
cax_a = fig.add_axes([pos_a.x0, pos_a.y0 - 0.055 * _k4, pos_a.width * 0.6,
                      0.014 * _k4])
cb_a = fig.colorbar(im_a, cax=cax_a, orientation="horizontal",
                    ticks=[0.0, 0.5, 1.0])
cb_a.set_label(r"correlation $\rho$", fontsize=hs.FS_LABEL)
cb_a.ax.tick_params(labelsize=hs.FS_TICK)
cb_a.outline.set_visible(False)

pos_b = ax_b.get_position()
cax_b = fig.add_axes([pos_b.x1 + 0.012, pos_b.y0, 0.010, pos_b.height * 0.7])
if FIG4B_COL_NORMALIZE:
    # Linear 0-1: each glomerulus's response relative to its own peak, not
    # an absolute ΔF/F -- pairs with `col_normalize` in the draw call above.
    cb_b = fig.colorbar(
        plt.cm.ScalarMappable(
            norm=plt.Normalize(vmin=FIG4B_NORM_VMIN, vmax=FIG4B_NORM_VMAX),
            cmap=hs.sequential_cmap(under=FIG4_EMPTY, bad=FIG4_EMPTY)),
        cax=cax_b)
    cb_b.set_label("resp., rel. to glom. peak", fontsize=hs.FS_LABEL)
    cb_b.set_ticks([0.0, 0.5, 1.0])
    cb_b.set_ticklabels(["0", "0.5", "1"])
else:
    cb_b = fig.colorbar(
        plt.cm.ScalarMappable(norm=LogNorm(vmin=FIG4_VMIN, vmax=FIG4_VMAX),
                              cmap=hs.sequential_cmap(under=FIG4_EMPTY,
                                                      bad=FIG4_EMPTY)),
        cax=cax_b)
    cb_b.set_label("ΔF/F", fontsize=hs.FS_LABEL)
    cb_b.set_ticks([1, 10, 100])
    cb_b.set_ticklabels(["1", "10", "100"])
cb_b.ax.minorticks_off()
cb_b.ax.tick_params(labelsize=hs.FS_TICK)

# --- bottom row: C (two clusters x four subjects) and D (class x cluster) ---
gs_bot = fig.add_gridspec(
    1, 3, width_ratios=[1, 1, 1.15], wspace=0.30,
    left=0.06, right=0.93, top=_y4(0.46), bottom=_y4(0.06))

# C's rows are B's rows: Fig 3/6's exact 59-odorant order (CLASS_ORDER_FIG3
# blocks, WITHIN_CLASS_ORDER within them), not the cluster-sorted rows `bic`
# carries. Only the rows move: columns stay one cluster's glomeruli per
# subject. C also takes B's peak-normalisation and raw floor, so a glomerulus
# reads identically in both and B's colorbar is C's too.
perm_c, _, _, class_bounds_c = similarity.class_sorted(
    d, CLASS_ORDER_FIG3, WITHIN_CLASS_ORDER)
assert FIG4B_ROW_ORDER == "odorant" and np.array_equal(perm_c, bic_b.rows), \
    "C's rows must be B's rows (FIG4B_ROW_ORDER='odorant', same class/within order)"
bic_c = dataclasses.replace(bic,
                            rows=perm_c,
                            row_boundaries=np.array([], dtype=int),
                            row_edges=np.array([0, len(perm_c)]))
FIG4C_LINE_COLOR, FIG4C_LINE_LW, FIG4C_LINE_ALPHA = "white", 0.6, 0.8

im_c = None
for g, k in enumerate(FIG4C_CLUSTERS):
    inner = GridSpecFromSubplotSpec(1, len(d.subjects), subplot_spec=gs_bot[0, g],
                                    wspace=0.18)
    axes_c = [fig.add_subplot(inner[0, j]) for j in range(len(d.subjects))]
    for j, subj in enumerate(d.subjects):
        im_c = P.clustering.per_subject_heatmap(
            d, axes_c[j], subj, cluster=k, conc=FIG4_CONC_A, order=bic_c,
            col_normalize=FIG4B_COL_NORMALIZE,
            raw_floor=FIG4B_RAW_FLOOR if FIG4B_COL_NORMALIZE else None,
            vmin=FIG4B_NORM_VMIN if FIG4B_COL_NORMALIZE else FIG4_VMIN,
            vmax=FIG4B_NORM_VMAX if FIG4B_COL_NORMALIZE else FIG4_VMAX,
            cmap=FIG4_CMAP, empty_color=FIG4_EMPTY, cluster_labels=False,
            ylabel=("odorants (by chemical class)" if j == 0 else None))
        for b in class_bounds_c:               # line between chemical classes
            axes_c[j].axhline(b - 0.5, color=FIG4C_LINE_COLOR,
                              lw=FIG4C_LINE_LW, alpha=FIG4C_LINE_ALPHA)
    mid = (axes_c[0].get_position().x0 + axes_c[-1].get_position().x1) / 2
    fig.text(mid, axes_c[0].get_position().y1 + 0.045 * _k4, f"cluster {k}",
             ha="center", va="bottom", fontsize=hs.FS_LABEL,
             color=hs.cluster_color(k))
    fig.text(mid, axes_c[0].get_position().y0 - 0.032 * _k4, "glomeruli",
             ha="center", va="top", fontsize=hs.FS_TICK)

# D — the class x cluster summary. Each concentration is normalised within itself
# (every cluster column / its own top class, separately for low and high), so low
# and high each get their own colour scale and colorbar.
inner_d = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_bot[0, 2], hspace=0.30)
ax_d_low = fig.add_subplot(inner_d[0, 0])
ax_d_high = fig.add_subplot(inner_d[1, 0])
summary_kw_4 = dict(class_order=FIG4_CHEM_ORDER, cmap=FIG4_CMAP,
                    empty_color=FIG4_EMPTY, boundary_color=FIG4_BOUNDARY,
                    boundary_alpha=FIG4_BOUNDARY_ALPHA)
mesh_d_low = P.clustering.cluster_summary(d, ax_d_low, conc="low", xlabel=None,
                                          **summary_kw_4)
mesh_d_high = P.clustering.cluster_summary(d, ax_d_high, conc="high",
                                           **summary_kw_4)

pos_dl = ax_d_low.get_position()
d4_label = "responding glomeruli, relative to cluster's top class"
for ax_d, mesh_d in ((ax_d_low, mesh_d_low), (ax_d_high, mesh_d_high)):
    pos_d = ax_d.get_position()
    cax_d = fig.add_axes([pos_d.x1 + 0.012, pos_d.y0, 0.010, pos_d.height])
    cb_d = fig.colorbar(mesh_d, cax=cax_d)
    cb_d.set_label(d4_label, fontsize=hs.FS_LABEL)
    cb_d.ax.tick_params(labelsize=hs.FS_TICK)

# --- E: where each glomerulus goes when the recruited non-responders are folded in ---
# Glomeruli silent at low but responding at high have their high column folded into the
# low matrix and NMF is run again. That fit is cached by `dev/fold_nonresponders.py` (seeded,
# same settings as the original clustering) and only READ here -- never refit at draw time.
# Left: the original clusters 1-7 (grey = glomeruli not in the original clustering, i.e. the
# folded-in ones). Right: the new clusters, numbered 1-7 by best overlap with the originals.
fig4_flows = folding.flows(d)
gs_e = fig.add_gridspec(1, 1, left=0.10, right=0.90, top=_y0_4 - 0.01, bottom=0.02)
ax_e = fig.add_subplot(gs_e[0, 0])
P.sankey.flow_diagram(ax_e, fig4_flows, left_title="original clustering",
                      right_title="with non-responders folded in")
print(f"E: {fig4_flows.n_matched} original + {fig4_flows.n_added} folded-in glomeruli; "
      f"share of each original cluster that stays in its matched new cluster: "
      f"{', '.join(f'{c + 1}: {s:.0%}' for c, s in enumerate(fig4_flows.stayed))}")

hs.panel_label(ax_a, "A", x=-0.22, y=1.12)
hs.panel_label(ax_b, "B", x=-0.035, y=1.06)
fig.text(0.02, _y4(0.48), "C", fontsize=hs.FS_LABEL + 4, weight="bold", va="top")
fig.text(pos_dl.x0 - 0.055, _y4(0.48), "D", fontsize=hs.FS_LABEL + 4,
         weight="bold", va="top")
fig.text(0.02, _y0_4 - 0.01, "E", fontsize=hs.FS_LABEL + 4, weight="bold", va="top")

save(fig, "fig4_tuning_by_cluster")
plt.show()
# export_figure("fig4_tuning_by_cluster")   # -> notebooks/outputs/fig4_tuning_by_cluster.{png,pdf,svg}
