# %%
# Membership only: every cluster member is one solid colour, and no response
# value is drawn -- so there is no concentration to choose and no colourbar.
FIG5_MEMBER_COLOR = "red"
FIG5_NONRESPONDER = "hollow"        # or "context"/"requested", or a dict of scatter kwargs
FIG5_CLUSTER_EDGE = "black"         # see below

# This section has no heatmap panel to key off (see §1), so it gets its own
# scale rather than reading FIG_SCALE directly, solved for MAIN_TARGET_W.
# `_fig5_k` is how much smaller that makes the page than it was under the
# old bare-FIG_SCALE sizing -- marker size (area, so k**2) and linewidth
# (linear, so k) are rescaled by it too, or a page this much smaller would
# leave 7 giant, overlapping dots per row instead of 7 legible clusters.
FIG5_SCALE = MAIN_TARGET_W / 16.0
_fig5_k = FIG5_SCALE / FIG_SCALE
FIG5_CLUSTER_SIZE = 62 * _fig5_k ** 2   # the cluster's own points, drawn on top
FIG5_CLUSTER_LW = 0.6 * _fig5_k
FIG5_CONTEXT_SIZE = 16 * _fig5_k ** 2   # "hollow"'s own base size (NONRESPONDER_MARKERS), scaled
FIG5_CONTEXT_LW = 0.4 * _fig5_k

fig = plt.figure(figsize=(16.0 * FIG5_SCALE, 3.6 * FIG5_SCALE), dpi=FIG_DPI)
gs5 = fig.add_gridspec(1, cl.N_CLUSTERS, wspace=0.06,
                       left=0.05, right=0.97, top=0.82, bottom=0.06)
xlim5, ylim5 = P.clustering.bulb_limits(d)

for k in range(1, cl.N_CLUSTERS + 1):
    ax = fig.add_subplot(gs5[0, k - 1])
    P.clustering.topography(
        d, ax, k, color=FIG5_MEMBER_COLOR,
        size=FIG5_CLUSTER_SIZE, edgecolor=FIG5_CLUSTER_EDGE,
        linewidth=FIG5_CLUSTER_LW, nonresponder_marker=FIG5_NONRESPONDER,
        context_size=FIG5_CONTEXT_SIZE, context_linewidth=FIG5_CONTEXT_LW,
        title=f"Cluster {k}")
    ax.set_xlim(*xlim5)
    ax.set_ylim(*ylim5)
    if k == 1:
        hs.panel_label(ax, "A", x=-0.30, y=1.12)

save(fig, "fig5_spatial_by_cluster")
plt.show()
# export_figure("fig5_spatial_by_cluster")   # -> notebooks/outputs/fig5_spatial_by_cluster.{png,pdf,svg}

# %%
for k in range(1, cl.N_CLUSTERS + 1):
    mask = d.cluster_mask(k)
    xy = d.xy[mask]
    spread = float(np.sqrt(((xy - xy.mean(0)) ** 2).sum(1).mean()))
    print(f"cluster {k}: n = {int(mask.sum()):3d}, "
          f"rms spread {spread:6.1f} px about its centroid")
