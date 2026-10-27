# %%
FIG7_METRIC = "pct_err"      # "pct_err" (default) | "pct_var" | "dq2"
FIG7_ALPHA = 0.05            # BH-FDR level; the black boxes
FIG7_CLIP = 0.1              # p above this is the palest cell — colour only
FIG7_NULL = "broken_pairing"  # the published null; see the markdown above
FIG7_MIN_CLASS = 15          # glomeruli a chemical class needs to be tested
FIG7_BOX_COLOR = "#2b9743"   # the box on significant cells (green)
FIG7_P_CMAP = "OrRd"         # A keeps the original sequential ramp (hs.SEQ is magma now)
FIG7_XLABEL = "cluster (predicted at high)"
FIG7_YLABEL = "added non-source cluster\nprediction improvement"

fig7_clusters = transfer.pairwise(m, "clusters", alpha_fdr=FIG7_ALPHA,
                                  null=FIG7_NULL)
fig7_classes = transfer.pairwise(m, "chem_class", alpha_fdr=FIG7_ALPHA,
                                 null=FIG7_NULL)

for t in (fig7_clusters, fig7_classes):
    s = t.summary(FIG7_ALPHA)
    print(f"{s['grouping']:>11}: {s['n_significant']} of {s['n_pairs']} "
          f"ordered pairs at q < {FIG7_ALPHA} "
          f"({s['n_one_way']} one-way), N = {s['n_perm']:,}")
print()
print(P.transfer.edge_table(fig7_clusters, alpha=FIG7_ALPHA,
                            metric=FIG7_METRIC))

# %%
# One row: A and B, the cluster analysis. The chemical-class control (what
# were C and D) is no longer drawn; it is still computed, printed and checked
# in the cells below.
# Scale the whole figure by one factor, solved so the rendered width is
# MAIN_TARGET_W.
FIG7_SCALE = MAIN_TARGET_W / 13.0
FIG7_W, FIG7_H = 13.0 * FIG7_SCALE, 6.0 * FIG7_SCALE
fig = plt.figure(figsize=(FIG7_W, FIG7_H), dpi=FIG_DPI)

gs7 = fig.add_gridspec(1, 2, left=0.09, right=0.90, top=0.82, bottom=0.17,
                       wspace=0.42)

# Only the significant squares carry a number, and they are boxed in green.
ax_p = fig.add_subplot(gs7[0, 0])
P.transfer.pvalue_matrix(ax_p, fig7_clusters, alpha=FIG7_ALPHA, clip=FIG7_CLIP,
                         xlabel=FIG7_XLABEL, ylabel=FIG7_YLABEL, cmap=FIG7_P_CMAP,
                         box_color=FIG7_BOX_COLOR, text_significant_only=True)
hs.panel_label(ax_p, "A", x=-0.26, y=1.12)

ax_m = fig.add_subplot(gs7[0, 1])
P.transfer.improvement_matrix(ax_m, fig7_clusters, metric=FIG7_METRIC,
                              alpha=FIG7_ALPHA, xlabel=FIG7_XLABEL, ylabel="",
                              box_color=FIG7_BOX_COLOR,
                              text_significant_only=True)
hs.panel_label(ax_m, "B", x=-0.16, y=1.12)
P.transfer.significance_legend(ax_m, fig7_clusters, alpha=FIG7_ALPHA,
                               loc=(1.02, -0.16), box_color=FIG7_BOX_COLOR)

fig.text(0.5, 0.98, "Motif transfer between glomerular groups "
                    f"(broken-pairing null, N = {fig7_clusters.n_perm:,})",
         ha="center", va="top", fontsize=hs.FS_TITLE + 1)

save(fig, "fig7_motif_transfer")
plt.show()
# export_figure("fig7_motif_transfer")   # -> notebooks/outputs/fig7_motif_transfer.{png,pdf,svg}

# %%
# The control's actual claim, checked rather than asserted: each significant
# cluster edge is looked up against the class edge between the two clusters'
# dominant preferred classes. Seven of eight should match. If this ever prints
# fewer, the two analyses have drifted apart — which is exactly what the
# control exists to detect.
fig7_match = transfer.matching_edges(fig7_clusters, fig7_classes, m,
                                     alpha=FIG7_ALPHA)
print(f"{'cluster edge':>16}  {'ΔQ²':>8}  {'via':>26}  chemical-class counterpart")
for row in fig7_match:
    e, c = row["cluster_edge"], row["counterpart"]
    via = f"{row['via'][0]} → {row['via'][1]}"
    tail = (f"{c['dq2']:+.4f}  (q = {c['q']:.4f})" if c else "— none —")
    print(f"{e['source'] + ' → ' + e['target']:>16}  {e['dq2']:+8.4f}  "
          f"{via:>26}  {tail}")
n_match = sum(r["counterpart"] is not None for r in fig7_match)
print(f"\n{n_match} of {len(fig7_match)} cluster edges have a class counterpart")

# %%
# Per-group intra baselines — the denominator of both percentage readings, and
# the reason `pct_var` is not the default: cluster 5 is predicted from its own
# mean at Q² = 0.046 and cluster 2 at 0.438, so the same absolute increment
# would read ten times larger on one than on the other.
for t in (fig7_clusters, fig7_classes):
    print(f"{t.grouping.id:>11} intra Q²: " +
          "  ".join(f"{lab}={q:.3f}"
                    for lab, q in zip(t.labels, t.q2_intra)))
print(f"\nchemical classes tested: {list(fig7_classes.labels)}")
print(f"glomeruli excluded (class below n = {FIG7_MIN_CLASS}): "
      f"{fig7_classes.grouping.dropped} of {m.n_gloms}")
