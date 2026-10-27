"""S2: two separate figures. A = lab's heatmaps (class-ordered rows, linear magma), B = Fig 5 styling."""
from pathlib import Path

ROOT = Path(r"C:\Users\u6068031\Documents\highconc-merge")
FINAL = ROOT / "final"
PANEL = FINAL / "highconc" / "panels" / "clustering.py"


def rd(p):
    return p.read_bytes().decode("utf-8").replace("\r\n", "\n")


def edit(text, edits):
    for old, new in edits:
        assert text.count(old) == 1, (text.count(old), old[:70])
        text = text.replace(old, new)
    return text


# ---- panel: per_subject_heatmap(norm=...) and per_subject_topography(context_alpha=...) ----
t = edit(rd(PANEL), [
    # heatmap signature (anchored on the col_normalize/raw_floor lines I added earlier)
    ('''    col_normalize: bool = False,
    raw_floor: float | None = None,
    boundaries: bool = True,
    boundary_color: str = "#1f3fff",
    boundary_lw: float = 0.7,
    boundary_alpha: float = 0.45,
    cluster_labels: bool = True,
    label_fontsize: float = 7.5,''',
     '''    col_normalize: bool = False,
    raw_floor: float | None = None,
    norm=None,
    boundaries: bool = True,
    boundary_color: str = "#1f3fff",
    boundary_lw: float = 0.7,
    boundary_alpha: float = 0.45,
    cluster_labels: bool = True,
    label_fontsize: float = 7.5,'''),
    ('''    cluster_labels : bool
        Write each block's cluster number above it, in that cluster's colour.
        With `cluster` set there is only one block, so this writes one number.
''',
     '''    norm : matplotlib Normalize, optional
        The colour normalisation, used as given (`vmin`/`vmax` are then
        ignored). Omitted, it is log with `vmin`/`vmax`, or linear 0-1 with
        `col_normalize`. Pass a linear `Normalize` to draw raw dF/F on a
        linear scale, as Fig S2A does.
    cluster_labels : bool
        Write each block's cluster number above it, in that cluster's colour.
        With `cluster` set there is only one block, so this writes one number.
'''),
    ('''    norm = (plt.Normalize(vmin=vmin, vmax=vmax) if col_normalize
            else LogNorm(vmin=vmin, vmax=vmax))
    im = ax.imshow(data, aspect="auto", cmap=heat_cmap, norm=norm,''',
     '''    if norm is None:
        norm = (plt.Normalize(vmin=vmin, vmax=vmax) if col_normalize
                else LogNorm(vmin=vmin, vmax=vmax))
    im = ax.imshow(data, aspect="auto", cmap=heat_cmap, norm=norm,'''),
    # topography
    ('''    context_size: float = 14,
    context_linewidth: float = 0.4,
    limits=None,''',
     '''    context_size: float = 14,
    context_linewidth: float = 0.4,
    context_alpha: float | None = None,
    limits=None,'''),
    ('''        ax.scatter(x, y, facecolors="none", edgecolors=context_color,
                   s=context_size, linewidths=context_linewidth, zorder=1)''',
     '''        ax.scatter(x, y, facecolors="none", edgecolors=context_color,
                   s=context_size, linewidths=context_linewidth,
                   alpha=context_alpha, zorder=1)'''),
])
PANEL.write_bytes(t.encode("utf-8"))

# ---- notebook cells ---------------------------------------------------------
cells = '''# %%
# Own constants — this section inherits no bare name from anywhere.
# Two separate figures: A (per-subject heatmaps) and B (per-subject maps).
SUPP_SUBJ_DPI = 200
SUPP_SUBJ_LABEL_FS = 7
SUPP_SUBJ_SCALE = SUPP_TARGET_W / 12.5
_supp_subj_k = SUPP_SUBJ_SCALE / (HEATMAP_IN / 2.352)   # vs. this section's old scale
SUPP_SUBJ_TITLE_PAD = 22 * _supp_subj_k   # lifts 'S1'..'S4' clear of the 1-7 cluster row
# A: raw dF/F on a linear scale, odorant rows in Fig 3 / Fig 4's chemical-class order.
SUPP_SUBJ_HEAT_VMIN, SUPP_SUBJ_HEAT_VMAX = 0.0, 100.0
SUPP_SUBJ_CMAP = "magma"
SUPP_SUBJ_ODOR_ORDER = CLASS_ORDER_FIG3        # "matt", the same blocks as Fig 3 and Fig 4C
SUPP_SUBJ_ODOR_WITHIN = WITHIN_CLASS_ORDER     # "matt", within each block

# Shared odorant row sort: chemical-class blocks in the Fig 3 order, so these
# panels line up with the Fig 3 matrices. The columns stay sorted by cluster.
# The row boundaries are cleared -- they were odorant-cluster lines, which no
# longer correspond to the rows.
_perm, *_ = similarity.class_sorted(d, SUPP_SUBJ_ODOR_ORDER, SUPP_SUBJ_ODOR_WITHIN)
supp_subj_order = dataclasses.replace(
    P.clustering.bicluster_order(d),
    rows=_perm,
    row_boundaries=np.array([], dtype=int),
    row_edges=np.array([0, len(_perm)]),
)
supp_subj_norm = Normalize(SUPP_SUBJ_HEAT_VMIN, SUPP_SUBJ_HEAT_VMAX)

# %%
# --- Figure S2A: per-subject heatmaps ---------------------------------------
# Canvas in the original 12.5in-wide units (3.6in tall: the old A block plus its
# title room), scaled so the width is exactly SUPP_TARGET_W.
SUPP_SUBJ_W, SUPP_SUBJ_H = 12.5 * SUPP_SUBJ_SCALE, 3.6 * SUPP_SUBJ_SCALE
fig = plt.figure(figsize=(SUPP_SUBJ_W, SUPP_SUBJ_H), dpi=SUPP_SUBJ_DPI)

gs_heat = fig.add_gridspec(1, 4, left=0.11, right=0.97, top=0.767, bottom=0.097,
                           wspace=0.14)
ax_heat = [fig.add_subplot(gs_heat[0, i]) for i in range(4)]
for ax, s in zip(ax_heat, d.subjects):
    P.clustering.per_subject_heatmap(d, ax, s, conc="low",
                                     order=supp_subj_order,
                                     cmap=SUPP_SUBJ_CMAP,
                                     norm=supp_subj_norm,
                                     boundary_color="white")
    # The subject title and the coloured cluster numbers both live above the
    # axes and were landing on the same line. The numbers are written by the
    # panel just above the data; the title is lifted clear of them here.
    ax.set_title(ax.get_title(), fontsize=7, pad=SUPP_SUBJ_TITLE_PAD)
ax_heat[0].set_ylabel("odorants\\n(by chemical class)", fontsize=5)
ax_heat[0].yaxis.set_label_coords(-0.16, 0.5)   # clear of the class strip

# Chemical-class strip on the leftmost heatmap only, in the plotted row order.
_row_cls = d.chem_class[supp_subj_order.rows]
_cls_name = {c: cl.class_name(c) for c in dict.fromkeys(_row_cls)}
hs.class_strip(ax_heat[0], [_cls_name[c] for c in _row_cls],
               colors={_cls_name[c]: cl.class_color(c) for c in _cls_name},
               orientation="y", labels=False)

# The colourbar is shortened to 70% of the panel height and centred on it. It
# sits left of the y-label, which sits left of the class strip.
p0 = ax_heat[0].get_position()
_cb_h = p0.height * 0.70
cax = fig.add_axes([0.030, p0.y0 + (p0.height - _cb_h) / 2, 0.012, _cb_h])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=supp_subj_norm,
                                        cmap=SUPP_SUBJ_CMAP),
                  cax=cax, extend="max")
cax.yaxis.set_ticks_position("left")
cax.yaxis.set_label_position("left")
cb.set_label("\\u0394F/F", fontsize=10)
cb.set_ticks(np.linspace(SUPP_SUBJ_HEAT_VMIN, SUPP_SUBJ_HEAT_VMAX, 5))
cb.ax.tick_params(labelsize=8)
cb.outline.set_visible(False)

save(fig, "supp_per_subject_heatmaps", dpi=SUPP_SUBJ_DPI)
plt.show()
# export_figure("supp_per_subject_heatmaps")   # -> notebooks/outputs/supp_per_subject_heatmaps.{png,pdf,svg}

# %%
# --- Figure S2B: where each mouse's members of each cluster sit -------------
# Drawn exactly like Fig 5: the same page scale (so each map is Fig 5's size),
# the same marker sizes, black outline on the members and the same pale hollow
# context. Restated here, not read from §5, so each section runs on its own.
# Only the fill differs: it changes with the cluster, as Fig 5's red does not.
SUPP_SUBJ_TOPO_SCALE = MAIN_TARGET_W / 16.0
_supp_topo_k = SUPP_SUBJ_TOPO_SCALE / FIG_SCALE
SUPP_SUBJ_CLUSTER_EDGE = "black"
SUPP_SUBJ_CLUSTER_SIZE = 62 * _supp_topo_k ** 2
SUPP_SUBJ_CLUSTER_LW = 0.6 * _supp_topo_k
SUPP_SUBJ_CONTEXT_COLOR = P.clustering.CONTEXT_COLOR   # "darkgray", as Fig 5's hollow dots
SUPP_SUBJ_CONTEXT_ALPHA = 0.5
SUPP_SUBJ_CONTEXT_SIZE = 16 * _supp_topo_k ** 2
SUPP_SUBJ_CONTEXT_LW = 0.4 * _supp_topo_k
SUPP_SUBJ_TOPO_H = 5.2      # 4 mice, each map about as tall as Fig 5's, plus the title row

supp_subj_limits = P.clustering.bulb_limits(d)      # shared spatial frame

fig = plt.figure(figsize=(MAIN_TARGET_W, SUPP_SUBJ_TOPO_H), dpi=SUPP_SUBJ_DPI)
gs_topo = fig.add_gridspec(4, cl.N_CLUSTERS, left=0.07, right=0.97, top=0.92,
                           bottom=0.02, wspace=0.06, hspace=0.06)
for si, s in enumerate(d.subjects):
    for j, number in enumerate(range(1, cl.N_CLUSTERS + 1)):
        P.clustering.per_subject_topography(
            d, fig.add_subplot(gs_topo[si, j]), s, number,
            limits=supp_subj_limits,
            size=SUPP_SUBJ_CLUSTER_SIZE, edgecolor=SUPP_SUBJ_CLUSTER_EDGE,
            linewidth=SUPP_SUBJ_CLUSTER_LW,
            context_color=SUPP_SUBJ_CONTEXT_COLOR,
            context_alpha=SUPP_SUBJ_CONTEXT_ALPHA,
            context_size=SUPP_SUBJ_CONTEXT_SIZE,
            context_linewidth=SUPP_SUBJ_CONTEXT_LW,
            empty_fontsize=8,
            title=(si == 0),                     # cluster number on row 1 only
            ylabel=(f"S{s}" if j == 0 else None))

save(fig, "supp_per_subject_topography", dpi=SUPP_SUBJ_DPI)
plt.show()
# export_figure("supp_per_subject_topography")   # -> notebooks/outputs/supp_per_subject_topography.{png,pdf,svg}
'''
out = FINAL / "notes" / "sections" / "09_figS2_code.py"
out.write_text(cells, encoding="utf-8", newline="\n")
print("wrote panel patch and", out, len(cells.splitlines()), "lines")
