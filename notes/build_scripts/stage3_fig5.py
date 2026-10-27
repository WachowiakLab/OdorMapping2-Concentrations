"""Stage 3: new section 5 = the paper's Fig 5B (class x cluster, low beside high), after Fig 4."""
import re
from pathlib import Path

P = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\notebooks\figures.py")
text = P.read_bytes().decode("utf-8").replace("\r\n", "\n")
cells = [c.rstrip("\n") for c in re.split(r"(?m)^(?=# %%)", text) if c.startswith("# %%")]

MARK = ("# > **[Caption not yet written for the merged figure — see "
        "notes/merge_decisions.md]**")
md_human = f'''# %% [markdown] tags=["human"]
# ## 5. Figure 5 — which chemical classes each glomerular group responds to
#
{MARK}
#
# Panel B of the paper's Figure 5. For each of the seven glomerular groups, how many of its
# glomeruli respond to each chemical class, scaled so the group's top class is 1 — at low
# concentration (left) and at high concentration (right). Panel A (the per-odorant tuning
# bars) is drawn separately and is not built here.'''
md_llm = f'''# %% [markdown] tags=["llm"]
# ## 5. Figure 5 — chemical-class composition of each NMF cluster, low vs high
#
{MARK}
#
# Paper Fig 5B. `class_cluster_glom_composition(d, conc, order)`: for each (class, cluster) cell the
# number of distinct glomeruli of that cluster with a nonzero response to at least one odorant of
# that class, divided by the cluster's own maximum over classes, separately for each concentration
# (columns do not sum to 1). Class order is `CLASS_ORDER_FIG3`. Panel A is drawn elsewhere.'''

code = '''# %%
# Own constants — this section inherits no bare name from anywhere.
FIG5_CHEM_ORDER = CLASS_ORDER_FIG3     # same broad chemical order as Figs 3, 4 and 7 ("matt")
FIG5_CMAP = "magma"
FIG5_EMPTY = "gainsboro"
FIG5_BOUNDARY, FIG5_BOUNDARY_ALPHA = "blue", 0.5
FIG5_H = 2.9                           # inches; the width is MAIN_TARGET_W

# B's numbers, one row per chemical class and one column per cluster, computed once and reused
# by the drawing calls and the CSV export, so the saved numbers are exactly what the panel draws.
# Each concentration is normalised within itself (every cluster column / its own top class).
_fig5_low, _fig5_letters = P.clustering.class_cluster_glom_composition(d, "low", FIG5_CHEM_ORDER)
_fig5_high, _ = P.clustering.class_cluster_glom_composition(d, "high", FIG5_CHEM_ORDER)
fig5_table = pd.concat(
    {conc: pd.DataFrame(matrix, index=_fig5_letters,
                        columns=[f"cluster_{k}" for k in range(1, 8)])
     for conc, matrix in (("low", _fig5_low), ("high", _fig5_high))},
    names=["conc", "chem_class"])
if SAVE:
    fig5_csv = OUT / "fig5_class_by_cluster_B_glom_composition.csv"
    fig5_table.to_csv(fig5_csv)
    print(f"saved: {fig5_csv.relative_to(ROOT)}")

# %%
# --- Figure 5B: low concentration beside high concentration ----------------------
fig = plt.figure(figsize=(MAIN_TARGET_W, FIG5_H), dpi=FIG_DPI)
gs5 = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.045], wspace=0.55,
                       left=0.11, right=0.92, top=0.84, bottom=0.17)
summary_kw_5 = dict(class_order=FIG5_CHEM_ORDER, cmap=FIG5_CMAP, empty_color=FIG5_EMPTY,
                    boundary_color=FIG5_BOUNDARY, boundary_alpha=FIG5_BOUNDARY_ALPHA,
                    xlabel=None, label_fontsize=hs.FS_LABEL, tick_fontsize=hs.FS_TICK)
ax_low = fig.add_subplot(gs5[0, 0])
ax_high = fig.add_subplot(gs5[0, 1])
mesh_low = P.clustering.cluster_summary(d, ax_low, conc="low", **summary_kw_5)
mesh_high = P.clustering.cluster_summary(d, ax_high, conc="high", **summary_kw_5)
for ax_s, title in ((ax_low, "low concentration"), (ax_high, "high concentration")):
    ax_s.set_title(title, fontsize=hs.FS_LABEL)
    for lbl in ax_s.get_xticklabels():                      # cluster numbers, in cluster colours
        lbl.set_color(hs.cluster_color(int(lbl.get_text())))
    for lbl, letter in zip(ax_s.get_yticklabels(), _fig5_letters):   # class names, in class colours
        lbl.set_color(cl.class_color(letter))

cax5 = fig.add_subplot(gs5[0, 2])
cb5 = fig.colorbar(mesh_high, cax=cax5)
cb5.set_label("responding glomeruli,\\nrelative to cluster's top class", fontsize=hs.FS_LABEL)
cb5.ax.tick_params(labelsize=hs.FS_TICK)
hs.panel_label(ax_low, "B", x=-0.42, y=1.10)

save(fig, "fig5_class_by_cluster")
plt.show()
# export_figure("fig5_class_by_cluster")   # -> notebooks/outputs/fig5_class_by_cluster.{png,pdf,svg}'''

new_cells = [md_human, md_llm] + [c for c in code.split("\n\n# %%\n")]
new_cells = [new_cells[0], new_cells[1], new_cells[2]] + ["# %%\n" + new_cells[3]]
i4 = max(i for i, c in enumerate(cells) if 'save(fig, "fig4_tuning_by_cluster")' in c)
cells[i4 + 1:i4 + 1] = new_cells
P.write_bytes(("\n\n".join(cells) + "\n").encode("utf-8"))
print("inserted", len(new_cells), "cells after cell", i4)
