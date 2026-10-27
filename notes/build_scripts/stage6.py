"""Stage 6: layout fixes found by rendering (Fig 7, S1, S2, Fig 6)."""
from pathlib import Path

P = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\notebooks\figures.py")
t = P.read_bytes().decode("utf-8").replace("\r\n", "\n")


def sub(old, new):
    global t
    assert t.count(old) == 1, (t.count(old), old[:70])
    t = t.replace(old, new)


# ---- Fig 7: the paper's words, and a W panel that fits ----
sub('''        observed_floor=fig7_obs_floor, **({} if j == 0 else {"ylabel": ""}))''',
    '''        observed_floor=fig7_obs_floor, label_observed="measured",
        label_predicted="predicted", **({} if j == 0 else {"ylabel": ""}))''')
sub('''P.pls.weight_blocks(ax_d, W6_reord, bnd6, centers6, norm=norm6,
                    title="W (low → high concentration mapping)")
hs.panel_label(ax_d, "D", x=-0.11, y=1.05)''',
    '''P.pls.weight_blocks(ax_d, W6_reord, bnd6, centers6, norm=norm6,
                    title="W (low → high)")
ax_d.title.set_fontsize(HEATMAP_TITLE_FONTSIZE)
ax_d.xaxis.label.set_fontsize(HEATMAP_TITLE_FONTSIZE)
ax_d.yaxis.label.set_fontsize(HEATMAP_TITLE_FONTSIZE)
ax_d.tick_params(labelsize=hs.FS_TICK)
hs.panel_label(ax_d, "D", x=-0.30, y=1.10)''')

# ---- S1: no title (the paper has none; it ran into the B letter), letters clear of the labels ----
sub('''    ax_var.set_title(f"{conc} concentration — chosen k = {chosen}",
                     fontsize=hs.FS_PANEL, loc="left")
''', '')
sub('''    for ax, letter in zip((ax_var, ax_slope, ax_ari), letters):
        hs.panel_label(ax, letter)''',
    '''    for ax, letter, (lx, ly) in zip((ax_var, ax_slope, ax_ari), letters,
                                    ((-0.20, 1.06), (-0.22, 1.06), (-0.09, 1.16))):
        hs.panel_label(ax, letter, x=lx, y=ly)''')

# ---- S2: class names on the left of the first heatmap, colourbar on the right (as in the paper) ----
sub("gs_heat = fig.add_gridspec(1, 4, left=0.11, right=0.97, top=0.767, bottom=0.097,",
    "gs_heat = fig.add_gridspec(1, 4, left=0.13, right=0.90, top=0.767, bottom=0.14,")
sub('''ax_heat[0].set_ylabel("odorants\\n(by chemical class)", fontsize=5)
ax_heat[0].yaxis.set_label_coords(-0.16, 0.5)   # clear of the class strip
''', '''fig.text(0.515, 0.03, "glomerulus #", ha="center", va="bottom", fontsize=hs.FS_LABEL)
''')
sub('''# Chemical-class strip on the leftmost heatmap only, in the plotted row order.''',
    '''# Chemical-class strip, with the class names, left of the first heatmap, in the plotted row order.''')
sub('''               orientation="y", labels=False)''',
    '''               orientation="y", side="left", labels=True, fontsize=CLASS_STRIP_FONTSIZE)''')
sub('''# The colourbar is shortened to 70% of the panel height and centred on it. It
# sits left of the y-label, which sits left of the class strip.
p0 = ax_heat[0].get_position()
_cb_h = p0.height * 0.70
cax = fig.add_axes([0.030, p0.y0 + (p0.height - _cb_h) / 2, 0.012, _cb_h])''',
    '''# The colourbar is shortened to 70% of the panel height and centred on it, right of the last heatmap.
p0 = ax_heat[-1].get_position()
_cb_h = p0.height * 0.70
cax = fig.add_axes([p0.x1 + 0.012, p0.y0 + (p0.height - _cb_h) / 2, 0.012, _cb_h])''')
sub('''cax.yaxis.set_ticks_position("left")
cax.yaxis.set_label_position("left")
cb.set_label("\\u0394F/F", fontsize=10)''', '''cb.set_label("\\u0394F (log)", fontsize=hs.FS_LABEL)''')

# ---- Fig 6: cluster numbers in cluster colours as the headings, as in the paper ----
sub('''        title=f"Cluster {k}")''', '''        title=str(k), title_fontsize=hs.FS_LABEL + 2)
    ax.title.set_color(hs.cluster_color(k))''')

P.write_bytes(t.encode("utf-8"))
print("stage 6 applied")
