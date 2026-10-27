"""Stage 2: Fig 4 = A, B, C, and the Sankey as D (bottom right). Class x cluster summary leaves (-> Fig 5)."""
from pathlib import Path

P = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\notebooks\figures.py")
t = P.read_bytes().decode("utf-8").replace("\r\n", "\n")


def sub(old, new):
    global t
    assert t.count(old) == 1, (t.count(old), old[:70])
    t = t.replace(old, new)


# --- constants cell: the class x cluster numbers now belong to Fig 5 ---
a = t.index("# Panel D's underlying numbers, one row per chemical class")
b = t.index("bm = gsim.block_means()")
t = t[:a] + t[b:]

# --- undo the third-row remap (E is gone; D is the Sankey in D's old slot) ---
sub('''FIG4_SCALE = MAIN_TARGET_W / 14.0
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

''', '''FIG4_SCALE = MAIN_TARGET_W / 14.0
FIG4_W, FIG4_H = 16.0 * FIG4_SCALE, 11.0 * FIG4_SCALE
''')
sub("left=0.06, right=0.93, top=_y4(0.95), bottom=_y4(0.60))", "left=0.06, right=0.93, top=0.95, bottom=0.60)")
sub('''cax_a = fig.add_axes([pos_a.x0, pos_a.y0 - 0.055 * _k4, pos_a.width * 0.6,
                      0.014 * _k4])''',
    "cax_a = fig.add_axes([pos_a.x0, pos_a.y0 - 0.055, pos_a.width * 0.6, 0.014])")
sub("left=0.06, right=0.93, top=_y4(0.46), bottom=_y4(0.06))", "left=0.06, right=0.93, top=0.46, bottom=0.06)")
sub('fig.text(mid, axes_c[0].get_position().y1 + 0.045 * _k4, f"cluster {k}",',
    'fig.text(mid, axes_c[0].get_position().y1 + 0.045, f"cluster {k}",')
sub('fig.text(mid, axes_c[0].get_position().y0 - 0.032 * _k4, "glomeruli",',
    'fig.text(mid, axes_c[0].get_position().y0 - 0.032, "glomeruli",')

# --- D (class summary) + E (Sankey, third row) -> D = Sankey in the bottom-right slot ---
a = t.index("# D — the class x cluster summary. Each concentration is normalised within itself")
b = t.index('save(fig, "fig4_tuning_by_cluster")')
new_d = '''# D — where each glomerulus goes when the recruited non-responders are folded in.
# Glomeruli silent at low but responding at high have their high column folded into the
# low matrix and NMF is run again. That fit is cached by `dev/fold_nonresponders.py` (seeded,
# same settings as the original clustering) and only READ here -- never refit at draw time.
# Left: the original clusters 1-7 (grey = the recruited glomeruli, which were not in the
# original clustering). Right: the new clusters, numbered 1-7 by best overlap with the originals.
fig4_flows = folding.flows(d)
ax_d = fig.add_subplot(gs_bot[0, 2])
P.sankey.flow_diagram(ax_d, fig4_flows, left_title="Low-Conc.\\nOnly",
                      right_title="Low-Conc. +\\nHigh-Conc. recruited",
                      unmatched_label="recruited\\nat high")
print(f"D: {fig4_flows.n_matched} original + {fig4_flows.n_added} recruited glomeruli; "
      f"share of each original cluster that stays in its matched new cluster: "
      f"{', '.join(f'{c + 1}: {s:.0%}' for c, s in enumerate(fig4_flows.stayed))}")

hs.panel_label(ax_a, "A", x=-0.22, y=1.12)
hs.panel_label(ax_b, "B", x=-0.035, y=1.06)
pos_d = ax_d.get_position()
fig.text(0.02, 0.48, "C", fontsize=hs.FS_LABEL + 4, weight="bold", va="top")
fig.text(pos_d.x0 - 0.055, 0.48, "D", fontsize=hs.FS_LABEL + 4,
         weight="bold", va="top")

'''
t = t[:a] + new_d + t[b:]
P.write_bytes(t.encode("utf-8"))
print("stage 2 applied; leftover tokens:", {k: t.count(k) for k in ("_y4", "_k4", "_y0_4", "FIG4_E_H", "ax_e", "ax_d_low")})
