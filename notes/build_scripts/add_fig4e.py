"""Add panel E (the folded-in flow diagram) to Fig 4: third row, A-D unchanged in size."""
from pathlib import Path

F = Path(r"C:\Users\u6068031\Documents\highconc-merge\final")

FIG4_EDITS = [
    ('''FIG4_SCALE = MAIN_TARGET_W / 14.0
FIG4_W, FIG4_H = 16.0 * FIG4_SCALE, 11.0 * FIG4_SCALE
''',
     '''FIG4_SCALE = MAIN_TARGET_W / 14.0
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

'''),
    ('''left=0.06, right=0.93, top=0.95, bottom=0.60)''',
     '''left=0.06, right=0.93, top=_y4(0.95), bottom=_y4(0.60))'''),
    ('''cax_a = fig.add_axes([pos_a.x0, pos_a.y0 - 0.055, pos_a.width * 0.6, 0.014])''',
     '''cax_a = fig.add_axes([pos_a.x0, pos_a.y0 - 0.055 * _k4, pos_a.width * 0.6,
                      0.014 * _k4])'''),
    ('''left=0.06, right=0.93, top=0.46, bottom=0.06)''',
     '''left=0.06, right=0.93, top=_y4(0.46), bottom=_y4(0.06))'''),
    ('''fig.text(mid, axes_c[0].get_position().y1 + 0.045, f"cluster {k}",''',
     '''fig.text(mid, axes_c[0].get_position().y1 + 0.045 * _k4, f"cluster {k}",'''),
    ('''fig.text(mid, axes_c[0].get_position().y0 - 0.032, "glomeruli",''',
     '''fig.text(mid, axes_c[0].get_position().y0 - 0.032 * _k4, "glomeruli",'''),
    ('''hs.panel_label(ax_a, "A", x=-0.22, y=1.12)
hs.panel_label(ax_b, "B", x=-0.035, y=1.06)
fig.text(0.02, 0.48, "C", fontsize=hs.FS_LABEL + 4, weight="bold", va="top")
fig.text(pos_dl.x0 - 0.055, 0.48, "D", fontsize=hs.FS_LABEL + 4,
         weight="bold", va="top")
''',
     '''# --- E: where each glomerulus goes when the recruited non-responders are folded in ---
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
'''),
]


def apply(path, edits):
    t = path.read_bytes().decode("utf-8").replace("\r\n", "\n")
    for old, new in edits:
        assert t.count(old) == 1, (path.name, t.count(old), old[:70])
        t = t.replace(old, new)
    path.write_bytes(t.encode("utf-8"))
    print("edited", path.relative_to(F))


apply(F / "notes" / "sections" / "04_fig4_code.py", FIG4_EDITS)
apply(F / "notebooks" / "figures.py", FIG4_EDITS + [(
    "    dimsweep,\n    flow,\n    load,", "    dimsweep,\n    flow,\n    folding,\n    load,")])
