"""Stage 4b: Fig 3 in the paper's outline: A MACCS, B most-sensitive, C low, D high / E recruitment, F 2x2 spectra."""
from pathlib import Path

P = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\notebooks\figures.py")
t = P.read_bytes().decode("utf-8").replace("\r\n", "\n")


def sub(old, new):
    global t
    assert t.count(old) == 1, (t.count(old), old[:70])
    t = t.replace(old, new)


# constants cell: letters now follow the paper, so the draw-order switch is gone
sub('''FIG3_SPLIT_AROMATICS = True    # split "R" into 3 display subclasses on panels A-D
FIG3_PANEL_ORDER = "DABC"      # left-to-right draw order for the top row -- a
                               # panel's letter stays tied to its content
                               # (D is still the MACCS matrix) regardless of
                               # where it's drawn
SPECTRUM_SCALE = "per-panel"   # "per-panel" | "shared"
EG_GLOMERULI = similarity.EG_GLOMERULI     # Fig 6B-D must reuse these
# Fourth example glomerulus (panel H), by dataset glomerulus id like EG_GLOMERULI.
# Kept out of EG_GLOMERULI, which other code imports.
EXTRA_GLOM = 155''',
    '''FIG3_SPLIT_AROMATICS = True    # split "R" into 3 display subclasses on panels A-D
SPECTRUM_SCALE = "per-panel"   # "per-panel" | "shared"
EG_GLOMERULI = similarity.EG_GLOMERULI     # panel F's first three example glomeruli
# Fourth example glomerulus (panel F), by dataset glomerulus id like EG_GLOMERULI.
# Kept out of EG_GLOMERULI, which other code imports.
EXTRA_GLOM = 155''')
sub('# Panel I\'s matrix is computed here rather than inherited from §7, so the two',
    '# Panel E\'s matrix is computed here rather than inherited from §7, so the two')
sub('print(f"panel I: rows built', 'print(f"panel E: rows built')

# figure cell
a = t.index("# --- Figure 3 ---------------------------------------------------------------\nfig = plt.figure(figsize=(14.0, 11.6)")
b = t.index('save(fig, "fig3_correlational_structure")')
new_fig = '''# --- Figure 3 ---------------------------------------------------------------
# The paper's outline: row 1 is the four odorant x odorant matrices (A the chemical MACCS
# matrix, then the three response matrices); row 2 is the recruitment matrix (E, left) and
# the four example spectra in a 2 x 2 block (F, right).
fig = plt.figure(figsize=(14.0, 10.6), dpi=200)
gs = fig.add_gridspec(2, 4, height_ratios=[1.00, 1.20], hspace=0.40, wspace=0.30,
                      left=0.07, right=0.91, top=0.965, bottom=0.045)

top_panels = (("A", mats[3], "MACCS cosine"),
              ("B", mats[0], "'most sensitive'"),
              ("C", mats[1], "'low' (1x)"),
              ("D", mats[2], "'high' (25x)"))
axes_top4 = []
for j, (L, mat, title) in enumerate(top_panels):
    ax = fig.add_subplot(gs[0, j])
    # The x-axis chemical-class strip, with its family names, is drawn on all
    # four panels (they share the same class order).
    im = P.correlation_matrix(ax, mat, title=title, boundaries=False,
                              strip_labels=True)
    hs.panel_label(ax, L, x=-0.16)
    axes_top4.append(ax)
fig.canvas.draw()                          # aspect="equal" fixes the box size only now
pos_d = axes_top4[-1].get_position()
cax = fig.add_axes([pos_d.x1 + 0.012, pos_d.y0, 0.010, pos_d.height])
P.similarity_colorbar(fig, im, cax)

# E: the recruitment matrix, square, on the left half of row 2.
axE = fig.add_subplot(gs[1, 0:2])
axE.set_box_aspect(1)
imE = P.recruitment.matrix(axE, rmx3, "observed",
                           vmax=P.recruitment.shared_vmax(rmx3), ylabel=True,
                           title="binary recruitment matrix (observed)")
hs.panel_label(axE, "E", x=-0.16, y=1.05)
fig.canvas.draw()
pos_e = axE.get_position()
caxE = fig.add_axes([pos_e.x1 + 0.012, pos_e.y0, 0.010, pos_e.height])
cbE = fig.colorbar(imE, cax=caxE, extend="max")
cbE.set_label("recruitment probability", fontsize=hs.FS_LABEL)
cbE.ax.tick_params(labelsize=hs.FS_TICK)

# F: one 2 x 2 block of example spectra, filled row by row (173, 211 / 97, 155).
gs_F = gs[1, 2:4].subgridspec(2, 2, hspace=0.85, wspace=0.42)
LIM = P.shared_spectrum_limit(spectra) if SPECTRUM_SCALE == "shared" else None
for k, s in enumerate(spectra):
    ax = fig.add_subplot(gs_F[k // 2, k % 2])
    P.response_spectrum(ax, s, limit=LIM,
                        ylabel="Resp. amplitude" if k % 2 == 0 else "")
    ax.set_title(f"glomerulus {s.glom} · cluster {s.cluster}",
                 fontsize=hs.FS_LABEL)
    if k == 0:
        hs.panel_label(ax, "F", x=-0.30)

'''
t = t[:a] + new_fig + t[b:]
P.write_bytes(t.encode("utf-8"))
print("Fig 3 restructured; leftover FIG3_PANEL_ORDER:", t.count("FIG3_PANEL_ORDER"), "| titles_by_letter:", t.count("titles_by_letter"))
