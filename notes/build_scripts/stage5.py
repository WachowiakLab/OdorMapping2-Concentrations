"""Stage 5: Fig 7 layout, Fig 8 letters, S1 insets, S2/S3 split, S4 colour, variable renames."""
import re
from pathlib import Path

P = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\notebooks\figures.py")
text = P.read_bytes().decode("utf-8").replace("\r\n", "\n")
cells = [c.rstrip("\n") for c in re.split(r"(?m)^(?=# %%)", text) if c.startswith("# %%")]


def sec_list(cs):
    out, cur = [], 0
    for c in cs:
        m = re.search(r"(?m)^# ## (\d+)\.", c) if c.startswith("# %% [markdown]") else None
        if m:
            cur = int(m.group(1))
        out.append(cur)
    return out


def sub_cell(i, old, new):
    assert cells[i].count(old) == 1, (i, cells[i].count(old), old[:70])
    cells[i] = cells[i].replace(old, new)


secs = sec_list(cells)
idx = {s: [i for i, x in enumerate(secs) if x == s] for s in set(secs)}
code = lambda s: [i for i in idx[s] if not cells[i].startswith("# %% [markdown]")]

# ---------------- lowercase variable renames: section 8 first (fig7_ -> fig8_), then 7 (fig6_ -> fig7_) ----------
for i in idx[8]:
    cells[i] = re.sub(r"\bfig7_", "fig8_", cells[i])
for i in idx[7]:
    cells[i] = re.sub(r"\bfig6_", "fig7_", cells[i])

# ---------------- Fig 7 (section 7): B side by side + D --------------------------------------------------------
c1, c2 = code(7)
sub_cell(c1, 'print(f"6C: glomerulus {s.glom}', 'print(f"7B: glomerulus {s.glom}')
sub_cell(c1, 'print(f"6C: smallest nonzero', 'print(f"7B: smallest nonzero')
sub_cell(c1, 'print(f"6D: W is', 'print(f"7D: W is')
sub_cell(c1, 'fig7_pls_C_prediction_spectra.csv', 'fig7_pls_B_prediction_spectra.csv')
sub_cell(c1, "FIG7_GLOMERULI = (211, 97)      # C's two example glomeruli (dataset glomerulus ids)",
         "FIG7_GLOMERULI = (211, 97)      # B's two example glomeruli (dataset glomerulus ids)")
sub_cell(c1, "# 6C's own numbers", "# 7B's own numbers")
a = cells[c2].index("# One row, in inches: C")
b = cells[c2].index("# --- D: the full W, every weight")
new_layout = '''# One row, in inches: B (two example glomeruli, side by side) | D (W, square, + colourbar).
# Panels A (the X -> W -> Y schematic) and C (the three weight-matrix cartoons) of the paper's
# Figure 7 are drawn separately; the letters here are the paper's. The width is MAIN_TARGET_W
# exactly; the right margin is whatever is left, so changing any other width cannot silently
# change the page width.
FIG7_MARGIN_L, FIG7_MARGIN_TOP, FIG7_MARGIN_BOT = 0.58, 0.55, 0.60
FIG7_B_W, FIG7_B_H, FIG7_B_GAP = 1.40, 1.05, 0.50   # each spectrum, and the gap between them
FIG7_COL_GAP = 0.45                                  # between B and D
FIG7_D_SIDE = 1.75                                   # D is square
FIG7_CBAR_PAD, FIG7_CBAR_W = 0.30, 0.12
FIG7_MARGIN_R = (MAIN_TARGET_W - FIG7_MARGIN_L - 2 * FIG7_B_W - FIG7_B_GAP
                 - FIG7_COL_GAP - FIG7_D_SIDE - FIG7_CBAR_PAD - FIG7_CBAR_W)
assert FIG7_MARGIN_R > 0.30, f"Fig 7 row does not fit: right margin {FIG7_MARGIN_R:.2f} in"
FIG7_W = MAIN_TARGET_W
FIG7_H = FIG7_MARGIN_TOP + FIG7_D_SIDE + FIG7_MARGIN_BOT


def x6(inches):
    return inches / FIG7_W


def y6(inches):
    return inches / FIG7_H


fig = plt.figure(figsize=(FIG7_W, FIG7_H), dpi=FIG_DPI)
b_bottom = FIG7_MARGIN_BOT + (FIG7_D_SIDE - FIG7_B_H) / 2     # B centred on D's height

# --- B: two example glomeruli, experimental high up / predicted high down ----
for j, s in enumerate(fig7_specs):
    x0 = FIG7_MARGIN_L + j * (FIG7_B_W + FIG7_B_GAP)
    ax = fig.add_axes([x6(x0), y6(b_bottom), x6(FIG7_B_W), y6(FIG7_B_H)])
    P.similarity.prediction_spectrum(
        ax, s, limit=FIG7_LIMIT, strip_labels=True, annotate_r=False,
        observed_floor=fig7_obs_floor, **({} if j == 0 else {"ylabel": ""}))
    # r goes in the title, as in Fig 3's F, not inside the axes
    ax.set_title(f"glomerulus {s.glom}  (cluster {s.cluster})    r = {s.r:.2f}",
                 fontsize=HEATMAP_TITLE_FONTSIZE)
    if j == 0:
        hs.panel_label(ax, "B", x=-0.30, y=1.18)

'''
cells[c2] = cells[c2][:a] + new_layout + cells[c2][b:]
sub_cell(c2, "d_left = FIG7_MARGIN_L + FIG7_C_W + FIG7_COL_GAP",
         "d_left = FIG7_MARGIN_L + 2 * FIG7_B_W + FIG7_B_GAP + FIG7_COL_GAP")

# ---------------- Fig 8 (section 8): paper letters B (p-values) and C (effect size) ---------------------------
f8 = [i for i in code(8) if 'hs.panel_label(ax_p, "A"' in cells[i]]
assert len(f8) == 1
sub_cell(f8[0], 'hs.panel_label(ax_p, "A", x=-0.26, y=1.12)', 'hs.panel_label(ax_p, "B", x=-0.26, y=1.12)')
sub_cell(f8[0], 'hs.panel_label(ax_m, "B", x=-0.16, y=1.12)', 'hs.panel_label(ax_m, "C", x=-0.16, y=1.12)')
sub_cell(f8[0], "# One row: A and B, the cluster analysis.", "# One row: B and C (the paper's letters; its A is the procedure schematic, drawn elsewhere).")

# ---------------- S1 (section 9): no consensus insets ----------------------------------------------------------
s1 = code(9)[0]
sub_cell(s1, '''    inner = GridSpecFromSubplotSpec(len(consensus_ks), 1,
                                    subplot_spec=gs_bot[0, 5], hspace=0.30)
    for i, k in enumerate(consensus_ks):
        P.dimchoice.consensus(consensus[k], fig.add_subplot(inner[i, 0]))

''', '''    # The consensus-matrix insets that used to sit beside C are omitted (a reviewer
    # asked what they show, and the paper's Figure S1 can do without them).

''')
sub_cell(s1, "consensus shown at {consensus_ks} | panels", "consensus matrices not drawn | panels")

# ---------------- S2 / S3 split (section 10) -------------------------------------------------------------------
c_const, c_a, c_b = code(10)
assert "Figure S2B" in cells[c_b]
B = cells[c_b]
B = B.replace("# --- Figure S2B: where each mouse's members of each cluster sit -----",
              "# --- Figure S3: where each mouse's members of each cluster sit -----")
B = B.replace("SUPP_SUBJ_TOPO_SCALE = MAIN_TARGET_W / 16.0",
              "SUPP_TOPO_DPI = 200\nSUPP_SUBJ_TOPO_SCALE = MAIN_TARGET_W / 16.0")
B = B.replace("dpi=SUPP_SUBJ_DPI)\ngs_topo", "dpi=SUPP_TOPO_DPI)\ngs_topo")
B = B.replace('save(fig, "supp_per_subject_topography", dpi=SUPP_SUBJ_DPI)', 'save(fig, "supp_per_subject_topography", dpi=SUPP_TOPO_DPI)')
assert "SUPP_SUBJ_DPI" not in B, "S3 would still inherit a name from S2's constants cell"
B = re.sub(r"(?m)^# %%\n", "# %%\n# Own constants -- this section inherits no bare name from anywhere (S3).\n", B, count=1)
MARK = ("# > **[Caption not yet written for the merged figure — see "
        "notes/merge_decisions.md]**")
md_h = f'''# %% [markdown] tags=["human"]
# ## 11. Figure S3 — where each mouse's clusters sit on the bulb
#
{MARK}
#
# For each of the four mice (rows) and each of the seven glomerular groups (columns), where that
# mouse's members of the group lie on the dorsal bulb, over every glomerulus drawn as a hollow grey
# circle. Coloured by cluster.'''
md_l = f'''# %% [markdown] tags=["llm"]
# ## 11. Figure S3 — per-subject topography of the seven clusters
#
{MARK}
#
# Paper Fig S3. `per_subject_topography(d, ax, subject, cluster)` over subjects x clusters, drawn exactly
# like Fig 6 (marker sizes, outlines, hollow context); the fill is the cluster's own colour.'''
del cells[c_b]
pos = c_a + 1
cells[pos:pos] = [md_h, md_l, B]
# S2's heading text
for i in idx[10]:
    if cells[i].startswith("# %% [markdown]") and "## 10. Figure S2" in cells[i]:
        cells[i] = cells[i].replace("# ## 10. Figure S2 — the same structure within each mouse",
                                    "# ## 10. Figure S2 — the same structure within each mouse (heatmaps)")

# ---------------- S4 (section 12): pink-to-purple -------------------------------------------------------------
s4 = [i for i, c in enumerate(cells) if 'SUPP_PAR_C_CMAP = "Purples"' in c]
assert len(s4) == 1
sub_cell(s4[0], 'SUPP_PAR_C_CMAP = "Purples"          # hexbin colormap',
         'SUPP_PAR_C_CMAP = "RdPu"             # hexbin colormap: pink to purple, as in the paper')

P.write_bytes(("\n\n".join(cells) + "\n").encode("utf-8"))
print("stage 5 applied")
