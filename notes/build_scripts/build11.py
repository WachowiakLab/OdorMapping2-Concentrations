"""S4: only panel C (parity of responding pairs), Purples, dotted line at the observed minimum."""
from pathlib import Path

FINAL = Path(r"C:\Users\u6068031\Documents\highconc-merge\final")


def edit_file(path, edits):
    t = path.read_bytes().decode("utf-8").replace("\r\n", "\n")
    for old, new in edits:
        assert t.count(old) == 1, (path.name, t.count(old), old[:70])
        t = t.replace(old, new)
    path.write_bytes(t.encode("utf-8"))


edit_file(FINAL / "highconc" / "panels" / "pls.py", [
    ('''                      cbar=False, title="responding pairs", title_fontsize=11,
                      note_fontsize=9):
    """Per-pair parity over responding pairs (observed > 0). Returns the hexbin.''',
     '''                      cbar=False, title="responding pairs", title_fontsize=11,
                      note_fontsize=9, cmap=None, floor=None, floor_color="0.35"):
    """Per-pair parity over responding pairs (observed > 0). Returns the hexbin.'''),
    ('''    changes how much ink is on the page.
    """
    obs, pred = _pls.responding_pairs(m)''',
     '''    changes how much ink is on the page.

    `cmap` is the hexbin colormap (default `highconc.style.SEQ`). `floor`, if
    given, draws a dotted horizontal line at that predicted value -- meant for
    `highconc.pls.observed_high_floor`, the smallest nonzero response ever
    measured, so a prediction under the line is one below anything observed.
    """
    obs, pred = _pls.responding_pairs(m)'''),
    ('''    hb = ax.hexbin(obs, pred, gridsize=gridsize, cmap=hs.SEQ, mincnt=1,''',
     '''    hb = ax.hexbin(obs, pred, gridsize=gridsize,
                   cmap=hs.SEQ if cmap is None else cmap, mincnt=1,'''),
    ('''        hs.despine(ax)

    if annotate:
        q2 = _pls.baseline_q2(m)''',
     '''        hs.despine(ax)

    if floor is not None:
        ax.axhline(floor, color=floor_color, ls=":", lw=1.0, zorder=3)

    if annotate:
        q2 = _pls.baseline_q2(m)'''),
])

cell = '''# %%
# Own constants — this section inherits no bare name from anywhere.
# Only panel C of the old three-panel figure remains: the parity of every
# responding glomerulus x odorant pair, observed vs held-out predicted.
SUPP_PAR_C_CMAP = "Purples"          # hexbin colormap
SUPP_PAR_W, SUPP_PAR_H, SUPP_PAR_DPI = SUPP_TARGET_W / 2, 2.4, 200

# The smallest nonzero observed high response in the whole dataset (log1p
# space, which is this panel's space): the dotted horizontal line. A predicted
# value under it is lower than any response the experiment ever measured.
supp_par_floor = pls.observed_high_floor(m, space="log")

fig, ax = plt.subplots(1, 1, figsize=(SUPP_PAR_W, SUPP_PAR_H),
                       dpi=SUPP_PAR_DPI, constrained_layout=True)
# Bare call = the pinned framing at its defaults, which is the draft's panel.
supp_par_hb = P.pls.responding_parity(
    ax, m, cbar=True, title="all responding pairs", title_fontsize=11.5,
    cmap=SUPP_PAR_C_CMAP, floor=supp_par_floor)

print(f"responding pairs   Q2 = {pls.baseline_q2(m):.3f}  "
      f"(n = {pls.responding_pairs(m)[0].size:,}); "
      "the draft page labels this same value R2")
print(f"dotted line at the smallest nonzero observed response: "
      f"{supp_par_floor:.4f} (log1p dF/F)")

save(fig, "supp_parity", dpi=SUPP_PAR_DPI)
plt.show()
# export_figure("supp_parity")   # -> notebooks/outputs/supp_parity.{png,pdf,svg}
'''
out = FINAL / "notes" / "sections" / "11_figS4_code.py"
out.write_text(cell, encoding="utf-8", newline="\n")
print("wrote", out, len(cell.splitlines()), "lines")
