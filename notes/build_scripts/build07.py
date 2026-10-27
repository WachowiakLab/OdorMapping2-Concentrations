"""Fig 7: only A and B; green significance boxes; values only on significant cells; new axis labels."""
from pathlib import Path

ROOT = Path(r"C:\Users\u6068031\Documents\highconc-merge")
FINAL = ROOT / "final"
HERE = Path(__file__).parent


def rd(p):
    return p.read_bytes().decode("utf-8").replace("\r\n", "\n")


def edit(text, edits):
    for old, new in edits:
        assert text.count(old) == 1, (text.count(old), old[:70])
        text = text.replace(old, new)
    return text


# ---- panel: box_color + text_significant_only (baseline panel file as the base) ----
P_ = FINAL / "highconc" / "panels" / "transfer.py"
t = edit(rd(P_), [
    ('''              box_lw=2.0, fontsize=None, labels=None, rotation=0,
              xlabel=XLABEL, ylabel=YLABEL, diagonal=hs.GAINSBORO):''',
     '''              box_lw=2.0, fontsize=None, labels=None, rotation=0,
              xlabel=XLABEL, ylabel=YLABEL, diagonal=hs.GAINSBORO,
              box_color=hs.BOUNDARY, text_significant_only=False):'''),
    ('''    The diagonal is drawn as an explicit grey square rather than left blank: a''',
     '''    `box_color` is the colour of the significance box. `text_significant_only`
    prints a cell's number only where the pair is significant (BH q < `alpha`),
    leaving every other cell blank -- the panel then reads as the claim alone.

    The diagonal is drawn as an explicit grey square rather than left blank: a'''),
    ('''            if text is not None:
                ax.text(T, F, text(F, T), ha="center", va="center",''',
     '''            if text is not None and (sig[F, T] or not text_significant_only):
                ax.text(T, F, text(F, T), ha="center", va="center",'''),
    ('''                                       edgecolor=hs.BOUNDARY, lw=box_lw,
                                       zorder=5))''',
     '''                                       edgecolor=box_color, lw=box_lw,
                                       zorder=5))'''),
    ('''    The black box, not the colour, is the claim: it marks BH q < `alpha` across''',
     '''    The box, not the colour, is the claim: it marks BH q < `alpha` across'''),
    ('''def significance_legend(ax, t, *, alpha=0.05, loc="upper right"):
    """One proxy handle explaining the black box. Cheap, and it removes the
    only piece of the panel a reader has to be told about in prose."""
    handle = Rectangle((0, 0), 1, 1, fill=False, edgecolor=hs.BOUNDARY, lw=2.0)''',
     '''def significance_legend(ax, t, *, alpha=0.05, loc="upper right",
                        box_color=hs.BOUNDARY):
    """One proxy handle explaining the significance box. Cheap, and it removes
    the only piece of the panel a reader has to be told about in prose."""
    handle = Rectangle((0, 0), 1, 1, fill=False, edgecolor=box_color, lw=2.0)'''),
])
P_.write_bytes(t.encode("utf-8"))

# ---- notebook cell ----------------------------------------------------------
src = (HERE / "sec" / "07_main_fig_code.py").read_text(encoding="utf-8")


def sub(old, new):
    global src
    assert src.count(old) == 1, (src.count(old), old[:70])
    src = src.replace(old, new)


sub('FIG7_MIN_CLASS = 15          # glomeruli a chemical class needs to be tested\n',
    'FIG7_MIN_CLASS = 15          # glomeruli a chemical class needs to be tested\n'
    'FIG7_BOX_COLOR = "#2b9743"   # the box on significant cells (green)\n'
    'FIG7_XLABEL = "cluster (predicted at high)"\n'
    'FIG7_YLABEL = "added non-source cluster\\nprediction improvement"\n')

a = src.index("# Scale the whole figure by one factor, solved so the rendered width is")
b = src.index("fig.text(0.5, 0.98")
new_block = '''# One row: A and B, the cluster analysis. The chemical-class control (what
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
                         xlabel=FIG7_XLABEL, ylabel=FIG7_YLABEL,
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

'''
src = src[:a] + new_block + src[b:]
sub('''# export_figure("fig7_motif_transfer")   # -> notebooks/outputs/fig7_motif_transfer.{png,pdf,svg}
# export_figure("fig7_motif_transfer", formats=("pdf",), outdir="~/Desktop")''',
    '''# export_figure("fig7_motif_transfer")   # -> notebooks/outputs/fig7_motif_transfer.{png,pdf,svg}''')
out = FINAL / "notes" / "sections" / "07_fig7_code.py"
out.write_text(src, encoding="utf-8", newline="\n")
print("wrote panel + ", out, len(src.splitlines()), "lines")
