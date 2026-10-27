"""Fig 5: add `color` (solid membership marker) to topography(), and assemble the §5 cell."""
from pathlib import Path

HERE = Path(__file__).parent
PANEL = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\highconc\panels\clustering.py")
OUT = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\notes\sections\05_fig5_code.py")


def patch(path, edits, crlf_ok=True):
    raw = path.read_bytes()
    crlf = b"\r\n" in raw
    t = raw.decode("utf-8").replace("\r\n", "\n")
    for old, new in edits:
        assert t.count(old) == 1, (t.count(old), old[:70])
        t = t.replace(old, new)
    path.write_bytes((t.replace("\n", "\r\n") if crlf else t).encode("utf-8"))
    return t


# ---- panel: topography(color=...) ------------------------------------------
patch(PANEL, [
    ('''    context: bool = True,
    nonresponder_marker="context",
    context_color: str | None = None,''',
     '''    context: bool = True,
    color: str | None = None,
    nonresponder_marker="context",
    context_color: str | None = None,'''),
    ('''    context : bool
        Draw the non-cluster glomeruli. Off gives the cluster alone.
    nonresponder_marker : str or dict''',
     '''    context : bool
        Draw the non-cluster glomeruli. Off gives the cluster alone.
    color : str, optional
        Draw every cluster member as one solid colour, so the map shows
        *membership* and nothing else -- no response value is encoded, and
        `conc`, `vmin`, `vmax` and `cmap` are ignored. Off by default (points
        coloured by mean dF/F as above). The returned collection then has no
        colour scale, so there is nothing to hand to `fig.colorbar`.
    nonresponder_marker : str or dict'''),
    ('''    points = ax.scatter(x[mask], y[mask], c=values[mask], cmap=cmap,
                        vmin=vmin, vmax=vmax, s=size, edgecolors=edgecolor,
                        linewidths=linewidth, zorder=2)''',
     '''    if color is not None:
        points = ax.scatter(x[mask], y[mask], color=color, s=size,
                            edgecolors=edgecolor, linewidths=linewidth,
                            zorder=2)
    else:
        points = ax.scatter(x[mask], y[mask], c=values[mask], cmap=cmap,
                            vmin=vmin, vmax=vmax, s=size, edgecolors=edgecolor,
                            linewidths=linewidth, zorder=2)'''),
])

# ---- notebook cell ---------------------------------------------------------
src = (HERE / "sec" / "05_main_fig_code.py").read_text(encoding="utf-8")


def sub(old, new):
    global src
    assert src.count(old) == 1, (src.count(old), old[:70])
    src = src.replace(old, new)


sub('''# %%
FIG5_CONC = "low"
FIG5_NONRESPONDER''',
    '''# %%
# Membership only: every cluster member is one solid colour, and no response
# value is drawn -- so there is no concentration to choose and no colourbar.
FIG5_MEMBER_COLOR = "red"
FIG5_NONRESPONDER''')
sub('''# Restated here for the same reason §4 restates its own: one section, one name.
FIG5_TOPO_VMIN = P.clustering.TOPO_VMIN
FIG5_TOPO_VMAX = P.clustering.TOPO_VMAX

''', '')
sub('left=0.05, right=0.90, top=0.82, bottom=0.06)',
    'left=0.05, right=0.97, top=0.82, bottom=0.06)')
sub('''    pts5 = P.clustering.topography(
        d, ax, k, conc=FIG5_CONC, vmin=FIG5_TOPO_VMIN, vmax=FIG5_TOPO_VMAX,
        size=FIG5_CLUSTER_SIZE, edgecolor=FIG5_CLUSTER_EDGE,
        linewidth=FIG5_CLUSTER_LW, nonresponder_marker=FIG5_NONRESPONDER,
        context_size=FIG5_CONTEXT_SIZE, context_linewidth=FIG5_CONTEXT_LW,
        title=f"Cluster {k}",
        ylabel=({"low": "low", "high": "high"}[FIG5_CONC] if k == 1 else None))''',
    '''    P.clustering.topography(
        d, ax, k, color=FIG5_MEMBER_COLOR,
        size=FIG5_CLUSTER_SIZE, edgecolor=FIG5_CLUSTER_EDGE,
        linewidth=FIG5_CLUSTER_LW, nonresponder_marker=FIG5_NONRESPONDER,
        context_size=FIG5_CONTEXT_SIZE, context_linewidth=FIG5_CONTEXT_LW,
        title=f"Cluster {k}")''')
sub('''cax5 = fig.add_axes([0.915, 0.18, 0.008, 0.56])
cb5 = fig.colorbar(pts5, cax=cax5)
cb5.set_label("mean ΔF/F", fontsize=hs.FS_LABEL)
cb5.ax.tick_params(labelsize=hs.FS_TICK)

''', '')
sub('''# export_figure("fig5_spatial_by_cluster")   # -> notebooks/outputs/fig5_spatial_by_cluster.{png,pdf,svg}
# export_figure("fig5_spatial_by_cluster", formats=("pdf",), outdir="~/Desktop")''',
    '''# export_figure("fig5_spatial_by_cluster")   # -> notebooks/outputs/fig5_spatial_by_cluster.{png,pdf,svg}''')

OUT.write_text(src, encoding="utf-8", newline="\n")
print("wrote", OUT, len(src.splitlines()), "lines")
