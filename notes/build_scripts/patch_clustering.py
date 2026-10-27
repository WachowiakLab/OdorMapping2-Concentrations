"""final/highconc/panels/clustering.py = main's file + col_normalize/raw_floor on per_subject_heatmap."""
from pathlib import Path

SRC = Path(r"C:\Users\u6068031\Documents\highconc-merge\highconc-figures-lab-main\highconc\panels\clustering.py")
DST = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\highconc\panels\clustering.py")
BASE = Path(r"C:\Users\u6068031\Documents\highconc-merge\github\original_codebase\highconc\panels\clustering.py")

raw = SRC.read_bytes()
crlf = b"\r\n" in raw
text = raw.decode("utf-8").replace("\r\n", "\n")
base_crlf = b"\r\n" in BASE.read_bytes()


def sub(old, new):
    global text
    assert text.count(old) == 1, (text.count(old), old[:70])
    text = text.replace(old, new)


# signature (per_subject_heatmap only: anchored on its own unique lines)
sub('''    vmin: float = VMIN,
    vmax: float = VMAX,
    cmap: str = CMAP_NAME,
    empty_color: str = "gainsboro",
    boundaries: bool = True,
    boundary_color: str = "#1f3fff",
    boundary_lw: float = 0.7,
    boundary_alpha: float = 0.45,
    cluster_labels: bool = True,
    label_fontsize: float = 7.5,
    title: str | bool = True,''',
    '''    vmin: float = VMIN,
    vmax: float = VMAX,
    cmap: str = CMAP_NAME,
    empty_color: str = "gainsboro",
    col_normalize: bool = False,
    raw_floor: float | None = None,
    boundaries: bool = True,
    boundary_color: str = "#1f3fff",
    boundary_lw: float = 0.7,
    boundary_alpha: float = 0.45,
    cluster_labels: bool = True,
    label_fontsize: float = 7.5,
    title: str | bool = True,''')

# docstring
sub('''    cluster_labels : bool
        Write each block's cluster number above it, in that cluster's colour.
        With `cluster` set there is only one block, so this writes one number.
''',
    '''    col_normalize, raw_floor
        Same meaning as in `biclustered`: peak-normalise each glomerulus to
        [0, 1] by its own maximum over all odorants, on the linear raw dF/F
        values, and draw on a LINEAR scale (pass ``vmin=0.0, vmax=1.0``);
        `raw_floor` greys out cells whose raw value was below it. Because this
        panel always carries all 59 odorant rows, a glomerulus normalises to
        exactly the same values here as in `biclustered`. Off by default (raw
        dF/F on a log scale), so the other figures using this panel are
        unchanged.
    cluster_labels : bool
        Write each block's cluster number above it, in that cluster's colour.
        With `cluster` set there is only one block, so this writes one number.
''')

# logic
sub('''    heat_cmap = plt.get_cmap(cmap).copy()
    heat_cmap.set_under(empty_color)
    heat_cmap.set_bad(empty_color)      # LogNorm masks zeros
    im = ax.imshow(data, aspect="auto", cmap=heat_cmap,
                   norm=LogNorm(vmin=vmin, vmax=vmax),
                   interpolation="nearest")''',
    '''    if col_normalize:
        peak = data.max(0, keepdims=True)
        peak = np.where(peak > 0, peak, np.inf)
        normed = data / peak
        if raw_floor is not None:
            normed = np.where(data >= raw_floor, normed, np.nan)
        data = normed

    heat_cmap = plt.get_cmap(cmap).copy()
    heat_cmap.set_under(empty_color)
    heat_cmap.set_bad(empty_color)      # LogNorm masks zeros; raw_floor uses this too
    norm = (plt.Normalize(vmin=vmin, vmax=vmax) if col_normalize
            else LogNorm(vmin=vmin, vmax=vmax))
    im = ax.imshow(data, aspect="auto", cmap=heat_cmap, norm=norm,
                   interpolation="nearest")''')

out = text.replace("\n", "\r\n") if base_crlf else text
DST.write_bytes(out.encode("utf-8"))
print("main crlf:", crlf, "| base crlf:", base_crlf, "| wrote", DST, len(text.splitlines()), "lines")
