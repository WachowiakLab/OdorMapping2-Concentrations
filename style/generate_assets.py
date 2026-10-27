"""Generate the tool-agnostic style assets for collaborators.

Emits, into palette/:
  highconc_palette.csv          named colors -> hex / RGB255 / RGB01
  seq_magma_256.csv             256-step sequential map (activation)
  div_teal_gainsboro_orrd_256.csv   256-step house diverging map (weights)
  highconc_swatches.ase         Adobe swatch exchange (Illustrator/Photoshop)
  swatches.png/.pdf             visual reference sheet

Run:  python generate_assets.py
"""
import csv
import struct
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb

mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

OUT = Path(__file__).parent / "palette"
OUT.mkdir(exist_ok=True)

# --- the house palette, transcribed from figures/*/generate.py --------------
DIV_ANCHORS = ["#01665e", "#5ab4ac", "gainsboro", "#fc8d59", "#b30000"]
DIV = LinearSegmentedColormap.from_list("teal_gainsboro_orrd", DIV_ANCHORS)
SEQ = plt.cm.magma

NAMED = [
    # name,                 hex,        role
    ("div_teal_dark",       "#01665e", "Diverging map, negative extreme (mixing / loss)"),
    ("div_teal_light",      "#5ab4ac", "Diverging map, negative midtone"),
    ("neutral_gainsboro",   "#dcdcdc", "Diverging midpoint; also the NOT-MEASURED grey (zero is not missing)"),
    ("div_orange",          "#fc8d59", "Diverging map, positive midtone"),
    ("div_red_dark",        "#b30000", "Diverging map, positive extreme; also POINT_COLOR"),
    ("cond_most_sensitive", "#e69f00", "Condition 'most-sensitive' (Okabe-Ito amber)"),
    ("cond_low",            "#0072b2", "Condition 'low' / 1x (same value as accent_blue)"),
    ("cond_high",           "#d55e00", "Condition 'high' / 25x (same value as accent_vermilion)"),
    ("accent_blue",         "#0072b2", "Okabe-Ito blue - Q2 / primary metric series"),
    ("accent_vermilion",    "#d55e00", "Okabe-Ito vermilion - overfit / canonical comparison"),
    ("accent_green",        "#009e73", "Okabe-Ito green - 'wins' series"),
    ("accent_grey",         "#888888", "'loses' / neutral series"),
    ("shuffled_grey",       "#bbbbbb", "Shuffled / null distributions"),
    ("pca_grey",            "#777777", "PCA reference series (vs NMF)"),
    ("nmf_purple",          "#9b5bd0", "NMF series"),
    ("identity_grey",       "#737373", "Identity / unity lines (matplotlib '0.45')"),
]

CLUSTERS = [(f"cluster_{c}", to_hex(plt.cm.tab10(i))) for i, c in enumerate("ABCDEFG")]


def write_named():
    p = OUT / "highconc_palette.csv"
    with p.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "hex", "R255", "G255", "B255", "R01", "G01", "B01", "role"])
        for name, hx, role in NAMED:
            r, g, b = to_rgb(hx)
            w.writerow([name, hx, round(r*255), round(g*255), round(b*255),
                        f"{r:.4f}", f"{g:.4f}", f"{b:.4f}", role])
        for name, hx in CLUSTERS:
            r, g, b = to_rgb(hx)
            w.writerow([name, hx, round(r*255), round(g*255), round(b*255),
                        f"{r:.4f}", f"{g:.4f}", f"{b:.4f}",
                        "Glomerulus/odorant cluster identity (matplotlib tab10)"])
    print("wrote", p.name)


def write_cmap(cmap, stem):
    p = OUT / f"{stem}.csv"
    with p.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["i", "hex", "R255", "G255", "B255", "R01", "G01", "B01"])
        for i in range(256):
            r, g, b, _ = cmap(i / 255)
            w.writerow([i, to_hex((r, g, b)), round(r*255), round(g*255), round(b*255),
                        f"{r:.6f}", f"{g:.6f}", f"{b:.6f}"])
    print("wrote", p.name)


def _ase_block(name, rgb):
    nb = (name + "\0").encode("utf-16-be")
    body = struct.pack(">H", len(name) + 1) + nb + b"RGB " + \
        struct.pack(">fff", *rgb) + struct.pack(">H", 2)   # 2 = normal colour
    return struct.pack(">HI", 0x0001, len(body)) + body


def write_ase():
    entries = [(n, to_rgb(h)) for n, h, _ in NAMED] + [(n, to_rgb(h)) for n, h in CLUSTERS]
    blocks = b"".join(_ase_block(n, rgb) for n, rgb in entries)
    p = OUT / "highconc_swatches.ase"
    p.write_bytes(b"ASEF" + struct.pack(">HHI", 1, 0, len(entries)) + blocks)
    print("wrote", p.name, f"({len(entries)} swatches)")


def sheet():
    n_named = len(NAMED)
    ncol = 7
    nrow = -(-n_named // ncol)                      # ceil
    fig = plt.figure(figsize=(9.5, 2.0 + 1.05 * nrow + 3.0))
    gs = fig.add_gridspec(5, 1, height_ratios=[1.05 * nrow, 0.62, 0.62, 1.15, 1.15],
                          hspace=0.95)

    def swatch_row(ax, items, title, label_fs=6.5, y0=0.30, h=0.42, w=0.09):
        ax.set_axis_off()
        ax.set_title(title, fontsize=11, loc="left", fontweight="bold")
        step = 1.0 / max(len(items), 1)
        for i, (nm, hx) in enumerate(items):
            ax.add_patch(plt.Rectangle((i * step, y0), w, h, facecolor=hx,
                                       edgecolor="0.4", lw=.5,
                                       transform=ax.transAxes, clip_on=False))
            ax.text(i * step + w / 2, y0 - 0.14, nm, fontsize=label_fs, ha="center",
                    va="top", transform=ax.transAxes)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # --- named colours, wrapped over as many rows as it takes ---------------
    ax = fig.add_subplot(gs[0]); ax.set_axis_off()
    ax.set_title("Named colours", fontsize=11, loc="left", fontweight="bold")
    row_h = 1.0 / nrow
    for i, (name, hx, _) in enumerate(NAMED):
        col, row = i % ncol, i // ncol
        y0 = 1.0 - (row + 1) * row_h + row_h * 0.42
        ax.add_patch(plt.Rectangle((col * 0.143, y0), 0.105, row_h * 0.36,
                                   facecolor=hx, edgecolor="0.4", lw=.5,
                                   transform=ax.transAxes, clip_on=False))
        ax.text(col * 0.143, y0 - row_h * 0.06, f"{name}\n{hx}", fontsize=6.5,
                va="top", transform=ax.transAxes)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # --- clusters -----------------------------------------------------------
    swatch_row(fig.add_subplot(gs[1]),
               [(c[-1], hx) for c, hx in CLUSTERS],
               "Cluster identity — tab10, clusters A–G (fixed across the paper)",
               label_fs=8)

    # --- conditions ---------------------------------------------------------
    swatch_row(fig.add_subplot(gs[2]),
               [("'most-sensitive'", "#e69f00"), ("'low'  (1x)", "#0072b2"),
                ("'high'  (25x)", "#d55e00")],
               "Concentration conditions — the paper's primary categorical scheme",
               label_fs=7.5, w=0.13)

    # --- colormaps ----------------------------------------------------------
    grad = np.linspace(0, 1, 256).reshape(1, -1)
    for k, (cm, lab) in enumerate([
            (SEQ, "Sequential — magma — activation, |values|, densities"),
            (DIV, "Diverging — teal_gainsboro_orrd — signed weights, differences")]):
        ax = fig.add_subplot(gs[3 + k])
        ax.imshow(grad, aspect="auto", cmap=cm, interpolation="nearest")
        ax.set_title(lab, fontsize=11, loc="left", fontweight="bold")
        ax.set_yticks([])
        ax.set_xticks([0, 127, 255])
        ax.set_xticklabels(["0.0", "0.5", "1.0"] if k == 0 else ["−max", "0", "+max"],
                           fontsize=8)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        if k == 1:
            for x, hx in zip(np.linspace(0, 255, 5), DIV_ANCHORS):
                ax.annotate(hx, (x, 0.5), xytext=(0, -26), textcoords="offset points",
                            ha="center", fontsize=6.5)

    fig.suptitle("highconcmss figure style — palette reference", fontsize=13,
                 fontweight="bold", x=0.02, ha="left", y=0.995)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"swatches.{ext}", dpi=200, bbox_inches="tight",
                    facecolor="white" if ext == "png" else "none",
                    transparent=(ext != "png"))
    plt.close(fig)
    print("wrote swatches.png / swatches.pdf")


if __name__ == "__main__":
    write_named()
    write_cmap(SEQ, "seq_magma_256")
    write_cmap(DIV, "div_teal_gainsboro_orrd_256")
    write_ase()
    sheet()
