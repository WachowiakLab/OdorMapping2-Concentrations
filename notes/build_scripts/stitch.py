"""Stitch final/notebooks/figures.py: main's notebook as the skeleton, code cells replaced per section
from final/notes/sections/*.py, section 10 (S3) removed, section 1 edited. Every edit asserts."""
import re
from pathlib import Path

ROOT = Path(r"C:\Users\u6068031\Documents\highconc-merge")
MAIN_FIG = ROOT / "highconc-figures-lab-main" / "notebooks" / "figures.py"
SEC = ROOT / "final" / "notes" / "sections"
OUT = ROOT / "final" / "notebooks" / "figures.py"
SCRATCH = Path(__file__).parent


def rd(p):
    return Path(p).read_bytes().decode("utf-8").replace("\r\n", "\n")


def cells_of(text):
    return [c.rstrip("\n") for c in re.split(r"(?m)^(?=# %%)", text) if c.startswith("# %%")]


def sub(text, old, new, label):
    assert text.count(old) == 1, (label, text.count(old), old[:70])
    return text.replace(old, new)


parts = cells_of(rd(MAIN_FIG))
assert len(parts) == 54, len(parts)

# ---------------- section 1: edit the setup cells (indices 4-8) ----------------
setup = "\n\n".join(parts[4:9])
setup = sub(setup, "from pathlib import Path\n", "import dataclasses\nfrom pathlib import Path\n", "dataclasses")
setup = sub(setup, 'mpl.use("Agg")\n', "", "Agg")
setup = sub(setup, "import os\n", "", "os")
setup = sub(setup, '''#: A second, SVG-only copy of every figure, kept separate from OUT so it can
#: be handed off (e.g. to a collaborator) without the PNG/PDF alongside it.
JOEL_FIGURES = ROOT / "joel_figures"
JOEL_FIGURES.mkdir(parents=True, exist_ok=True)

''', "", "JOEL const")
setup = sub(setup, '''you do not have to do anything to get them. It also drops a second,
    SVG-only copy into `joel_figures/` at the repo root. `export_figure()` is
    for when you want them somewhere else, in fewer formats, or at a different
    DPI.''', '''you do not have to do anything to get them. `export_figure()` is for when
    you want them somewhere else, in fewer formats, or at a different DPI.''', "JOEL doc")
setup = sub(setup, '''        _write(fig, stem, JOEL_FIGURES, ("svg",), dpi or FIG_DPI)
        print(f"saved: {stem}.{{png,pdf,svg}} (+ svg -> {JOEL_FIGURES.name}/)")''',
            '''        print(f"saved: {stem}.{{png,pdf,svg}}")''', "JOEL write")
assert "JOEL" not in setup and "joel" not in setup.lower()
setup_cells = cells_of(setup + "\n")
assert len(setup_cells) == 5, len(setup_cells)

# ---------------- section 2: main's LIVE Fig 2 (from its figures_human.py), tidied ----------------
hum = cells_of(rd(ROOT / "highconc-figures-lab-main" / "notebooks" / "figures_human.py"))
fig2 = [c for c in hum if not c.startswith("# %% [markdown]")
        and ("FIG2_ORDER" in c or 'save(fig, "fig2_response_matrix")' in c)]
assert len(fig2) == 2, len(fig2)          # the constants cell and the figure cell, in that order
assert "FIG2_ORDER" in fig2[0] and 'save(fig, "fig2_response_matrix")' in fig2[1]
fig2_text = "\n\n".join(fig2)
fig2_text = sub(fig2_text, "# cax = fig.add_axes([0.84, 0.15, 0.03, 0.7])\n", "", "dead cax")
fig2_text = sub(fig2_text, '''# export_figure("fig2_response_matrix")   # -> notebooks/outputs/fig2_response_matrix.{png,pdf,svg}
# export_figure("fig2_response_matrix", formats=("pdf",), outdir="~/Desktop")''', "", "fig2 export")
fig2_cells_new = cells_of(fig2_text.rstrip("\n") + "\n")

# ---------------- per-section replacement code cells ----------------
def sec(name):
    return cells_of(rd(SEC / name))


new = {
    3: sec("03_fig3_code.py"), 4: sec("04_fig4_code.py"), 5: sec("05_fig5_code.py"),
    6: sec("06_fig6_code.py"), 7: sec("07_fig7_code.py"), 8: sec("08_figS1_code.py"),
    9: sec("09_figS2_code.py"), 11: sec("11_figS4_code.py"),
}

# main's code-cell indices per section (from the cell listing)
old_code = {2: [11], 3: [14, 15, 16, 19, 20], 4: [23, 24], 5: [27, 28], 6: [31, 32, 33],
            7: [36, 37, 38, 39], 8: [42], 9: [45], 11: [51]}
new[2] = fig2_cells_new

# sanity: section 3 and 7 and 4, 5 must map one-to-one in count
for s in (3, 4, 5, 7, 8, 11):
    assert len(old_code[s]) == len(new[s]), (s, len(old_code[s]), len(new[s]))
assert len(new[6]) == 2 and len(old_code[6]) == 3          # Fig 6: third cell (parity check) dropped
assert len(new[9]) == 3 and len(old_code[9]) == 1          # S2: one cell becomes three
print("Fig 2 cells:", len(new[2]), "(old", len(old_code[2]), ")")

MARK = ("# > **[Caption not yet updated for the merged figure — see "
        "notes/merge_decisions.md]**\n#")
CHANGED_SECTIONS = {3, 4, 5, 6, 7, 8, 9, 11}


def mark_md(c):
    """Insert the caption-todo marker right after the heading line of a markdown cell."""
    lines = c.split("\n")
    # lines[0] is the "# %% [markdown] tags=..." header, lines[1] is the "# ## N. ..." heading
    for i, l in enumerate(lines[1:], start=1):
        if l.startswith("# ## ") or l.startswith("# ### "):
            lines.insert(i + 1, "#")
            lines.insert(i + 2, MARK)
            return "\n".join(lines)
    return c


out_cells = []
sec_of = {}
cur = 0
for i, c in enumerate(parts):
    md = c.startswith("# %% [markdown]")
    m = re.search(r"(?m)^# ## (\d+)\.", c) if md else None
    if m:
        cur = int(m.group(1))
    sec_of[i] = cur

emitted_sec = set()
for i, c in enumerate(parts):
    s = sec_of[i]
    md = c.startswith("# %% [markdown]")
    if 4 <= i <= 8:                      # section 1 code
        if i == 4:
            out_cells.extend(setup_cells)
        continue
    if s == 10:                          # S3 removed
        continue
    if md:
        out_cells.append(mark_md(c) if (s in CHANGED_SECTIONS and ("# ## " in c or "# ### " in c)) else c)
        continue
    # code cell of a replaced section: cell-for-cell, so captions between code cells keep their place
    if s in old_code and i in old_code[s]:
        k = old_code[s].index(i)
        emitted_sec.add(s)
        if k < len(new[s]):
            if k == len(old_code[s]) - 1:        # last old cell: also take any extra new cells
                out_cells.extend(new[s][k:])
            else:
                out_cells.append(new[s][k])
        continue                                 # (an old cell with no new counterpart is dropped)
    out_cells.append(c)                  # untouched (nothing falls here today)

text = "\n\n".join(out_cells) + "\n"
OUT.write_bytes(text.encode("utf-8"))
print("wrote", OUT, "|", len(out_cells), "cells,", text.count("\n"), "lines")
print("sections emitted:", sorted(emitted_sec))
