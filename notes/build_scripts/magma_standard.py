"""Make magma the house sequential standard across style/ (text sources). Every edit asserts."""
from pathlib import Path

S = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\style")


def rd(p):
    return p.read_bytes().decode("utf-8").replace("\r\n", "\n")


def wr(p, t):
    p.write_bytes(t.encode("utf-8"))
    print("updated", p.relative_to(S))


def edit(p, pairs):
    t = rd(p)
    for o, n in pairs:
        assert t.count(o) == 1, (p.name, t.count(o), o[:80])
        t = t.replace(o, n)
    wr(p, t)


# ---------------- STYLE_GUIDE.md ----------------
p = S / "STYLE_GUIDE.md"
t = rd(p)
# (a) remove my status note
a = t.index("> **Status note (2026-10).**")
b = t.index("> is an open decision.\n") + len("> is an open decision.\n\n")
t = t[:a] + t[b:]
# (b) scope paragraph
a = t.index("**Scope.** These conventions were derived")
b = t.index("keep the acquisition colormap — see §10.\n") + len("keep the acquisition colormap — see §10.\n")
t = t[:a] + (
    "**Scope.** These conventions were derived from the analysis figures in\n"
    "`figures/` (the NMF/PLS work — Fig 4 and supplements) and the manuscript figure\n"
    "set. **The standard sequential map is magma** (§2.1): black at zero, pale yellow\n"
    "at the ceiling, which reads on the dark ground the paper's response matrices\n"
    "use. The earlier drafts of Figs 1–3, 5 and 6 used a black-background hot\n"
    "colormap for ΔF; magma replaces it with a perceptually even, colourblind-safe\n"
    "map that sits on the same ground. (An earlier version of this guide set OrRd on\n"
    "white as the standard; the repo's `highconc.style.SEQ` and the paper are magma\n"
    "now.) One documented exception: raw photometric bulb activity maps keep the\n"
    "acquisition colormap — see §10.\n") + t[b:]
wr(p, t)
edit(p, [
    ("| `seq_orrd_256.csv` | The sequential map, 256 steps.", "| `seq_magma_256.csv` | The sequential map (magma), 256 steps."),
    ("3. **Sequential data → OrRd. Signed data", "3. **Sequential data → magma. Signed data"),
    ("- Glomeruli as filled circles on a fixed spatial frame, OrRd, generous marker", "- Glomeruli as filled circles on a fixed spatial frame, magma, generous marker"),
    ("a *measurement*: it belongs at the pale end of the sequential ramp.", "a *measurement*: it belongs at the floor of the sequential ramp (black in magma)."),
    ("rather than OrRd. They are photometric images", "rather than magma. They are photometric images"),
])
t = rd(p)
a = t.index("**OrRd**, full range.")
b = t.index('brewer.pal(9, "OrRd")`).') + len('brewer.pal(9, "OrRd")`).')
t = t[:a] + (
    "**magma** (`plt.cm.magma`), full range. Used for ΔF/F, |values|, densities,\n"
    "response matrices — anything that runs from \"nothing\" to \"a lot\" with no\n"
    "meaningful zero-crossing.\n\n"
    "Black at the floor, pale yellow at the ceiling, through purple and orange-pink.\n"
    "Because zero is black, **missing data must be drawn explicitly** in `gainsboro`,\n"
    "never left to the axes background (§2.5). Load `palette/seq_magma_256.csv` if\n"
    "your tool doesn't ship magma (R has it via `viridisLite::magma(256)`; Fiji ships\n"
    "the `mpl-magma` LUT).\n\n"
    "One panel keeps a different ramp on purpose: Figure 8's p-value matrix is\n"
    "OrRd, passed to that panel explicitly (`FIG8_P_CMAP` in the notebook). OrRd on\n"
    "white, the earlier standard, is retired for everything else.") + t[b:]
wr(p, t)
# §10 bullet about the drafts
t = rd(p)
a = t.index("- **The Fig 1–3/5/6 drafts are on a different colormap.**")
b = t.index("- **Raw bulb activity maps are a deliberate exception.**")
t = t[:a] + (
    "- **The older Fig 1–3/5/6 drafts used black-background hot** (black → red →\n"
    "  yellow) for ΔF heatmaps, ranked response maps and correlation matrices. They\n"
    "  are on magma now, which sits on the same dark ground. Because zero is black,\n"
    "  the zero/missing distinction of §2.5 has to be made explicitly — missing\n"
    "  cells in opaque `gainsboro` — where the old black ground hid both.\n") + t[b:]
wr(p, t)

# ---------------- CHECKLIST.md ----------------
edit(S / "CHECKLIST.md", [
    ("Sequential quantity → **OrRd**;", "Sequential quantity → **magma**;"),
    ("**Zeros are pale, not grey**", "**Zeros are the ramp's floor (black), not grey**"),
])

# ---------------- style/README.md ----------------
p = S / "README.md"
t = rd(p)
a = t.index("> **Status note (2026-10).**")
b = t.index("> is an open decision.\n") + len("> is an open decision.\n\n")
t = t[:a] + t[b:]
wr(p, t)
edit(p, [("| `palette/seq_orrd_256.csv` | Sequential map (OrRd), 256 steps", "| `palette/seq_magma_256.csv` | Sequential map (magma), 256 steps")])

# ---------------- generate_assets.py ----------------
edit(S / "generate_assets.py", [
    ("  seq_orrd_256.csv              256-step sequential map (activation)", "  seq_magma_256.csv             256-step sequential map (activation)"),
    ("SEQ = plt.cm.OrRd", "SEQ = plt.cm.magma"),
    ('"Sequential — OrRd — activation, |values|, densities"', '"Sequential — magma — activation, |values|, densities"'),
    ('write_cmap(SEQ, "seq_orrd_256")', 'write_cmap(SEQ, "seq_magma_256")'),
])

# ---------------- for-collaborators ----------------
edit(S / "for-collaborators" / "build_prompt.py", [
    ("SEQ = plt.cm.OrRd                                     # sequential: magnitudes", "SEQ = plt.cm.magma                                    # sequential: magnitudes"),
    ("→ `SEQ` (OrRd).", "→ `SEQ` (magma)."),
    ("pale end of the ramp; `GAINSBORO` grey means NOT MEASURED.", "floor of the ramp (black in magma); `GAINSBORO` grey means NOT MEASURED."),
])
p = S / "for-collaborators" / "claude-skill" / "highconc-figure" / "SKILL.md"
t = rd(p)
a = t.index("The Fig 1–3/5/6 drafts use a black-background hot colormap")
b = t.index("- Raw photometric bulb activity maps are the documented exception")
t = t[:a] + (
    "The older Fig 1–3/5/6 drafts used a black-background hot colormap (black → red →\n"
    "yellow). **The standard is now magma**, which sits on the same dark ground. When\n"
    "converting one:\n\n"
    "- Rebuild from the data, not the image.\n"
    "- The black ground used to hide both zeros and gaps. With magma, zero is the\n"
    "  ramp's black floor (`vmin=0`) and a missing cell must be drawn explicitly in\n"
    "  `gainsboro`. Ask if the data does not distinguish them.\n") + t[b:]
wr(p, t)
edit(p, [("pale end of the ramp; grey means *not measured*.", "floor of the ramp (black in magma); grey means *not measured*.")])
print("done")
