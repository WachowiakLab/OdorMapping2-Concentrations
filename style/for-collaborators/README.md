# Start here

You have data — probably a CSV — and you need a figure that matches the rest of
the manuscript. Pick the row that describes you.

| If you… | Use | Setup |
|---|---|---|
| use **any AI chatbot** (ChatGPT, Claude, Gemini, Copilot) | **`FIGURE_PROMPT.md`** | none — paste it in with your CSV |
| write **Python** yourself | **`highconc_style.py`** | drop the file next to your script |
| use **Claude Code** | **the skill** | `./install-claude-skill.sh` |
| use **MATLAB, R, Prism, Igor** | `../STYLE_GUIDE.md` + `../palette/*.csv` | none |

All four produce the same thing, because all four come from the same source:
`highconc_style.py` is the reference implementation, and the prompt is generated
from it. They cannot drift apart.

---

## The chatbot route (no setup at all)

1. Open `FIGURE_PROMPT.md`, copy the block between the `````text` fences.
2. Paste it into ChatGPT / Claude / whatever you use.
3. Attach or paste your CSV and say what you want plotted.

Two things to tell it up front, because it cannot see them and will otherwise
guess wrong:

- **What the columns mean** — units, and whether values are raw or transformed.
- **Whether this panel will sit next to another one showing the same quantity.**
  If it will, they must share a colour scale, and it needs the limits.

It will look at your data, tell you which figure type fits, write a script,
export PNG + PDF + SVG, and check that the text survived as editable text.

## The Python route

```python
import highconc_style as hs
hs.apply()

fig, ax = plt.subplots(figsize=(hs.WIDTH_1COL, 2.8))
hs.parity(ax, df.observed, df.predicted)     # square, unity line, in-panel r
hs.panel_label(ax, "A")
hs.save(fig, "my_panel")                     # -> .png / .pdf / .svg
hs.check_export("my_panel.svg")              # -> PASS / FAIL
```

Helpers: `parity`, `heatmap`, `topomap`, `vector_matrix`, `ecdf`, `box`,
`rank_align`, `class_strip`, `despine`, `panel_label`, `identity_line`,
`sequential_cmap`, `save`, `check_export`. Palette: `SEQ`, `DIV`,
`CLUSTER_COLORS`, `CONDITION_COLORS`, `condition_color`, `GAINSBORO`,
`IDENTITY`, `C_PRIMARY`, `C_OVERFIT`, `C_WIN`, `C_LOSE`, `C_SHUFFLED`.

Anything split by concentration goes through the condition colours:

```python
ax.plot(x, y, color=hs.condition_color("25x"))   # accepts 'high', '25x', 'low', '1x', ...
hs.ecdf(ax, {"low": low_vals, "high": high_vals})   # ordered + coloured for you
hs.box(ax, {"low": low_c, "high": high_c}, log=True)
```

Check it works in your environment — this renders a demo of every helper and
verifies the export:

```bash
python highconc_style.py
```

If it prints `resolved sans-serif font: DejaVu Sans` rather than `Arial`, install
Arial or expect your text to differ from everyone else's.

## The Claude Code route

```bash
./install-claude-skill.sh
```

Then restart Claude Code and just ask — "make a figure from this CSV for the
paper". The installer copies the style guide, the checklist, the palette files
and the module in beside the skill, so it works from any directory.

Project-local instead of global: `./install-claude-skill.sh .claude/skills`

## The other-software route

`../STYLE_GUIDE.md` states the target in tool-neutral terms; `../palette/*.csv`
gives every colour as hex and RGB, and both colormaps as 256-row tables you can
load straight into MATLAB, R, or ImageJ. `../CHECKLIST.md` is the pre-handover
tick list. §9 of `../style_reference.ipynb` has a matplotlib → MATLAB → R
equivalence table.

---

## Whatever route you take

**Send three things, not one:** the PDF or SVG, the script or steps that made
it, and the numbers behind it. A figure can be rebuilt from data when a reviewer
asks for a change. It cannot be rebuilt from a PDF.

**Do not hand-edit the exported file** to fix something. If a panel is wrong,
fix it upstream and re-export — otherwise the fix disappears the next time the
figure is regenerated.

**If the export check says FAIL**, stop and fix it. Text that arrives in
Illustrator as outlines or pixels cannot be corrected at assembly, and that is
the single most common reason a panel has to be redone.

## Files here

| | |
|---|---|
| `FIGURE_PROMPT.md` | The paste-anywhere prompt. Generated — do not hand-edit. |
| `highconc_style.py` | The reference implementation. Source of truth for every value. |
| `claude-skill/` | The Claude Code skill definition. |
| `install-claude-skill.sh` | Installs the skill with current style files bundled. |
| `build_prompt.py` | Regenerates `FIGURE_PROMPT.md` from `highconc_style.py`. |

Changed a style value? Edit `highconc_style.py`, then run `python
build_prompt.py`, then update `../STYLE_GUIDE.md` to match.
