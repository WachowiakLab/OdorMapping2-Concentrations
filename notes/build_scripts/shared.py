"""Shared modules + tooling into final/: main's files, plus the agreed adjustments. Every edit asserts."""
import shutil
from pathlib import Path

ROOT = Path(r"C:\Users\u6068031\Documents\highconc-merge")
M = ROOT / "highconc-figures-lab-main"
F = ROOT / "final"
log = []


def rd(p):
    return p.read_bytes().decode("utf-8").replace("\r\n", "\n")


def wr(p, text):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(text.replace("\r\n", "\n").encode("utf-8"))


def edit(text, edits, label):
    for old, new in edits:
        assert text.count(old) == 1, (label, text.count(old), old[:70])
        text = text.replace(old, new)
    return text


# 1. main's files, taken as they are (copy == main; lab changed nothing or is corrupt)
for rel in ["highconc/similarity.py", "highconc/clusters.py", "highconc/flow.py",
            "highconc/recruitment.py", "highconc/panels/recruitment.py",
            "highconc/panels/__init__.py", "dev/build_notebooks.py"]:
    wr(F / rel, rd(M / rel))
    log.append(f"main's {rel}")

# data needed by Fig 2 (identical in main and copy): binary-safe copy
for p in sorted((M / "data" / "burton_et_al").glob("*")):
    dst = F / "data" / "burton_et_al" / p.name
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(p, dst)
    log.append(f"data/burton_et_al/{p.name}")

# 2. style.py: main's (SEQ = magma, class_strip side=) + lab's class_strip offset=
t = edit(rd(M / "highconc" / "style.py"), [
    ('''                side="right"):''', '''                side="right", offset=-0.5):'''),
    ('''    Returns {category: colour} so the same mapping can be reused on sibling''',
     '''    `offset` is where cell 0 starts in data coordinates. The default -0.5 fits
    `imshow`/`hs.heatmap` (cell i centred on i); pass `offset=0` for
    `sns.heatmap` or `pcolormesh` (cell i spans i..i+1).

    Returns {category: colour} so the same mapping can be reused on sibling'''),
    ('''ax.add_patch(plt.Rectangle((a - .5, -pad - height), b - a, height,''',
     '''ax.add_patch(plt.Rectangle((a + offset, -pad - height), b - a, height,'''),
    ('''ax.text((a + b - 1) / 2, -pad - height - .012, name,''',
     '''ax.text((a + b) / 2 + offset, -pad - height - .012, name,'''),
    ('''ax.add_patch(plt.Rectangle((1 + pad, a - .5), height, b - a,''',
     '''ax.add_patch(plt.Rectangle((1 + pad, a + offset), height, b - a,'''),
    ('''ax.text(1 + pad + height + .008, (a + b - 1) / 2, name,''',
     '''ax.text(1 + pad + height + .008, (a + b) / 2 + offset, name,'''),
    ('''ax.add_patch(plt.Rectangle((-pad - height, a - .5), height, b - a,''',
     '''ax.add_patch(plt.Rectangle((-pad - height, a + offset), height, b - a,'''),
    ('''ax.text(-pad - height - .008, (a + b - 1) / 2, name,''',
     '''ax.text(-pad - height - .008, (a + b) / 2 + offset, name,'''),
], "style.py")
wr(F / "highconc" / "style.py", t)
log.append("style.py = main's + lab's class_strip offset")

# 3. Fig 4B: its class strip sits on a seaborn heatmap (row i spans i..i+1) -> offset=0
p = F / "highconc" / "panels" / "clustering.py"
wr(p, edit(rd(p), [(
    '''colors=strip_colors, orientation="y", side="left",''',
    '''colors=strip_colors, orientation="y", side="left", offset=0,''')], "clustering.py"))
log.append("panels/clustering.py: class_strip offset=0 on the seaborn heatmap")

# 4. Fig 7 A pinned to the original (baseline) OrRd ramp
p = F / "highconc" / "panels" / "transfer.py"
wr(p, edit(rd(p), [
    ('''def pvalue_matrix(ax, t, *, alpha=0.05, clip=0.1, values=True, cbar=True,
                  cax=None, title=None, **kw):''',
     '''def pvalue_matrix(ax, t, *, alpha=0.05, clip=0.1, values=True, cbar=True,
                  cax=None, title=None, cmap=None, **kw):'''),
    ('''    cmap = hs.SEQ
    sig = t.significant(alpha)''',
     '''    cmap = hs.SEQ if cmap is None else cmap
    sig = t.significant(alpha)'''),
], "panels/transfer.py"))
p = F / "notes" / "sections" / "07_fig7_code.py"
wr(p, edit(rd(p), [
    ('''FIG7_XLABEL = "cluster (predicted at high)"''',
     '''FIG7_P_CMAP = "OrRd"         # A keeps the original sequential ramp (hs.SEQ is magma now)
FIG7_XLABEL = "cluster (predicted at high)"'''),
    ('''                         xlabel=FIG7_XLABEL, ylabel=FIG7_YLABEL,
                         box_color=FIG7_BOX_COLOR, text_significant_only=True)''',
     '''                         xlabel=FIG7_XLABEL, ylabel=FIG7_YLABEL, cmap=FIG7_P_CMAP,
                         box_color=FIG7_BOX_COLOR, text_significant_only=True)'''),
], "07 section"))
log.append("Fig 7 A: cmap option + FIG7_P_CMAP = OrRd")

# 5. pyproject: private dependencies over HTTPS (standard form, no `git@` user)
p = F / "pyproject.toml"
t = rd(p)
assert t.count("git+ssh://git@github.com/") == 2
wr(p, t.replace("git+ssh://git@github.com/", "git+https://github.com/"))
log.append("pyproject.toml: ssh -> https")
p = F / "README.md"
wr(p, edit(rd(p), [(
    "Needs Python ≥ 3.11 and SSH access to the two private workhorse repos.",
    "Needs Python ≥ 3.11 and HTTPS access to the two private workhorse repos "
    "(a GitHub account with access to CastroLab, signed in through a personal "
    "access token or the GitHub CLI credential helper).")], "README.md"))
log.append("README.md: SSH -> HTTPS wording")

# 6. .claude/settings.json: baseline's exact commands + PowerShell mirrors; no where.py, no wide prefixes
cmds = ['python notebooks/figures.py',
        'python dev/concordance.py --check',
        'python dev/concordance.py --status',
        'python dev/concordance.py --panels',
        'python dev/build_notebooks.py',
        'python dev/build_notebooks.py --check',
        'python dev/build_site.py',
        'python -m nbconvert --to notebook --execute:*',
        'python -m ipykernel install --user --name highconc:*',
        'python -c \\"import highconc; print(highconc.__file__)\\"']
rows = [f'      "Bash({c})"' for c in cmds] + [f'      "PowerShell({c})"' for c in cmds]
settings = ('{\n  "permissions": {\n    "allow": [\n' + ",\n".join(rows) + "\n    ]\n  }\n}\n")
import json
json.loads(settings)   # must be valid JSON
wr(F / ".claude" / "settings.json", settings)
log.append(".claude/settings.json: baseline commands + PowerShell mirrors")

print("\n".join(log))
