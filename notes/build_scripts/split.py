"""Dump per-section code (and prose) of every notebook source into scratchpad/sec/."""
import re, sys
from pathlib import Path

ROOT = Path(r"C:\Users\u6068031\Documents\highconc-merge")
OUT = Path(__file__).parent / "sec"
VERS = {"base": "github/original_codebase", "lab": "highconc-figures-lab",
        "main": "highconc-figures-lab-main", "copy": "highconc-figures-lab-main - Copy"}
SRCS = {"fig": "figures.py", "hum": "figures_human.py", "llm": "figures_LLM.py"}


def read(p):
    b = p.read_bytes()
    try:
        t = b.decode("utf-8")
    except UnicodeDecodeError:
        t = b.decode("cp1252")
    return t.replace("\r\n", "\n")


def split(text):
    cells = [c for c in re.split(r"(?m)^(?=# %%)", text) if c.startswith("# %%")]
    sec, out = "0", {}
    for c in cells:
        md = c.startswith("# %% [markdown]")
        m = re.search(r"(?m)^# ## (\d+)\.", c) if md else None
        if m:
            sec = m.group(1)
        out.setdefault(sec, {"code": [], "md": []})["md" if md else "code"].append(c.rstrip("\n"))
    return out


OUT.mkdir(exist_ok=True)
for v, d in VERS.items():
    for s, fn in SRCS.items():
        p = ROOT / d / "notebooks" / fn
        if not p.exists():
            continue
        for sec, parts in split(read(p)).items():
            (OUT / f"{sec.zfill(2)}_{v}_{s}_code.py").write_text("\n\n".join(parts["code"]) + "\n", encoding="utf-8", newline="\n")
            (OUT / f"{sec.zfill(2)}_{v}_{s}_md.md").write_text("\n\n".join(parts["md"]) + "\n", encoding="utf-8", newline="\n")
print(sorted({f.name.split("_")[0] for f in OUT.iterdir()}))
