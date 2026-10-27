"""Helpers for the caption pass: replace a markdown cell of figures.py, found by tag + heading."""
import re
from pathlib import Path

FIG = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\notebooks\figures.py")


def load():
    t = FIG.read_bytes().decode("utf-8").replace("\r\n", "\n")
    return [c.rstrip("\n") for c in re.split(r"(?m)^(?=# %%)", t) if c.startswith("# %%")]


def save(cells):
    FIG.write_bytes(("\n\n".join(cells) + "\n").encode("utf-8"))


def md(tag, text):
    """text = plain markdown; becomes a commented cell with the right header."""
    body = "\n".join(("# " + l).rstrip() if l.strip() else "#" for l in text.strip("\n").split("\n"))
    return f'# %% [markdown] tags=["{tag}"]\n{body}'


def replace_md(cells, tag, heading_prefix, new_text):
    """Replace the markdown cell with this tag whose first non-header line starts with heading_prefix."""
    hits = [i for i, c in enumerate(cells)
            if c.startswith(f'# %% [markdown] tags=["{tag}"]')
            and c.split("\n")[1].startswith(heading_prefix)]
    assert len(hits) == 1, (tag, heading_prefix, hits)
    cells[hits[0]] = md(tag, new_text)
    return hits[0]
