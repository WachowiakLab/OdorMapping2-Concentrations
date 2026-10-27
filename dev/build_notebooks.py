"""Build the two reader-facing notebooks from the one source.

`notebooks/figures.py` is the **source of truth**. It carries every figure's
code once, and each caption twice: a `human`-tagged markdown cell and an
`llm`-tagged one. This script splits it into two notebooks that differ only in
prose:

    notebooks/figures_human.{py,ipynb}   friendly: the claim, then how it was made
    notebooks/figures_LLM.{py,ipynb}     dense: provenance, parameters, exact values

**Why a generator and not two notebooks.** The whole value of having two is
that they draw the *same figures*. Two hand-maintained copies drift — someone
fixes a colour scale in one and not the other, and now the lab and the record
disagree about what the figure is. Here the code cells are copied verbatim into
both, so they cannot differ. Only prose can.

    python dev/build_notebooks.py            # rebuild both
    python dev/build_notebooks.py --check    # are they current? (exit 1 if not)

`--check` is the one to run before shipping: it rebuilds in memory and compares,
so a source edit that was never rebuilt fails loudly instead of leaving a stale
notebook in the clone someone else pulls.

Tagging rules, all of them:

* markdown tagged `human`  -> human notebook only
* markdown tagged `llm`    -> LLM notebook only
* markdown with no tag     -> both
* code                     -> both, always, byte for byte

Exit codes: 0 fine, 1 out of date (--check only), 2 source unreadable.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SOURCE = REPO / "notebooks" / "figures.py"
VARIANTS = {"human": REPO / "notebooks" / "figures_human.py",
            "llm": REPO / "notebooks" / "figures_LLM.py"}

BANNER = {
    "human": "the friendly one — what each figure shows, and how it was made",
    "llm": "the exact one — provenance, parameters and values, stated in full",
}

_CELL = re.compile(r"^# %%")
_TAGS = re.compile(r'tags=\[([^\]]*)\]')


def parse_cells(text: str) -> list[dict]:
    """Split percent-format source into cells with their kind and tags."""
    cells: list[dict] = []
    cur: dict | None = None
    for line in text.split("\n"):
        if _CELL.match(line):
            if cur is not None:
                cells.append(cur)
            m = _TAGS.search(line)
            tags = ([t.strip().strip('"\'') for t in m.group(1).split(",") if t.strip()]
                    if m else [])
            cur = {"markdown": "[markdown]" in line, "tags": tags, "body": []}
        elif cur is not None:
            cur["body"].append(line)
    if cur is not None:
        cells.append(cur)
    return cells


def render(cells: list[dict], variant: str) -> str:
    """The percent-format text of one variant."""
    other = {"human", "llm"} - {variant}
    out = [
        "# %% [markdown]",
        f"# <!-- GENERATED FILE — do not edit. Source: notebooks/figures.py.",
        f"#      Rebuild with `python dev/build_notebooks.py`. -->",
        f"#",
        f"# *This is the **{variant}** edition of the figure notebook — "
        f"{BANNER[variant]}.*",
        f"# *Same code, same figures, different words: the other edition is "
        f"`figures_{'human' if variant == 'llm' else 'LLM'}`.*",
        f"# *Edits go in `notebooks/figures.py`, never here — this file is "
        f"overwritten.*",
        "",
    ]
    for c in cells:
        if c["markdown"] and (set(c["tags"]) & other):
            continue
        out.append("# %% [markdown]" if c["markdown"] else "# %%")
        out.extend(c["body"])
    return "\n".join(out).rstrip("\n") + "\n"


def build(check: bool) -> int:
    if not SOURCE.exists():
        print(f"no source at {SOURCE}", file=sys.stderr)
        return 2
    cells = parse_cells(SOURCE.read_text(encoding="utf-8"))
    n_code = sum(1 for c in cells if not c["markdown"])
    stale = []

    for variant, path in VARIANTS.items():
        text = render(cells, variant)
        n_md = sum(1 for c in cells
                   if c["markdown"] and not (set(c["tags"]) & ({"human", "llm"} - {variant})))
        if check:
            current = path.read_text(encoding="utf-8") if path.exists() else None
            mark = "ok" if current == text else "STALE"
            if current != text:
                stale.append(path.name)
            print(f"  {path.name:24} {n_md:3d} markdown + {n_code:3d} code cells   {mark}")
            continue

        path.write_text(text, encoding="utf-8")
        # `sys.executable -m jupytext`, not a bare `jupytext`: the console
        # script is only on PATH if the env is activated, and this is run from
        # notebooks and CI as often as from an activated shell.
        proc = subprocess.run(
            [sys.executable, "-m", "jupytext", "--to", "ipynb",
             "-o", str(path.with_suffix(".ipynb")), str(path)],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            print(proc.stderr.strip(), file=sys.stderr)
            raise SystemExit(
                f"jupytext failed for {path.name}. Install it into this "
                f"interpreter: {sys.executable} -m pip install jupytext"
            )
        print(f"  {path.name:24} {n_md:3d} markdown + {n_code:3d} code cells   written"
              f"  (+ {path.with_suffix('.ipynb').name})")

    if check and stale:
        print(f"\n{len(stale)} notebook(s) out of date: {', '.join(stale)}")
        print("Run `python dev/build_notebooks.py` and commit the result.")
        return 1
    if check:
        print("\nboth notebooks are current.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="verify the built notebooks match the source")
    args = ap.parse_args()
    print(f"source: {SOURCE.relative_to(REPO)}")
    return build(args.check)


if __name__ == "__main__":
    raise SystemExit(main())
