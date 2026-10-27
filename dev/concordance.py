"""Concordance harness — does the notebook produce exactly the agreed figures?

The refactor's own harness asked "does this code draw the same picture as
before?" That question is answered and closed. This one asks something else,
and needs a different instrument:

    Does the notebook emit EXACTLY the figures in dev/CONCORDANCE.yaml —
    no more, no fewer?

That is the actual specification now. "These are the only figures we need to
show" is a claim about a *set*, so the harness checks a set. An extra render is
scope creep; a missing one is unfinished work; and both are silent failures
without something that looks.

    python dev/concordance.py --check      # notebook outputs vs the manifest
    python dev/concordance.py --status     # per-figure progress table
    python dev/concordance.py --panels     # every panel and its status

--check reads `notebooks/outputs/`, so run the notebook first. It does NOT run
the notebook itself: rendering is minutes, and a harness that silently triggers
it is the same trap `dimsweep.load_cached` exists to avoid.

Exit codes: 0 clean, 1 mismatch (extra or missing), 2 manifest unreadable,
3 outputs directory absent.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent.parent
MANIFEST = Path(__file__).resolve().parent / "CONCORDANCE.yaml"
OUTPUTS = REPO / "notebooks" / "outputs"

#: A figure is "shippable" once every panel is done or deliberately skipped.
OPEN_STATUSES = {"new", "gap", "restyle"}


def load_manifest() -> dict:
    if not MANIFEST.exists():
        sys.exit(f"no manifest at {MANIFEST}")
    try:
        return yaml.safe_load(MANIFEST.read_text())
    except yaml.YAMLError as exc:
        print(f"manifest is not valid YAML: {exc}", file=sys.stderr)
        sys.exit(2)


def figures(man: dict) -> list[dict]:
    """Every figure entry, main then supplemental, in manifest order."""
    return list(man.get("main") or []) + list(man.get("supplemental") or [])


def expected_stems(man: dict) -> set[str]:
    return {f["stem"] for f in figures(man)}


def rendered_stems() -> set[str]:
    """Distinct stems in notebooks/outputs/, ignoring the format suffix."""
    if not OUTPUTS.is_dir():
        print(f"no outputs at {OUTPUTS} — run the notebook first", file=sys.stderr)
        sys.exit(3)
    return {p.stem for p in OUTPUTS.iterdir()
            if p.suffix.lower() in {".png", ".pdf", ".svg"}}


def panel_counts(fig: dict) -> tuple[int, int]:
    """(open panels, total panels) for one figure, variants included."""
    rows = list(fig.get("panels") or []) + list(fig.get("variants") or [])
    total = len(rows)
    still_open = sum(1 for r in rows if r.get("status") in OPEN_STATUSES)
    return still_open, total


def cmd_check(man: dict) -> int:
    want = expected_stems(man)
    have = rendered_stems()
    missing = sorted(want - have)
    extra = sorted(have - want)

    retired = {r["stem"] for r in (man.get("retired") or [])}
    undecided = {r["stem"] for r in (man.get("undecided") or [])}
    # Earlier stems that a current figure absorbed. Every replacement has
    # landed and every source section was removed, so these are gone too.
    superseded = {s for f in figures(man) for s in (f.get("supersedes") or [])}

    scope_creep = [s for s in extra
                   if s not in retired and s not in undecided and s not in superseded]
    regressions = [s for s in extra if s in retired or s in superseded]
    open_questions = [s for s in extra if s in undecided]

    print(f"manifest: {len(want)} figures | rendered: {len(have)} stems")

    if missing:
        print(f"\nMISSING ({len(missing)}) — in the manifest, not rendered:")
        for s in missing:
            fig = next(f for f in figures(man) if f["stem"] == s)
            op, tot = panel_counts(fig)
            print(f"    {s:34} {fig['id']:18} {op}/{tot} panels still open")
    if regressions:
        print(f"\nREGRESSION ({len(regressions)}) — retired renders that "
              f"something is drawing again:")
        for s in regressions:
            owner = next((f["id"] for f in figures(man)
                          if s in (f.get("supersedes") or [])), None)
            where = f"-> absorbed by {owner}" if owner else "-> retired outright"
            print(f"    {s:34} {where}")
        print("  A notebook section came back, or an old file survived in "
              "notebooks/outputs/.")
    if open_questions:
        print(f"\nUNDECIDED ({len(open_questions)}) — referenced in the MS but "
              f"with no settled home. Ask before shipping or retiring:")
        for s in open_questions:
            why = next(r.get("why", "") for r in man["undecided"] if r["stem"] == s)
            print(f"    {s:34} {why}")
    if scope_creep:
        print(f"\nUNACCOUNTED ({len(scope_creep)}) — rendered, and in none of "
              f"the manifest, retirement, supersedes or undecided lists:")
        for s in scope_creep:
            print(f"    {s}")
        print("  Either add it to the manifest with a reason, or stop drawing it.")

    if not (missing or scope_creep or regressions):
        if open_questions:
            print("\nCONCORDANCE CLEAN except for the undecided stems above.")
            return 1
        print("\nCONCORDANCE CLEAN — the notebook draws exactly the agreed set.")
        return 0
    return 1


def cmd_status(man: dict) -> int:
    print(f"{'figure':22} {'stem':34} {'stage':7} {'panels':>10}  status")
    print("-" * 92)
    for key in ("main", "supplemental"):
        rows = man.get(key) or []
        if rows:
            print(f"[{key}]")
        for fig in rows:
            op, tot = panel_counts(fig)
            verdict = "READY" if op == 0 else f"{op} open"
            print(f"  {fig['id']:20} {fig['stem']:34} {fig.get('stage',''):7} "
                  f"{tot - op:>4}/{tot:<5}  {verdict}")
    retired = man.get("retired") or []
    print(f"\n{len(retired)} earlier renders retired.")
    return 0


def cmd_panels(man: dict) -> int:
    for fig in figures(man):
        print(f"\n{fig['id']}  —  {fig.get('title', fig.get('source',''))}")
        for row in list(fig.get("panels") or []) + list(fig.get("variants") or []):
            name = row.get("panel", row.get("id", "?"))
            print(f"    {str(name):8} {row.get('status',''):9} {row.get('what','')}")
        if fig.get("open_question"):
            print(f"    ! blocks on open question: {fig['open_question']}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--check", action="store_true",
                   help="notebook outputs vs the manifest (the real test)")
    g.add_argument("--status", action="store_true", help="per-figure progress")
    g.add_argument("--panels", action="store_true", help="every panel and status")
    args = ap.parse_args()

    man = load_manifest()
    if args.check:
        return cmd_check(man)
    if args.status:
        return cmd_status(man)
    return cmd_panels(man)


if __name__ == "__main__":
    raise SystemExit(main())
