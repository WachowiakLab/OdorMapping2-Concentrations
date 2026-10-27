"""Fig 3's 4th example glomerulus: dataset id 155, set directly (no subject lookup)."""
from pathlib import Path

F = Path(r"C:\Users\u6068031\Documents\highconc-merge\final")
OLD = '''# Fourth example glomerulus (panel H). Not part of EG_GLOMERULI, so Fig 6 is unaffected.
# Position within one mouse, counted from 0 among the glomeruli that survive the
# `nonresponders_dropped` filter.
EXTRA_SUBJECT = 1
EXTRA_GLOM_IN_SUBJECT = 8
EXTRA_GLOM = int(np.flatnonzero(d.subject == EXTRA_SUBJECT)[EXTRA_GLOM_IN_SUBJECT])
'''
NEW = '''# Fourth example glomerulus (panel H), by dataset glomerulus id like EG_GLOMERULI.
# Kept out of EG_GLOMERULI, which other code imports.
EXTRA_GLOM = 155
'''
for rel in ("notebooks/figures.py", "notes/sections/03_fig3_code.py"):
    p = F / rel
    t = p.read_bytes().decode("utf-8").replace("\r\n", "\n")
    assert t.count(OLD) == 1, (rel, t.count(OLD))
    p.write_bytes(t.replace(OLD, NEW).encode("utf-8"))
    print("edited", rel)

log = F / "notes" / "merge_decisions.md"
t = log.read_bytes().decode("utf-8")
old = "- 4th example glomerulus from lab (`EXTRA_GLOM`, subject 1, #8): spectra are E-H,"
assert t.count(old) == 1
t = t.replace(old, "- 4th example glomerulus: **dataset id 155** (`EXTRA_GLOM = 155`; user's choice after lab's\n"
                   "  \"subject 1, #8\" lookup raised IndexError on real data). Spectra are E-H,")
log.write_bytes(t.encode("utf-8"))
print("log updated")
