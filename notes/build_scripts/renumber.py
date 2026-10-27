"""Stage 1: renumber sections/constants/stems to the paper's figure numbers.

Old section -> new section / figure:
  5 bulb maps   -> 6   (FIG5_ -> FIG6_, fig5_spatial_by_cluster -> fig6_spatial_by_cluster)
  6 PLS         -> 7   (FIG6_ -> FIG7_, fig6_pls -> fig7_pls)
  7 transfer    -> 8   (FIG7_ -> FIG8_, fig7_motif_transfer -> fig8_motif_transfer)
  8 S1          -> 9
  9 S2 (A + B)  -> 10 (S2 heatmaps) + 11 (S3 topography, split out)
 11 S4          -> 12
New section 5 (class x cluster) is added in stage 3.
"""
import re
from pathlib import Path

F = Path(r"C:\Users\u6068031\Documents\highconc-merge\final")
P = F / "notebooks" / "figures.py"
text = P.read_bytes().decode("utf-8").replace("\r\n", "\n")
cells = [c.rstrip("\n") for c in re.split(r"(?m)^(?=# %%)", text) if c.startswith("# %%")]

# current section of every cell
sec, cur = [], 0
for c in cells:
    m = re.search(r"(?m)^# ## (\d+)\.", c) if c.startswith("# %% [markdown]") else None
    if m:
        cur = int(m.group(1))
    sec.append(cur)

NEW = {5: 6, 6: 7, 7: 8, 8: 9, 9: 10, 11: 12}
PREFIX = {5: ("FIG5_", "FIG6_"), 6: ("FIG6_", "FIG7_"), 7: ("FIG7_", "FIG8_")}
STEM = {5: ("fig5_spatial_by_cluster", "fig6_spatial_by_cluster"),
        6: ("fig6_pls", "fig7_pls"), 7: ("fig7_motif_transfer", "fig8_motif_transfer")}
FIGNAME = {5: ("Figure 5", "Figure 6"), 6: ("Figure 6", "Figure 7"), 7: ("Figure 7", "Figure 8")}

out = []
for c, s in zip(cells, sec):
    if s in NEW:
        md = c.startswith("# %% [markdown]")
        if s in PREFIX:
            a, b = PREFIX[s]
            c = c.replace(a, b)                                  # constants (and mentions in prose)
        if s in STEM:
            a, b = STEM[s]
            c = c.replace(a, b)
        if md:
            # heading numbers: "## 5. Figure 5 —" -> "## 6. Figure 6 —"
            c = re.sub(r"(?m)^# ## %d\. " % s, "# ## %d. " % NEW[s], c)
            if s in FIGNAME:
                a, b = FIGNAME[s]
                c = re.sub(r"(?m)^(# ## \d+\. )%s " % re.escape(a), r"\g<1>%s " % b, c)
            if s in (8, 9, 11):                                   # supplements: figure names unchanged
                pass
    out.append(c)

# comments that point at another section by number
text2 = "\n\n".join(out) + "\n"
for old, new in [("rather than inherited from §6, so the two", "rather than inherited from §7, so the two")]:
    if old in text2:
        text2 = text2.replace(old, new)
P.write_bytes(text2.encode("utf-8"))

# report
cells2 = [c for c in re.split(r"(?m)^(?=# %%)", text2) if c.startswith("# %%")]
heads = [re.search(r"(?m)^# ## (\d+\. .*)$", c).group(1)[:70] for c in cells2
         if c.startswith("# %% [markdown]") and re.search(r"(?m)^# ## \d+\.", c)]
print("\n".join(heads))
print("leftover old tokens:", {t: text2.count(t) for t in ("FIG5_", "fig5_spatial", "fig6_pls", "fig7_motif")})
