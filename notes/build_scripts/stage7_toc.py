"""Stage 7: the intro's contents tables (human + llm) in the paper's numbering."""
import re
from pathlib import Path

P = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\notebooks\figures.py")
t = P.read_bytes().decode("utf-8").replace("\r\n", "\n")

# ---- human table ----
a = t.index("# | § | Figure | In one line |")
b = t.index("# | 11 | **Fig S4** | Predicted vs measured, at three levels of averaging |")
b = t.index("\n", b) + 1
human = """# | § | Figure | In one line |
# |---|---|---|
# | 1 | setup | Loads the data and fits the model. Run this first. |
# | 2 | **Fig 2** | Odorant responses across Burton's own glomeruli, at three concentrations |
# | 3 | **Fig 3** | Odorant–odorant similarity, the recruitment matrix, and four example glomeruli |
# | 4 | **Fig 4** | The seven glomerular groups, and how they carry over when high-only glomeruli are folded in |
# | 5 | **Fig 5** | Which chemical classes each group responds to |
# | 6 | **Fig 6** | Where those groups sit on the bulb |
# | 7 | **Fig 7** | A linear model predicting the low → high change |
# | 8 | **Fig 8** | Directional transfer of tuning between groups |
# | 9 | **Fig S1** | How the number of groups was chosen |
# | 10 | **Fig S2** | The same structure within each mouse |
# | 11 | **Fig S3** | Where each mouse's groups sit on the bulb |
# | 12 | **Fig S4** | Predicted vs measured, over every responding pair |
"""
t = t[:a] + human + t[b:]
old = "# Section numbers are the paper's figure numbers, so §3 draws Figure 3."
assert t.count(old) == 1
t = t.replace(old, "# Section numbers are the paper's figure numbers, so §3 draws Figure 3; the\n"
                   "# supplements follow, §9 to §12 being Figures S1 to S4. Panels the paper has but this\n"
                   "# repo does not draw (schematics, photographs, Fig 2B-E, Fig 5A) are left to the figure's own\n"
                   "# panel letters, so a letter can start at B or D.")

# ---- llm table ----
a = t.index("# | § | Figure | Stem | Panels |")
b = t.index("`fig4_tuning_by_cluster`", a)   # sanity: still inside the old table
b = t.index("# | 11 | **Figure S4**", a)
b = t.index("\n", b) + 1
llm = """# | § | Figure | Stem | Panels |
# |---|---|---|---|
# | 1 | Setup — style, paths, shared data and model fit | — | — |
# | 2 | **Figure 2** — odorant x glomerulus response, Burton's own dataset, three concentrations | `fig2_response_matrix` | `panels.similarity` |
# | 3 | **Figure 3** — correlational structure vs. concentration, recruitment, example spectra; and its threshold supplement | `fig3_correlational_structure`, `fig3_supp_threshold_robustness` | `panels.similarity`, `panels.recruitment` |
# | 4 | **Figure 4** — NMF clusters, per-subject detail, and the folded-in Sankey | `fig4_tuning_by_cluster` | `panels.clustering`, `panels.sankey` (analysis: `folding`) |
# | 5 | **Figure 5** — chemical-class composition of each cluster, low vs high | `fig5_class_by_cluster` | `panels.clustering` |
# | 6 | **Figure 6** — spatial mapping of glomeruli by cluster | `fig6_spatial_by_cluster` | `panels.clustering` |
# | 7 | **Figure 7** — PLS-R: predicted spectra and the W matrix | `fig7_pls` | `panels.similarity`, `panels.pls` |
# | 8 | **Figure 8** — motif transfer between glomerular groups | `fig8_motif_transfer` | `panels.transfer` |
# | 9 | **Figure S1** — NMF dimensionality choice (low concentration) | `supp_nmf_dim_choice` | `panels.dimchoice` |
# | 10 | **Figure S2** — per-subject response heatmaps | `supp_per_subject_heatmaps` | `panels.clustering` |
# | 11 | **Figure S3** — per-subject topography of the clusters | `supp_per_subject_topography` | `panels.clustering` |
# | 12 | **Figure S4** — parity over every responding pair | `supp_parity` | `panels.pls` |
"""
t = t[:a] + llm + t[b:]
P.write_bytes(t.encode("utf-8"))
print("contents tables rewritten")
