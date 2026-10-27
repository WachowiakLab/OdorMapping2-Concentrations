"""S2A back to the log scale; drop the now-unused `norm` option from per_subject_heatmap."""
from pathlib import Path

FINAL = Path(r"C:\Users\u6068031\Documents\highconc-merge\final")


def edit_file(path, edits):
    t = path.read_bytes().decode("utf-8").replace("\r\n", "\n")
    for old, new in edits:
        assert t.count(old) == 1, (path.name, t.count(old), old[:70])
        t = t.replace(old, new)
    path.write_bytes(t.encode("utf-8"))


# section
edit_file(FINAL / "notes" / "sections" / "09_figS2_code.py", [
    ('''# A: raw dF/F on a linear scale, odorant rows in Fig 3 / Fig 4's chemical-class order.
SUPP_SUBJ_HEAT_VMIN, SUPP_SUBJ_HEAT_VMAX = 0.0, 100.0''',
     '''# A: raw dF/F on a log scale, odorant rows in Fig 3 / Fig 4's chemical-class order.
SUPP_SUBJ_HEAT_VMIN, SUPP_SUBJ_HEAT_VMAX = 1.0, 140.0'''),
    ('''supp_subj_norm = plt.Normalize(SUPP_SUBJ_HEAT_VMIN, SUPP_SUBJ_HEAT_VMAX)
''', ''),
    ('''                                     cmap=SUPP_SUBJ_CMAP,
                                     norm=supp_subj_norm,
                                     boundary_color="white")''',
     '''                                     cmap=SUPP_SUBJ_CMAP,
                                     vmin=SUPP_SUBJ_HEAT_VMIN,
                                     vmax=SUPP_SUBJ_HEAT_VMAX,
                                     boundary_color="white")'''),
    ('''cb = fig.colorbar(plt.cm.ScalarMappable(norm=supp_subj_norm,
                                        cmap=SUPP_SUBJ_CMAP),
                  cax=cax, extend="max")''',
     '''cb = fig.colorbar(plt.cm.ScalarMappable(
    norm=LogNorm(SUPP_SUBJ_HEAT_VMIN, SUPP_SUBJ_HEAT_VMAX),
    cmap=SUPP_SUBJ_CMAP), cax=cax)'''),
    ('''cb.set_ticks(np.linspace(SUPP_SUBJ_HEAT_VMIN, SUPP_SUBJ_HEAT_VMAX, 5))
cb.ax.tick_params(labelsize=8)
cb.outline.set_visible(False)''',
     '''cb.set_ticks([1, 10, 100])
cb.set_ticklabels(["1", "10", "100"])
cb.ax.tick_params(labelsize=8)
cb.ax.minorticks_off()
cb.outline.set_visible(False)'''),
])

# log
edit_file(FINAL / "notes" / "merge_decisions.md", [
    ('''boundaries, class strip on the first heatmap, y label "odorants (by chemical class)". User: ordering''',
     '''boundaries, class strip on the first heatmap, y label "odorants (by chemical class)". User: ordering'''),
    ('''not lab's hard-coded `("matt", "table")`. Scale (linear 0-100 magma, not peak-normalised like Fig 4C) was
  my lean and NOT confirmed. Width `SUPP_TARGET_W` (5.5 in), fonts as main (7 pt title, 5 pt y label).''',
     '''not lab's hard-coded `("matt", "table")`. **Scale: user chose LOG (dF/F 1-140, ticks 1/10/100), not lab's
  linear 0-100.** Colormap stays magma (lab's; matches Fig 4) — user did not comment on the colormap.
  Width `SUPP_TARGET_W` (5.5 in), fonts as main (7 pt title, 5 pt y label).'''),
    ('''- Panel patches in `final/highconc/panels/clustering.py`: `per_subject_heatmap(norm=None)` (clean version of
  lab's half-done arg; default behaviour unchanged) and `per_subject_topography(context_alpha=None)`.''',
     '''- Panel patch in `final/highconc/panels/clustering.py`: `per_subject_topography(context_alpha=None)` only.
  (An earlier `per_subject_heatmap(norm=...)` patch for lab's linear scale was removed once log was chosen.)'''),
])
print("ok")
