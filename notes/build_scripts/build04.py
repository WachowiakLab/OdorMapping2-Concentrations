"""Assemble final's Fig 4 code section from main's figures.py + the user's choices."""
from pathlib import Path

HERE = Path(__file__).parent
src = (HERE / "sec" / "04_main_fig_code.py").read_text(encoding="utf-8")
OUT = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\notes\sections\04_fig4_code.py")


def sub(old, new):
    global src
    assert src.count(old) == 1, (src.count(old), old[:70])
    src = src.replace(old, new)


# B and C share the 'matt' class order AND the 'matt' within-class order (== §3's)
sub('''bic = P.clustering.bicluster_order(
    d, chem_order=(FIG4_CHEM_ORDER if FIG4B_SORT_BY_CLASS else None))
bic_b = P.clustering.bicluster_order(
    d, chem_order=FIG4_CHEM_ORDER, row_order=FIG4B_ROW_ORDER)''',
    '''bic = P.clustering.bicluster_order(
    d, chem_order=(FIG4_CHEM_ORDER if FIG4B_SORT_BY_CLASS else None),
    within_order=WITHIN_CLASS_ORDER)
bic_b = P.clustering.bicluster_order(
    d, chem_order=FIG4_CHEM_ORDER, within_order=WITHIN_CLASS_ORDER,
    row_order=FIG4B_ROW_ORDER)''')

# the normalisation comment now covers C as well
sub('# Peak-normalise each glomerulus (column) to [0, 1] by its own max response,',
    '# Panels B and C both: peak-normalise each glomerulus (column) to [0, 1] by its\n'
    '# own max response,')

# C: B's exact row order (asserted, not recomputed), and the same normalisation
sub('''# C's rows are Fig 3/6's exact 59-odorant order (CLASS_ORDER_FIG3 blocks,
# WITHIN_CLASS_ORDER within them), not the cluster-sorted rows `bic` carries.
# Only the rows move: columns stay one cluster's glomeruli per subject.
perm_c, _, _, class_bounds_c = similarity.class_sorted(
    d, CLASS_ORDER_FIG3, WITHIN_CLASS_ORDER)
bic_c = dataclasses.replace(bic,
                            rows=perm_c,''',
    '''# C's rows are B's rows: Fig 3/6's exact 59-odorant order (CLASS_ORDER_FIG3
# blocks, WITHIN_CLASS_ORDER within them), not the cluster-sorted rows `bic`
# carries. Only the rows move: columns stay one cluster's glomeruli per
# subject. C also takes B's peak-normalisation and raw floor, so a glomerulus
# reads identically in both and B's colorbar is C's too.
perm_c, _, _, class_bounds_c = similarity.class_sorted(
    d, CLASS_ORDER_FIG3, WITHIN_CLASS_ORDER)
assert FIG4B_ROW_ORDER == "odorant" and np.array_equal(perm_c, bic_b.rows), \\
    "C's rows must be B's rows (FIG4B_ROW_ORDER='odorant', same class/within order)"
bic_c = dataclasses.replace(bic,
                            rows=perm_c,''')
sub('''            d, axes_c[j], subj, cluster=k, conc=FIG4_CONC_A, order=bic_c,
            vmin=FIG4_VMIN, vmax=FIG4_VMAX, cmap=FIG4_CMAP,
            empty_color=FIG4_EMPTY, cluster_labels=False,''',
    '''            d, axes_c[j], subj, cluster=k, conc=FIG4_CONC_A, order=bic_c,
            col_normalize=FIG4B_COL_NORMALIZE,
            raw_floor=FIG4B_RAW_FLOOR if FIG4B_COL_NORMALIZE else None,
            vmin=FIG4B_NORM_VMIN if FIG4B_COL_NORMALIZE else FIG4_VMIN,
            vmax=FIG4B_NORM_VMAX if FIG4B_COL_NORMALIZE else FIG4_VMAX,
            cmap=FIG4_CMAP, empty_color=FIG4_EMPTY, cluster_labels=False,''')

# D: normalised within each concentration (class_cluster_glom_composition divides each
# cluster column by its own max, separately for low and high) -- say so in the cell.
sub('# D — the class x cluster summary, low and high each with its own colour scale.',
    '# D — the class x cluster summary. Each concentration is normalised within itself\n'
    '# (every cluster column / its own top class, separately for low and high), so low\n'
    '# and high each get their own colour scale and colorbar.')

# dead commented export examples
sub('''# export_figure("fig4_tuning_by_cluster")   # -> notebooks/outputs/fig4_tuning_by_cluster.{png,pdf,svg}
# export_figure("fig4_tuning_by_cluster", formats=("pdf",), outdir="~/Desktop")''',
    '''# export_figure("fig4_tuning_by_cluster")   # -> notebooks/outputs/fig4_tuning_by_cluster.{png,pdf,svg}''')

OUT.write_text(src, encoding="utf-8", newline="\n")
print("wrote", OUT, len(src.splitlines()), "lines")
