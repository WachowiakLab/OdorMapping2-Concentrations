"""Assemble final's Fig 3 code section: Copy's code + the user's chosen changes."""
from pathlib import Path

HERE = Path(__file__).parent
src = (HERE / "sec" / "03_copy_fig_code.py").read_text(encoding="utf-8")
OUT = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\notes\sections\03_fig3_code.py")


def sub(old, new):
    global src
    assert src.count(old) == 1, (src.count(old), old[:60])
    src = src.replace(old, new)


# order 'matt', 'matt'
sub('WITHIN_CLASS_ORDER = "burton"  # "table" | "matt" (pyrazine only) | "burton" (full order)',
    'WITHIN_CLASS_ORDER = "matt"    # "table" | "matt" (pyrazine only) | "burton" (full order)')

# 4th example glomerulus (from lab), kept out of EG_GLOMERULI so Fig 6 is unaffected
sub('EG_GLOMERULI = similarity.EG_GLOMERULI     # Fig 6B-D must reuse these\n',
    'EG_GLOMERULI = similarity.EG_GLOMERULI     # Fig 6B-D must reuse these\n'
    '# Fourth example glomerulus (panel H). Not part of EG_GLOMERULI, so Fig 6 is unaffected.\n'
    '# Position within one mouse, counted from 0 among the glomeruli that survive the\n'
    '# `nonresponders_dropped` filter.\n'
    'EXTRA_SUBJECT = 1\n'
    'EXTRA_GLOM_IN_SUBJECT = 8\n'
    'EXTRA_GLOM = int(np.flatnonzero(d.subject == EXTRA_SUBJECT)[EXTRA_GLOM_IN_SUBJECT])\n')
sub('spectra = similarity.eg_spectra(d, EG_GLOMERULI, order=CLASS_ORDER_FIG3,\n'
    '                                within=WITHIN_CLASS_ORDER)',
    'spectra = similarity.eg_spectra(d, [*EG_GLOMERULI, EXTRA_GLOM],\n'
    '                                order=CLASS_ORDER_FIG3, within=WITHIN_CLASS_ORDER)')

# recruitment matrix moves H -> I
sub("# Panel H's matrix is computed here", "# Panel I's matrix is computed here")
sub('print(f"panel H: rows built', 'print(f"panel I: rows built')
sub('hs.panel_label(axH, "H", x=-0.16, y=1.05)', 'hs.panel_label(axH, "I", x=-0.16, y=1.05)')

# chemical-class labels on every matrix, not only A
sub('    # The x-axis chemical-class strip is drawn on every panel here (all four\n'
    '    # share the same class order), but its text names only on A -- repeating\n'
    '    # "acids, aldehydes, ..." four times across one row is pure clutter.\n',
    '    # The x-axis chemical-class strip, with its family names, is drawn on all\n'
    '    # four panels (they share the same class order).\n')
sub('boundaries=False, strip_labels=(L == "A"))', 'boundaries=False, strip_labels=True)')

# E-H: four spectra
sub('# E/F/G keep their original 3-way split, nested inside the outer row\'s full\n'
    '# span rather than 3 of the top row\'s 4 columns -- so adding panel A doesn\'t\n'
    '# also squeeze E/F/G.\n'
    'gs_efg = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1, :], wspace=0.30)',
    '# E-H: one column per spectrum, nested inside the outer row\'s full span rather\n'
    '# than the top row\'s 4 columns.\n'
    'gs_efg = GridSpecFromSubplotSpec(1, len(spectra), subplot_spec=gs[1, :],\n'
    '                                 wspace=0.30)')
sub('for j, (s, L) in enumerate(zip(spectra, "EFG")):', 'for j, (s, L) in enumerate(zip(spectra, "EFGH")):')

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(src, encoding="utf-8", newline="\n")
print("wrote", OUT, len(src.splitlines()), "lines")
