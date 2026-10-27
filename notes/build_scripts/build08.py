"""S1: main's cell, 1x block only (panels A-C)."""
from pathlib import Path

HERE = Path(__file__).parent
src = (HERE / "sec" / "08_main_fig_code.py").read_text(encoding="utf-8")
OUT = Path(r"C:\Users\u6068031\Documents\highconc-merge\final\notes\sections\08_figS1_code.py")


def sub(old, new):
    global src
    assert src.count(old) == 1, (src.count(old), old[:70])
    src = src.replace(old, new)


sub('''    ("1x",  FIGURES / "nmf_dim_choice_1x" / "data" / "cache",  [7, 9], 7, "ABC"),
    ("25x", FIGURES / "nmf_dim_choice_25x" / "data" / "cache", [5, 7], 5, "DEF"),
]''',
    '''    ("1x",  FIGURES / "nmf_dim_choice_1x" / "data" / "cache",  [7, 9], 7, "ABC"),
    # The 25x block (panels D-F) is no longer drawn. Its sweep cache is still on
    # disk under figures/nmf_dim_choice_25x/ and is simply not read here.
]''')
sub('''# Scaled from the original 9.0 x 12.5in hand layout so the width is exactly
# SUPP_TARGET_W, proportions held.
SUPP_DIM_SCALE = SUPP_TARGET_W / 9.0
SUPP_DIM_W, SUPP_DIM_H, SUPP_DIM_DPI = 9.0 * SUPP_DIM_SCALE, 12.5 * SUPP_DIM_SCALE, 300''',
    '''# One block now, so the 12.5in-tall two-block hand layout is cut to 5.5in; the
# whole thing is scaled so the width is exactly SUPP_TARGET_W, proportions held.
SUPP_DIM_SCALE = SUPP_TARGET_W / 9.0
SUPP_DIM_W, SUPP_DIM_H, SUPP_DIM_DPI = 9.0 * SUPP_DIM_SCALE, 5.5 * SUPP_DIM_SCALE, 300''')
sub('''# One block per concentration. Inside a block the layout is the retired''',
    '''# One block (the 1x concentration). Inside a block the layout is the retired''')
sub('gs_rows = fig.add_gridspec(2, 1, hspace=0.32)',
    'gs_rows = fig.add_gridspec(len(SUPP_DIM_ROWS), 1, hspace=0.32)')
sub('''# export_figure("supp_nmf_dim_choice")   # -> notebooks/outputs/supp_nmf_dim_choice.{png,pdf,svg}
# export_figure("supp_nmf_dim_choice", formats=("pdf",), outdir="~/Desktop")''',
    '''# export_figure("supp_nmf_dim_choice")   # -> notebooks/outputs/supp_nmf_dim_choice.{png,pdf,svg}''')

OUT.write_text(src, encoding="utf-8", newline="\n")
print("wrote", OUT, len(src.splitlines()), "lines")
