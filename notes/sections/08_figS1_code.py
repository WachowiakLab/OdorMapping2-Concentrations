# %%
# Own constants — this section inherits no bare name from anywhere.
SUPP_DIM_K_VALUES = list(range(1, 21))   # the sweep on disk; changing this invalidates it
SUPP_DIM_N_SEEDS, SUPP_DIM_N_SHUFFLES = 50, 20
SUPP_DIM_ROWS = [
    # (conc, cache dir, consensus k's, chosen k, panel letters)
    ("1x",  FIGURES / "nmf_dim_choice_1x" / "data" / "cache",  [7, 9], 7, "ABC"),
    # The 25x block (panels D-F) is no longer drawn. Its sweep cache is still on
    # disk under figures/nmf_dim_choice_25x/ and is simply not read here.
]
# One block now, so the 12.5in-tall two-block hand layout is cut to 5.5in; the
# whole thing is scaled so the width is exactly SUPP_TARGET_W, proportions held.
SUPP_DIM_SCALE = SUPP_TARGET_W / 9.0
SUPP_DIM_W, SUPP_DIM_H, SUPP_DIM_DPI = 9.0 * SUPP_DIM_SCALE, 5.5 * SUPP_DIM_SCALE, 300

fig = plt.figure(figsize=(SUPP_DIM_W, SUPP_DIM_H), dpi=SUPP_DIM_DPI)
# One block (the 1x concentration). Inside a block the layout is the retired
# per-concentration figure's: the
# two line panels on top, then a bottom strip holding the ARI curve and the
# stacked consensus matrices, whose height matches the original 7 x 1.75-inch
# rectangle.
gs_rows = fig.add_gridspec(len(SUPP_DIM_ROWS), 1, hspace=0.32)

for row, (conc, cache_dir, consensus_ks, chosen, letters) in enumerate(SUPP_DIM_ROWS):
    sweep, consensus = dimsweep.load_cached(cache_dir, conc=conc, ks=consensus_ks)

    gs_block = gs_rows[row].subgridspec(2, 1, height_ratios=[2, 1], hspace=0.40)
    gs_top = gs_block[0].subgridspec(1, 2, wspace=0.40)
    gs_bot = gs_block[1].subgridspec(1, 6, width_ratios=[0.65] * 5 + [1.0],
                                     wspace=0.20)

    ax_var = fig.add_subplot(gs_top[0, 0])
    P.dimchoice.variance_explained(sweep, ax_var, n_nmf_seeds=SUPP_DIM_N_SEEDS,
                                   n_shuffles=SUPP_DIM_N_SHUFFLES)
    ax_slope = fig.add_subplot(gs_top[0, 1])
    P.dimchoice.variance_slope(sweep, ax_slope)
    # The ARI curve stops before the consensus column so it does not run under
    # the matrices.
    ax_ari = fig.add_subplot(gs_bot[0, 0:5])
    P.dimchoice.stability(sweep, ax_ari, n_nmf_seeds=SUPP_DIM_N_SEEDS)

    inner = GridSpecFromSubplotSpec(len(consensus_ks), 1,
                                    subplot_spec=gs_bot[0, 5], hspace=0.30)
    for i, k in enumerate(consensus_ks):
        P.dimchoice.consensus(consensus[k], fig.add_subplot(inner[i, 0]))

    for ax, letter in zip((ax_var, ax_slope, ax_ari), letters):
        hs.panel_label(ax, letter)
    ax_var.set_title(f"{conc} concentration — chosen k = {chosen}",
                     fontsize=hs.FS_PANEL, loc="left")

    print(f"{conc}: k = {sweep.k_values.tolist()} | chosen k = {chosen} | "
          f"consensus shown at {consensus_ks} | panels {'/'.join(letters)}")

save(fig, "supp_nmf_dim_choice", dpi=SUPP_DIM_DPI)
plt.show()
# export_figure("supp_nmf_dim_choice")   # -> notebooks/outputs/supp_nmf_dim_choice.{png,pdf,svg}
