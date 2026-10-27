# %%
# Own constants — this section inherits no bare name from anywhere.
# Only panel C of the old three-panel figure remains: the parity of every
# responding glomerulus x odorant pair, observed vs held-out predicted.
SUPP_PAR_C_CMAP = "Purples"          # hexbin colormap
SUPP_PAR_W, SUPP_PAR_H, SUPP_PAR_DPI = SUPP_TARGET_W / 2, 2.4, 200

# The smallest nonzero observed high response in the whole dataset (log1p
# space, which is this panel's space): the dotted horizontal line. A predicted
# value under it is lower than any response the experiment ever measured.
supp_par_floor = pls.observed_high_floor(m, space="log")

fig, ax = plt.subplots(1, 1, figsize=(SUPP_PAR_W, SUPP_PAR_H),
                       dpi=SUPP_PAR_DPI, constrained_layout=True)
# Bare call = the pinned framing at its defaults, which is the draft's panel.
supp_par_hb = P.pls.responding_parity(
    ax, m, cbar=True, title="all responding pairs", title_fontsize=11.5,
    cmap=SUPP_PAR_C_CMAP, floor=supp_par_floor)

print(f"responding pairs   Q2 = {pls.baseline_q2(m):.3f}  "
      f"(n = {pls.responding_pairs(m)[0].size:,}); "
      "the draft page labels this same value R2")
print(f"dotted line at the smallest nonzero observed response: "
      f"{supp_par_floor:.4f} (log1p dF/F)")

save(fig, "supp_parity", dpi=SUPP_PAR_DPI)
plt.show()
# export_figure("supp_parity")   # -> notebooks/outputs/supp_parity.{png,pdf,svg}
