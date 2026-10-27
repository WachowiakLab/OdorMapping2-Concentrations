import sys
sys.path.insert(0, __file__.rsplit("\\", 1)[0])
from cap_lib import load, save, md


def put(cells, tag, prefix, text):
    hits = [i for i, c in enumerate(cells) if c.startswith(f'# %% [markdown] tags=["{tag}"]')
            and c.split("\n")[1].startswith(prefix)]
    assert len(hits) == 1, (tag, prefix, hits)
    cells[hits[0]] = md(tag, text)


S8_HUMAN = r"""
## 8. Figure 8 — directional transfer of tuning between groups

**The test.** Predict a group's high-concentration responses from its own low-concentration average — that is the
baseline. Then add a second group's low-concentration average and see whether the prediction improves. Repeat for
every ordered pair, in both directions. Prediction is scored on held-out odorants, so an improvement is not simply
the result of adding a parameter.

- **B** — the permutation *p* for each pair, against 10,000 shuffles. A **green box** marks the pairs that survive
  Benjamini–Hochberg correction at q < 0.05, and only those boxed squares carry a number.
- **C** — the size of each improvement, as a percentage of the error the target's own average left behind. Most
  cells are near zero; again only the boxed pairs carry a number.

Rows are the **source** group whose low-concentration average is added ("added non-source cluster"); columns are the
**target** group being predicted at high concentration. The matrix is not symmetric and is not meant to be.

**A control that is no longer drawn.** The same test can group glomeruli by the chemical class each one prefers
instead of by the NMF clusters, which does not use the clustering at all. That analysis is still run, and the
cells after the figure print the cluster edges next to their chemical-class counterparts (seven of the eight
significant cluster results have one) — but the paper's Figure 8 shows only the cluster version, so only that is
drawn.

**The explanatory header in the paper (panel A) is drawn elsewhere.**
"""

S8_LLM = r"""
## 8. Figure 8 — motif transfer between glomerular groups

**What it shows.** The direct test that replaces reading the off-diagonal weights of `W`. For every *ordered* pair
of glomerular groups `F → T`: does adding group **F**'s mean low-concentration response to the prediction of group
**T**'s glomeruli at high improve held-out performance, beyond what **T**'s own low-concentration mean already
achieves?

    intra    = predict T's high from  [T's low mean]
    intra+F  = predict T's high from  [T's low mean, F's low mean]
    ΔQ²      = Q²(intra+F) − Q²(intra)        # out-of-fold, pooled over T

A positive ΔQ² says a tuning motif carried by F at low concentration reappears *inside* T at high. That is the
manuscript's claim of structured retuning, stated as a prediction rather than as an inspection of weights.

**Why not just read W.** Because `W`'s off-diagonal weights cannot identify which groups inherit whose tuning: the
diagonal and off-diagonal predictors are 94% redundant (`cov(X)` is not diagonal), so the same held-out performance
can be assembled many ways. The pairwise test avoids that by asking whether the *prediction* improves, which is a
property of performance rather than of a particular parameterisation.

**Panels (the paper's letters; its A is the procedure schematic, drawn elsewhere).**

- **B — permutation p per pair**, clipped at `FIG8_CLIP = 0.1`, on the original sequential ramp (OrRd, passed as
  `cmap`; `hs.SEQ` is magma elsewhere). Colour separates the transfers that happened from the ones that did not. The
  **box** is the claim: BH q < `FIG8_ALPHA` (0.05) across all ordered off-diagonal pairs (42 for the clusters), drawn
  in `FIG8_BOX_COLOR` (`#2b9743`, green). **Only the boxed squares carry a number** (`text_significant_only=True`); it
  is the unclipped p, and `<1e-4` means the 10,000 draws cannot resolve it further. Axes: x "cluster (predicted at
  high)", y "added non-source cluster / prediction improvement"; rows are the source, columns the target.
- **C — the size of the transfer**, as a percentage of the error the target's own mean left behind
  (`FIG8_METRIC = "pct_err"`), on the diverging ramp with zero at grey, same boxes and same numbers-only-where-boxed
  rule. Most cells are grey. It has no y label of its own (B's applies). A legend in the same green explains the box.

**The chemical-class control is computed but not drawn.** The same test groups each glomerulus by the class evoking
its largest mean low-concentration response, keeping the nine classes with at least 15 glomeruli (`FIG8_MIN_CLASS`;
372 of 376, the ester-preferring group, n = 4, is excluded). It never saw the NMF factorisation, so if the transfers
survive it they are a property of the population and not of the clustering. The cells after the figure print the
cluster edges beside their class counterparts — **seven of the eight cluster edges have one** — and the per-group
intra baselines, but the paper's Figure 8 has no class panels, so none is drawn (the earlier two-row layout with
its C and D is gone). 72 ordered pairs are tested for classes.

**One group is predicted worse than its own mean.** Furan-preferring glomeruli have an intra baseline of
Q² = −0.084: predicting them from their own low-concentration group mean is *worse* than predicting the mean of the
data, so any percentage of that baseline runs against a negative denominator. It is a real property of the group,
and a second reason `FIG8_METRIC` defaults to `"pct_err"`: `"pct_var"` divides by that baseline and would flip the
sign of every furan cell while inflating its magnitude. (It was most visible in the class row; see the printed
baselines.)

**Read the axes, not the shape.** The matrix is *not* symmetric and is not meant to be: 11 of the 13 class-level
transfers are one-way. The exception is worth reporting rather than smoothing over — **amines ↔ thiazoles is
reciprocal** at class resolution, whereas NMF cluster 5 is a strict sink. The amine → thiazole increment is the
smallest significant transfer in either analysis (+0.013).

**The null, and that it is a choice.** Significance comes from 10,000 *broken-pairing* permutations: the odorant
order of Y is shuffled against X, which preserves the marginal distribution of high responses and all correlation
structure within X and within Y, while destroying the low-to-high correspondence. This is the null the manuscript
publishes. A size-matched control run on 2026-08-20 leaves only three of the eight cluster edges standing and removes
the cluster-5 sink example — a real caveat, recorded in the methods audit, and deliberately not what this figure
shows (`FIG8_NULL`).

**Nothing here recomputes a permutation.** Both nulls are cached in `data/pairwise/` and ship with the clone;
`transfer.pairwise` reads them and recomputes only the intra-group baselines, which are seconds. Rebuilding a cache
is an explicit `transfer.compute_null(...)` call — 11 minutes for the clusters, 16 for the classes serially — and
`pairwise` raises rather than silently starting one.

**Layout and colours.** One row (B, C) scaled to `MAIN_TARGET_W` (7.0 in); the box colour, p-panel ramp and labels
are constants of the cell (`FIG8_BOX_COLOR`, `FIG8_P_CMAP = "OrRd"`, `FIG8_XLABEL`, `FIG8_YLABEL`), not hard-coded in
`panels.transfer`, which takes `box_color`, `cmap` and `text_significant_only`.
"""

cells = load()
put(cells, "human", "# ## 8. Figure 8", "# " + S8_HUMAN.strip("\n"))
put(cells, "llm", "# ## 8. Figure 8", "# " + S8_LLM.strip("\n"))
save(cells)
print("S8 captions replaced")
