"""Fold the recruited non-responders into the low matrix, recluster once, and cache the labels.

Fig 4E (the Sankey) compares each glomerulus's original low-concentration cluster with the
cluster it lands in when the non-responders that are *recruited* at high concentration are
folded into the low matrix and NMF is run again. That reclustering is a fit, and the notebook
never refits at draw time -- so it is done here, once, seeded, and the labels are cached:

    data/nonresponders_folded/glom_clusters.csv    subject, gloms, new_label, recruited
    data/nonresponders_folded/meta.json            seed, settings and counts of this run

    python dev/fold_nonresponders.py               # (re)build the cache
    python dev/fold_nonresponders.py --check       # print what a rebuild would give; write nothing

What "folded in" means. A glomerulus is a *recruited non-responder* if it is silent at low
concentration (zero for every odorant) but responds to at least one odorant at high. Its HIGH
response column replaces its all-zero LOW column, so it can be clustered with everything else.
Glomeruli that are flat after folding (no response at either concentration) cannot be clustered
and are left out, exactly as the original analysis dropped them.

Seeding. The clustering uses `glom_explorer.assign_cluster_indices` with its own defaults (the
package's NMF seed and iteration cap, 7 components) -- the same settings that made the original
"nonresponders_dropped" clustering. Passing `params={"n_components": 7}` instead would silently
discard the seed and the iteration cap and give different labels on every run.

Needs `glom_explorer` with `load_processed_data` and `assign_cluster_indices`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "data" / "nonresponders_folded"
N_CLUSTERS = 7


def fold(low_raw: pd.DataFrame, high_raw: pd.DataFrame):
    """Return ``(folded_low, recruited)``: the folded matrix and the recruited-glomerulus mask."""
    silent_at_low = (low_raw == 0).all(axis=0)
    responds_at_high = (high_raw != 0).any(axis=0)
    recruited = (silent_at_low & responds_at_high).to_numpy()
    folded = low_raw.copy()
    folded.loc[:, recruited] = high_raw.loc[:, recruited]
    return folded, recruited


def build() -> tuple[pd.DataFrame, dict]:
    import glom_explorer as gx

    _, low_raw, high_raw = gx.load_processed_data()
    folded, recruited = fold(low_raw, high_raw)

    flat = (folded.std() == 0).to_numpy()
    cluster_on = folded.loc[:, ~flat]
    # params=None -> the package's seeded defaults (n_clusters, nmf_seed, nmf_max_iter)
    _, col_labels = gx.assign_cluster_indices(cluster_on, method="nmf")

    cols = cluster_on.columns
    table = pd.DataFrame({
        "subject": np.asarray(cols.get_level_values("subject")).astype(str),
        "gloms": np.asarray(cols.get_level_values("gloms")).astype(str),
        "new_label": np.asarray(col_labels, dtype=int),
        "recruited": recruited[~flat],
    })
    defaults = getattr(getattr(gx, "config", None), "ANALYSIS_DEFAULTS", {}) or {}
    meta = {
        "n_clusters": int(len(set(col_labels))),
        "nmf_seed": defaults.get("nmf_seed"),
        "nmf_max_iter": defaults.get("nmf_max_iter"),
        "n_glomeruli_total": int(low_raw.shape[1]),
        "n_recruited_nonresponders": int(recruited.sum()),
        "n_flat_left_out": int(flat.sum()),
        "n_clustered": int(table.shape[0]),
        "rule": "silent at low (all zero) and responding at high: HIGH column replaces the LOW column",
        "cluster_sizes_raw_label": {int(k): int(v) for k, v in
                                    table["new_label"].value_counts().sort_index().items()},
        "numpy": np.__version__,
    }
    try:
        import sklearn
        meta["scikit_learn"] = sklearn.__version__
    except Exception:
        pass
    return table, meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="print the result; write nothing")
    args = ap.parse_args()
    table, meta = build()
    print(json.dumps(meta, indent=2))
    if args.check:
        return 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT_DIR / "glom_clusters.csv", index=False)
    (OUT_DIR / "meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {OUT_DIR.relative_to(REPO)}/glom_clusters.csv ({len(table)} rows) and meta.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
