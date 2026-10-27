"""Non-responders folded in: how the clusters carry over (Fig 4E).

Analysis only -- numbers in, numbers out, nothing drawn and nothing refitted. The reclustering
itself is a fit and was done once by `dev/fold_nonresponders.py`, which cached each glomerulus's
new cluster label in `data/nonresponders_folded/glom_clusters.csv`. This module reads that cache,
lines it up with the original clustering, and counts where each glomerulus went.

Matching. A glomerulus is identified by ``(subject, gloms)``. Those in the original
("nonresponders_dropped") clustering carry their house number 1-7 (`Dataset.cluster_number`).
The rest are the *recruited non-responders* that folding added: they have no original cluster, so
they are counted in an extra "not in the original clustering" row.

Numbering the new clusters. A fresh NMF run hands out arbitrary labels 0-6. They are renamed 1-7
by best one-to-one overlap with the original clusters (Hungarian assignment on the matched
glomeruli), so new cluster k is the one that most resembles original cluster k and a cluster that
simply carries over reads as a straight band. This is a naming choice, not a result: the counts are
the same under any numbering.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import clusters as _clusters

__all__ = ["FoldedFlows", "CACHE_PATH", "load_cache", "flows"]

CACHE_PATH = Path(__file__).resolve().parent.parent / "data" / "nonresponders_folded" / "glom_clusters.csv"


def load_cache(path: str | Path | None = None) -> pd.DataFrame:
    """The cached new cluster label per glomerulus. Raises if it is missing.

    It is never recomputed here: a silent refit could give different labels from the ones the
    manuscript shows. Rebuild it deliberately with ``python dev/fold_nonresponders.py``.
    """
    path = Path(path) if path is not None else CACHE_PATH
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing. The folded reclustering is cached, not recomputed at draw time; "
            "rebuild it with `python dev/fold_nonresponders.py`.")
    table = pd.read_csv(path, dtype={"subject": str, "gloms": str})
    need = {"subject", "gloms", "new_label", "recruited"}
    if not need <= set(table.columns):
        raise ValueError(f"{path} lacks columns {sorted(need - set(table.columns))}")
    return table


@dataclass(frozen=True)
class FoldedFlows:
    """Where each glomerulus went when the recruited non-responders were folded in.

    Attributes
    ----------
    counts : ndarray, shape (8, 7)
        ``counts[i, j]`` = glomeruli from original cluster ``i + 1`` that are in new cluster
        ``j + 1``. Row 7 is the glomeruli that were not in the original clustering at all.
    label_to_number : dict[int, int]
        Raw NMF label -> the 1-7 number it was given by best overlap with the original clusters.
    """

    counts: np.ndarray
    label_to_number: dict

    @property
    def n_clusters(self) -> int:
        return self.counts.shape[1]

    @property
    def n_total(self) -> int:
        return int(self.counts.sum())

    @property
    def n_matched(self) -> int:
        """Glomeruli that were in the original clustering."""
        return int(self.counts[:-1].sum())

    @property
    def n_added(self) -> int:
        """Glomeruli added by folding -- the recruited non-responders."""
        return int(self.counts[-1].sum())

    @property
    def original_sizes(self) -> np.ndarray:
        return self.counts[:-1].sum(axis=1)

    @property
    def new_sizes(self) -> np.ndarray:
        return self.counts.sum(axis=0)

    @property
    def stayed(self) -> np.ndarray:
        """Per original cluster: the fraction that stays in its best-matched new cluster."""
        k = self.n_clusters
        return np.diag(self.counts[:k]) / np.maximum(self.original_sizes, 1)

    def __repr__(self) -> str:
        return (f"FoldedFlows({self.n_matched} original + {self.n_added} added = "
                f"{self.n_total} glomeruli, {self.n_clusters} clusters)")


def flows(d, table: pd.DataFrame | None = None) -> FoldedFlows:
    """Count the flow from the original clusters to the folded-in clusters.

    `d` is the `highconc.data.Dataset` of the original ("nonresponders_dropped") clustering;
    `table` is `load_cache()` (loaded here if omitted).
    """
    from scipy.optimize import linear_sum_assignment

    table = load_cache() if table is None else table
    cols = d.low.columns
    key_old = list(zip(np.asarray(cols.get_level_values("subject")).astype(str),
                       np.asarray(cols.get_level_values("gloms")).astype(str)))
    if len(set(key_old)) != len(key_old):
        raise ValueError("duplicate (subject, gloms) in the original clustering")
    key_new = list(zip(table["subject"].astype(str), table["gloms"].astype(str)))
    if len(set(key_new)) != len(key_new):
        raise ValueError("duplicate (subject, gloms) in the folded-clustering cache")
    missing = set(key_old) - set(key_new)
    if missing:
        raise ValueError(
            f"{len(missing)} glomeruli of the original clustering are absent from the cache "
            "(was it built from the same data?)")

    k = _clusters.N_CLUSTERS
    old_number = dict(zip(key_old, np.asarray(d.cluster_number, dtype=int)))
    src = np.array([old_number.get(key, 0) for key in key_new])          # 0 = not in the original
    labels = np.sort(table["new_label"].unique())
    if len(labels) != k:
        raise ValueError(f"expected {k} new clusters, the cache has {len(labels)}")
    col = {int(lab): i for i, lab in enumerate(labels)}
    new_col = table["new_label"].map(col).to_numpy()

    overlap = np.zeros((k, k), dtype=int)                                # original x new, matched only
    for s, c in zip(src, new_col):
        if s:
            overlap[s - 1, c] += 1
    rows, cols_ = linear_sum_assignment(-overlap)                        # original i <-> new column
    number_of_col = {int(c): int(r) + 1 for r, c in zip(rows, cols_)}

    counts = np.zeros((k + 1, k), dtype=int)
    for s, c in zip(src, new_col):
        counts[(s - 1) if s else k, number_of_col[int(c)] - 1] += 1
    return FoldedFlows(counts=counts,
                       label_to_number={int(lab): number_of_col[col[int(lab)]] for lab in labels})
