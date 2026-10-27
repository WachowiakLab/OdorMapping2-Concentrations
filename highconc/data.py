"""One call, one object: the data every panel draws from.

    from highconc.data import load
    d = load()
    d.low.shape, d.high.shape        # (59, 376) each
    d.cluster_index                  # per-glomerulus component index 0-6
    d.cluster_number                 # the same gloms as Joel's 1-7
    d.xy                             # (376, 2) dorsal-bulb coordinates

Panels take a `Dataset` plus an `ax`; they never load anything themselves.
That is what lets the notebook load once and draw twelve figures.

The underlying data ships inside the `glom-explorer` package, so there is no
external path to configure. `load()` caches per variant — calling it in every
notebook cell is free.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd

from . import clusters as _clusters

#: The clustering variant the manuscript uses.
DEFAULT_VARIANT = "nonresponders_dropped"


@dataclass(frozen=True)
class Dataset:
    """Paired low/high-concentration glomerular responses, plus their labels.

    Attributes
    ----------
    low, high : DataFrame, 59 odorants x 376 glomeruli
        Response matrices at 1x and 25x. Identical row and column order, so a
        glomerulus means the same thing in both.
    odorants : DataFrame
        One row per odorant: Odorant, SMILES, Name, Group, Group_num. `Group`
        is the chemical class (C acids, D aldehydes, E esters, F furans,
        K ketones, L alcohols, N amines, P pyrazines, R aromatics,
        S thiazoles).
    cluster_index : ndarray, shape (376,)
        Our NMF component index 0-6 for each glomerulus.
    cluster_number : ndarray, shape (376,)
        The same assignment in Joel's 1-7 numbering. Prefer this for display.
    subject : ndarray, shape (376,)
        Mouse each glomerulus came from.
    xy : ndarray, shape (376, 2)
        Dorsal-bulb coordinates.
    variant : str
        Which clustering variant was loaded.
    """

    low: pd.DataFrame
    high: pd.DataFrame
    odorants: pd.DataFrame
    cluster_index: np.ndarray
    cluster_number: np.ndarray
    subject: np.ndarray
    xy: np.ndarray
    variant: str

    @property
    def n_gloms(self) -> int:
        return self.low.shape[1]

    @property
    def n_odorants(self) -> int:
        return self.low.shape[0]

    @property
    def odorant_cluster(self) -> np.ndarray:
        """Bi-cluster row label per odorant, in row order.

        Pairs with `cluster_index`: odorant cluster k is the block of odorants
        that glomerular component k responds to, which is what puts the
        bi-clustered heatmap's highlighted blocks on the diagonal.
        """
        return np.asarray(
            self.low.index.get_level_values("odorant_cluster"), dtype=int
        )

    @property
    def chem_class(self) -> np.ndarray:
        """Chemical-class letter per odorant, in row order."""
        return self.odorants["Group"].to_numpy()

    @property
    def subjects(self) -> list:
        """The four mice, sorted — the order per-subject figures lay them out in."""
        return sorted(np.unique(self.subject).tolist())

    def subject_mask(self, subject) -> np.ndarray:
        """Boolean column mask for one mouse."""
        mask = self.subject == subject
        if not mask.any():
            raise ValueError(
                f"no glomeruli for subject {subject!r}; have {self.subjects}"
            )
        return mask

    def cluster_mask(self, joel: int) -> np.ndarray:
        """Boolean column mask for one of Joel's clusters, 1-7."""
        return self.cluster_number == _clusters._check(
            joel, 1, 7, "Joel cluster number"
        )

    def __repr__(self) -> str:  # keep notebook output short
        return (f"Dataset(variant={self.variant!r}, "
                f"{self.n_odorants} odorants x {self.n_gloms} gloms, "
                f"{_clusters.N_CLUSTERS} clusters)")


@lru_cache(maxsize=4)
def load(variant: str = DEFAULT_VARIANT) -> Dataset:
    """Load one clustering variant. Cached, so repeated calls are free."""
    import glom_explorer as gx

    low, high = gx.load_clustered_variant(variant)
    odorants = gx.load_odorants()

    cols = low.columns
    cluster_index = np.asarray(cols.get_level_values("parent_cluster"), dtype=int)
    cluster_number = np.array(
        [_clusters.to_joel(i) for i in cluster_index], dtype=int
    )
    subject = np.asarray(cols.get_level_values("subject"))
    xy = np.column_stack([
        np.asarray(cols.get_level_values("x"), dtype=float),
        np.asarray(cols.get_level_values("y"), dtype=float),
    ])

    if not low.columns.equals(high.columns):
        raise ValueError("low and high disagree on glomerulus order")
    if not low.index.equals(high.index):
        raise ValueError("low and high disagree on odorant order")

    return Dataset(
        low=low, high=high, odorants=odorants,
        cluster_index=cluster_index, cluster_number=cluster_number,
        subject=subject, xy=xy, variant=variant,
    )
