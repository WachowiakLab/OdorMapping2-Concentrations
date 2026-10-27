"""The NMF dimensionality sweep — compute, cache, and the containers between.

Promoted from ``figures/_shared/nmf_dim_sweep.py`` so the notebook can reach it
without a ``sys.path`` insert. Two halves:

*Compute* (:func:`run_dim_sweep`, :func:`run_consensus_at_k`) refits NMF at
k = 1..20 across 50 seeds and is **slow — minutes per concentration**.

*Cache* (:func:`load_sweep`, :func:`load_consensus`, :func:`load_cached`) reads
the ``.npz`` files that compute already wrote, in milliseconds.

The notebook must only ever call the cache half::

    from highconc import dimsweep

    sweep, consensus = dimsweep.load_cached(
        "figures/nmf_dim_choice_1x/data/cache", conc="1x", ks=(7, 9)
    )

:func:`load_cached` raises if a cache is missing rather than quietly falling
back to computing, because a notebook that silently recomputes the sweep looks
identical to one that has hung. Regenerating a cache is a deliberate act: run
the figure's ``generate.py`` with ``USE_CACHE = False``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

__all__ = [
    "DimSweepResult", "ConsensusResult",
    "N_NMF_SEEDS", "N_SHUFFLES", "K_VALUES",
    "run_dim_sweep", "run_consensus_at_k",
    "save_sweep", "load_sweep", "save_consensus", "load_consensus",
    "sweep_cache_path", "consensus_cache_path", "load_cached",
]

#: Sweep parameters the cached manuscript figures were computed with. Changing
#: one of these means the caches on disk no longer describe what you asked for,
#: so recompute rather than mixing them.
K_VALUES = list(range(1, 21))
N_NMF_SEEDS = 50
N_SHUFFLES = 20


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------
@dataclass
class DimSweepResult:
    """Per-k sweep results for variance-explained / slope / stability.

    The ``_std`` fields are standard deviations across seeds (NMF) or across
    permutations (shuffled PCA). The panels draw standard *error* bands, so
    they divide by ``sqrt(n)`` — which is why the seed counts are panel
    arguments and not baked into this container: the cache files predate them.
    """

    k_values: np.ndarray                # shape (K,)
    pca_var_explained: np.ndarray       # shape (K,)
    pca_shuffled_var_mean: np.ndarray   # shape (K,)
    pca_shuffled_var_std: np.ndarray    # shape (K,)
    nmf_var_explained_mean: np.ndarray  # shape (K,)
    nmf_var_explained_std: np.ndarray   # shape (K,)
    nmf_ari_mean: np.ndarray            # shape (K,)
    nmf_ari_std: np.ndarray             # shape (K,)


@dataclass
class ConsensusResult:
    """Consensus matrix at a single k: how often each pair co-clusters."""

    k: int
    consensus: np.ndarray               # shape (n_samples, n_samples)


# ---------------------------------------------------------------------------
# Compute — slow. Only run to regenerate a cache.
# ---------------------------------------------------------------------------
def run_dim_sweep(
    X: np.ndarray,
    k_values: Sequence[int] = tuple(K_VALUES),
    n_nmf_seeds: int = N_NMF_SEEDS,
    n_shuffle: int = N_SHUFFLES,
    rng: np.random.Generator | None = None,  # unused; glom-explorer owns its RNG
) -> DimSweepResult:
    """Sweep PCA, shuffled-PCA, and NMF across ``k_values``.

    Variance-explained curves come from one call to
    ``gx.run_reconstruction_analysis``. NMF clustering stability (adjusted Rand
    index across initialisations) iterates ``gx.consensus_clustering_stability``
    over k — that loop is the slow part.
    """
    import glom_explorer as gx

    K = np.asarray(list(k_values))
    max_k = int(K.max())

    # gx returns a dict of DataFrames each shaped (max_k, n_iterations). PCA is
    # deterministic so its iteration columns coincide; NMF and shuffled vary.
    recon = gx.run_reconstruction_analysis(
        X, max_components=max_k, n_iterations=n_shuffle
    )
    pca_arr = recon["PCA"].values
    shuf_arr = recon["PCA (Shuffled)"].values
    nmf_arr = recon["NMF"].values

    # gx indexes k from 1, so k maps to row k-1.
    pca_ve = pca_arr[K - 1, 0]
    pca_shuf_mean = shuf_arr[K - 1, :].mean(axis=1)
    pca_shuf_std = shuf_arr[K - 1, :].std(axis=1)
    nmf_ve_mean = nmf_arr[K - 1, :].mean(axis=1)
    nmf_ve_std = nmf_arr[K - 1, :].std(axis=1)

    nmf_ari_mean = np.zeros_like(K, dtype=float)
    nmf_ari_std = np.zeros_like(K, dtype=float)
    for i, k in enumerate(K):
        if k < 2:
            nmf_ari_mean[i] = 1.0  # one cluster is trivially stable
            nmf_ari_std[i] = 0.0
            continue
        out = gx.consensus_clustering_stability(
            X, n_components=int(k), n_runs=n_nmf_seeds
        )
        nmf_ari_mean[i] = float(out["mean_ari_to_consensus"])
        nmf_ari_std[i] = float(out["std_ari_to_consensus"])

    return DimSweepResult(
        k_values=K,
        pca_var_explained=pca_ve,
        pca_shuffled_var_mean=pca_shuf_mean,
        pca_shuffled_var_std=pca_shuf_std,
        nmf_var_explained_mean=nmf_ve_mean,
        nmf_var_explained_std=nmf_ve_std,
        nmf_ari_mean=nmf_ari_mean,
        nmf_ari_std=nmf_ari_std,
    )


def run_consensus_at_k(
    X: np.ndarray,
    k: int,
    n_seeds: int = N_NMF_SEEDS,
    rng: np.random.Generator | None = None,  # unused
) -> ConsensusResult:
    """Consensus matrix at a single k. Slow; cache the result."""
    import glom_explorer as gx

    out = gx.consensus_clustering_stability(
        X, n_components=int(k), n_runs=n_seeds
    )
    return ConsensusResult(k=int(k), consensus=np.asarray(out["consensus_matrix"]))


# ---------------------------------------------------------------------------
# Cache — fast. This is the path the notebook takes.
# ---------------------------------------------------------------------------
def save_sweep(result: DimSweepResult, path) -> None:
    np.savez(
        path,
        k_values=result.k_values,
        pca_var_explained=result.pca_var_explained,
        pca_shuffled_var_mean=result.pca_shuffled_var_mean,
        pca_shuffled_var_std=result.pca_shuffled_var_std,
        nmf_var_explained_mean=result.nmf_var_explained_mean,
        nmf_var_explained_std=result.nmf_var_explained_std,
        nmf_ari_mean=result.nmf_ari_mean,
        nmf_ari_std=result.nmf_ari_std,
    )


def load_sweep(path) -> DimSweepResult:
    d = np.load(path)
    return DimSweepResult(**{k: d[k] for k in d.files})


def save_consensus(result: ConsensusResult, path) -> None:
    np.savez(path, k=result.k, consensus=result.consensus)


def load_consensus(path) -> ConsensusResult:
    d = np.load(path)
    return ConsensusResult(k=int(d["k"]), consensus=d["consensus"])


def sweep_cache_path(cache_dir, conc: str) -> Path:
    """Where the sweep cache for a concentration lives. ``conc`` is '1x'/'25x'."""
    return Path(cache_dir) / f"sweep_{conc}.npz"


def consensus_cache_path(cache_dir, conc: str, k: int) -> Path:
    return Path(cache_dir) / f"consensus_{conc}_k{int(k)}.npz"


def load_cached(cache_dir, *, conc: str, ks: Sequence[int]):
    """Read a concentration's cached sweep and consensus matrices.

    Parameters
    ----------
    cache_dir : path
        The figure folder's ``data/cache``. The caches stay with the figure
        that owns them; nothing copies them into the package.
    conc : {"1x", "25x"}
        Which concentration's caches to read.
    ks : sequence of int
        The k values whose consensus matrices the figure shows — ``(7, 9)`` at
        1x, ``(5, 7)`` at 25x.

    Returns
    -------
    (DimSweepResult, dict[int, ConsensusResult])

    Raises
    ------
    FileNotFoundError
        If any cache is absent. Deliberate: recomputing takes minutes, and a
        notebook that starts a sweep behind your back is indistinguishable from
        one that has frozen. Regenerate with the figure's ``generate.py`` and
        ``USE_CACHE = False``.
    """
    cache_dir = Path(cache_dir)
    sweep_path = sweep_cache_path(cache_dir, conc)
    if not sweep_path.exists():
        raise FileNotFoundError(
            f"no sweep cache at {sweep_path}. Recomputing takes minutes, so "
            f"this is not done implicitly — run figures/nmf_dim_choice_{conc}/"
            f"generate.py with USE_CACHE = False to rebuild it."
        )
    consensus = {}
    for k in ks:
        path = consensus_cache_path(cache_dir, conc, k)
        if not path.exists():
            raise FileNotFoundError(
                f"no consensus cache at {path}. Rebuild it the same way as "
                f"the sweep cache — see the message above."
            )
        consensus[int(k)] = load_consensus(path)
    return load_sweep(sweep_path), consensus
