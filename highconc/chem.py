"""Chemical similarity between the 59 odorants — Figure 3D.

MACCS structural keys from the SMILES already in the odorant table, compared by
cosine similarity. The result is a 59 x 59 matrix that never changes, so it is
**computed once and cached to disk**:

    data/chem/maccs_cosine.npz

`maccs_cosine()` reads that cache. RDKit is needed only to *rebuild* it, which
is why it lives in an optional extra:

    pip install -e ".[chem]"      # only if you are rebuilding the cache
    python -m highconc.chem --rebuild

Drawing Figure 3D needs neither RDKit nor a network. That is deliberate: the
clustering half of this bundle is Joel's priority, and it must install and run
without a chemistry toolkit.

The cache stores the odorant identifiers alongside the matrix and
`maccs_cosine()` checks them against the dataset it is handed, so a cache built
from a different odorant panel raises instead of silently misaligning 59 rows.

Note for whoever reconciles this with the manuscript: the MS says "106?
structural features". RDKit's MACCS implementation is **167 bits** (166 keys
plus an unused bit 0). We report what we compute; the discrepancy is Matt's
text to fix, not something to paper over by picking a different fingerprint.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

#: 167-bit RDKit MACCS. Recorded in the cache so a future fingerprint change is
#: visible rather than inferred.
FINGERPRINT = "rdkit-MACCS-167"

#: Where the cache lives, relative to the repo root.
CACHE = Path(__file__).resolve().parent.parent / "data" / "chem" / "maccs_cosine.npz"


def _keys(d) -> np.ndarray:
    """The odorant identifiers a cache is keyed on. SMILES, not names: names
    get edited for display, SMILES do not."""
    return np.asarray([str(s) for s in d.odorants["SMILES"].values])


def fingerprints(d) -> np.ndarray:
    """(59, 167) binary MACCS matrix. **Requires RDKit.**

    Raises on any SMILES RDKit cannot parse. A silently dropped odorant would
    shift every row below it, which is the one failure mode this matrix has.
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import MACCSkeys
    except ImportError as e:                                # pragma: no cover
        raise ImportError(
            "RDKit is needed to rebuild the MACCS cache but not to draw "
            "Figure 3D. Install it with `pip install -e \".[chem]\"`, or use "
            "the cached matrix that ships in data/chem/."
        ) from e

    RDLogger.DisableLog("rdApp.*")
    smiles = _keys(d)
    bits, bad = [], []
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            bad.append((i, smi))
            continue
        bits.append(np.array(MACCSkeys.GenMACCSKeys(mol), dtype=np.uint8))
    if bad:
        raise ValueError(f"RDKit could not parse {len(bad)} SMILES: {bad}")
    return np.vstack(bits)


def cosine(F: np.ndarray) -> np.ndarray:
    """Cosine similarity between binary fingerprint rows.

    Cosine rather than Tanimoto because Figure 3D shares a 0-1 colour scale
    with the two *correlation* matrices beside it, and cosine is the same
    normalised-inner-product quantity they are. Tanimoto is the more usual
    cheminformatics choice and would be defensible; it is not the one that
    makes the three panels comparable.
    """
    F = np.asarray(F, float)
    norms = np.linalg.norm(F, axis=1, keepdims=True)
    if not np.all(norms > 0):
        raise ValueError("a fingerprint has no bits set; cannot normalise")
    Fn = F / norms
    return Fn @ Fn.T


def rebuild(d, path: Path | str = CACHE, *, verbose: bool = True) -> np.ndarray:
    """Recompute the matrix from SMILES and write the cache. Needs RDKit."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    F = fingerprints(d)
    M = cosine(F)
    np.savez_compressed(path, M=M, smiles=_keys(d), n_bits=F.shape[1],
                        fingerprint=FINGERPRINT, variant=d.variant)
    if verbose:
        off = M[~np.eye(len(M), dtype=bool)]
        print(f"wrote {path}  ({len(M)} odorants, {F.shape[1]} bits, "
              f"off-diagonal mean {off.mean():.3f} "
              f"min {off.min():.3f} max {off.max():.3f})")
    return M


def maccs_cosine(d, path: Path | str = CACHE, *, allow_rebuild: bool = False
                 ) -> np.ndarray:
    """The cached 59 x 59 MACCS cosine-similarity matrix.

    Raises if the cache is missing rather than quietly recomputing it — the
    same rule as `dimsweep.load_cached` and `pls.load_weights`. Pass
    ``allow_rebuild=True`` to build it on the spot (needs RDKit).
    """
    path = Path(path)
    if not path.exists():
        if allow_rebuild:
            return rebuild(d, path)
        raise FileNotFoundError(
            f"no MACCS cache at {path}. It ships with the repo; if you are "
            "working from a partial checkout, rebuild it with "
            "`python -m highconc.chem --rebuild` (needs the `chem` extra).")
    with np.load(path, allow_pickle=False) as z:
        M, smiles = z["M"], z["smiles"]
    want = _keys(d)
    if list(smiles) != list(want):
        raise ValueError(
            f"the cache at {path} was built for a different odorant set "
            f"({len(smiles)} SMILES) than the dataset given ({len(want)}). "
            "Rebuild it rather than reordering the matrix by hand.")
    return M


def _main(argv=None):                                       # pragma: no cover
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rebuild", action="store_true",
                    help="recompute from SMILES and overwrite the cache")
    ap.add_argument("--path", default=str(CACHE))
    args = ap.parse_args(argv)

    from .data import load
    d = load()
    if args.rebuild:
        rebuild(d, args.path)
    else:
        M = maccs_cosine(d, args.path)
        off = M[~np.eye(len(M), dtype=bool)]
        print(f"{args.path}: {M.shape}, off-diagonal mean {off.mean():.3f} "
              f"min {off.min():.3f} max {off.max():.3f}")


if __name__ == "__main__":                                  # pragma: no cover
    _main()
