"""Canonical cluster identity and numbering.

Two numberings for the same seven NMF clusters exist in the wild:

* **Ours** — the raw NMF component index 0-6, historically displayed as the
  letters A-G.
* **Joel's** — 1-7, the numbering used in the manuscript figures and in every
  conversation with the Utah lab.

Everything downstream of this module uses **1-7**. A-G survives only as a
legacy mapping, so that old panel code and old figure exports can still be
read.

The mapping is *exact*, not inferred: the seven clusters have seven distinct
sizes, so matching sizes determines the correspondence uniquely.

    >>> from highconc import clusters as cl
    >>> cl.to_joel(0)            # our component index -> Joel's number
    1
    >>> cl.to_index(2)           # Joel's number -> our component index
    4
    >>> cl.LABELS[cl.to_index(2)]
    '2'
"""

from __future__ import annotations

import numpy as np

# Joel's cluster n (1-based) is our component index JOEL_TO_INDEX[n - 1].
JOEL_TO_INDEX: tuple[int, ...] = (0, 4, 2, 5, 3, 1, 6)

# Our component index i is Joel's cluster INDEX_TO_JOEL[i].
INDEX_TO_JOEL: tuple[int, ...] = tuple(
    JOEL_TO_INDEX.index(i) + 1 for i in range(7)
)

# Glomerulus counts per cluster, keyed by our component index 0-6. These are
# what pin the mapping down.
SIZES_BY_INDEX: tuple[int, ...] = (45, 78, 40, 30, 41, 75, 67)

# The same counts read out in Joel's 1-7 order.
SIZES_BY_JOEL: tuple[int, ...] = tuple(
    SIZES_BY_INDEX[JOEL_TO_INDEX[n - 1]] for n in range(1, 8)
)

#: Display labels in *our* component-index order. This is what panels use.
LABELS: tuple[str, ...] = tuple(str(INDEX_TO_JOEL[i]) for i in range(7))

#: Display labels in Joel's 1-7 order — the order figures lay clusters out in.
#: Use this for axis ticks; use LABELS when the axis is in component order.
LABELS_BY_JOEL: tuple[str, ...] = tuple(str(n) for n in range(1, 8))

#: The retired A-G labels, also in our component-index order. Present so old
#: exports remain readable; do not use in new panels.
LEGACY_LABELS: tuple[str, ...] = tuple("ABCDEFG")

N_CLUSTERS = 7


def to_joel(index: int) -> int:
    """Our component index 0-6 -> Joel's cluster number 1-7."""
    return INDEX_TO_JOEL[_check(index, 0, 6, "component index")]


def to_index(joel: int) -> int:
    """Joel's cluster number 1-7 -> our component index 0-6."""
    return JOEL_TO_INDEX[_check(joel, 1, 7, "Joel cluster number") - 1]


def label(index: int) -> str:
    """Display label ('1'-'7') for one of our component indices."""
    return LABELS[_check(index, 0, 6, "component index")]


def legacy_label(index: int) -> str:
    """Retired A-G label for one of our component indices."""
    return LEGACY_LABELS[_check(index, 0, 6, "component index")]


def joel_order() -> list[int]:
    """Our component indices, ordered so cluster 1 comes first, 7 last.

    Use this to reorder rows/columns of a matrix that is in component order
    into the order the manuscript figures display.
    """
    return [JOEL_TO_INDEX[n - 1] for n in range(1, 8)]


def reorder(a, axis: int = 0):
    """Reorder an array along `axis` from component order into Joel's 1-7 order."""
    return np.take(np.asarray(a), joel_order(), axis=axis)


# ---------------------------------------------------------------------------
# Chemical classes
# ---------------------------------------------------------------------------
# Not clusters, but the same kind of thing: a canonical label set with two
# competing orderings in circulation. It lives here so a panel needs one import
# for "what are the categories called and in what order do they go".

#: Class letter -> name, as used on axis labels.
CHEM_CLASS_NAMES: dict[str, str] = {
    "C": "acids", "D": "aldehydes", "E": "esters", "F": "furans",
    "K": "ketones", "L": "alcohols", "N": "amines", "P": "pyrazines",
    "R": "aromatics", "S": "thiazoles",
}

#: Our odorant-table order — alphabetical by letter, which is how the odorant
#: table is laid out.
CLASS_ORDER_TABLE: tuple[str, ...] = ("C", "D", "E", "F", "K", "L", "N", "P",
                                      "R", "S")

#: Best estimate of Matt's manuscript ordering. Reconstructed, not supplied:
#: it lifts agreement with the published block layout from 0.45 to 0.77. Ship
#: both and let the Utah lab pick — this is a call-site argument, not a fact we
#: have to establish.
CLASS_ORDER_MATT: tuple[str, ...] = ("C", "D", "E", "K", "L", "R", "F", "N",
                                     "P", "S")

#: Matt's handwritten note on Figure 3 (2026-08-23): "consider reordering with
#: esters right after acids". CLASS_ORDER_MATT with E moved ahead of D. Shipped
#: so the request is one argument away, not a code edit; it is NOT the default,
#: because every existing export uses the order actually drawn in his PDF.
CLASS_ORDER_MATT_ESTERS_FIRST: tuple[str, ...] = ("C", "E", "D", "K", "L", "R",
                                                  "F", "N", "P", "S")

CLASS_ORDERS: dict[str, tuple[str, ...]] = {
    "table": CLASS_ORDER_TABLE,
    "matt": CLASS_ORDER_MATT,
    "matt_esters_first": CLASS_ORDER_MATT_ESTERS_FIRST,
}

#: Colour per chemical class, fixed by LETTER so that reordering the blocks
#: never recolours them. Keyed in table order over tab10 — ten classes, ten
#: colours, no collisions. Matt's PDF uses a different palette; matching it
#: exactly would mean freezing an unpublished choice, and the letter->colour
#: binding is the property that actually matters (Fig 3B, 3C, 3D and 3E-G must
#: all say "amine" in the same colour).
#: Written out as hex rather than indexed off tab10 at import time, so this
#: module stays free of matplotlib — cluster and class IDENTITY should not
#: depend on a plotting library. The values ARE tab10, in table order.
CHEM_CLASS_COLORS: dict[str, str] = {
    "C": "#1f77b4",   # acids
    "D": "#ff7f0e",   # aldehydes
    "E": "#2ca02c",   # esters
    "F": "#d62728",   # furans
    "K": "#9467bd",   # ketones
    "L": "#8c564b",   # alcohols
    "N": "#e377c2",   # amines
    "P": "#7f7f7f",   # pyrazines
    "R": "#bcbd22",   # aromatics
    "S": "#17becf",   # thiazoles
}


#: Figure 3 panels A-D only: the ten "R" (aromatics) odorants broken into
#: three display subclasses along the odorant axis. Nothing else in the
#: package ever produces these labels -- Fig 4/5/6's `chem_class` and the
#: 3E/F/G spectra still see plain "R" -- so extending CHEM_CLASS_NAMES/COLORS
#: with them is safe: no other panel can accidentally encounter one. Keys are
#: two/three letters so they cannot collide with a real odorant code like
#: "R3" or "R10" (`similarity.split_block` matches on those codes).
#: Grouping and colours supplied directly (2026-09-02); the split follows
#: `data/burton_et_al/odorant_order.txt`'s own within-block sequence for the
#: R block: first 1 odorant Ar. Ketone, next 7 Ar. Ketone/Ester, last 2
#: Phenol/Methoxy.
AROMATIC_SUBCLASS_NAMES: dict[str, str] = {
    "RK": "Ar. Ketone", "RKE": "Ar. Ketone/Ester", "RP": "Phenol/Methoxy",
}
AROMATIC_SUBCLASS_COLORS: dict[str, str] = {
    "RK": "#a32e49", "RKE": "#00bbb4", "RP": "#d59537",
}
AROMATIC_SUBCLASS_BY_CODE: dict[str, str] = {
    "R27": "RK",
    "R11": "RKE", "R12": "RKE", "R21": "RKE", "R31": "RKE",
    "R3": "RKE", "R29": "RKE", "R40": "RKE",
    "R10": "RP", "R23": "RP",
}
CHEM_CLASS_NAMES.update(AROMATIC_SUBCLASS_NAMES)
CHEM_CLASS_COLORS.update(AROMATIC_SUBCLASS_COLORS)


#: Within a chemical class, whose odorant order? Ours is the odorant table's.
#: Matt's figure PDF agrees with it everywhere EXCEPT the pyrazine block, where
#: he puts the two pyridines last and moves P19 forward. Recovered exactly, not
#: guessed: reading his Fig 3E/F/G bar geometry out of the PDF and matching it
#: against all 376 glomeruli reproduces every bar to within 0.008 pt, and the
#: pyrazine block is the only place the two orders disagree. Keyed by the
#: odorant CODE (the leading token of the `Odorant` column), so it survives any
#: renaming of the display names.
WITHIN_BLOCK_MATT: dict[str, tuple[str, ...]] = {
    "P": ("P1", "P2", "P4", "P16", "P19", "P3", "P5", "P18", "P12", "P13"),
}

#: The full 59-odorant order supplied in data/burton_et_al/odorant_order.txt
#: (2026-09-02), covering every class rather than only the pyrazine reshuffle
#: WITHIN_BLOCK_MATT makes. Its block sequence is CLASS_ORDER_MATT's, so pair
#: this with order="matt".
WITHIN_BLOCK_BURTON: dict[str, tuple[str, ...]] = {
    "C": ("C7", "C1", "C2", "C3", "C6"),
    "D": ("D15", "D16", "D1", "D2", "D8", "D14"),
    "E": ("E1", "E2", "E9", "E4", "E14"),
    "K": ("K3", "K1", "K15", "K2", "K11"),
    "L": ("L2", "L8"),
    "R": ("R27", "R11", "R12", "R21", "R31", "R3", "R29", "R40", "R10", "R23"),
    "F": ("F5", "F6"),
    "N": ("N6", "N11", "N15", "N9", "N7", "N1", "N3"),
    "P": ("P1", "P2", "P4", "P16", "P19", "P3", "P5", "P18", "P12", "P13"),
    "S": ("S1", "S2", "S3", "S5", "S7", "S13", "S17"),
}

WITHIN_BLOCK_ORDERS = {"table": {}, "matt": WITHIN_BLOCK_MATT,
                       "burton": WITHIN_BLOCK_BURTON}


def class_color(letter: str):
    """Colour for a chemical-class letter. Raises on an unknown class rather
    than inventing a colour — a silent grey block is worse than a traceback."""
    try:
        return CHEM_CLASS_COLORS[letter]
    except KeyError:
        raise KeyError(
            f"unknown chemical class {letter!r}; "
            f"expected one of {sorted(CHEM_CLASS_COLORS)}") from None


def class_order(order="table", present=None) -> list[str]:
    """Chemical-class letters in a chosen order.

    `order` is ``"table"``, ``"matt"``, or an explicit sequence of letters.
    `present` optionally restricts the result to the classes actually in the
    data, keeping the requested order and appending anything unrecognised so
    nothing is silently dropped.
    """
    if isinstance(order, str):
        try:
            seq = list(CLASS_ORDERS[order.lower()])
        except KeyError:
            raise ValueError(
                f"unknown class order {order!r}; "
                f"expected one of {sorted(CLASS_ORDERS)} or a sequence"
            ) from None
    else:
        seq = list(order)
    if present is None:
        return seq
    present = list(dict.fromkeys(present))
    kept = [c for c in seq if c in present]
    return kept + [c for c in present if c not in kept]


def class_name(letter: str) -> str:
    """Display name for a chemical-class letter; the letter itself if unknown."""
    return CHEM_CLASS_NAMES.get(letter, letter)


def fraction_responding(dataset, conc: str, perm=None) -> np.ndarray:
    """Per odorant and cluster, the fraction of the cluster's glomeruli that respond.

    Returns ``(n_odorants, 7)``: column ``n - 1`` is cluster ``n`` (1-7). A
    glomerulus responds to an odorant when its response is above zero, at the
    concentration `conc` (``"low"`` or ``"high"``). `perm` reorders the odorant
    rows, e.g. the permutation from `highconc.similarity.class_sorted`.
    """
    if conc not in ("low", "high"):
        raise ValueError(f"conc must be 'low' or 'high', got {conc!r}")
    resp = getattr(dataset, conc).to_numpy(float) > 0
    out = np.zeros((resp.shape[0], N_CLUSTERS))
    for n in range(1, N_CLUSTERS + 1):
        out[:, n - 1] = resp[:, dataset.cluster_mask(n)].mean(axis=1)
    return out if perm is None else out[np.asarray(perm)]


def _check(value: int, lo: int, hi: int, what: str) -> int:
    value = int(value)
    if not lo <= value <= hi:
        raise ValueError(f"{what} must be in [{lo}, {hi}], got {value}")
    return value


def _self_check() -> None:
    """The mapping is only exact while the sizes stay distinct. Assert both."""
    assert sorted(JOEL_TO_INDEX) == list(range(7)), "mapping is not a permutation"
    assert len(set(SIZES_BY_INDEX)) == 7, "sizes are no longer distinct"
    assert all(to_index(to_joel(i)) == i for i in range(7)), "round-trip broken"
    assert SIZES_BY_JOEL == (45, 41, 40, 75, 30, 78, 67), "sizes drifted"


_self_check()
