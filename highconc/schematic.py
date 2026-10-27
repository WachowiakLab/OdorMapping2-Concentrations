"""Synthetic data for the Y = X · W schematic.

**Nothing in this module is measured.** It builds a small cartoon that is
internally consistent -- the Y it returns really is ``X @ W`` -- so the
schematic illustrates the *operation* without implying a result. It is kept
apart from `highconc.pls` for exactly that reason: no figure that makes a
claim should be able to reach a generator of made-up numbers by accident.

The cartoon is seeded, so the schematic is reproducible; change `seed` and you
get a different arbitrary example, not a different answer.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

#: The seed the manuscript's schematic was drawn with. Changing it changes the
#: picture, not the point -- but the exported figure would no longer match.
SEED = 7

#: Per-glomerulus gain range for W's diagonal. Chosen so every glomerulus
#: clearly amplifies (>1) while still differing from its neighbours.
GAIN_RANGE = (1.3, 2.6)

#: Fraction of off-diagonal cells carrying a cross-glomerular term, and the
#: scale of those terms. Sparse and weak relative to the diagonal, which is the
#: qualitative claim the real W supports (see `highconc.pls.block_energy`).
OFF_DIAGONAL_RATE = 0.16
OFF_DIAGONAL_SCALE = 0.45


@dataclass(frozen=True)
class Cartoon:
    """One synthetic X, W, Y triple plus cartoon map positions.

    Attributes
    ----------
    x, y : ndarray, (n,)
        Low- and high-concentration activation vectors. ``y == clip(x @ W, 0)``.
    W : ndarray, (n, n)
        Diagonal gains plus a few sparse signed off-diagonal terms.
    px, py : ndarray, (n,)
        Positions for the cartoon glomerular maps, in a unit square.
    callout : tuple[int, int]
        ``(row, col)`` of the strongest clearly-off-diagonal weight
        (``|i - j| >= 2``) -- where the "across-glomerular combination" arrow
        should point. Picked from the data so the arrow always lands on a cell
        that actually carries a visible term.
    """

    x: np.ndarray
    W: np.ndarray
    y: np.ndarray
    px: np.ndarray
    py: np.ndarray
    callout: tuple

    @property
    def n(self) -> int:
        return self.x.size

    @property
    def vmax(self) -> float:
        """One activation scale for X, Y and both maps -- they must be
        comparable, since the whole point is that Y is bigger than X."""
        return float(max(self.x.max(), self.y.max()))

    def diagonal_callout(self, k: int = 3) -> tuple:
        """``(row, col)`` of a diagonal cell for the "gain" arrow, clamped to
        fit however many glomeruli are drawn."""
        i = min(k, self.n - 1)
        return (i, i)


@lru_cache(maxsize=4)
def make(n: int = 10, seed: int = SEED) -> Cartoon:
    """Build the seeded cartoon for `n` glomeruli."""
    rng = np.random.default_rng(seed)
    x = np.clip(rng.normal(0.45, 0.28, n), 0.03, 1.0)
    W = np.diag(rng.uniform(*GAIN_RANGE, n))
    offmask = rng.random((n, n)) < OFF_DIAGONAL_RATE
    np.fill_diagonal(offmask, False)
    W[offmask] = rng.normal(0.0, OFF_DIAGONAL_SCALE, int(offmask.sum()))
    y = np.clip(x @ W, 0, None)
    px, py = rng.uniform(0, 1, n), rng.uniform(0, 1, n)

    far = offmask & (np.abs(np.subtract.outer(np.arange(n), np.arange(n))) >= 2)
    masked = np.where(far, np.abs(W), -1.0)
    oi, oj = np.unravel_index(int(np.argmax(masked)), W.shape)
    return Cartoon(x=x, W=W, y=y, px=px, py=py, callout=(int(oi), int(oj)))


# ---------------------------------------------------------------------------
# The concentration-flow cartoon
#
# Four groups of points in a 2-D "tuning space", each moving coherently in its
# own direction, plus one recruited glomerulus per group that starts silent at
# the origin and joins the motion. It is the idealized version of what
# `highconc.flow` measures in the real data.
#
# The cartoon's groups are lettered A-D, NOT numbered. Real clusters are Joel's
# 1-7 everywhere in this codebase, so a letter here can only mean "illustrative"
# -- which is the point. Do not renumber them to match.
# ---------------------------------------------------------------------------

#: Seed for the flow cartoon. As with `SEED`, changing it changes the picture,
#: not the point -- but the exported figure would no longer match.
FLOW_SEED = 3

#: Per group: centre in tuning space, coherent motion direction, palette index.
#: The palette index picks the same tab10 entries the real clusters use, so the
#: cartoon and the data figure look like they belong to one story.
FLOW_GROUPS = (
    dict(name="cluster A", centre=(-2.4, 1.3), direction=(1.5, 1.4), palette=0),
    dict(name="cluster B", centre=(2.2, 1.6), direction=(1.3, -1.5), palette=2),
    dict(name="cluster C", centre=(-2.0, -1.8), direction=(-1.4, -1.2), palette=4),
    dict(name="cluster D", centre=(2.0, -1.6), direction=(1.6, -0.4), palette=5),
)

#: Points per group, positional scatter, angular jitter (radians) and the
#: magnitude range each point's step is drawn from. Enough spread to look like
#: data, tight enough that the coherence is unmistakable.
FLOW_N = 11
FLOW_SPREAD = 0.55
FLOW_JITTER = 0.18
FLOW_MAGNITUDE = (0.8, 1.15)

#: Where a recruited glomerulus starts, as a fraction of its group's centre.
#: Near the origin -- silent at low concentration is the whole idea -- but not
#: exactly on it, or four arrows would overlap into one.
FLOW_RECRUIT_START = 0.18


@dataclass(frozen=True)
class FlowGroup:
    """One cartoon cluster: its points, their steps, and its recruited glom."""

    name: str
    palette: int
    points: np.ndarray
    steps: np.ndarray
    centre: np.ndarray
    direction: np.ndarray
    recruit_start: np.ndarray
    recruit_step: np.ndarray

    @property
    def label_at(self) -> np.ndarray:
        """Where the group's name goes -- half a step out from the centre, so
        it sits inside the flow rather than on top of the points."""
        return self.centre + self.direction * 0.5


@lru_cache(maxsize=2)
def make_flow(seed: int = FLOW_SEED) -> tuple:
    """Build the four-group concentration-flow cartoon."""
    rng = np.random.default_rng(seed)
    out = []
    for g in FLOW_GROUPS:
        centre = np.array(g["centre"], float)
        direction = np.array(g["direction"], float)
        points = rng.normal(centre, FLOW_SPREAD, size=(FLOW_N, 2))
        angle = np.arctan2(direction[1], direction[0]) + rng.normal(
            0, FLOW_JITTER, len(points))
        magnitude = np.hypot(*direction) * rng.uniform(*FLOW_MAGNITUDE,
                                                       len(points))
        steps = np.c_[np.cos(angle), np.sin(angle)] * magnitude[:, None]
        start = centre * FLOW_RECRUIT_START
        out.append(FlowGroup(name=g["name"], palette=g["palette"],
                             points=points, steps=steps, centre=centre,
                             direction=direction, recruit_start=start,
                             recruit_step=(centre + direction) - start))
    return tuple(out)


# ---------------------------------------------------------------------------
# The retuning-hypothesis cartoons -- Figure 6F
#
# Three 7x7 sketches of what W could have looked like, drawn as a hypothesis
# ladder: independent per-glomerulus gain (diagonal only), unstructured
# retuning (diagonal plus scattered off-diagonal terms), and structured
# retuning (diagonal plus coherent blocks). The real W, in 6G, is the third.
#
# All three share ONE diagonal, `RETUNING_GAINS`. That is the whole design: the
# panels differ only in what is off the diagonal, so a reader comparing them is
# comparing the hypothesis and not an incidental change of gain. Nothing here
# is measured -- these are cartoons, and they live beside the other cartoons
# for the same reason.
# ---------------------------------------------------------------------------

#: Diagonal gains, shared by all three cartoons. Seven values spanning the
#: colour ramp, deliberately not monotonic, so the diagonal reads as "each
#: glomerulus has its own gain" rather than as a gradient.
RETUNING_GAINS = (1.00, 0.30, 0.95, 0.55, 0.45, 0.90, 0.55)

#: Off-diagonal terms for the "unstructured retuning" cartoon, as
#: ``(row, col, value)``. Scattered on purpose: no two share a row-block or a
#: column-block, and one is negative, so the panel says "mixing, with no
#: organisation" and cannot be read as a faint block structure.
RETUNING_SCATTER = (
    (1, 4, 0.28), (4, 5, -0.55), (5, 3, 0.85),
    (6, 3, 0.55), (6, 5, 0.55), (6, 1, 0.30),
)

#: Off-diagonal terms for the "structured retuning" cartoon: two coherent
#: blocks. A negative 2x2 at rows 0-1 x cols 2-3 and a positive 3x3 at rows
#: 4-6 x cols 4-6 -- the second is the shape the real W actually has (see
#: `highconc.pls.block_energy`), which is why the cartoon ladder ends here.
RETUNING_BLOCKS = (
    dict(rows=(0, 2), cols=(2, 4), value=-0.62),
    dict(rows=(4, 7), cols=(4, 7), value=0.48),
)

#: Extra within-block texture, ``(row, col, value)``, so the structured blocks
#: are not flat rectangles of one colour -- real blocks are dense, not uniform.
RETUNING_BLOCK_TEXTURE = (
    (0, 1, 0.45), (1, 0, 0.55), (1, 1, 0.30),
    (4, 5, 0.55), (4, 6, 0.30), (5, 4, 0.55),
    (5, 6, 0.45), (6, 4, 0.55), (6, 5, 0.45),
)


@dataclass(frozen=True)
class Retuning:
    """One 6F cartoon: a name and the 7x7 matrix that illustrates it."""

    name: str
    W: np.ndarray

    @property
    def n(self) -> int:
        return self.W.shape[0]


def retuning_cartoons(gains=RETUNING_GAINS) -> tuple:
    """The three Fig 6F cartoons, in the order they are drawn.

    ``(independent gain control, unstructured retuning, structured retuning)``.
    Every entry is on one 0-1-ish scale, so the three panels share a diverging
    colour range and can be compared cell for cell.
    """
    g = np.asarray(gains, float)
    n = g.size

    independent = np.diag(g)

    unstructured = np.diag(g)
    for i, j, v in RETUNING_SCATTER:
        unstructured[i, j] = v

    structured = np.diag(g)
    for b in RETUNING_BLOCKS:
        r0, r1 = b["rows"]
        c0, c1 = b["cols"]
        structured[r0:r1, c0:c1] = b["value"]
    for i, j, v in RETUNING_BLOCK_TEXTURE:
        structured[i, j] = v
    np.fill_diagonal(structured, g)     # the shared diagonal always wins

    names = ("independent glomerular\ngain control", "unstructured retuning",
             "structured retuning")
    mats = (independent, unstructured, structured)
    if any(M.shape != (n, n) for M in mats):
        raise ValueError("retuning cartoons must all be n x n")
    return tuple(Retuning(name=nm, W=M) for nm, M in zip(names, mats))
