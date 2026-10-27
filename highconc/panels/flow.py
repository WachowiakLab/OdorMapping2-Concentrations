"""Panels for the concentration-flow figures, reading `highconc.flow`.

Together they make one claim two ways: glomeruli in the same NMF cluster move
in the same direction through tuning space as concentration rises. `flow_3d`
shows it (per-glomerulus arrows, per-cluster spread cones, cluster-mean
trajectories); `angle_hist`, `angle_by_cluster` and `coherence_matrix` measure
it (each glomerulus's angle to its own cluster's mean direction, against the
angle to everyone else's).

Clusters are Joel's 1-7, and colour comes from `highconc.style.cluster_color`,
so a cluster is the same colour here as in every other figure.
"""

from __future__ import annotations

import numpy as np

from .. import style as hs

#: Colour of the shared low-concentration anchor and the chance references.
ANCHOR = "0.25"
CHANCE = "0.4"

#: Warm end of the diverging map, reused for the "own centroid" series so the
#: histogram and the matrix agree on what "aligned" looks like.
C_OWN = "#b30000"
C_OTHER = "0.6"


def _ortho_basis(u):
    """Two unit vectors spanning the plane perpendicular to `u`."""
    u = u / np.linalg.norm(u)
    a = np.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(u, a)
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(u, e1)


def _cone(ax, apex, u, length, half_angle, color, alpha=0.13, n=32):
    """Translucent cone: apex at `apex`, axis `u`, given length and half-angle."""
    e1, e2 = _ortho_basis(u)
    phi = np.linspace(0, 2 * np.pi, n)
    P, T = np.meshgrid(phi, np.array([0.0, length]))
    R = T * np.tan(half_angle)
    coords = [apex[i] + u[i] * T + R * (np.cos(P) * e1[i] + np.sin(P) * e2[i])
              for i in range(3)]
    ax.plot_surface(*coords, color=color, alpha=alpha, linewidth=0, shade=False)


def flow_3d(ax, low, high, cluster_number, flows, *, faint_arrows=True,
            cones=True, origin=True, arrow_alpha=0.18, arrow_lw=0.5,
            mean_lw=3.0, spread_alpha=0.13):
    """Concentration trajectories in a 3-PC tuning space, on a 3-D `ax`.

    Three layers, and they are meant to be read in this order: the faint
    per-glomerulus arrows are the data, the shaded cone is the SPREAD of
    directions within a cluster (from `flows[i]["half_angle"]` -- not a
    confidence interval), and the bold arrow is the cluster mean.

    `flows` comes from `highconc.flow.cluster_flow`. Cluster labels are
    deliberately not drawn: in a 3-D projection matplotlib cannot place them
    without collisions, and they are added by hand at assembly.
    """
    low, high = np.asarray(low, float), np.asarray(high, float)
    cluster_number = np.asarray(cluster_number, int)
    dv = high - low

    if origin:
        ax.scatter([0], [0], [0], marker="o", s=26, color=ANCHOR,
                   depthshade=False, zorder=10)
    if faint_arrows:
        for num in np.unique(cluster_number):
            m = cluster_number == num
            if m.sum() < 2:
                continue
            ax.quiver(low[m, 0], low[m, 1], low[m, 2],
                      dv[m, 0], dv[m, 1], dv[m, 2],
                      color=hs.cluster_color(int(num)), alpha=arrow_alpha,
                      lw=arrow_lw, arrow_length_ratio=0.05)
    for f in flows:
        color = hs.cluster_color(f["number"])
        if cones:
            _cone(ax, f["origin"], f["unit"], f["length"],
                  np.radians(f["half_angle"]), color, alpha=spread_alpha)
        ax.scatter(*f["origin"], s=40, color=color, edgecolors="white",
                   linewidths=0.8, depthshade=False)
        ax.quiver(*f["origin"], *f["mean"], color=color, lw=mean_lw,
                  arrow_length_ratio=0.2)
    return ax


def style_3d(ax, labels=("PC 1", "PC 2", "PC 3"), *, label_fontsize=9,
             tick_fontsize=7, pane="white", pane_alpha=0.4, edge="0.85"):
    """Pale panes and small ticks for a 3-D axis. The default 3-D styling is
    far too heavy next to the 2-D panels it sits beside."""
    for setter, text in zip((ax.set_xlabel, ax.set_ylabel, ax.set_zlabel), labels):
        setter(text, fontsize=label_fontsize, labelpad=1)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(pane)
        axis.pane.set_alpha(pane_alpha)
        axis.pane.set_edgecolor(edge)
    ax.tick_params(labelsize=tick_fontsize)
    return ax


def data_box_aspect(ax, points):
    """Equal-unit framing, so the 3-D view matches a picked camera.

    matplotlib's default cube stretches each axis to the box, which silently
    changes every angle in the picture -- in a figure whose whole claim is
    about directions, that is not cosmetic.
    """
    p = np.vstack(points)
    ax.set_box_aspect(tuple(p.max(0) - p.min(0)))
    return ax


def angle_hist(ax, own, other, *, bins=61, own_color=C_OWN,
               other_color=C_OTHER, alpha=0.6, medians=True,
               xlabel="angle to cluster-mean direction (°)", ylabel="density",
               legend=True, legend_fontsize=8.5):
    """Angle to own vs. other cluster-mean directions, as densities.

    Densities, not counts: there are six times as many "other" pairs as "own"
    ones, so raw counts would show the sampling and not the effect.
    """
    edges = np.linspace(0, 180, bins)
    ax.hist(other, bins=edges, density=True, color=other_color, alpha=alpha,
            label=f"to other centroids (med {np.median(other):.0f}°)")
    ax.hist(own, bins=edges, density=True, color=own_color, alpha=alpha,
            label=f"to own centroid (med {np.median(own):.0f}°)")
    if medians:
        ax.axvline(np.median(other), color=CHANCE, ls="--", lw=1)
        ax.axvline(np.median(own), color=own_color, ls="--", lw=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if legend:
        ax.legend(fontsize=legend_fontsize, frameon=False, loc="upper right")
    hs.despine(ax)
    return ax


def angle_by_cluster(ax, groups, numbers, *, reference=None, widths=0.62,
                     alpha=0.85, xlabel="cluster",
                     ylabel="angle to own centroid (°)", legend=True,
                     legend_fontsize=8.5):
    """One box per cluster: how tightly its glomeruli track its own mean.

    `reference` draws the to-other median as a chance line -- without it a
    reader has no way to tell a tight cluster from a loose one, since the
    absolute angle depends on the dimensionality of the embedding.
    """
    bp = ax.boxplot(list(groups), patch_artist=True, showfliers=False,
                    widths=widths, medianprops=dict(color="black"))
    for patch, num in zip(bp["boxes"], numbers):
        patch.set_facecolor(hs.cluster_color(int(num)))
        patch.set_alpha(alpha)
    if reference is not None:
        ax.axhline(reference, color=CHANCE, ls="--", lw=1.2,
                   label=f"to-other median ({reference:.0f}°)")
    ax.set_xticks(range(1, len(numbers) + 1))
    ax.set_xticklabels([str(int(n)) for n in numbers])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if legend and reference is not None:
        ax.legend(fontsize=legend_fontsize, frameon=False, loc="lower right")
    hs.despine(ax)
    return ax


def coherence_matrix(ax, M, labels, *, chance, vmin=0.0, vmax=1.0, cmap=None,
                     split=None, annotate=True, annotate_fontsize=7,
                     line_color="0.2", line_width=1.4,
                     xlabel="aligned to cluster-mean (centroid)",
                     ylabel="glomeruli from cluster"):
    """The 7x7 coherence / cross-coherence matrix, pivoted at chance.

    The diverging scale is centred on `chance` -- the across-cluster mean cosine
    -- not on zero. Zero is not the null here: two random directions in three
    dimensions are not orthogonal on average once both are drawn from the same
    anisotropic cloud, so a zero-pivoted map would read almost every cell as
    "aligned". Warm = above chance, gainsboro = chance, teal = below.

    Pass `M` and `labels` from `highconc.flow.CentroidCoherence.reorder`, and
    `split` to rule off the coherence super-groups. Returns the image.
    """
    from matplotlib.colors import TwoSlopeNorm

    M = np.asarray(M, float)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=chance, vmax=vmax)
    im = ax.imshow(M, cmap=hs.DIV if cmap is None else cmap, norm=norm)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.tick_params(length=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if annotate:
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                v = M[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=annotate_fontsize,
                        color="white" if (v < 0.5 * chance or v > 0.90) else "0.15")
    if split:
        ax.axvline(split - 0.5, color=line_color, lw=line_width)
        ax.axhline(split - 0.5, color=line_color, lw=line_width)
    return im
