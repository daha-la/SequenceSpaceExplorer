"""2D region geometry for cluster overlays in the visualizer.

Given the projected points of one cluster, produce filled-polygon *rings* that
outline the area it occupies on the 2D map. Two shapes, same output format so
the visualizer renders them through one code path:

  concave_hull_polygons - an alpha shape that hugs the points (falls back to the
                          convex hull when the alpha shape degenerates).
  kde_polygons          - a smooth density contour enclosing a chosen fraction
                          of the cluster's points.

Each returns a list of rings; a ring is a closed (M, 2) float array. Multiple
rings mean a disconnected or holed region. All functions are defensive: too few
points, collinear points, or a numerical failure returns [] (draw nothing)
rather than raising, because they run inside the figure callback.
"""

from collections import defaultdict

import numpy as np


def _convex_ring(pts):
    from scipy.spatial import ConvexHull, QhullError

    pts = np.asarray(pts, float)
    if len(pts) < 3:
        return []
    try:
        hull = ConvexHull(pts)
    except (QhullError, ValueError):
        return []
    ring = pts[hull.vertices]
    return [np.vstack([ring, ring[:1]])]


def _stitch_rings(boundary_edges, pts):
    """Walk undirected boundary edges into closed vertex rings."""
    adj = defaultdict(list)
    for a, b in boundary_edges:
        adj[a].append(b)
        adj[b].append(a)

    remaining = set(map(lambda e: tuple(sorted(e)), boundary_edges))
    rings = []
    while remaining:
        a, b = remaining.pop()
        ring = [a, b]
        prev, cur = a, b
        while cur != ring[0]:
            nxts = [n for n in adj[cur]
                    if n != prev and tuple(sorted((cur, n))) in remaining]
            if not nxts:
                break
            nxt = nxts[0]
            remaining.discard(tuple(sorted((cur, nxt))))
            ring.append(nxt)
            prev, cur = cur, nxt
        if len(ring) >= 3:
            coords = pts[ring]
            if not np.allclose(coords[0], coords[-1]):
                coords = np.vstack([coords, coords[:1]])
            rings.append(coords)
    return rings


def concave_hull_polygons(points, tightness=0.5):
    """Alpha-shape rings for 2D points.

    `tightness` in [0, 1]: 0 keeps every Delaunay triangle (== convex hull),
    1 keeps only the most compact ones (tight, concave boundary that follows
    gaps between sub-groups).
    """
    from scipy.spatial import Delaunay, QhullError

    pts = np.asarray(points, float)
    if len(pts) < 4:
        return _convex_ring(pts)
    try:
        tri = Delaunay(pts)
    except (QhullError, ValueError):
        return _convex_ring(pts)

    ia, ib, ic = tri.simplices.T
    a = np.linalg.norm(pts[ia] - pts[ib], axis=1)
    b = np.linalg.norm(pts[ib] - pts[ic], axis=1)
    c = np.linalg.norm(pts[ic] - pts[ia], axis=1)
    s = (a + b + c) / 2.0
    area = np.sqrt(np.maximum(s * (s - a) * (s - b) * (s - c), 1e-12))
    circumradius = (a * b * c) / (4.0 * area)

    # tightness -> keep triangles below a circumradius percentile.
    q = float(np.clip(1.0 - 0.9 * tightness, 0.1, 1.0))
    threshold = np.quantile(circumradius, q)
    keep = circumradius <= threshold

    edge_count = defaultdict(int)
    for simplex, k in zip(tri.simplices, keep):
        if not k:
            continue
        for e in ((simplex[0], simplex[1]),
                  (simplex[1], simplex[2]),
                  (simplex[2], simplex[0])):
            edge_count[tuple(sorted(e))] += 1
    boundary = [e for e, n in edge_count.items() if n == 1]

    rings = _stitch_rings(boundary, pts)
    return rings if rings else _convex_ring(pts)


def _contour_rings(xs, ys, zz, level):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure()
    try:
        cs = plt.contour(xs, ys, zz, levels=[level])
        segs = cs.allsegs[0] if cs.allsegs else []
    finally:
        plt.close(fig)

    rings = []
    for seg in segs:
        seg = np.asarray(seg, float)
        if len(seg) >= 3:
            if not np.allclose(seg[0], seg[-1]):
                seg = np.vstack([seg, seg[:1]])
            rings.append(seg)
    return rings


def kde_polygons(points, coverage=0.8, grid=72, max_fit=1200, seed=0):
    """Density-contour rings enclosing ~`coverage` of the cluster's points.

    A Gaussian KDE is evaluated on a padded grid; the contour is taken at the
    density level that forms the highest-density region containing `coverage`
    fraction of the points (outliers fall outside, unlike a hull). Large
    clusters are subsampled to `max_fit` points before fitting - the density
    estimate is visually identical but much cheaper.
    """
    from scipy.stats import gaussian_kde

    pts = np.asarray(points, float)
    if len(pts) < 5:
        return []
    fit = pts
    if len(pts) > max_fit:
        fit = pts[np.random.default_rng(seed).choice(len(pts), max_fit, replace=False)]
    try:
        kde = gaussian_kde(fit.T)
    except (np.linalg.LinAlgError, ValueError):
        return []

    lo = pts.min(axis=0)
    hi = pts.max(axis=0)
    pad = (hi - lo) * 0.15 + 1e-6
    xs = np.linspace(lo[0] - pad[0], hi[0] + pad[0], grid)
    ys = np.linspace(lo[1] - pad[1], hi[1] + pad[1], grid)
    xx, yy = np.meshgrid(xs, ys)
    zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)

    coverage = float(np.clip(coverage, 0.05, 0.99))
    level = float(np.quantile(kde(pts.T), 1.0 - coverage))
    return _contour_rings(xs, ys, zz, level)


def cluster_region_rings(points, shape, *, tightness=0.5, coverage=0.8):
    """Dispatch to the requested region shape. Always returns a list of rings."""
    if shape == "kde":
        return kde_polygons(points, coverage=coverage)
    return concave_hull_polygons(points, tightness=tightness)
