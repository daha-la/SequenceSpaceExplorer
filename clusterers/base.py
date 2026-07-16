"""Clusterer contract: a (reduced) embedding matrix -> a cluster label per row.

A clusterer is a plug-in axis like embedders and reducers. It owns the `label`
fragment that names its output columns (`<tag>_<label>_cluster`, `_representative`,
`_dist_to_center`), so results from different methods coexist in one datafile.

CONTRACT
--------
Subclass Clusterer, set `name` (--clusterer key) and `label` (column fragment).
Implement:
  cluster(X, args) -> (labels (N,), metadata dict)

`labels` is an integer array; -1 marks an unclustered / noise point (HDBSCAN).
The orchestrator turns labels into columns and derives centroids, medoids
(representatives), and per-point distance-to-centre, so a clusterer only has to
assign labels and report a little metadata.
"""

import numpy as np


class Clusterer:
    name: str = ""
    label: str = ""

    def cluster(self, X: np.ndarray, args) -> "tuple[np.ndarray, dict]":
        raise NotImplementedError


def silhouette(X, labels, *, sample_size: int = 2000, seed: int = 0) -> float:
    """Silhouette score, ignoring noise (-1) points and sampling for large N.

    Returns NaN when fewer than two real clusters remain, so callers can print
    it without special-casing. silhouette_score is O(N^2) in memory, so on big
    matrices we score a random sample (deterministic via `seed`).
    """
    from sklearn.metrics import silhouette_score

    labels = np.asarray(labels)
    mask = labels != -1
    if mask.sum() < 2 or len(set(labels[mask])) < 2:
        return float("nan")
    Xs, ls = X[mask], labels[mask]
    kwargs = {}
    if Xs.shape[0] > sample_size:
        kwargs = {"sample_size": sample_size, "random_state": seed}
    return float(silhouette_score(Xs, ls, **kwargs))
