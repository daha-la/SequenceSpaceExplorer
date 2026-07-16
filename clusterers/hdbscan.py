"""HDBSCAN clusterer (density-based; discovers cluster count, flags outliers).

Unlike k-means, HDBSCAN does not take a cluster count: it finds variable-density
clusters and labels points that fit none as noise (-1). It is sensitive to
dimensionality, so cluster on a PCA reduction (the sse_cluster default), not the
full embedding. Tune granularity with --min-cluster-size (the smallest group
called a cluster) and optionally --min-samples (higher = more conservative, more
noise). Uses scikit-learn's HDBSCAN, so no extra dependency.
"""

import numpy as np

from .base import Clusterer, silhouette


class HDBSCANClusterer(Clusterer):
    name = "hdbscan"
    label = "hdbscan"

    def cluster(self, X, args):
        from sklearn.cluster import HDBSCAN

        # copy=True: never mutate the caller's matrix in place (and it pins the
        # behaviour sklearn changes the default of in 1.10).
        hdb = HDBSCAN(min_cluster_size=args.min_cluster_size,
                      min_samples=args.min_samples, copy=True).fit(X)
        labels = hdb.labels_
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = int((labels == -1).sum())
        sil = silhouette(X, labels)
        print(f"  HDBSCAN: {n_clusters} clusters, {n_noise} noise points "
              f"(min_cluster_size={args.min_cluster_size}, "
              f"min_samples={args.min_samples}), silhouette={sil:.3f}")
        if n_clusters == 0:
            print("    note: no clusters found - try a smaller --min-cluster-size "
                  "or --min-samples.")
        return labels, {"n_clusters": n_clusters, "n_noise": n_noise,
                        "min_cluster_size": args.min_cluster_size,
                        "min_samples": args.min_samples, "silhouette": sil}
