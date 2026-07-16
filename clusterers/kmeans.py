"""K-means clusterer with silhouette-based selection of k.

k-means partitions every point into exactly one of k spherical clusters, so each
cluster has a well-defined centroid - which is what the orchestrator snaps a
representative (medoid) sequence to. Pass --k to fix k, or let the tool sweep
--k-min..--k-max and keep the k with the best silhouette score.
"""

import numpy as np

from .base import Clusterer, silhouette


class KMeansClusterer(Clusterer):
    name = "kmeans"
    label = "kmeans"

    def cluster(self, X, args):
        from sklearn.cluster import KMeans

        if args.k is not None:
            if args.k < 2:
                raise ValueError("--k must be at least 2.")
            km = KMeans(n_clusters=args.k, n_init=10, random_state=0).fit(X)
            sil = silhouette(X, km.labels_)
            print(f"  k-means: k={args.k} (manual), silhouette={sil:.3f}")
            return km.labels_, {"k": args.k, "selected": "manual",
                                "silhouette": sil, "inertia": float(km.inertia_)}

        k_min = max(2, args.k_min)
        k_max = min(args.k_max, X.shape[0] - 1)
        if k_max < k_min:
            raise ValueError(f"k range is empty: k_min={k_min}, k_max={k_max} "
                             f"for {X.shape[0]} points.")

        best = None
        sweep = []
        print(f"  k-means: sweeping k={k_min}..{k_max} by silhouette ...")
        for k in range(k_min, k_max + 1):
            km = KMeans(n_clusters=k, n_init=10, random_state=0).fit(X)
            sil = silhouette(X, km.labels_)
            sweep.append((k, sil))
            print(f"    k={k}: silhouette={sil:.3f}")
            if best is None or (not np.isnan(sil) and sil > best[1]):
                best = (k, sil, km)

        if best is None:
            raise ValueError("k-means sweep produced no valid clustering.")
        k, sil, km = best
        print(f"  k-means: selected k={k} (silhouette={sil:.3f})")
        return km.labels_, {"k": k, "selected": "auto", "silhouette": sil,
                            "inertia": float(km.inertia_),
                            "sweep": [(int(a), float(b)) for a, b in sweep]}
