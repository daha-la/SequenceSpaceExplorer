"""UMAP reducer.

Same contract as the PCA reducer (matrix -> coordinates + metadata + figure).
UMAP has no explained-variance analogue, so the diagnostic figure is a 2D
projection scatter (the first two components) rather than a scree plot.
"""

import numpy as np

from .base import Reducer
from ..common import abort


class UMAPReducer(Reducer):
    name = "umap"
    label = "UMAP"

    def reduce(self, X, n_components, args):
        try:
            import umap
        except ImportError:
            abort("umap-learn is required for --reducer umap "
                  "(pip install umap-learn).")

        k = max(1, min(n_components, X.shape[1]))
        n_neighbors = getattr(args, "umap_neighbors", 15)
        min_dist = getattr(args, "umap_min_dist", 0.1)
        metric = getattr(args, "umap_metric", "euclidean")
        # n_neighbors must be < n_samples.
        n_neighbors = min(n_neighbors, max(2, X.shape[0] - 1))

        reducer = umap.UMAP(n_components=k, n_neighbors=n_neighbors,
                            min_dist=min_dist, metric=metric, random_state=0)
        coords = reducer.fit_transform(X.astype(np.float32))

        print(f"  UMAP: {k} components, n_neighbors={n_neighbors}, "
              f"min_dist={min_dist}, metric={metric}")
        meta = {"k": k, "n_neighbors": n_neighbors, "min_dist": min_dist,
                "metric": metric,
                "_coords2d": coords[:, :2] if k >= 2 else None}
        return coords, meta

    def figure(self, metadata, path):
        c = metadata.get("_coords2d")
        if c is None:
            return
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(c[:, 0], c[:, 1], s=5, alpha=0.5, color="steelblue",
                   edgecolors="none")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title("UMAP projection")
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
