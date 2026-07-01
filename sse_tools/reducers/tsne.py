"""t-SNE reducer.

Same contract as the PCA / UMAP reducers. t-SNE has no explained-variance
analogue, so the diagnostic figure is a 2D projection scatter.

Note: t-SNE is intended for 2D (occasionally 3D) visualization; high
n_components is slow and rarely meaningful. scikit-learn's exact/barnes-hut
solver caps at n_components <= 3 for barnes_hut, so this reducer keeps the
default low and lets the user override.
"""

import numpy as np

from .base import Reducer
from ..common import abort


class TSNEReducer(Reducer):
    name = "tsne"
    label = "TSNE"

    def reduce(self, X, n_components, args):
        from sklearn.manifold import TSNE

        k = max(1, min(n_components, X.shape[1]))
        perplexity = getattr(args, "tsne_perplexity", 30.0)
        # perplexity must be < n_samples (sklearn requires perplexity < n_samples).
        perplexity = float(min(perplexity, max(5.0, (X.shape[0] - 1) / 3.0)))
        # barnes_hut (the fast default) only supports k <= 3.
        method = "barnes_hut" if k <= 3 else "exact"

        tsne = TSNE(n_components=k, perplexity=perplexity, method=method,
                    init="pca", random_state=0)
        coords = tsne.fit_transform(X.astype(np.float32))

        print(f"  t-SNE: {k} components, perplexity={perplexity:g}, method={method}")
        meta = {"k": k, "perplexity": perplexity, "method": method,
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
        ax.set_xlabel("t-SNE 1")
        ax.set_ylabel("t-SNE 2")
        ax.set_title("t-SNE projection")
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
