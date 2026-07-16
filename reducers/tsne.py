"""t-SNE reducer.

Same contract as the PCA / UMAP reducers. t-SNE has no explained-variance
analogue, so the diagnostic figure is a 2D projection scatter.

Note: t-SNE is intended for 2D (occasionally 3D) visualization; high
n_components is slow and rarely meaningful. scikit-learn's exact/barnes-hut
solver caps at n_components <= 3 for barnes_hut, so this reducer keeps the
default low and lets the user override.

t-SNE reads its input only through pairwise distances, so the many low-variance
directions of a wide embedding contribute noise to every distance it uses. The
standard remedy (van der Maaten's own recommendation) is to project onto the
top ~50 principal components first: it denoises the neighbor structure and cuts
the cost of the neighbor search. That is --tsne-pca, applied here before t-SNE
runs. It is a pre-processing step, not the reduction; t-SNE still produces the
final coordinates. Note that TSNE(init="pca") is unrelated -- that only seeds
the low-dimensional layout.
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

        X, pca_meta = self._pre_reduce(X, getattr(args, "tsne_pca", 50), k)

        tsne = TSNE(n_components=k, perplexity=perplexity, method=method,
                    init="pca", random_state=0)
        coords = tsne.fit_transform(X.astype(np.float32))

        print(f"  t-SNE: {k} components, perplexity={perplexity:g}, method={method}")
        meta = {"k": k, "perplexity": perplexity, "method": method,
                "_coords2d": coords[:, :2] if k >= 2 else None}
        meta.update(pca_meta)
        return coords, meta

    def _pre_reduce(self, X, pca_dims, k):
        """Project onto the top `pca_dims` PCs before t-SNE. 0/negative disables."""
        if pca_dims is None or pca_dims <= 0:
            return X, {"pca_dims": None}
        # Nothing to gain once the request meets or exceeds the input width, and
        # a PCA can never return more components than samples or features.
        limit = min(X.shape[0], X.shape[1])
        if pca_dims >= X.shape[1]:
            print(f"  t-SNE PCA pre-reduction: skipped "
                  f"({X.shape[1]} dims <= --tsne-pca {pca_dims}).")
            return X, {"pca_dims": None}
        d = max(k, min(pca_dims, limit))

        from sklearn.decomposition import PCA
        pca = PCA(n_components=d, random_state=0)
        reduced = pca.fit_transform(X.astype(np.float32))
        retained = float(pca.explained_variance_ratio_.sum())
        print(f"  t-SNE PCA pre-reduction: {X.shape[1]} -> {d} dims, "
              f"{retained * 100:.1f}% of variance retained.")
        return reduced, {"pca_dims": d, "pca_variance_retained": retained}

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
        d = metadata.get("pca_dims")
        retained = metadata.get("pca_variance_retained")
        sub = f"perplexity={metadata['perplexity']:g}"
        if d:
            sub += f", PCA pre-reduction to {d} dims ({retained * 100:.1f}% variance)"
        else:
            sub += ", no PCA pre-reduction"
        # Axes are unitless and inter-cluster distance is not meaningful in t-SNE;
        # the caption keeps that attached to the figure.
        ax.text(0.5, -0.11, sub, transform=ax.transAxes, ha="center",
                fontsize=8, color="dimgray")
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
