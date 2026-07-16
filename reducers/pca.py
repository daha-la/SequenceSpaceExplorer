"""PCA reducer."""

import numpy as np

from .base import Reducer


class PCAReducer(Reducer):
    name = "pca"
    label = "PC"

    def reduce(self, X, n_components, args):
        from sklearn.decomposition import PCA
        k = min(n_components, X.shape[0], X.shape[1])
        if k < n_components:
            print(f"  note: {n_components} requested, {k} possible "
                  f"({X.shape[0]} sequences, {X.shape[1]} dims).")
        pca = PCA(n_components=k, random_state=0)
        coords = pca.fit_transform(X.astype(np.float32))
        meta = {"k": k,
                "explained_variance_ratio": pca.explained_variance_ratio_.tolist()}
        ratios = pca.explained_variance_ratio_
        print("  explained variance:",
              ", ".join(f"PC{i+1} {r*100:.1f}%" for i, r in enumerate(ratios[:5]))
              + (" ..." if len(ratios) > 5 else ""))
        return coords, meta

    def figure(self, metadata, path):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        ratios = np.array(metadata["explained_variance_ratio"])
        cum = np.cumsum(ratios)
        pcs = range(1, len(ratios) + 1)
        fig, ax = plt.subplots(figsize=(max(8, len(ratios)), 5))
        ax.bar(pcs, ratios * 100, color="steelblue", edgecolor="black",
               label="Per-PC variance")
        ax.plot(pcs, cum * 100, marker="o", color="darkred", linewidth=1.5,
                label="Cumulative variance")
        ax.set_xlabel("Principal Component")
        ax.set_ylabel("Explained Variance (%)")
        ax.set_xticks(list(pcs))
        ax.legend()
        ax.set_ylim(0, min(cum[-1] * 100 * 1.1, 105))
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
