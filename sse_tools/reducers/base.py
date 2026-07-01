"""Reducer contract: an embedding matrix -> low-dimensional coordinates.

A reducer is the second plug-in axis (the first is embedders). It owns the
component label that names the coordinate columns (PCA -> "PC", UMAP -> "UMAP"),
so a column is named `<embedder-tag>_<LABEL><n>` and coordinate systems separate
cleanly in the datafile and the visualizer.

CONTRACT
--------
Subclass Reducer, set `name` (--reducer key) and `label` (component prefix).
Implement:
  reduce(X, n_components, args) -> (coords (N,k), metadata dict)
  figure(metadata, path)        -> optional diagnostic plot (may be a no-op)
"""

import numpy as np


class Reducer:
    name: str = ""
    label: str = ""

    def reduce(self, X: np.ndarray, n_components: int, args):
        raise NotImplementedError

    def figure(self, metadata: dict, path):
        pass
