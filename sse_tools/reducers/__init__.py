"""Reducer registry — the explicit list of available --reducer values.

To add a reducer: write `your_method.py` here (subclass Reducer from base.py),
import it, add one instance to REGISTRY.
"""

from .base import Reducer
from .pca import PCAReducer
from .umap import UMAPReducer
from .tsne import TSNEReducer

REGISTRY = {
    "pca": PCAReducer(),
    "umap": UMAPReducer(),
    "tsne": TSNEReducer(),
}

__all__ = ["Reducer", "REGISTRY"]
