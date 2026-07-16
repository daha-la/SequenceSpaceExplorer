"""Clusterer registry — the explicit list of available --clusterer values.

To add a clusterer: write `your_method.py` here (subclass Clusterer from
base.py), import it, add one instance to REGISTRY.
"""

from .base import Clusterer, silhouette
from .kmeans import KMeansClusterer
from .hdbscan import HDBSCANClusterer

REGISTRY = {
    "kmeans": KMeansClusterer(),
    "hdbscan": HDBSCANClusterer(),
}

__all__ = ["Clusterer", "REGISTRY", "silhouette"]
