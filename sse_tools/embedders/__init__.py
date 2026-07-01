"""Embedder registry — the explicit list of available --embedder values.

To add an embedder:
  1. Write `your_model.py` here, subclassing Embedder (base.py). Copy esmc.py
     for a sequence model, or prostt5.py for a structure model.
  2. Import it and add one instance to REGISTRY.
"""

from .base import Embedder, EmbedContext, resolve_device, warn_if_slow
from .esmc import EsmcEmbedder
from .prostt5 import ProstT5Embedder
from .saprot import SaProtEmbedder

REGISTRY = {
    "esmc": EsmcEmbedder(),
    "prostt5": ProstT5Embedder(),
    "saprot": SaProtEmbedder(),
}

__all__ = ["Embedder", "EmbedContext", "REGISTRY", "resolve_device", "warn_if_slow"]
