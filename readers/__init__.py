"""Readers — one module per source format, plus the explicit registry.

To add a source format:
  1. Write `your_source.py` here, following the contract in base.py (copy
     foldseek.py if your source is not already tabular).
  2. Import its reader function below.
  3. Add one line to REGISTRY under a short --source tag.
Nothing else in the codebase needs to change, and no user-facing script changes.

Shared utilities (abort, reserved names, datafile I/O) live in sse_tools/common.py.
"""

from .base import ReaderResult
from .em import read_em
from .fasta import read_fasta
from .foldseek import read_foldseek

REGISTRY = {
    "em": read_em,
    "fasta": read_fasta,
    "fs": read_foldseek,
}

__all__ = ["ReaderResult", "REGISTRY"]
