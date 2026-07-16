"""Taxonomy-resolution strategy registry.

Each strategy module resolves NCBI taxIds for a datafile's rows using
whatever source-specific information is available (spec §5's per-source
split); lineage expansion once a taxId is known is shared, in base.py.
Mirrors sse_tools/readers/ (creation) and embedders//reducers/
(coordinates): one file per strategy, an explicit registry, no detection
magic hidden inside the caller.
"""
from . import em, foldseek

REGISTRY = {
    em.NAME: em,
    foldseek.NAME: foldseek,
}

# detect() order matters: foldseek's detect() is a positive check (specific
# source columns present, spec §5.2); em has no positive signal of its own
# and is the universal fallback, so it is tried last.
_DETECT_ORDER = [foldseek.NAME, em.NAME]


def detect_strategy(df, types) -> str:
    """Pick a strategy by datafile shape. Falls back to 'em' (needs only the
    id column) if nothing in _DETECT_ORDER reports a positive match.
    """
    for name in _DETECT_ORDER:
        if REGISTRY[name].detect(df, types):
            return name
    return em.NAME
