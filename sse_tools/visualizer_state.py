"""State helpers for the SSE Dash visualizer.

The visualizer opens an entry, reads one typed .sse.tsv datafile, and derives
UI state from the Type row. This module contains no Dash callbacks.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

from .common import (
    COL_ID,
    COL_QUERY,
    COL_SEQ,
    TYPE_COORDINATE,
    TYPE_ID,
    TYPE_LABEL,
    SSEError,
    abort,
    id_column,
    read_datafile,
    resolve_entry_path,
)

BOOL_STRINGS = {"true", "false", "yes", "no", "1", "0"}
COMMA_THRESHOLD = 0.30
COMMA_THRESHOLD_LOW = 0.05
HIGH_CARD_UNIQUE_RATIO = 0.80
HIGH_CARD_MEAN_LEN = 30
NUMERIC_COL_THRESHOLD = 0.80
MAX_CAT_UNIQUE = 200

NAME_KEYWORDS = {"name", "label", "alias", "accession", "gene"}

# Nesting-detection guardrails for the high-cardinality categorical rescue pass.
_NEST_MIN_COVERAGE = 0.40
_NEST_MIN_CLEAN_FRAC = 0.98
_NEST_MIN_PAIR_OVERLAP = 0.40


@dataclass(frozen=True)
class EntryContext:
    stem: str
    entry_dir: Path
    datafile_path: Path
    logs_dir: Path
    figures_dir: Path
    structures_dir: Path
    msa_cache_dir: Path
    layers_path: Path
    jobs_path: Path
    manifest_path: Path


@dataclass
class VisualizerState:
    entry: EntryContext
    df: pd.DataFrame
    types: Dict[str, str]
    id_col: str
    label_cols: List[str]
    coord_cols: List[str]
    col_meta: Dict[str, dict]
    query_ids: List[str]
    name_cols: List[str]
    coord_systems: Dict[str, List[str]]
    x_col: str | None
    y_col: str | None
    warning: str = ""


def resolve_entry(arg: str | None, *, entries_dir=None) -> EntryContext:
    """Resolve an entry stem, entry directory, or .sse.tsv path, then set up
    this entry's working directories. Path resolution itself is shared with
    every other ENTRY-taking script via common.resolve_entry_path; this adds
    only what's specific to the visualizer (figures/structures/msa_cache
    dirs, and creating them).
    """
    if not arg:
        abort("No entry supplied. Usage: python scripts/sse_visualizer.py <entry-stem|entry-dir|datafile.sse.tsv>")
    datafile_path = resolve_entry_path(arg, entries_dir=entries_dir)
    entry_dir = datafile_path.parent
    stem = datafile_path.name[:-len(".sse.tsv")]

    logs_dir = entry_dir / "logs"
    figures_dir = entry_dir / "figures"
    structures_dir = entry_dir / "structures"
    msa_cache_dir = entry_dir / "msa_cache"
    logs_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    structures_dir.mkdir(parents=True, exist_ok=True)
    msa_cache_dir.mkdir(parents=True, exist_ok=True)

    return EntryContext(
        stem=stem,
        entry_dir=entry_dir,
        datafile_path=datafile_path,
        logs_dir=logs_dir,
        figures_dir=figures_dir,
        structures_dir=structures_dir,
        msa_cache_dir=msa_cache_dir,
        layers_path=logs_dir / "layers.json",
        jobs_path=logs_dir / "jobs.json",
        manifest_path=logs_dir / f"{stem}.sse.manifest.json",
    )


def _is_boolean(s: pd.Series) -> bool:
    if pd.api.types.is_bool_dtype(s):
        return True
    non_null = s.replace("", pd.NA).dropna().astype(str).str.strip().str.lower().unique()
    return len(non_null) > 0 and set(non_null).issubset(BOOL_STRINGS)


def _is_numeric(s: pd.Series) -> bool:
    if _is_boolean(s):
        return False
    non_null = s.replace("", pd.NA).dropna()
    if non_null.empty:
        return False
    return pd.to_numeric(non_null, errors="coerce").notna().mean() >= NUMERIC_COL_THRESHOLD


def _is_date_like(s: pd.Series) -> bool:
    non_null = s.replace("", pd.NA).dropna().astype(str)
    if non_null.empty:
        return False
    return non_null.str.match(r"^\d{4}-\d{2}-\d{2}").mean() > 0.6


def _is_numeric_list(s: pd.Series) -> bool:
    non_null = s.replace("", pd.NA).dropna().astype(str).head(50)
    if non_null.empty:
        return False
    return non_null.apply(
        lambda v: all(
            t.strip().lstrip("-").replace(".", "", 1).isdigit()
            for t in v.split(",") if t.strip()
        )
    ).all()


def classify_label_column(col: str, s: pd.Series) -> str:
    """Return continuous/boolean/categorical/tag_split/skip for label columns."""
    clean = s.replace("", pd.NA)
    n_non_null = clean.notna().sum()
    n_unique = clean.nunique(dropna=True)
    n_total = len(s)

    if n_non_null == 0:
        return "skip"
    if n_unique == 1:
        # Keep query as boolean even if all False, because it is a structural UI column.
        if col == COL_QUERY and _is_boolean(s):
            return "boolean"
        return "skip"
    if _is_numeric(s):
        return "continuous"
    if _is_boolean(s):
        return "boolean"

    non_null_str = clean.dropna().astype(str)
    mean_len = non_null_str.str.len().mean()
    unique_ratio = n_unique / max(n_total, 1)
    comma_frac = non_null_str.str.contains(",").mean()

    if _is_date_like(s):
        return "skip"
    if comma_frac >= COMMA_THRESHOLD and _is_numeric_list(s):
        return "skip"
    if unique_ratio >= HIGH_CARD_UNIQUE_RATIO and mean_len >= HIGH_CARD_MEAN_LEN:
        return "skip"
    if n_unique > MAX_CAT_UNIQUE:
        return "skip"
    if comma_frac >= COMMA_THRESHOLD:
        return "tag_split"
    if comma_frac >= COMMA_THRESHOLD_LOW:
        all_tags: set[str] = set()
        for v in non_null_str:
            all_tags.update(t.strip() for t in v.split(",") if t.strip())
        if len(all_tags) < n_unique:
            return "tag_split"
    return "categorical"


# ---------------------------------------------------------------------------
# High-cardinality categorical rescue
# ---------------------------------------------------------------------------
# classify_label_column skips anything above MAX_CAT_UNIQUE because a flat
# dropdown with hundreds of entries is unusable. But some of those columns
# (e.g. genus, species) are legitimate categories -- they are only large, not
# free text -- and become usable once a *related* column is filtered, because
# narrow_cat_options (in the visualizer) then shrinks them to a handful of
# live values. So an oversized-but-legitimate column is "rescued" back to
# categorical only when it has a nesting partner: another categorical column
# it nests under, or that nests under it. A column with no such partner would
# just be unfilterable clutter and stays skipped. This never assumes an order
# between columns and never inspects column names (taxonomy ranks are not
# special-cased; they are rescued because they happen to nest cleanly).

def _is_high_card_categorical(s: pd.Series) -> bool:
    """True if a column is a genuine categorical that was only skipped for
    being above MAX_CAT_UNIQUE -- not numeric, boolean, date, numeric-list, or
    free text. Unusable as a flat dropdown alone; usable once narrowed.
    """
    clean = s.replace("", pd.NA)
    nn = clean.notna().sum()
    nu = clean.nunique(dropna=True)
    nt = len(s)
    if nn == 0 or nu <= 1 or nt == 0:
        return False
    if nu <= MAX_CAT_UNIQUE:
        return False  # would already be a normal categorical
    if _is_numeric(s) or _is_boolean(s) or _is_date_like(s):
        return False
    if _is_numeric_list(s):
        return False
    non_null = clean.dropna().astype(str)
    mean_len = non_null.str.len().mean()
    unique_ratio = nu / nt
    if unique_ratio >= HIGH_CARD_UNIQUE_RATIO and mean_len >= HIGH_CARD_MEAN_LEN:
        return False  # free-text guard
    return True


def _nests(df: pd.DataFrame, child: str, parent: str) -> bool:
    """True if `child` is a strict refinement of `parent`: child has more
    unique values, and >= 98% of child values map to exactly one parent value
    (a functional dependency child -> parent). Order-agnostic.
    """
    n = len(df)
    child_s = df[child].replace("", pd.NA)
    parent_s = df[parent].replace("", pd.NA)
    if child_s.nunique(dropna=True) <= parent_s.nunique(dropna=True):
        return False
    sub = pd.DataFrame({"child": child_s, "parent": parent_s}).dropna()
    if len(sub) < _NEST_MIN_PAIR_OVERLAP * n:
        return False
    parents_per_child = sub.groupby("child")["parent"].nunique()
    return (parents_per_child <= 1).mean() >= _NEST_MIN_CLEAN_FRAC


def has_nesting_partner(df: pd.DataFrame, col: str, partner_cols: List[str]) -> bool:
    """True if `col` nests with at least one column in `partner_cols`, either
    as parent or as child.
    """
    if df.empty or col not in df.columns:
        return False
    s = df[col].replace("", pd.NA)
    if s.notna().mean() < _NEST_MIN_COVERAGE:
        return False
    for other in partner_cols:
        if other == col or other not in df.columns:
            continue
        other_s = df[other].replace("", pd.NA)
        if other_s.notna().mean() < _NEST_MIN_COVERAGE:
            continue
        if _nests(df, col, other) or _nests(df, other, col):
            return True
    return False


def _rescue_high_card_categoricals(df: pd.DataFrame, types: Dict[str, str], meta: Dict[str, dict]) -> None:
    """Promote skipped label columns back to categorical when they nest with
    an already-categorical/tag_split column. Iterates so a chain like
    family -> genus -> species is fully recovered: genus rescues via family,
    then species rescues via the now-categorical genus.
    """
    label_cols = [c for c, t in types.items() if t == TYPE_LABEL]
    candidates = [c for c in label_cols
                  if meta[c]["type"] == "skip" and _is_high_card_categorical(df[c])]
    partner_pool = [c for c in label_cols if meta[c]["type"] in ("categorical", "tag_split")]
    changed = True
    while changed and candidates:
        changed = False
        still_skipped = []
        for col in candidates:
            if has_nesting_partner(df, col, partner_pool):
                meta[col]["type"] = "categorical"
                partner_pool.append(col)
                changed = True
            else:
                still_skipped.append(col)
        candidates = still_skipped


def build_col_meta(df: pd.DataFrame, types: Dict[str, str]) -> Dict[str, dict]:
    meta: Dict[str, dict] = {}
    for col in df.columns:
        t = types[col]
        if t == TYPE_ID:
            meta[col] = {"type": "id", "tags": [], "override": None, "sse_type": t}
        elif t == TYPE_COORDINATE:
            meta[col] = {"type": "coordinate", "tags": [], "override": "coordinate", "sse_type": t}
        else:
            col_type = classify_label_column(col, df[col])
            tags: list[str] = []
            if col_type == "tag_split":
                all_tags: set[str] = set()
                for val in df[col].replace("", pd.NA).dropna().astype(str):
                    all_tags.update(x.strip() for x in val.split(",") if x.strip())
                tags = sorted(all_tags)
            meta[col] = {"type": col_type, "tags": tags, "override": None, "sse_type": t}
    _rescue_high_card_categoricals(df, types, meta)
    return meta


def boolean_mask(s: pd.Series) -> pd.Series:
    return (
        s.astype(str).str.strip().str.lower()
        .str.replace(r"\.0$", "", regex=True)
        .map({"true": True, "1": True, "yes": True,
              "false": False, "0": False, "no": False})
    )


def query_ids_from_df(df: pd.DataFrame) -> List[str]:
    if COL_QUERY not in df.columns or COL_ID not in df.columns:
        return []
    q = boolean_mask(df[COL_QUERY]).fillna(False)
    return df.loc[q, COL_ID].astype(str).tolist()


def detect_name_cols(df: pd.DataFrame, id_col: str, label_cols: Iterable[str]) -> List[str]:
    out = []
    for c in label_cols:
        if c == id_col:
            continue
        if any(kw in c.lower() for kw in NAME_KEYWORDS):
            if df[c].dtype == object and df[c].nunique(dropna=True) <= MAX_CAT_UNIQUE:
                out.append(c)
    return out


def coordinate_system_key(col: str) -> str:
    """Prefix up to the trailing integer, e.g. esmc_PC10 -> esmc_PC."""
    m = re.match(r"^(.*?)(\d+)$", col)
    return m.group(1) if m else col


def coordinate_axis_number(col: str) -> int:
    m = re.search(r"(\d+)$", col)
    return int(m.group(1)) if m else 0


def group_coordinate_systems(coord_cols: Iterable[str]) -> Dict[str, List[str]]:
    systems: Dict[str, List[str]] = {}
    for c in coord_cols:
        systems.setdefault(coordinate_system_key(c), []).append(c)
    for k in list(systems):
        systems[k] = sorted(systems[k], key=lambda c: (coordinate_axis_number(c), c))
    return dict(sorted(systems.items()))


def default_axes(coord_systems: Dict[str, List[str]]) -> Tuple[str | None, str | None]:
    for cols in coord_systems.values():
        if len(cols) >= 2:
            return cols[0], cols[1]
    cols = [c for group in coord_systems.values() for c in group]
    if len(cols) == 1:
        return cols[0], cols[0]
    return None, None


def load_visualizer_state(entry: EntryContext) -> VisualizerState:
    df, types = read_datafile(entry.datafile_path)
    id_col = id_column(types)
    if id_col != COL_ID:
        # The spec standardizes the output name as id, but tolerate future variants.
        pass
    label_cols = [c for c, t in types.items() if t == TYPE_LABEL]
    coord_cols = [c for c, t in types.items() if t == TYPE_COORDINATE]
    col_meta = build_col_meta(df, types)
    coord_systems = group_coordinate_systems(coord_cols)
    x_col, y_col = default_axes(coord_systems)
    warning = ""
    if not coord_cols:
        warning = "No coordinate columns found yet. Run sse_coordinates.py, then reload the datafile."
    return VisualizerState(
        entry=entry,
        df=df,
        types=types,
        id_col=id_col,
        label_cols=label_cols,
        coord_cols=coord_cols,
        col_meta=col_meta,
        query_ids=query_ids_from_df(df),
        name_cols=detect_name_cols(df, id_col, label_cols),
        coord_systems=coord_systems,
        x_col=x_col,
        y_col=y_col,
        warning=warning,
    )


def axis_range(df: pd.DataFrame, x_col: str | None, y_col: str | None, pad_fraction: float = 0.05):
    if not x_col or not y_col or x_col not in df.columns or y_col not in df.columns:
        return None
    xv = pd.to_numeric(df[x_col], errors="coerce")
    yv = pd.to_numeric(df[y_col], errors="coerce")
    m = xv.notna() & yv.notna()
    if not m.any():
        return None
    xmin, xmax = float(xv[m].min()), float(xv[m].max())
    ymin, ymax = float(yv[m].min()), float(yv[m].max())
    xpad = (xmax - xmin) * pad_fraction if not np.isclose(xmin, xmax) else 1.0
    ypad = (ymax - ymin) * pad_fraction if not np.isclose(ymin, ymax) else 1.0
    return {"xrange": [xmin - xpad, xmax + xpad], "yrange": [ymin - ypad, ymax + ypad]}
