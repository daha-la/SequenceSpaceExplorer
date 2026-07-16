"""Distances from SSE query sequences in a cached ESM-C embedding matrix."""

from pathlib import Path
import re

import numpy as np
import pandas as pd

from .common import COL_ID, COL_QUERY, abort


def query_ids(df: pd.DataFrame, requested=None) -> list[str]:
    """Resolve requested IDs, or all rows marked by initialization as queries."""
    ids = set(df[COL_ID].astype(str))
    if requested:
        selected = list(dict.fromkeys(map(str, requested)))
        missing = [value for value in selected if value not in ids]
        if missing:
            abort(f"query ID(s) not found in datafile: {missing}")
        return selected

    if COL_QUERY not in df.columns:
        abort(f"datafile has no {COL_QUERY!r} column; pass --query-id explicitly.")
    marked = df[COL_QUERY].astype(str).str.strip().str.lower().isin(
        {"true", "1", "yes", "y"}
    )
    selected = df.loc[marked, COL_ID].astype(str).tolist()
    if not selected:
        abort("datafile has no rows marked as queries; pass --query-id explicitly.")
    return selected


def read_embedding_matrix(path) -> tuple[pd.DataFrame, list[str]]:
    """Read and validate the ID + numeric-dimension embedding cache contract."""
    path = Path(path)
    if not path.exists():
        abort(f"embedding cache not found: {path}")
    matrix = pd.read_csv(path, sep="\t", dtype={"ID": str})
    if "ID" not in matrix.columns:
        abort(f"embedding cache has no 'ID' column: {path}")
    if matrix["ID"].duplicated().any():
        dupes = matrix.loc[matrix["ID"].duplicated(keep=False), "ID"].unique()[:10]
        abort(f"embedding cache contains duplicate IDs: {list(dupes)}")
    dimensions = [column for column in matrix.columns if column != "ID"]
    if not dimensions:
        abort(f"embedding cache contains no embedding dimensions: {path}")
    numeric = matrix[dimensions].apply(pd.to_numeric, errors="coerce")
    bad = numeric.isna().any(axis=1)
    if bad.any():
        ids = matrix.loc[bad, "ID"].astype(str).tolist()[:10]
        abort(f"embedding cache has missing or non-numeric values for ID(s): {ids}")
    matrix[dimensions] = numeric
    return matrix, dimensions


def safe_id(value: str) -> str:
    """Make a stable, readable query-ID fragment for an SSE column name."""
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return cleaned or "query"


def distance_columns(matrix: pd.DataFrame, dimensions: list[str], queries,
                     embedding_tag: str) -> pd.DataFrame:
    """Calculate Euclidean distance from every cached vector to each query."""
    indexed = matrix.set_index("ID", drop=False)
    missing = [query for query in queries if query not in indexed.index]
    if missing:
        abort(f"query ID(s) missing from embedding cache: {missing}")

    names = [f"{embedding_tag}_distance_to_{safe_id(query)}" for query in queries]
    if len(names) != len(set(names)):
        abort("query IDs produce colliding distance column names; select them "
              "individually with --query-id.")

    values = matrix[dimensions].to_numpy(dtype=np.float64)
    out = pd.DataFrame({COL_ID: matrix["ID"].astype(str)})
    for query, name in zip(queries, names):
        reference = indexed.loc[query, dimensions].to_numpy(dtype=np.float64)
        out[name] = np.linalg.norm(values - reference, axis=1)
    return out
