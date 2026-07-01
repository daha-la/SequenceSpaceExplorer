"""EnzymeMiner / generic tabular reader."""

import pandas as pd

from .base import ReaderResult
from ..common import abort

DEFAULT_ID_COL = "Accession"
DEFAULT_SEQ_COL = "Sequence"


def read_em(path, args) -> ReaderResult:
    """EnzymeMiner / generic tabular TSV.

    ID and sequence come from named columns (--id_col / --seq_col, defaulting to
    Accession / Sequence). Queries are auto-derived from `Closest query`: a row
    is a query if its ID appears as a *value* in that column (spec §8.1).
    `--query` matches the raw --id_col value.
    """
    id_col = getattr(args, "id_col", DEFAULT_ID_COL)
    seq_col = getattr(args, "seq_col", DEFAULT_SEQ_COL)
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                     encoding="utf-8-sig")
    df.columns = [c.strip() for c in df.columns]

    if id_col not in df.columns:
        abort(f"--source em: ID column {id_col!r} not found. "
              f"Columns present: {list(df.columns)}. Use --id_col to set it.")
    if seq_col not in df.columns:
        # No sequence information at all -> structural abort (spec §6.1.3).
        abort(f"--source em: sequence column {seq_col!r} not found (no sequence "
              f"information). Columns present: {list(df.columns)}. Use --seq_col.")

    auto_query = None
    if "Closest query" in df.columns:
        query_ids = {v.strip() for v in df["Closest query"] if v and v.strip()}
        auto_query = df[id_col].isin(query_ids)

    return ReaderResult(table=df, id_col=id_col, seq_col=seq_col, source="em",
                        match_col=id_col, match_label=f"{id_col} (ID column)",
                        auto_query=auto_query,
                        notes={"id_col": id_col, "seq_col": seq_col})
