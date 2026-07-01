"""The reader contract: what a reader returns (ReaderResult).

Shared utilities (abort, reserved column names) live in sse_tools/common.py;
readers import those from there. This module holds only the reader-specific
contract, so it has no dependency on the rest of the toolchain beyond pandas.

WHAT A READER RETURNS
---------------------
A reader is a callable `reader(path, args) -> ReaderResult` where `path` is the
resolved input file and `args` is the parsed CLI namespace (a reader may read
the flags relevant to it, e.g. args.id_col for tabular sources). It returns a
ReaderResult with:

    table       pd.DataFrame   resolved table. Read every column as str
                               (dtype=str, keep_default_na=False). Must contain
                               id_col and seq_col, both already resolved (IDs
                               final and intended-unique; sequences raw strings).
                               Extra columns become `label` columns. Columns whose
                               name starts with "_" are internal, dropped on write.
    id_col      str            column holding the final unique ID.
    seq_col     str            column holding the amino-acid sequence.
    source      str            short tag for the manifest ("em", "fs", "fasta").
    match_col   str            column --query is tested against (RAW pre-resolution
                               names). Defaults to id_col.
    match_label str            human name of the match field for --query errors.
                               Defaults to match_col.
    auto_query  pd.Series|None boolean mask of source-derived queries aligned to
                               `table`, or None if the source has no native query
                               concept. Overridden entirely when --query is given.
    notes       dict           source facts for the manifest header. May be empty.

The pipeline renames id_col -> "id" and seq_col -> "Sequence" on output, so a
reader may use any internal names. To add a source: write a reader module here,
then register it in __init__.py (one import + one REGISTRY line).
"""

from dataclasses import dataclass, field
from typing import Optional

import pandas as pd


@dataclass
class ReaderResult:
    table: pd.DataFrame
    id_col: str
    seq_col: str
    source: str
    match_col: Optional[str] = None
    match_label: Optional[str] = None
    auto_query: Optional[pd.Series] = None
    notes: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.match_col is None:
            self.match_col = self.id_col
        if self.match_label is None:
            self.match_label = self.match_col
