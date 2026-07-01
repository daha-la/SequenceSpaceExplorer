"""Foldseek webserver JSON reader — the worked example for a custom reader.

This is the template to copy when writing a reader for a non-tabular source:
it parses a nested structure, deduplicates, aggregates a tag-split column, and
resolves a query, then returns a flat table like any other reader.
"""

import json
from collections import OrderedDict

import pandas as pd

from .base import ReaderResult
from ..common import abort, COL_ID, COL_SEQ


def read_foldseek(path, args) -> ReaderResult:
    """Foldseek webserver JSON.

    The JSON is a one-element top list whose [0] holds `queries` (header +
    sequence of the search protein) and `results` (one block per database, each
    with an `alignments` dict of hit records). This reader:
      * collects one row per unique `target` (dedup across databases),
      * aggregates the databases each target hit into a `Databases` tag-split,
      * resolves the query: a hit with seqId==100 whose tSeq equals the query
        sequence IS the query (flagged, no separate row); if none exists, a
        synthetic row is added from queries[0] using the header token as ID
        (spec §5.2, §8.2).
    `--query` matches the raw `target` (== the resolved ID).
    """
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, list) or not data:
        abort("--source fs: unexpected JSON shape (expected a non-empty list).")
    block = data[0]
    if "results" not in block or "queries" not in block:
        abort("--source fs: JSON missing 'results'/'queries' (not Foldseek output?).")

    qrec = block["queries"][0]
    query_seq = (qrec.get("sequence") or "").replace("\n", "").strip().upper()
    query_token = (qrec.get("header") or "query").split(None, 1)[0]

    targets = OrderedDict()           # target -> row dict
    self_hit_targets = []
    for rblock in block["results"]:
        db = rblock.get("db") or ""
        for hit_list in rblock.get("alignments", {}).values():
            for hit in hit_list:
                tgt = hit.get("target")
                if tgt is None:
                    continue
                tseq = (hit.get("tSeq") or "").replace("\n", "").strip().upper()
                row = targets.get(tgt)
                if row is None:
                    row = {COL_ID: tgt, COL_SEQ: tseq, "_dbs": set(),
                           "taxName": hit.get("taxName") or "",
                           "taxId": str(hit.get("taxId") or ""),
                           "Description": hit.get("description") or ""}
                    targets[tgt] = row
                row["_dbs"].add(db)
                if hit.get("seqId") == 100 and tseq and tseq == query_seq:
                    self_hit_targets.append(tgt)

    rows = list(targets.values())
    notes = {"unique_targets": len(rows),
             "self_hit": self_hit_targets[0] if self_hit_targets else None}

    # Query resolution.
    if self_hit_targets:
        query_ids = set(self_hit_targets)            # existing rows are the query
    else:
        # No self-hit: add a synthetic query row from queries[0].
        if not query_seq:
            abort("--source fs: query has no sequence and no seqId=100 self-hit; "
                  "cannot create the reference row.")
        rows.append({COL_ID: query_token, COL_SEQ: query_seq, "_dbs": set(),
                     "taxName": "", "taxId": "",
                     "Description": qrec.get("header") or ""})
        query_ids = {query_token}
        notes["synthetic_query_id"] = query_token

    df = pd.DataFrame(rows)
    df["Databases"] = df["_dbs"].apply(
        lambda s: ", ".join(sorted(x for x in s if x)))
    df = df.drop(columns=["_dbs"])
    auto_query = df[COL_ID].isin(query_ids)

    return ReaderResult(table=df, id_col=COL_ID, seq_col=COL_SEQ, source="fs",
                        match_col=COL_ID, match_label="target",
                        auto_query=auto_query, notes=notes)
