#!/usr/bin/env python3
"""Append distances to query sequence(s) in cached ESM-C embedding space."""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sse_tools.common import (SSEError, TYPE_LABEL, merge_columns,
                              read_datafile, resolve_entry_path,
                              resolve_embedding_cache)
from sse_tools.esmc_distance import (distance_columns, query_ids,
                                     read_embedding_matrix)

TOOL_NAME = "sse_esmc_distance"
TOOL_VERSION = "1.0.0"


def build_parser():
    parser = argparse.ArgumentParser(
        description="Append Euclidean distances from cached ESM-C embeddings "
                    "to the entry's query sequence(s)."
    )
    parser.add_argument("entry", help="entry stem, entry directory, or .sse.tsv path")
    parser.add_argument("--embedding", help="embedding tag or path (auto if exactly one exists)")
    parser.add_argument("--raw", action="store_true",
                        help="use the raw embedding cache even if an L2-normalized one exists")
    parser.add_argument("--query-id", nargs="+", help="query ID(s); default: rows marked query=True")
    parser.add_argument("--force", action="store_true", help="overwrite existing distance columns")
    parser.add_argument("--entries-dir")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        datafile = resolve_entry_path(args.entry, args.entries_dir)
        df, _ = read_datafile(datafile)
        queries = query_ids(df, args.query_id)
        emb_path, normalized = resolve_embedding_cache(
            datafile.parent, args.embedding, prefer_normalized=not args.raw)
        print(f"using {'L2-normalized' if normalized else 'raw'} "
              f"embeddings: {emb_path}")
        matrix, dimensions = read_embedding_matrix(emb_path)
        tag = emb_path.name.removesuffix(".emb.tsv")
        result = distance_columns(matrix, dimensions, queries, tag)
        distance_names = [c for c in result.columns if c != "id"]

        stem = datafile.name.removesuffix(".sse.tsv")
        manifest = datafile.parent / "logs" / f"{stem}.sse.manifest.json"
        geometry = "l2" if normalized else "none"
        params = (f"embedding={tag};metric=euclidean;normalize={geometry};"
                  f"dimensions={len(dimensions)};queries={','.join(queries)}")
        notes = ("Euclidean distance in "
                 + ("L2-normalized " if normalized else "")
                 + "pooled ESM-C embedding space")
        merge_columns(
            datafile, result, {c: TYPE_LABEL for c in distance_names},
            manifest_path=manifest if manifest.exists() else None,
            provenance_source="sse", tool=TOOL_NAME, version=TOOL_VERSION,
            params=params, notes=notes,
            force=args.force,
        )
    except SSEError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"merged {len(distance_names)} distance column(s) for "
          f"{len(result)}/{len(df)} embedded rows into {datafile}")
    for query, column in zip(queries, distance_names):
        print(f"  {query}: {column}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
