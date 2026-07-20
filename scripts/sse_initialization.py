#!/usr/bin/env python3
"""
sse_initialization.py — create an SSE datafile (an "entry") from a source file.

Bootstrap-only creation tool. Turns one source file into
`entries/<stem>/<stem>.sse.tsv` plus the entry's subfolders and a provenance
manifest. Everything after creation (coordinates, RMSD, third-party labels) is
added later by separate tools via the merge contract — never by re-running this
script. See docs/SSE_datafile_spec.md.

This script owns the creation pipeline: initialise folders, run a reader, run
feature computation, assemble and write the datafile + manifest. The readers and
the feature computation live in the importable `sse_tools/` library; to support
a new source format you add a reader there (see sse_tools/readers/), and never
touch this script.
"""

import argparse
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

# --- make the top-level sse_tools/ package importable regardless of cwd --------
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sse_tools.common import (SSEError, abort, COL_ID, COL_SEQ, COL_QUERY,
                              TYPE_ID, TYPE_LABEL, write_datafile,
                              resolve_input_path, entry_dir,
                              ENTRIES_DIR, INITIAL_FILES_DIR,
                              column_coverage, make_column_entry, write_manifest)
from sse_tools.readers import REGISTRY, ReaderResult
from sse_tools.compute_seq_features import (compute_seq_features, is_usable_sequence,
                                            FEATURE_COLS, INT_FEATURES)

SCHEMA_VERSION = 1
TOOL_NAME = "sse_initialization"
TOOL_VERSION = "1.0.0"
DEFAULT_ID_COL = "Accession"
DEFAULT_SEQ_COL = "Sequence"


# ----------------------------------------------------------- creation pipeline

def _resolve_query_mask(result: ReaderResult, query_values) -> pd.Series:
    """User --query OVERRIDES source auto-detection (the user knows their data)."""
    df = result.table
    if query_values:
        qvals = set(query_values)
        match_vals = df[result.match_col].astype(str)
        missing = sorted(qvals - set(match_vals))
        if missing:
            abort(f"--query value(s) not found in {result.match_label}: "
                  f"{missing}. Nothing was written.")
        for v in sorted(qvals):
            n = int((match_vals == v).sum())
            if n > 1:
                print(f"  note: --query {v!r} matched {n} rows (all flagged).")
        return match_vals.isin(qvals).reset_index(drop=True)
    if result.auto_query is not None:
        return result.auto_query.reset_index(drop=True)
    return pd.Series([False] * len(df))


def build_entry(result: ReaderResult, query_values, report: dict) -> pd.DataFrame:
    """Source-agnostic pipeline -> final datafile DataFrame (without Type row)."""
    df = result.table.reset_index(drop=True).copy()
    id_col, seq_col = result.id_col, result.seq_col

    # 1. Uniqueness of resolved IDs (spec §5).
    dups = df[id_col][df[id_col].duplicated(keep=False)]
    if not dups.empty:
        sample = sorted(set(dups))[:10]
        abort(f"duplicate IDs after resolution ({dups.nunique()} distinct): "
              f"{sample}{' ...' if dups.nunique() > 10 else ''}. Nothing written.")

    # 2. Query membership BEFORE any drop (spec §6.1).
    df[COL_QUERY] = _resolve_query_mask(result, query_values).values

    # 3. Structural sequence check.
    if seq_col not in df.columns or df[seq_col].fillna("").str.strip().eq("").all():
        abort("no sequence information in source. Nothing written.")

    # 4. Per-row usability, split drop rule (spec §6.1.4).
    usable = df[seq_col].fillna("").apply(is_usable_sequence)
    bad_query = (~usable) & df[COL_QUERY]
    if bad_query.any():
        ids = df.loc[bad_query, id_col].tolist()
        abort(f"query row(s) have no usable sequence: {ids[:10]}"
              f"{' ...' if len(ids) > 10 else ''}. A reference protein cannot be "
              "dropped. Nothing written.")
    dropped = df.loc[~usable, id_col].tolist()
    if dropped:
        report["dropped_count"] = len(dropped)
        report["dropped_ids"] = dropped
        print(f"  dropped {len(dropped)} row(s) with unusable sequence "
              f"(reported in manifest).")
    df = df.loc[usable].reset_index(drop=True)

    # 5. Dangling Closest query check for EM (spec §6.1 note).
    if result.source == "em" and "Closest query" in df.columns:
        surviving = set(df[id_col])
        referenced = {v.strip() for v in df["Closest query"] if v and v.strip()}
        dangling = sorted(referenced - surviving)
        if dangling:
            abort("a surviving row's 'Closest query' references dropped/absent "
                  f"ID(s): {dangling[:10]}. A referenced query must have a usable "
                  "sequence. Nothing written.")

    # 6. Features (every surviving row has a usable sequence).
    feats = df[seq_col].apply(compute_seq_features).apply(pd.Series)
    for int_col in INT_FEATURES:
        if int_col in feats.columns:
            feats[int_col] = feats[int_col].astype("Int64")  # clean integer render
    df = pd.concat([df, feats], axis=1)

    # 7. Assemble final column order and reserved names.
    df = df.rename(columns={id_col: COL_ID, seq_col: COL_SEQ})
    df = df.drop(columns=[c for c in df.columns if c.startswith("_")],
                 errors="ignore")
    df[COL_QUERY] = df[COL_QUERY].map(lambda b: "True" if b else "False")

    reserved = [COL_ID, COL_SEQ, COL_QUERY]
    source_labels = [c for c in df.columns
                     if c not in reserved and c not in FEATURE_COLS]
    ordered = reserved + source_labels + [c for c in FEATURE_COLS if c in df.columns]
    df = df[ordered]

    report["rows_kept"] = len(df)
    report["query_count"] = int((df[COL_QUERY] == "True").sum())
    report["columns"] = ordered
    return df


def build_manifest(df: pd.DataFrame, types: dict, result: ReaderResult,
                   report: dict, source_file: Path) -> dict:
    total = len(df)
    cov = column_coverage(df)
    columns = [make_column_entry(c, types[c], provenance_source="sse",
                                 tool=TOOL_NAME, version=TOOL_VERSION,
                                 coverage=cov[c]) for c in df.columns]
    return {
        "schema_version": SCHEMA_VERSION,
        "source_file": source_file.name,
        "source_type": result.source,
        "creation_tool": TOOL_NAME,
        "creation_tool_version": TOOL_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "id_resolution": result.notes,
        "row_counts": {
            "kept": report.get("rows_kept", total),
            "dropped": report.get("dropped_count", 0),
            "dropped_reason": "unusable sequence" if report.get("dropped_count") else "",
            "queries": report.get("query_count", 0),
        },
        "dropped_ids": report.get("dropped_ids", []),
        "columns": columns,
    }


# --------------------------------------------------------------------- the CLI

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sse_initialization.py",
        description="Create an SSE datafile (entry) from a source file. "
                    "Bootstrap-only: never re-run to update an entry — use the "
                    "merge tools instead.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="examples:\n"
               "  sse_initialization.py my_enzymes.tsv --source em\n"
               "  sse_initialization.py search.json --source fs --name oleD_hits\n"
               "  sse_initialization.py seqs.fasta --source fasta --query "
               "'sp|P12345|OLED_STRAN OleD'\n"
               "  sse_initialization.py /path/to/data.tsv --source em "
               "--id_col ProteinID --seq_col Seq --force\n",
    )
    parser.add_argument("input_file",
                        help="source file: a bare name (looked up in "
                             "initial_files/) or an explicit/relative path.")
    parser.add_argument("--source", required=True, choices=sorted(REGISTRY),
                        help="source type / reader to use (mandatory).")
    parser.add_argument("--id_col", default=DEFAULT_ID_COL,
                        help=f"em/tabular only: ID column (default {DEFAULT_ID_COL!r}).")
    parser.add_argument("--seq_col", default=DEFAULT_SEQ_COL,
                        help=f"em/tabular only: sequence column (default {DEFAULT_SEQ_COL!r}).")
    parser.add_argument("--query", nargs="+", metavar="VALUE",
                        help="mark query sequences. Matches raw --id_col (em), the "
                             "full header (fasta), or raw target (fs). OVERRIDES "
                             "source auto-detection. Unmatched values abort.")
    parser.add_argument("--name", metavar="STEM",
                        help="override the entry stem (default: source filename stem).")
    parser.add_argument("--entries-dir", metavar="DIR", help="override the entries/ root.")
    parser.add_argument("--initial-files-dir", metavar="DIR",
                        help="override the initial_files/ lookup dir for bare names.")
    parser.add_argument("--force", action="store_true",
                        help="delete and rebuild an existing entry.")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    entries_dir = Path(args.entries_dir) if args.entries_dir else ENTRIES_DIR
    initial_dir = Path(args.initial_files_dir) if args.initial_files_dir else INITIAL_FILES_DIR

    try:
        src = resolve_input_path(args.input_file, initial_dir)
        stem = args.name or src.stem
        edir = entry_dir(stem, entries_dir)

        if edir.exists() and not args.force:
            abort(f"entry already exists: {edir}. Creation is bootstrap-only "
                  "— use --force to delete and rebuild, or a merge tool to add "
                  "columns.")

        reader = REGISTRY[args.source]
        print(f"reading {src} as --source {args.source} ...")
        result = reader(src, args)

        report = {}
        df = build_entry(result, args.query, report)
        types = {c: (TYPE_ID if c == COL_ID else TYPE_LABEL) for c in df.columns}
        manifest = build_manifest(df, types, result, report, src)

        # Atomic: build the whole entry in a temp dir, then move it into place.
        entries_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=entries_dir) as tmp:
            tmp = Path(tmp)
            for sub in ("external_data", "structures", "figures", "msa_cache", "logs"):
                (tmp / sub).mkdir()
            write_datafile(df, types, tmp / f"{stem}.sse.tsv")
            write_manifest(tmp / "logs" / f"{stem}.sse.manifest.json",
                           manifest, TOOL_NAME)
            if edir.exists():
                shutil.rmtree(edir)
            shutil.move(str(tmp), str(edir))

    except SSEError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    print(f"created {edir / (stem + '.sse.tsv')}")
    print(f"  rows: {report['rows_kept']}  queries: {report['query_count']}  "
          f"columns: {len(report['columns'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
