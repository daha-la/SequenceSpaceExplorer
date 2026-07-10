#!/usr/bin/env python3
"""Merge an external dataset's columns into an SSE entry's datafile.

Reads a TSV/CSV with its own id column (default: its first column) and
merges every other column (or a chosen subset) into the entry's .sse.tsv via
common.merge_columns -- the only write path onto a datafile. external_file
and --translator each accept a full/relative path, or a bare filename looked
up in the entry's external_data/ directory (common.resolve_input_path, the
same convention used elsewhere in the codebase). External ids are matched
against the datafile's id column by exact string equality; an optional
translator table remaps external ids onto datafile ids first, for cases
where the two sources don't share an id scheme.

External rows whose id has no match in the datafile are warned about and
dropped (unmatched datafile rows are already handled by merge_columns' left
join -- they're kept, with the new columns left blank). Column-name
collisions with the datafile follow the same --force convention as
fetch_taxonomy.py: merge_columns aborts and names the colliding column(s)
unless --force is passed.

Usage:
    python merge_external.py ENTRY external_data.csv
    python merge_external.py ENTRY external_data.tsv --id-col Accession \\
        --columns pI,MW --translator id_map.tsv --force
    # bare filenames resolve against entries/<stem>/external_data/:
    python merge_external.py ENTRY dummy_biophysical.csv
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sse_tools import common


def sniff_delimiter(path, explicit=None):
    if explicit:
        return explicit
    ext = Path(path).suffix.lower()
    if ext == ".tsv":
        return "\t"
    if ext == ".csv":
        return ","
    common.abort(f"cannot infer delimiter from extension {ext!r} on {path}; "
                 "pass --delimiter.")


def load_translator(path, delimiter=None):
    """Two-column id-translation table: column 0 = datafile id, column 1 =
    external id. Positional, not by header name -- the two sources rarely
    agree on what to call their id columns, so this avoids requiring
    specific header text. Every row must be non-empty in both columns.
    """
    delim = sniff_delimiter(path, delimiter)
    tdf = pd.read_csv(path, sep=delim, dtype=str, keep_default_na=False)
    if tdf.shape[1] < 2:
        common.abort(f"translator table must have >= 2 columns: {path}")
    tdf = tdf.iloc[:, :2]
    tdf.columns = ["_datafile_id", "_external_id"]
    blank = tdf[(tdf["_datafile_id"] == "") | (tdf["_external_id"] == "")]
    if not blank.empty:
        common.abort(f"translator table has {len(blank)} row(s) with a blank "
                     f"id in one column: {path}")
    dupes = tdf["_external_id"][tdf["_external_id"].duplicated()].tolist()
    if dupes:
        common.abort(f"translator table has duplicate external id(s): {sorted(set(dupes))[:20]}")
    return dict(zip(tdf["_external_id"], tdf["_datafile_id"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("entry", help="entry stem (in entries/) or path to a .sse.tsv")
    ap.add_argument("external_file",
                    help="TSV/CSV of external data to merge in -- a full/relative "
                         "path, or a bare filename looked up in the entry's "
                         "external_data/ directory")
    ap.add_argument("--id-col", default=None,
                    help="column in external_file to use as its id (default: first column)")
    ap.add_argument("--columns", default=None,
                    help="comma-separated external columns to merge (default: all "
                         "columns except the id column)")
    ap.add_argument("--translator", default=None,
                    help="tsv/csv mapping datafile id (col 1) -> external id (col 2). "
                         "Same path resolution as external_file: full path, or a bare "
                         "filename looked up in the entry's external_data/ directory")
    ap.add_argument("--type", choices=[common.TYPE_LABEL, common.TYPE_COORDINATE],
                    default=common.TYPE_LABEL,
                    help="Type token applied to every merged column (default: label)")
    ap.add_argument("--delimiter", default=None,
                    help="override delimiter auto-detection (applies to both "
                         "external_file and --translator)")
    ap.add_argument("--force", action="store_true",
                    help="overwrite colliding column(s) already on the datafile")
    ap.add_argument("--entries-dir", default=None)
    args = ap.parse_args()

    datafile_path = common.resolve_entry_path(args.entry, args.entries_dir)
    stem = datafile_path.stem.removesuffix(".sse")
    manifest_path = datafile_path.parent / "logs" / f"{stem}.sse.manifest.json"
    external_data_dir = datafile_path.parent / "external_data"

    df, types = common.read_datafile(datafile_path)
    id_col = common.id_column(types)

    external_path = common.resolve_input_path(args.external_file, external_data_dir)
    delim = sniff_delimiter(external_path, args.delimiter)
    ext = pd.read_csv(external_path, sep=delim, dtype=str, keep_default_na=False)
    if ext.empty:
        common.abort(f"external file has no data rows: {external_path}")

    ext_id_col = args.id_col or ext.columns[0]
    if ext_id_col not in ext.columns:
        common.abort(f"--id-col {ext_id_col!r} not found in {external_path}. "
                     f"Columns: {list(ext.columns)}")

    if args.columns:
        wanted = [c.strip() for c in args.columns.split(",") if c.strip()]
        missing = [c for c in wanted if c not in ext.columns]
        if missing:
            common.abort(f"--columns name(s) not found in {external_path}: {missing}")
        merge_cols = wanted
    else:
        merge_cols = [c for c in ext.columns if c != ext_id_col]
    if not merge_cols:
        common.abort("nothing to merge: no columns besides the id column.")

    ext = ext[[ext_id_col] + merge_cols].copy()

    # translate external ids onto datafile ids, or use them as-is
    if args.translator:
        translator_path = common.resolve_input_path(args.translator, external_data_dir)
        translator = load_translator(translator_path, args.delimiter)
        ext["_mapped_id"] = ext[ext_id_col].map(translator)
        untranslated = ext[ext["_mapped_id"].isna()]
        if not untranslated.empty:
            sample = untranslated[ext_id_col].head(20).tolist()
            print(f"warning: {len(untranslated)} external id(s) not found in "
                  f"--translator, dropped: {sample}"
                  f"{' ...' if len(untranslated) > 20 else ''}")
        ext = ext.dropna(subset=["_mapped_id"])
        ext[ext_id_col] = ext["_mapped_id"]
        ext = ext.drop(columns=["_mapped_id"])
    ext = ext.rename(columns={ext_id_col: id_col})

    # warn-and-drop external rows whose (possibly translated) id has no
    # match in the datafile -- exact string match, no normalization.
    datafile_ids = set(df[id_col].astype(str))
    ext[id_col] = ext[id_col].astype(str)
    unmatched = ext[~ext[id_col].isin(datafile_ids)]
    if not unmatched.empty:
        sample = unmatched[id_col].head(20).tolist()
        print(f"warning: {len(unmatched)} external row(s) have no matching "
              f"datafile id, dropped: {sample}{' ...' if len(unmatched) > 20 else ''}")
        ext = ext[ext[id_col].isin(datafile_ids)]
    if ext.empty:
        common.abort("no external rows matched a datafile id; nothing to merge.")

    new_types = {id_col: common.TYPE_ID, **{c: args.type for c in merge_cols}}
    source_name = external_path.name

    merged = common.merge_columns(
        datafile_path, ext, new_types,
        manifest_path=manifest_path if manifest_path.exists() else None,
        provenance_source="external", tool="merge_external", version="1.0",
        params=f"source={source_name},id_col={ext_id_col},"
               f"translator={bool(args.translator)},type={args.type}",
        notes=f"merged from {source_name}",
        id_col=id_col, force=args.force,
    )

    print(f"done; merged {len(merge_cols)} column(s) for {len(ext)} row(s) "
          f"into {datafile_path}")
    print(f"  columns: {merge_cols}")


if __name__ == "__main__":
    main()
