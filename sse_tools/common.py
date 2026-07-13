"""Shared utilities for the SSE tools.

Everything here is tool-agnostic: errors, reserved column names / Type tokens,
datafile writing, path resolution, and manifest I/O. Every script and library
module imports from here rather than reaching into another tool's module. This
is the bottom of the dependency graph — common.py imports nothing else from
sse_tools.
"""

import csv
import json
import os
import threading
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


# --------------------------------------------------------------------- errors

class SSEError(Exception):
    """Raised to stop a tool with a clear, user-facing message. The CLIs catch
    this, print it to stderr, and exit 1. A single project-wide error type:
    these tools print-and-exit, so per-tool subclasses aren't needed."""


def abort(msg: str):
    raise SSEError(msg)


# ------------------------------------------------- reserved names + Type tokens

# Reserved output column names (datafile spec §4).
COL_ID = "id"
COL_SEQ = "Sequence"
COL_QUERY = "query"

# Type-row tokens (datafile spec §3.1).
TYPE_ID = "id"
TYPE_LABEL = "label"
TYPE_COORDINATE = "coordinate"

ALLOWED_TYPE_TOKENS = {TYPE_ID, TYPE_LABEL, TYPE_COORDINATE}

# Process-local guard for read -> write -> os.replace datafile updates.
# This prevents lost updates when app worker threads finish concurrently.
_DATAFILE_WRITE_LOCK = threading.Lock()


def _read_header_columns(path) -> list:
    """Read only the physical header record using the project TSV dialect.

    This catches duplicate column names before pandas can silently mangle them.
    """
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        try:
            return next(csv.reader(fh, delimiter="\t", quotechar='"',
                                   doublequote=True))
        except StopIteration:
            return []


def validate_type_map(types: dict, columns, *, context: str) -> dict:
    """Return a normalized Type map after enforcing the SSE datafile contract."""
    cols = list(columns)
    if len(cols) != len(set(cols)):
        seen, dupes = set(), []
        for c in cols:
            if c in seen and c not in dupes:
                dupes.append(c)
            seen.add(c)
        abort(f"{context}: duplicate column name(s): {dupes}")

    missing = [c for c in cols if c not in types or str(types.get(c, "")).strip() == ""]
    if missing:
        abort(f"{context}: no Type token for column(s): {missing}")

    norm = {c: str(types[c]).strip().lower() for c in cols}
    invalid = {c: t for c, t in norm.items() if t not in ALLOWED_TYPE_TOKENS}
    if invalid:
        allowed = ", ".join(sorted(ALLOWED_TYPE_TOKENS))
        abort(f"{context}: invalid Type token(s) {invalid}; allowed: {allowed}")

    id_cols = [c for c, t in norm.items() if t == TYPE_ID]
    if len(id_cols) != 1:
        abort(f"{context}: expected exactly one '{TYPE_ID}' Type token, found {len(id_cols)}: {id_cols}")
    if norm[cols[0]] != TYPE_ID:
        abort(f"{context}: first column must carry the '{TYPE_ID}' Type token; got {norm[cols[0]]!r}")
    return norm


# ----------------------------------------------------------------------- paths

# common.py lives at <repo>/sse_tools/common.py, so the repo root is two up.
REPO_ROOT = Path(__file__).resolve().parent.parent
ENTRIES_DIR = REPO_ROOT / "entries"
INITIAL_FILES_DIR = REPO_ROOT / "initial_files"


def entry_dir(stem: str, entries_dir=None) -> Path:
    """Path to entries/<stem>/ (under a custom root if given)."""
    return (Path(entries_dir) if entries_dir else ENTRIES_DIR) / stem


def resolve_input_path(name: str, fallback_dir) -> Path:
    """Resolve an input given as either a path (absolute/relative, used as-is if
    it exists) or a bare name looked up in `fallback_dir`. Aborts if neither."""
    p = Path(name)
    if p.exists():
        return p
    cand = Path(fallback_dir) / name
    if cand.exists():
        return cand
    abort(f"input file not found: {name!r}. Looked at {p} and {cand}.")


def resolve_entry_path(arg: str, entries_dir=None) -> Path:
    """Resolve ENTRY to its datafile path. Shared by every script/tool that
    takes an ENTRY argument, so 'point me at an entry' means the same thing
    everywhere: sse_coordinates.py, fetch_taxonomy.py, merge_external.py, and
    (wrapped with its own working-directory setup) the visualizer's
    EntryContext. Previously each had its own slightly different copy of
    this; this is the one place it's implemented.

    ENTRY may be:
      - a direct path to a .sse.tsv file
      - an entry directory containing exactly one .sse.tsv file
      - a bare stem, looked up as <entries_dir>/<stem>/<stem>.sse.tsv, or
        (if that exact name isn't found) the entry directory's sole
        .sse.tsv file
    """
    if not arg:
        abort("No entry supplied.")
    root = Path(entries_dir) if entries_dir else ENTRIES_DIR
    p = Path(arg)

    if p.exists() and p.is_file():
        if not p.name.endswith(".sse.tsv"):
            abort(f"Expected an .sse.tsv datafile, got: {p}")
        return p

    if p.exists() and p.is_dir():
        candidates = sorted(p.glob("*.sse.tsv"))
        if len(candidates) != 1:
            abort(f"Entry directory must contain exactly one .sse.tsv file, "
                  f"found {len(candidates)}: {p}")
        return candidates[0]

    entry_d = root / arg
    if not entry_d.exists():
        abort(f"Entry not found: {entry_d}")
    candidate = entry_d / f"{arg}.sse.tsv"
    if candidate.exists():
        return candidate
    candidates = sorted(entry_d.glob("*.sse.tsv"))
    if len(candidates) == 1:
        return candidates[0]
    abort(f"No datafile found for entry {arg!r}. Expected {candidate}")


# ------------------------------------------------------------- datafile writing

def write_datafile(df: pd.DataFrame, types: dict, path) -> None:
    """Write `df` as an SSE datafile: header row, Type row, then data.

    `types` maps every column name to its Type token (id / label / coordinate).
    Lossless QUOTE_MINIMAL TSV, UTF-8, empty field = null. Reused by every tool
    that writes a datafile (creation, merge, embedding, ...).
    """
    norm_types = validate_type_map(types, df.columns, context="write_datafile")
    type_row = pd.DataFrame([{c: norm_types[c] for c in df.columns}])
    out = pd.concat([type_row, df], ignore_index=True)
    out.to_csv(path, sep="\t", index=False, na_rep="",
               quoting=csv.QUOTE_MINIMAL, encoding="utf-8")


# ------------------------------------------------------------- manifest (logs/)
# The manifest is a JSON record of column provenance (datafile spec §11). It has
# a header (facts about the whole datafile) and a `columns` list (one entry per
# column). Creation builds it; merge reads it, appends column entries, and writes
# it back. The visualizer never depends on it.

def column_coverage(df: pd.DataFrame) -> dict:
    """Per-column 'N/total' populated-cell count (empty string = unpopulated)."""
    total = len(df)
    return {c: f"{int(df[c].replace('', pd.NA).notna().sum())}/{total}"
            for c in df.columns}


def make_column_entry(name: str, col_type: str, *, provenance_source: str = "sse",
                      tool: str = "", version: str = "", params: str = "",
                      notes: str = "", coverage: str = "") -> dict:
    """One column's provenance entry, uniform schema (every key always present;
    blanks mean unknown). `provenance_source` is 'sse' for first-party tools or
    'external' for user-supplied merges."""
    return {
        "name": name,
        "type": col_type,
        "provenance_source": provenance_source,
        "tool": tool,
        "version": version,
        "params": params,
        "notes": notes,
        "coverage": coverage,
    }


def read_manifest(path) -> dict:
    """Load an existing manifest (used by merge). Aborts if missing."""
    p = Path(path)
    if not p.exists():
        abort(f"manifest not found: {p}")
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def write_manifest(path, manifest: dict, tool: str) -> None:
    """Write the manifest, stamping who touched it last. Called by every tool
    that writes the manifest; at creation the stamp equals the creation tool,
    after a merge it shows the merge tool."""
    manifest["last_modified_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["last_tool"] = tool
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)


# ------------------------------------------------------------- datafile reading

def read_datafile(path):
    """Read an SSE datafile. Returns (df, types):
      df    : the data rows (Type row peeled off), index reset, all columns str.
      types : {column -> Type token}.

    Aborts if the Type row is missing or violates the SSE contract: the first
    column's Type-row cell must be 'id', every token must be valid, and exactly
    one id column must exist.
    """
    header_cols = _read_header_columns(path)
    if not header_cols:
        abort(f"datafile is empty or missing a header row: {path}")
    if len(header_cols) != len(set(header_cols)):
        seen, dupes = set(), []
        for c in header_cols:
            if c in seen and c not in dupes:
                dupes.append(c)
            seen.add(c)
        abort(f"datafile has duplicate column name(s): {dupes}")

    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                      encoding="utf-8-sig")
    if raw.empty:
        abort(f"datafile missing Type row: {path}")
    if list(raw.columns) != header_cols:
        abort("datafile header changed during parsing; check for duplicate or malformed columns.")

    type_row = raw.iloc[0]
    tokens = {c: str(type_row[c]).strip().lower() for c in raw.columns}
    if tokens.get(raw.columns[0]) != TYPE_ID:
        abort(f"datafile missing a Type row: first Type cell must be '{TYPE_ID}', "
              f"got {tokens.get(raw.columns[0], '')!r}: {path}")
    tokens = validate_type_map(tokens, raw.columns, context="read_datafile")
    df = raw.iloc[1:].reset_index(drop=True)
    return df, tokens


def id_column(types: dict) -> str:
    """Name of the column whose Type token is 'id'."""
    norm = validate_type_map(types, list(types), context="id_column")
    for c, t in norm.items():
        if t == TYPE_ID:
            return c
    abort("no id column in datafile types")


# ------------------------------------------------------ additive column merge

def _format_duplicate_ids(values) -> str:
    shown = sorted(set(map(str, values)))[:20]
    extra = len(set(map(str, values))) - len(shown)
    msg = str(shown)
    if extra > 0:
        msg += f" (+{extra} more)"
    return msg


def merge_columns(datafile_path, new_df, new_types, *, manifest_path=None,
                  provenance_source="sse", tool="", version="", params="",
                  notes="", id_col=COL_ID, force=False, drop_columns=None):
    """Additive left-join of `new_df` (must contain `id_col` + new columns) into
    the datafile at `datafile_path`. New columns are typed by `new_types`.

    Mechanics shared by every column-adding tool: rows are never reordered,
    removed, or multiplied; unmatched datafile rows get empty cells; a name
    collision aborts unless `force` (then full-column replace). `drop_columns`
    removes an old column system inside the same atomic write, used for replacing
    a whole coordinate system without orphaning PC3..PC10-style columns.

    The read -> write temp -> os.replace critical section is protected by a
    process-local lock so concurrent app worker threads cannot lose each other's
    updates. Cross-process writers should still be avoided. Returns the merged
    DataFrame.
    """
    with _DATAFILE_WRITE_LOCK:
        df, types = read_datafile(datafile_path)
        if id_col not in df.columns:
            abort(f"merge_columns: datafile missing id column {id_col!r}")
        if id_col not in new_df.columns:
            abort(f"merge_columns: new data missing id column {id_col!r}")

        if len(new_df.columns) != len(set(new_df.columns)):
            seen, dupes = set(), []
            for c in new_df.columns:
                if c in seen and c not in dupes:
                    dupes.append(c)
                seen.add(c)
            abort(f"merge_columns: incoming data has duplicate column name(s): {dupes}")

        # A left-join is row-count preserving only when both sides are unique on
        # the join key. Abort before merge to prevent row fan-out.
        base_dup_mask = df[id_col].duplicated(keep=False)
        if base_dup_mask.any():
            dupes = df.loc[base_dup_mask, id_col].astype(str).tolist()
            abort(f"merge_columns: datafile contains duplicate {id_col!r} value(s): "
                  f"{_format_duplicate_ids(dupes)}")

        incoming_dup_mask = new_df[id_col].duplicated(keep=False)
        if incoming_dup_mask.any():
            dupes = new_df.loc[incoming_dup_mask, id_col].astype(str).tolist()
            abort(f"merge_columns: incoming data contains duplicate {id_col!r} value(s): "
                  f"{_format_duplicate_ids(dupes)}. Merge would duplicate datafile rows; "
                  "deduplicate or resolve conflicts first.")

        new_names = [c for c in new_df.columns if c != id_col]
        reserved_hit = [c for c in new_names if c in (COL_ID, COL_SEQ, COL_QUERY)]
        if reserved_hit and not force:
            abort(f"refusing to overwrite reserved column(s): {reserved_hit}. "
                  f"Use --force to overwrite.")

        drop_set = set(drop_columns or [])
        reserved_drop = [c for c in drop_set if c in (COL_ID, COL_SEQ, COL_QUERY)]
        if reserved_drop:
            abort(f"refusing to drop reserved column(s): {reserved_drop}.")
        if drop_set:
            existing_drop = [c for c in df.columns if c in drop_set]
            if existing_drop:
                df = df.drop(columns=existing_drop)
                types = {c: t for c, t in types.items() if c not in existing_drop}

        collisions = [c for c in new_names if c in df.columns]
        if collisions and not force:
            abort(f"column(s) already exist: {collisions}. Use --force to overwrite.")
        if collisions:
            df = df.drop(columns=collisions)
            types = {c: t for c, t in types.items() if c not in collisions}

        n_before = len(df)
        try:
            merged = df.merge(new_df, on=id_col, how="left", validate="one_to_one")
        except pd.errors.MergeError as e:
            abort(f"merge_columns: merge is not one-to-one on {id_col!r}: {e}")
        if len(merged) != n_before:
            abort(f"merge_columns: internal error: merge changed row count "
                  f"from {n_before} to {len(merged)}")

        out_types = {**types, **new_types}

        tmp = Path(str(datafile_path) + ".tmp")
        write_datafile(merged, out_types, tmp)
        os.replace(tmp, datafile_path)

        if manifest_path and Path(manifest_path).exists():
            manifest = read_manifest(manifest_path)
            cov = column_coverage(merged)
            entries = {e["name"]: e for e in manifest.get("columns", [])}

            # Remove provenance entries for columns dropped as part of a full-system
            # replace, then add/replace provenance for the incoming columns.
            for c in drop_set:
                entries.pop(c, None)
            for c in new_names:
                entries[c] = make_column_entry(
                    c, new_types.get(c, TYPE_LABEL), provenance_source=provenance_source,
                    tool=tool, version=version, params=params, notes=notes,
                    coverage=cov.get(c, ""))
            for c in merged.columns:                      # refresh coverage everywhere
                if c in entries:
                    entries[c]["coverage"] = cov.get(c, entries[c].get("coverage", ""))
            manifest["columns"] = [entries[c] for c in merged.columns if c in entries]
            write_manifest(manifest_path, manifest, tool)

        return merged
