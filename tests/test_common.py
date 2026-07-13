"""Contract tests for ``sse_tools.common`` — the datafile foundation.

Every tool in SSE reads and writes ``.sse.tsv`` files through this module, so a
silent regression here corrupts every entry. These tests pin the guarantees
written up in ``docs/SSE_datafile_spec.md``:

  - the Type-row contract on read/write   (spec §3)
  - reserved columns and Type tokens      (spec §4)
  - the additive, row-preserving merge    (spec §10, §12)

They are pure CPU logic: no models, no network, no GPU. Run them with::

    python -m pytest tests/ -v
"""
import pandas as pd
import pytest

from sse_tools import common


# --------------------------------------------------------------------- helpers

def write_df(path, data: dict, types: dict):
    """Write a datafile from column->list ``data`` and column->token ``types``."""
    common.write_datafile(pd.DataFrame(data), types, path)


def toy_datafile(directory, name="toy.sse.tsv"):
    """A minimal but valid 3-row datafile: id + Sequence + one label column."""
    path = directory / name
    write_df(
        path,
        {"id": ["A", "B", "C"],
         "Sequence": ["MKL", "MKV", "MKI"],
         "group": ["x", "x", "y"]},
        {"id": "id", "Sequence": "label", "group": "label"},
    )
    return path


# --------------------------------------------------------- read/write round trip

def test_write_read_roundtrip(tmp_path):
    path = toy_datafile(tmp_path)
    df, types = common.read_datafile(path)

    assert list(df.columns) == ["id", "Sequence", "group"]
    assert types == {"id": "id", "Sequence": "label", "group": "label"}
    assert list(df["id"]) == ["A", "B", "C"]
    assert len(df) == 3          # the Type row is peeled off, not counted as data


def test_read_returns_string_values(tmp_path):
    # Numeric-looking values survive as strings; the visualizer coerces later (§3.2).
    path = tmp_path / "nums.sse.tsv"
    write_df(path, {"id": ["A", "B"], "score": ["1", "2"]},
             {"id": "id", "score": "label"})
    df, _ = common.read_datafile(path)
    assert list(df["score"]) == ["1", "2"]


# ------------------------------------------------- write-side validation (§3.1)

def test_write_rejects_unknown_token(tmp_path):
    with pytest.raises(common.SSEError):
        write_df(tmp_path / "bad.sse.tsv",
                 {"id": ["A"], "x": ["1"]},
                 {"id": "id", "x": "labl"})          # typo'd token


def test_write_rejects_missing_token(tmp_path):
    with pytest.raises(common.SSEError):
        write_df(tmp_path / "bad.sse.tsv",
                 {"id": ["A"], "x": ["1"]},
                 {"id": "id"})                        # no token for x


def test_write_rejects_no_id_column(tmp_path):
    with pytest.raises(common.SSEError):
        write_df(tmp_path / "bad.sse.tsv",
                 {"a": ["A"], "b": ["1"]},
                 {"a": "label", "b": "label"})


def test_write_rejects_multiple_id_columns(tmp_path):
    with pytest.raises(common.SSEError):
        write_df(tmp_path / "bad.sse.tsv",
                 {"id": ["A"], "id2": ["B"]},
                 {"id": "id", "id2": "id"})


def test_write_rejects_first_column_not_id(tmp_path):
    # §3.1: the first column must carry the 'id' token.
    with pytest.raises(common.SSEError):
        write_df(tmp_path / "bad.sse.tsv",
                 {"label_first": ["x"], "id": ["A"]},
                 {"label_first": "label", "id": "id"})


# -------------------------------------------------- read-side validation (§3.2)

def test_read_rejects_empty_file(tmp_path):
    path = tmp_path / "empty.sse.tsv"
    path.write_text("", encoding="utf-8")
    with pytest.raises(common.SSEError):
        common.read_datafile(path)


def test_read_rejects_missing_type_row(tmp_path):
    # A plain TSV with no Type row: the first data row's id cell isn't "id".
    path = tmp_path / "notype.sse.tsv"
    path.write_text("id\tSequence\nA\tMKL\nB\tMKV\n", encoding="utf-8")
    with pytest.raises(common.SSEError):
        common.read_datafile(path)


def test_read_rejects_duplicate_headers(tmp_path):
    # Duplicate physical column names are caught before pandas can mangle them.
    path = tmp_path / "dupe.sse.tsv"
    path.write_text("id\tgroup\tgroup\nid\tlabel\tlabel\nA\tx\ty\n", encoding="utf-8")
    with pytest.raises(common.SSEError):
        common.read_datafile(path)


# ------------------------------------------------------ merge invariants (§10)

def test_merge_preserves_row_order_and_count(tmp_path):
    path = toy_datafile(tmp_path)
    # incoming deliberately in a different row order than the datafile
    incoming = pd.DataFrame({"id": ["C", "A", "B"], "score": ["3", "1", "2"]})
    common.merge_columns(path, incoming, {"score": "label"})

    df, _ = common.read_datafile(path)
    assert list(df["id"]) == ["A", "B", "C"]        # order unchanged
    assert list(df["score"]) == ["1", "2", "3"]     # matched by id, not by position
    assert len(df) == 3                             # no rows added or dropped


def test_merge_unmatched_datafile_rows_get_empty(tmp_path):
    path = toy_datafile(tmp_path)
    incoming = pd.DataFrame({"id": ["A"], "score": ["1"]})   # covers only 1 of 3
    common.merge_columns(path, incoming, {"score": "label"})

    df, _ = common.read_datafile(path)
    assert list(df["score"]) == ["1", "", ""]       # unmatched rows -> empty cell
    assert len(df) == 3


def test_merge_ignores_extra_incoming_ids(tmp_path):
    path = toy_datafile(tmp_path)
    incoming = pd.DataFrame({"id": ["A", "Z"], "score": ["1", "9"]})  # Z not in datafile
    common.merge_columns(path, incoming, {"score": "label"})

    df, _ = common.read_datafile(path)
    assert "Z" not in list(df["id"])                # extra incoming id doesn't add a row
    assert len(df) == 3


def test_merge_aborts_on_duplicate_incoming_ids(tmp_path):
    path = toy_datafile(tmp_path)
    incoming = pd.DataFrame({"id": ["A", "A"], "score": ["1", "2"]})
    with pytest.raises(common.SSEError):
        common.merge_columns(path, incoming, {"score": "label"})


def test_merge_aborts_on_duplicate_datafile_ids(tmp_path):
    # write_datafile does not enforce id uniqueness; the merge must.
    path = tmp_path / "dupids.sse.tsv"
    write_df(path, {"id": ["A", "A"], "Sequence": ["MKL", "MKV"]},
             {"id": "id", "Sequence": "label"})
    incoming = pd.DataFrame({"id": ["A"], "score": ["1"]})
    with pytest.raises(common.SSEError):
        common.merge_columns(path, incoming, {"score": "label"})


def test_merge_aborts_on_column_collision_without_force(tmp_path):
    path = toy_datafile(tmp_path)
    incoming = pd.DataFrame({"id": ["A", "B", "C"], "group": ["p", "q", "r"]})
    with pytest.raises(common.SSEError):
        common.merge_columns(path, incoming, {"group": "label"})


def test_merge_replaces_column_with_force(tmp_path):
    path = toy_datafile(tmp_path)
    incoming = pd.DataFrame({"id": ["A", "B", "C"], "group": ["p", "q", "r"]})
    common.merge_columns(path, incoming, {"group": "label"}, force=True)

    df, _ = common.read_datafile(path)
    assert list(df["group"]) == ["p", "q", "r"]     # fully replaced
    assert len(df) == 3


def test_merge_drop_columns_removes_orphans(tmp_path):
    # The coordinate-system replace case (§11A.3): shrink PC1..PC3 down to PC1..PC2.
    path = tmp_path / "coords.sse.tsv"
    write_df(path,
             {"id": ["A", "B"], "PC1": ["1", "2"], "PC2": ["3", "4"], "PC3": ["5", "6"]},
             {"id": "id", "PC1": "coordinate", "PC2": "coordinate", "PC3": "coordinate"})
    incoming = pd.DataFrame({"id": ["A", "B"], "PC1": ["9", "9"], "PC2": ["9", "9"]})
    common.merge_columns(path, incoming,
                         {"PC1": "coordinate", "PC2": "coordinate"},
                         force=True, drop_columns=["PC3"])

    df, _ = common.read_datafile(path)
    assert "PC3" not in df.columns                  # orphan dropped in the same write
    assert list(df["PC1"]) == ["9", "9"]


def test_merge_aborts_on_reserved_column(tmp_path):
    # Merging onto a reserved name without force must abort (§4).
    path = toy_datafile(tmp_path)
    incoming = pd.DataFrame({"id": ["A", "B", "C"], "Sequence": ["Q", "Q", "Q"]})
    with pytest.raises(common.SSEError):
        common.merge_columns(path, incoming, {"Sequence": "label"})


def test_merge_reserved_column_replaced_with_force(tmp_path):
    # §4: reserved columns follow the same force rule as any other column, so a
    # forced overwrite of Sequence/query/id succeeds and fully replaces the
    # column. Regression guard for the former abort-before-force bug (the abort
    # used to fire even with force, making this path unreachable).
    path = toy_datafile(tmp_path)
    incoming = pd.DataFrame({"id": ["A", "B", "C"], "Sequence": ["Q", "Q", "Q"]})
    common.merge_columns(path, incoming, {"Sequence": "label"}, force=True)

    df, _ = common.read_datafile(path)
    assert list(df["Sequence"]) == ["Q", "Q", "Q"]


def test_merge_aborts_when_incoming_missing_id(tmp_path):
    path = toy_datafile(tmp_path)
    incoming = pd.DataFrame({"score": ["1", "2", "3"]})   # no id column at all
    with pytest.raises(common.SSEError):
        common.merge_columns(path, incoming, {"score": "label"})


# ---------------------------------------------------------------- small helpers

def test_column_coverage_counts_nonempty(tmp_path):
    path = tmp_path / "cov.sse.tsv"
    write_df(path, {"id": ["A", "B", "C"], "score": ["1", "", "3"]},
             {"id": "id", "score": "label"})
    df, _ = common.read_datafile(path)
    cov = common.column_coverage(df)
    assert cov["id"] == "3/3"
    assert cov["score"] == "2/3"        # the empty cell doesn't count as populated


def test_id_column_returns_id_typed_name(tmp_path):
    _, types = common.read_datafile(toy_datafile(tmp_path))
    assert common.id_column(types) == "id"


# ------------------------------------------------- entry-path resolution (§2)

def test_resolve_entry_path_direct_file(tmp_path):
    path = toy_datafile(tmp_path, name="thing.sse.tsv")
    assert common.resolve_entry_path(str(path)) == path


def test_resolve_entry_path_bare_stem(tmp_path):
    entries = tmp_path / "entries"
    stem_dir = entries / "myentry"
    stem_dir.mkdir(parents=True)
    path = toy_datafile(stem_dir, name="myentry.sse.tsv")
    assert common.resolve_entry_path("myentry", entries_dir=entries) == path


def test_resolve_entry_path_missing_aborts(tmp_path):
    with pytest.raises(common.SSEError):
        common.resolve_entry_path("does_not_exist", entries_dir=tmp_path)


# ------------------------------------------------------- null / quoting / encoding (§3.3)

def test_na_like_strings_survive_as_text(tmp_path):
    # Biological tokens that look like nulls must NOT be coerced to missing:
    # read uses keep_default_na=False (§3.3).
    path = tmp_path / "na.sse.tsv"
    write_df(path, {"id": ["A", "B", "C", "D"],
                    "note": ["NA", "NaN", "null", "N/A"]},
             {"id": "id", "note": "label"})
    df, _ = common.read_datafile(path)
    assert list(df["note"]) == ["NA", "NaN", "null", "N/A"]


def test_lossless_quoting_round_trip(tmp_path):
    # Embedded quotes, commas, and newlines must survive write -> read unchanged
    # (§3.3: a different quoting dialect is a contract violation).
    path = tmp_path / "q.sse.tsv"
    values = ['has "quotes"', "a,b,c", "line1\nline2"]
    write_df(path, {"id": ["A", "B", "C"], "note": values},
             {"id": "id", "note": "label"})
    df, _ = common.read_datafile(path)
    assert list(df["note"]) == values


def test_read_tolerates_utf8_bom(tmp_path):
    # A BOM (added by Excel round-trips) is tolerated on read via utf-8-sig (§3.3).
    path = toy_datafile(tmp_path)
    path.write_bytes(b"\xef\xbb\xbf" + path.read_bytes())
    df, types = common.read_datafile(path)
    assert list(df["id"]) == ["A", "B", "C"]
    assert types["id"] == "id"


def test_write_produces_no_bom(tmp_path):
    # Write side is UTF-8 without a BOM (§3.3).
    path = toy_datafile(tmp_path)
    assert not path.read_bytes().startswith(b"\xef\xbb\xbf")


def test_read_normalizes_type_tokens_case_and_whitespace(tmp_path):
    # Type tokens are matched after strip + lowercase (§3.1/§3.2).
    path = tmp_path / "mixed.sse.tsv"
    path.write_text("id\tgroup\n ID \t Label \nA\tx\n", encoding="utf-8")
    df, types = common.read_datafile(path)
    assert types == {"id": "id", "group": "label"}
    assert list(df["id"]) == ["A"]


# --------------------------------------------------- more merge invariants (§4, §10.7)

def test_merge_drop_columns_refuses_reserved_even_with_force(tmp_path):
    # The documented asymmetry (§4): force-*overwrite* of a reserved column is
    # allowed, but force-*drop* via drop_columns is refused unconditionally.
    path = toy_datafile(tmp_path)
    incoming = pd.DataFrame({"id": ["A", "B", "C"], "score": ["1", "2", "3"]})
    with pytest.raises(common.SSEError):
        common.merge_columns(path, incoming, {"score": "label"},
                             force=True, drop_columns=["Sequence"])


def test_merge_rejects_invalid_type_token(tmp_path):
    # A merged column's Type token must be valid (§10.7).
    path = toy_datafile(tmp_path)
    incoming = pd.DataFrame({"id": ["A", "B", "C"], "x": ["1", "2", "3"]})
    with pytest.raises(common.SSEError):
        common.merge_columns(path, incoming, {"x": "banana"})


def test_merge_rejects_id_type_for_merged_column(tmp_path):
    # 'id' is never valid for a merged column (would make two id columns) (§10.7).
    path = toy_datafile(tmp_path)
    incoming = pd.DataFrame({"id": ["A", "B", "C"], "x": ["1", "2", "3"]})
    with pytest.raises(common.SSEError):
        common.merge_columns(path, incoming, {"x": "id"})


# ------------------------------------------------------------------- manifest (§11)

def test_read_manifest_aborts_when_missing(tmp_path):
    with pytest.raises(common.SSEError):
        common.read_manifest(tmp_path / "nope.json")


def test_merge_appends_manifest_provenance(tmp_path):
    # A merge with a manifest_path appends the new column's provenance entry and
    # re-stamps last_tool (§11.2).
    path = toy_datafile(tmp_path)
    manifest_path = tmp_path / "m.json"
    common.write_manifest(manifest_path, {"columns": [
        common.make_column_entry("id", "id"),
        common.make_column_entry("Sequence", "label"),
        common.make_column_entry("group", "label"),
    ]}, tool="test_setup")

    incoming = pd.DataFrame({"id": ["A", "B", "C"], "score": ["1", "2", "3"]})
    common.merge_columns(path, incoming, {"score": "label"},
                         manifest_path=manifest_path, tool="sse_merge",
                         provenance_source="external")

    manifest = common.read_manifest(manifest_path)
    entries = {e["name"]: e for e in manifest["columns"]}
    assert entries["score"]["provenance_source"] == "external"
    assert entries["score"]["coverage"] == "3/3"
    assert manifest["last_tool"] == "sse_merge"      # re-stamped on write
