"""Tests for the creation pipeline: ``sse_initialization.build_entry`` (§6).

``build_entry`` is the source-agnostic core of entry creation — it takes a
reader's ``ReaderResult`` and applies the ID-uniqueness check, query resolution
(auto + ``--query`` override), the split drop rule for unusable sequences, the
EM dangling-``Closest query`` check, and feature attachment. These drive it end
to end through the real FASTA/EM readers.
"""
from types import SimpleNamespace

import pytest

from sse_tools import common
from sse_tools.readers.fasta import read_fasta
from sse_tools.readers.em import read_em

import sse_initialization as ssi          # from scripts/, importable via conftest


def fasta_result(tmp_path, text):
    p = tmp_path / "s.fasta"
    p.write_text(text, encoding="utf-8")
    return read_fasta(p, SimpleNamespace())


def em_result(tmp_path, text):
    p = tmp_path / "e.tsv"
    p.write_text(text, encoding="utf-8")
    return read_em(p, SimpleNamespace())


# ------------------------------------------------------------- happy paths

def test_build_entry_fasta_end_to_end(tmp_path):
    result = fasta_result(tmp_path,
                          ">P1 alpha\nMKLVACDEFGHIKLMNPQRS\n"
                          ">P2 beta\nMKIVACDEFGHIKLMNPQRS\n")
    report = {}
    df = ssi.build_entry(result, None, report)

    # reserved columns first, then source labels, then features (§6.1 step 7)
    assert list(df.columns[:4]) == \
        [common.COL_ID, common.COL_SEQ, common.COL_QUERY, "Description"]
    assert list(df[common.COL_ID]) == ["P1", "P2"]
    assert list(df[common.COL_QUERY]) == ["False", "False"]   # fasta has no query concept
    assert "length" in df.columns and "GRAVY" in df.columns    # features attached
    assert report["rows_kept"] == 2
    assert report["query_count"] == 0


def test_build_entry_em_auto_query_and_closest_query_survives(tmp_path):
    result = em_result(tmp_path,
                       "Accession\tSequence\tClosest query\n"
                       "Q1\tMKLVACDEFGHIKLMNPQRS\t\n"
                       "H1\tMKIVACDEFGHIKLMNPQRS\tQ1\n")
    df = ssi.build_entry(result, None, {})
    flags = dict(zip(df[common.COL_ID], df[common.COL_QUERY]))
    assert flags["Q1"] == "True"          # referenced as a Closest query -> query
    assert flags["H1"] == "False"
    assert "Closest query" in df.columns  # kept as a source label (§8.1)


# --------------------------------------------------------- validation aborts

def test_build_entry_aborts_on_duplicate_ids(tmp_path):
    # Two records collapse to the same id with no fallback field (§5.3 -> §6.4).
    result = fasta_result(tmp_path, ">P1 one\nMKLV\n>P1 two\nMKIV\n")
    with pytest.raises(common.SSEError):
        ssi.build_entry(result, None, {})


def test_build_entry_drops_unusable_non_query_row(tmp_path):
    result = fasta_result(tmp_path,
                          ">P1 good\nMKLVACDEFGHIKLMNPQRS\n"
                          ">P2 bad\nMKLVMKLV*\n")           # '*' -> unusable
    report = {}
    df = ssi.build_entry(result, None, report)
    assert list(df[common.COL_ID]) == ["P1"]               # bad row dropped
    assert report["dropped_count"] == 1
    assert "P2" in report["dropped_ids"]


def test_build_entry_aborts_when_a_query_row_is_unusable(tmp_path):
    # A query may never be silently dropped -- an unusable query aborts (§6.1.4).
    result = fasta_result(tmp_path,
                          ">P1 good\nMKLVACDEFGHIKLMNPQRS\n"
                          ">P2 bad\nMKLVMKLV*\n")
    with pytest.raises(common.SSEError):
        ssi.build_entry(result, ["P2 bad"], {})            # mark the bad row as query


def test_build_entry_query_override_flags_selected_row(tmp_path):
    result = fasta_result(tmp_path,
                          ">P1 alpha\nMKLVACDEFGHIKLMNPQRS\n"
                          ">P2 beta\nMKIVACDEFGHIKLMNPQRS\n")
    df = ssi.build_entry(result, ["P1 alpha"], {})         # --query by full header
    flags = dict(zip(df[common.COL_ID], df[common.COL_QUERY]))
    assert flags["P1"] == "True"
    assert flags["P2"] == "False"


def test_build_entry_aborts_on_unmatched_query(tmp_path):
    result = fasta_result(tmp_path, ">P1 alpha\nMKLVACDEFGHIKLMNPQRS\n")
    with pytest.raises(common.SSEError):
        ssi.build_entry(result, ["no such header"], {})


def test_build_entry_aborts_on_dangling_closest_query(tmp_path):
    # A surviving row references a Closest query id that is no row at all (§6.1 note).
    result = em_result(tmp_path,
                       "Accession\tSequence\tClosest query\n"
                       "H1\tMKLVACDEFGHIKLMNPQRS\tZZZ\n")
    with pytest.raises(common.SSEError):
        ssi.build_entry(result, None, {})


def test_build_entry_query_override_replaces_auto_detection(tmp_path):
    # --query REPLACES the auto-detected query set entirely; it is not a union (§8.5).
    result = em_result(tmp_path,
                       "Accession\tSequence\tClosest query\n"
                       "Q1\tMKLVACDEFGHIKLMNPQRS\t\n"
                       "H1\tMKIVACDEFGHIKLMNPQRS\tQ1\n"
                       "H2\tMKLLACDEFGHIKLMNPQRS\tQ1\n")
    # auto-detection would flag Q1 (its id appears in Closest query); override with H2.
    df = ssi.build_entry(result, ["H2"], {})
    flags = dict(zip(df[common.COL_ID], df[common.COL_QUERY]))
    assert flags["H2"] == "True"
    assert flags["Q1"] == "False"          # auto-detected query replaced, not unioned
    assert flags["H1"] == "False"


def test_build_entry_aborts_when_no_sequence_information(tmp_path):
    # Sequence column present but entirely blank -> structural abort (§6.1.3).
    result = em_result(tmp_path, "Accession\tSequence\nA1\t\nA2\t\n")
    with pytest.raises(common.SSEError):
        ssi.build_entry(result, None, {})


def test_build_entry_drops_internal_underscore_columns(tmp_path):
    # Columns whose name starts with '_' are internal and dropped on write (§5.0/§6.1).
    result = fasta_result(tmp_path, ">P1 some description\nMKLVACDEFGHIKLMNPQRS\n")
    df = ssi.build_entry(result, None, {})
    assert not any(c.startswith("_") for c in df.columns)   # _full_header gone
    assert "Description" in df.columns                       # real source label kept


def test_build_entry_integer_features_render_without_dot_zero(tmp_path):
    # length/acidic_count/basic_count are Int64 -> plain integers, no ".0" (§9.1).
    result = fasta_result(tmp_path, ">P1 desc\nMKLVACDEFGHIKLMNPQRS\n")   # length 20
    df = ssi.build_entry(result, None, {})
    types = {c: (common.TYPE_ID if c == common.COL_ID else common.TYPE_LABEL)
             for c in df.columns}
    out = tmp_path / "P1.sse.tsv"
    common.write_datafile(df, types, out)
    df2, _ = common.read_datafile(out)
    assert df2.loc[0, "length"] == "20"       # not "20.0"
