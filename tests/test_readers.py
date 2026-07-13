"""Tests for the source readers (``sse_tools/readers/``, spec §5).

Readers turn a raw source file into a flat ``ReaderResult`` table. Pure parsing
logic, no models or network. The FASTA two-pass ID de-duplication (§5.3) is the
subtlest code in the project, so it gets the most coverage here.
"""
import json
from types import SimpleNamespace

import pytest

from sse_tools import common
from sse_tools.readers.fasta import read_fasta
from sse_tools.readers.em import read_em
from sse_tools.readers.foldseek import read_foldseek


def _write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


# ============================================================= FASTA (§5.3)

def test_fasta_known_db_tag_takes_accession(tmp_path):
    p = _write(tmp_path / "s.fasta",
               ">sp|P12345|OLED_STRAN Oleandomycin glucosyltransferase\nMKLVACDEF\n")
    r = read_fasta(p, SimpleNamespace())
    row = r.table.iloc[0]
    assert row[common.COL_ID] == "P12345"                 # field[1] of a db-tagged token
    assert row[common.COL_SEQ] == "MKLVACDEF"
    assert row["Description"] == "Oleandomycin glucosyltransferase"
    assert row["_full_header"] == \
        "sp|P12345|OLED_STRAN Oleandomycin glucosyltransferase"
    assert r.source == "fasta"


def test_fasta_plain_token_is_the_id(tmp_path):
    p = _write(tmp_path / "s.fasta", ">MyProtein some description\nMKLV\n")
    r = read_fasta(p, SimpleNamespace())
    assert r.table.iloc[0][common.COL_ID] == "MyProtein"
    assert r.table.iloc[0]["Description"] == "some description"


def test_fasta_collision_resolved_with_third_pipe_field(tmp_path):
    # Two records with the same accession P1 disambiguate on field[2] (db-tag branch).
    p = _write(tmp_path / "s.fasta",
               ">sp|P1|NAME_A first\nMKLV\n"
               ">sp|P1|NAME_B second\nMKIV\n")
    ids = list(read_fasta(p, SimpleNamespace()).table[common.COL_ID])
    assert ids == ["P1_NAME_A", "P1_NAME_B"]
    assert len(set(ids)) == 2                             # collision resolved


def test_fasta_collision_resolved_with_second_pipe_field(tmp_path):
    # Plain (non-db-tag) branch disambiguates on field[1].
    p = _write(tmp_path / "s.fasta",
               ">clone|A extra\nMKLV\n"
               ">clone|B extra\nMKIV\n")
    ids = list(read_fasta(p, SimpleNamespace()).table[common.COL_ID])
    assert ids == ["clone_A", "clone_B"]


def test_fasta_unresolvable_collision_left_to_pipeline(tmp_path):
    # No pipe-field to fall back on: the reader does NOT abort, it returns
    # duplicates and lets creation's uniqueness check catch them (§5.3).
    p = _write(tmp_path / "s.fasta", ">P1 one\nMKLV\n>P1 two\nMKIV\n")
    ids = list(read_fasta(p, SimpleNamespace()).table[common.COL_ID])
    assert ids == ["P1", "P1"]


def test_fasta_no_records_aborts(tmp_path):
    p = _write(tmp_path / "empty.fasta", "\n\n")          # no '>' headers at all
    with pytest.raises(common.SSEError):
        read_fasta(p, SimpleNamespace())


def test_fasta_has_no_auto_query_and_matches_full_header(tmp_path):
    p = _write(tmp_path / "s.fasta", ">P1 desc\nMKLV\n")
    r = read_fasta(p, SimpleNamespace())
    assert r.auto_query is None                           # FASTA has no query concept
    assert r.match_col == "_full_header"                  # --query matches the full header
    assert r.match_label == "full header"


# ================================================================= EM (§5.1)

def test_em_reads_default_columns(tmp_path):
    p = _write(tmp_path / "e.tsv", "Accession\tSequence\nA1\tMKLV\nA2\tMKIV\n")
    r = read_em(p, SimpleNamespace())
    assert (r.id_col, r.seq_col) == ("Accession", "Sequence")
    assert list(r.table["Accession"]) == ["A1", "A2"]
    assert r.auto_query is None


def test_em_honours_custom_column_names(tmp_path):
    p = _write(tmp_path / "e.tsv", "ProteinID\tSeq\nP1\tMKLV\n")
    r = read_em(p, SimpleNamespace(id_col="ProteinID", seq_col="Seq"))
    assert r.id_col == "ProteinID"
    assert list(r.table["ProteinID"]) == ["P1"]


def test_em_missing_id_column_aborts(tmp_path):
    p = _write(tmp_path / "e.tsv", "Name\tSequence\nA1\tMKLV\n")
    with pytest.raises(common.SSEError):
        read_em(p, SimpleNamespace())


def test_em_missing_sequence_column_aborts(tmp_path):
    p = _write(tmp_path / "e.tsv", "Accession\tName\nA1\tfoo\n")
    with pytest.raises(common.SSEError):
        read_em(p, SimpleNamespace())


def test_em_auto_query_from_closest_query(tmp_path):
    # A row is a query iff its id appears as a *value* in 'Closest query' (§8.1).
    p = _write(tmp_path / "e.tsv",
               "Accession\tSequence\tClosest query\n"
               "Q1\tMKLV\t\n"
               "H1\tMKIV\tQ1\n"
               "H2\tMKLL\tQ1\n")
    r = read_em(p, SimpleNamespace())
    mask = dict(zip(r.table["Accession"], r.auto_query))
    assert bool(mask["Q1"]) is True
    assert bool(mask["H1"]) is False
    assert bool(mask["H2"]) is False


def test_em_without_closest_query_has_no_auto_query(tmp_path):
    p = _write(tmp_path / "e.tsv", "Accession\tSequence\nA1\tMKLV\n")
    assert read_em(p, SimpleNamespace()).auto_query is None


# =========================================================== Foldseek (§5.2)

def _foldseek_json(path, results, header="job_A ref", seq="MKLVMKLV"):
    data = [{"queries": [{"header": header, "sequence": seq}],
             "results": [{"db": db, "alignments": {"0": hits}}
                         for db, hits in results]}]
    return _write(path, json.dumps(data))


def test_foldseek_dedups_targets_and_aggregates_databases(tmp_path):
    p = _foldseek_json(tmp_path / "fs.json", [
        ("afdb",   [{"target": "T1", "tSeq": "AAAA", "seqId": 50}]),
        ("pdb100", [{"target": "T1", "tSeq": "AAAA", "seqId": 50},
                    {"target": "T2", "tSeq": "CCCC", "seqId": 40}]),
    ])
    df = read_foldseek(p, SimpleNamespace()).table
    assert list(df[common.COL_ID]).count("T1") == 1       # one row per unique target
    t1 = df[df[common.COL_ID] == "T1"].iloc[0]
    assert t1["Databases"] == "afdb, pdb100"              # aggregated + sorted tag-split


def test_foldseek_adds_synthetic_query_when_no_self_hit(tmp_path):
    p = _foldseek_json(tmp_path / "fs.json", [
        ("afdb", [{"target": "T1", "tSeq": "AAAA", "seqId": 50}]),
    ])
    r = read_foldseek(p, SimpleNamespace())
    qmask = dict(zip(r.table[common.COL_ID], r.auto_query))
    assert "job_A" in qmask                                # synthetic query row from queries[0]
    assert bool(qmask["job_A"]) is True
    assert bool(qmask["T1"]) is False


def test_foldseek_self_hit_becomes_the_query(tmp_path):
    # A seqId==100 hit whose tSeq equals the query sequence IS the query: it is
    # flagged, and no separate synthetic row is added (§8.2).
    p = _foldseek_json(tmp_path / "fs.json", [
        ("afdb", [{"target": "T1", "tSeq": "MKLVMKLV", "seqId": 100},
                  {"target": "T2", "tSeq": "CCCC", "seqId": 40}]),
    ], seq="MKLVMKLV")
    r = read_foldseek(p, SimpleNamespace())
    ids = list(r.table[common.COL_ID])
    assert "job_A" not in ids                              # no synthetic duplicate row
    qmask = dict(zip(ids, r.auto_query))
    assert bool(qmask["T1"]) is True
    assert bool(qmask["T2"]) is False


def test_foldseek_bad_json_shape_aborts(tmp_path):
    p = _write(tmp_path / "fs.json", json.dumps({"not": "a list"}))
    with pytest.raises(common.SSEError):
        read_foldseek(p, SimpleNamespace())


def test_foldseek_keeps_numeric_suffix_targets_distinct(tmp_path):
    # `_N`-suffixed targets are distinct proteins and must never be collapsed by
    # suffix-normalization (§5.2).
    p = _foldseek_json(tmp_path / "fs.json", [
        ("afdb", [{"target": "A0A1", "tSeq": "AAAA", "seqId": 50},
                  {"target": "A0A1_2", "tSeq": "CCCC", "seqId": 50},
                  {"target": "A0A1_3", "tSeq": "GGGG", "seqId": 50}]),
    ])
    ids = list(read_foldseek(p, SimpleNamespace()).table[common.COL_ID])
    for t in ("A0A1", "A0A1_2", "A0A1_3"):
        assert t in ids                                   # kept verbatim, not merged
