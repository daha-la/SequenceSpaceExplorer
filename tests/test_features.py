"""Tests for ``sse_tools.compute_seq_features`` — the sequence feature set (§9).

The interesting cases are the guards: divide-by-zero-shaped ratios that must
come back ``None`` instead of raising, and the pKa/ProtParam group that must be
suppressed to ``None`` on ambiguous residues rather than silently producing a
number for a model that isn't defined there.
"""
import pytest

from sse_tools import compute_seq_features as features


# ------------------------------------------------------ usability policy (§6.2)

@pytest.mark.parametrize("seq,usable", [
    ("MKLV", True),
    ("mklv", True),          # lower-cased before checking
    ("  MKLV  ", True),      # surrounding whitespace stripped
    ("MKLVX", True),         # X is a tolerated ambiguity code
    ("MKLVU", True),         # U (selenocysteine) tolerated
    ("MKLVBZ", True),        # B, Z tolerated
    ("MKLV*", False),        # stop symbol
    ("MK-LV", False),        # alignment gap
    ("MKLVO", False),        # O (pyrrolysine) not tolerated by the tokenizer
    ("MKLV1", False),        # digit
    ("", False),             # empty
    ("   ", False),          # whitespace-only
])
def test_is_usable_sequence(seq, usable):
    assert features.is_usable_sequence(seq) is usable


# ------------------------------------------------------------- robust group (§9.1)

def test_feature_counts():
    f = features.compute_seq_features("DEKKRR")   # acidic D+E=2, basic K+K+R+R=4
    assert f["length"] == 6
    assert f["acidic_count"] == 2
    assert f["basic_count"] == 4


def test_ed_rk_ratio_is_none_without_basic_residues():
    # K+R == 0 must yield None, not a ZeroDivisionError.
    f = features.compute_seq_features("DEAA")
    assert f["basic_count"] == 0
    assert f["ED_RK_ratio"] is None


def test_ed_ik_ratio_is_none_without_i_or_k():
    # I+K == 0 must yield None (its own formula, distinct from ED_RK).
    f = features.compute_seq_features("DERR")
    assert f["ED_IK_ratio"] is None


# ---------------------------------------------------- conditional group (§9.1)

def test_conditional_features_are_none_on_ambiguous_residue():
    f = features.compute_seq_features("MKLVDEX")   # X -> pKa model undefined
    # robust group is still populated ...
    assert f["length"] == 7
    assert f["acidic_count"] == 2
    # ... but the whole conditional group is suppressed.
    for k in features.CONDITIONAL_FEATURES:
        assert f[k] is None


def test_conditional_features_present_for_standard_sequence():
    f = features.compute_seq_features("MKLVDEACDEFGHIK")   # all 20-standard
    for k in features.CONDITIONAL_FEATURES:
        assert f[k] is not None


def test_empty_sequence_is_all_null_not_an_error():
    f = features.compute_seq_features("")
    assert f["length"] == 0
    assert f["acidic_ratio"] is None      # L == 0 -> None, not ZeroDivisionError
    assert f["ED_RK_ratio"] is None
    for k in features.CONDITIONAL_FEATURES:
        assert f[k] is None
