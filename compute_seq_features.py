"""Sequence-derived features for the SSE datafile.

Imported by the init script to compute features at creation. Kept separate so
the feature set can grow over time (and so a future standalone "recompute
features on an existing datafile" tool can reuse it) without touching creation.

Two groups (datafile spec §9):
  ROBUST     count-based, defined for every usable sequence.
  CONDITIONAL pKa-/ProtParam-based, undefined for ambiguous residues -> None.
"""

from collections import Counter

from Bio.SeqUtils.ProtParam import ProteinAnalysis

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")      # the 20 ProtParam handles
ALLOWED_AMBIGUOUS = set("UBZX")                # tolerated by the ESM-C tokeniser
USABLE_AA = STANDARD_AA | ALLOWED_AMBIGUOUS

# Integer-valued count features (rendered without a trailing ".0").
INT_FEATURES = ["length", "acidic_count", "basic_count"]
ROBUST_FEATURES = INT_FEATURES + ["acidic_ratio", "basic_ratio",
                                  "ED_RK_ratio", "ED_IK_ratio"]
CONDITIONAL_FEATURES = ["net_charge_pH7", "MW", "pI", "aromaticity",
                        "instability_index", "GRAVY"]
FEATURE_COLS = ROBUST_FEATURES + CONDITIONAL_FEATURES


def is_usable_sequence(seq: str) -> bool:
    """Usable iff non-empty and only standard + tolerated-ambiguous residues."""
    s = (seq or "").strip().upper()
    return bool(s) and set(s) <= USABLE_AA


def compute_seq_features(seq: str) -> dict:
    """All features for one sequence. The CONDITIONAL group is None when the
    sequence contains ambiguous residues (the pKa model is undefined there)."""
    s = seq.strip().upper()
    L = len(s)
    c = Counter(s)
    acidic = c["D"] + c["E"]
    basic = c["K"] + c["R"]
    ik = c["I"] + c["K"]
    feats = {
        "length": L,
        "acidic_count": acidic,
        "basic_count": basic,
        "acidic_ratio": acidic / L if L else None,
        "basic_ratio": basic / L if L else None,
        "ED_RK_ratio": acidic / basic if basic else None,
        "ED_IK_ratio": acidic / ik if ik else None,
    }
    if L > 0 and set(s) <= STANDARD_AA:
        pa = ProteinAnalysis(s)
        feats.update({
            "net_charge_pH7": pa.charge_at_pH(7.0),
            "MW": pa.molecular_weight(),
            "pI": pa.isoelectric_point(),
            "aromaticity": pa.aromaticity(),
            "instability_index": pa.instability_index(),
            "GRAVY": pa.gravy(),
        })
    else:
        for k in CONDITIONAL_FEATURES:
            feats[k] = None
    return feats
