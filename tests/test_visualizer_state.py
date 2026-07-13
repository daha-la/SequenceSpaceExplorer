"""Tests for ``sse_tools.visualizer_state`` — column classification and the
other pure-logic the visualizer derives from a datafile (§14.3, §14.2, §11A.1).

This is all pure pandas logic (no Dash, no network), and it decides the entire
filter UI: whether a column becomes a continuous slider, a boolean toggle, a
categorical dropdown, a tag-split dropdown, or is hidden. A misclassification
silently makes a column unfilterable, so the rules are worth pinning down.
"""
import pandas as pd
import pytest

from sse_tools import common
from sse_tools import visualizer_state as vs


def cls(values, col="x"):
    return vs.classify_label_column(col, pd.Series(values))


# ------------------------------------------------ classify_label_column (§14.3)

def test_classify_numeric_is_continuous():
    assert cls(["1.5", "2.0", "3.1", "4.2", "0.9"]) == "continuous"


def test_classify_boolean_words_is_boolean():
    assert cls(["true", "false", "yes", "no", "true"]) == "boolean"


def test_classify_zero_one_is_boolean_not_continuous():
    # 1/0 columns are booleans, not continuous — the numeric test defers to the
    # boolean test first.
    assert cls(["1", "0", "1", "0", "1"]) == "boolean"


def test_classify_plain_text_is_categorical():
    assert cls(["red", "blue", "green", "red", "blue"]) == "categorical"


def test_classify_constant_is_skip():
    assert cls(["x", "x", "x"]) == "skip"


def test_classify_empty_is_skip():
    assert cls(["", "", ""]) == "skip"


def test_classify_date_like_is_skip():
    assert cls(["2020-01-01", "2021-05-05", "2019-12-31", "2022-03-03"]) == "skip"


def test_classify_comma_values_is_tag_split():
    assert cls(["a,b", "b,c", "a,c", "c,d"]) == "tag_split"


def test_classify_comma_separated_numbers_is_skip_not_tag_split():
    assert cls(["1,2", "3,4", "5,6", "7,8"]) == "skip"


def test_classify_long_unique_free_text_is_skip():
    vals = [f"a rather long descriptive sentence number {i} with more detail"
            for i in range(20)]
    assert cls(vals) == "skip"


def test_classify_high_cardinality_is_skip():
    assert cls([f"v{i:04d}" for i in range(250)]) == "skip"


def test_classify_query_stays_boolean_even_when_constant(tmp_path):
    # Structural exception: `query` stays boolean even if every value is False,
    # because it's a fixed UI column (§14.3).
    assert vs.classify_label_column(common.COL_QUERY, pd.Series(["False"] * 5)) \
        == "boolean"


# --------------------------------- high-cardinality categorical rescue (§14.3)

def _hierarchy(levels):
    """Build a clean nested hierarchy: 210 genera over 3 families, one/two rows
    each, optionally with 420 unique species. Returns (df, types)."""
    rows = []
    for g in range(210):
        fam = f"fam{g // 70}"                     # 3 families
        for rep in range(2):
            sp = g * 2 + rep
            row = {"id": f"r{sp:04d}", "family": fam, "genus": f"g{g:03d}"}
            if "species" in levels:
                row["species"] = f"sp{sp:04d}"    # 420 unique
            rows.append(row)
    df = pd.DataFrame(rows)
    types = {c: (common.TYPE_ID if c == "id" else common.TYPE_LABEL)
             for c in df.columns}
    return df, types


def test_rescue_promotes_high_card_column_that_nests():
    # genus has 210 unique values (> MAX_CAT_UNIQUE) so it is skipped alone, but
    # it nests cleanly under family (categorical) and is rescued to categorical.
    df, types = _hierarchy(["family", "genus"])
    meta = vs.build_col_meta(df, types)
    assert meta["family"]["type"] == "categorical"
    assert meta["genus"]["type"] == "categorical"     # rescued via family


def test_rescue_recovers_a_multi_level_chain():
    # family -> genus -> species: both oversized columns recover.
    df, types = _hierarchy(["family", "genus", "species"])
    meta = vs.build_col_meta(df, types)
    assert meta["genus"]["type"] == "categorical"
    assert meta["species"]["type"] == "categorical"


def test_high_card_column_without_partner_stays_skip():
    df = pd.DataFrame({"id": [f"r{i}" for i in range(250)],
                       "code": [f"c{i:04d}" for i in range(250)]})
    types = {"id": "id", "code": "label"}
    meta = vs.build_col_meta(df, types)
    assert meta["code"]["type"] == "skip"             # no nesting partner


def test_tag_split_meta_lists_individual_tags():
    df = pd.DataFrame({"id": ["A", "B", "C"],
                       "dbs": ["afdb,pdb", "pdb,cath", "afdb,cath"]})
    types = {"id": "id", "dbs": "label"}
    meta = vs.build_col_meta(df, types)
    assert meta["dbs"]["type"] == "tag_split"
    assert meta["dbs"]["tags"] == ["afdb", "cath", "pdb"]   # split, deduped, sorted


# ---------------------------------------- coordinate-system grouping (§11A.1)

def test_coordinate_system_key_strips_trailing_number():
    assert vs.coordinate_system_key("esmc600m_mean_PC10") == "esmc600m_mean_PC"


def test_group_coordinate_systems_sorts_axes_numerically():
    systems = vs.group_coordinate_systems(
        ["esmc_PC1", "esmc_PC10", "esmc_PC2", "saprot_UMAP1", "saprot_UMAP2"])
    assert systems["esmc_PC"] == ["esmc_PC1", "esmc_PC2", "esmc_PC10"]   # not lexical
    assert systems["saprot_UMAP"] == ["saprot_UMAP1", "saprot_UMAP2"]


def test_default_axes_picks_first_system_with_two_axes():
    systems = {"esmc_PC": ["esmc_PC1", "esmc_PC2", "esmc_PC10"]}
    assert vs.default_axes(systems) == ("esmc_PC1", "esmc_PC2")


def test_default_axes_none_when_no_coordinates():
    assert vs.default_axes({}) == (None, None)


# ------------------------------------------------ boolean / query helpers (§8)

def test_boolean_mask_maps_all_truthy_and_falsy_forms():
    m = vs.boolean_mask(pd.Series(["True", "1", "yes", "0.0", "no", "false"]))
    assert list(m) == [True, True, True, False, False, False]   # "0.0" -> 0 -> False


def test_boolean_mask_unknown_value_is_na():
    assert vs.boolean_mask(pd.Series(["maybe"])).isna().all()


def test_query_ids_from_df_returns_truthy_rows():
    df = pd.DataFrame({"id": ["A", "B", "C"], "query": ["True", "False", "1"]})
    assert vs.query_ids_from_df(df) == ["A", "C"]


def test_detect_name_cols_matches_keyword_columns():
    df = pd.DataFrame({"id": ["A", "B"], "gene_name": ["x", "y"], "misc": ["p", "q"]})
    assert vs.detect_name_cols(df, "id", ["gene_name", "misc"]) == ["gene_name"]


# --------------------------------------- NaN-aware axis ranges (§14.2)

def test_axis_range_ignores_rows_missing_a_coordinate():
    # The x=99 row has no y value, so it must not stretch the x range.
    df = pd.DataFrame({"x": ["1", "99", "3"], "y": ["1", "", "3"]})
    r = vs.axis_range(df, "x", "y")
    assert r["xrange"][1] < 10


def test_axis_range_none_when_no_complete_rows():
    df = pd.DataFrame({"x": ["", "", ""], "y": ["1", "2", "3"]})
    assert vs.axis_range(df, "x", "y") is None
