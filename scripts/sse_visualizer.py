"""SSE visualizer: Dash viewer for one entry-centered SSE datafile.

Run:
    python scripts/sse_visualizer.py <entry-stem|entry-dir|datafile.sse.tsv>

The app reads exactly one .sse.tsv datafile, using the Type row as the contract:
id / label / coordinate. Selected points can be exported to a selection cache on
disk; the Boltz-2 structure-prediction and RMSD workflow that consumes those
selections runs in the pipeline (scripts/sse_boltz.py), not here.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import threading
import re
import sys
import uuid
from pathlib import Path
from typing import Optional

# Allow running as `python scripts/sse_visualizer.py` from a source checkout.
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output, State, ctx, no_update, ALL

from sse_tools.common import COL_ID, COL_QUERY, COL_SEQ, SSEError, TYPE_COORDINATE, TYPE_ID, TYPE_LABEL
from sse_tools.visualizer_state import (
    BOOL_STRINGS,
    EntryContext,
    MAX_CAT_UNIQUE,
    VisualizerState,
    axis_range,
    boolean_mask,
    coordinate_system_key,
    group_coordinate_systems,
    load_visualizer_state,
    resolve_entry,
)
from sse_tools import layers as layer_store
from sse_tools import selections as selection_cache

PORT = 8051

COLORMAP_OPTIONS = [
    "Viridis", "Plasma", "Inferno", "Cividis", "Turbo", "Magma",
    "Blues", "Reds", "Greens", "Oranges", "Purples",
    "RdBu", "RdYlBu", "Spectral", "Coolwarm",
]
FIXED_COLOR_OPTIONS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "crimson", "forestgreen", "gold", "navy", "darkorange",
]
MARKER_SYMBOL_OPTIONS = [
    {"label": "● Circle", "value": "circle"},
    {"label": "■ Square", "value": "square"},
    {"label": "◆ Diamond", "value": "diamond"},
    {"label": "▲ Triangle up", "value": "triangle-up"},
    {"label": "+ Cross", "value": "cross"},
    {"label": "✕ X", "value": "x"},
    {"label": "★ Star", "value": "star"},
]
SELECTION_COLOR_OPTIONS = ["#e74c3c", "#f39c12", "#2ecc71", "#9b59b6", "#1abc9c", "#ffffff"]
EXPORT_DPI_OPTIONS = [
    {"label": "150 dpi", "value": 150},
    {"label": "300 dpi", "value": 300},
    {"label": "600 dpi", "value": 600},
]
EXPORT_FORMAT_OPTIONS = [{"label": "PNG", "value": "png"}, {"label": "SVG", "value": "svg"}, {"label": "PDF", "value": "pdf"}]
EXPORT_DEFAULT_WIDTH = 1200
EXPORT_DEFAULT_HEIGHT = 800

DEFAULT_ALPHA = 0.7
DEFAULT_POINT_SIZE = 6
DEFAULT_BG_SIZE = 4
DEFAULT_MARKER_SIZE = 14
DEFAULT_MARKER_ALPHA = 0.9
DEFAULT_COLOR_MODE = "fixed"
DEFAULT_FIXED_COLOR = "#1f77b4"
DEFAULT_COLORMAP = "Viridis"
DEFAULT_COLOR_RANGE = "subset"
DEFAULT_WF_POSITION = "top"
DEFAULT_MARKER_MODE = "top"
DEFAULT_SYMBOL = "circle"
DEFAULT_SELECTION_COLOR = "#e74c3c"

LABEL_STYLE = {"fontSize": "13px", "marginLeft": "6px", "cursor": "pointer"}
SECTION_STYLE = {"fontWeight": "600", "fontSize": "12px", "margin": "10px 0 6px 0", "color": "var(--text-muted)"}
CONTROL_WRAPPER = {"marginLeft": "22px", "marginTop": "6px", "marginBottom": "4px"}

ENTRY: EntryContext
_STATE: VisualizerState
_ann_df = pd.DataFrame()
_types: dict[str, str] = {}
_col_meta: dict[str, dict] = {}
_id_col = COL_ID
_x_col: Optional[str] = None
_y_col: Optional[str] = None
_query_ids: list[str] = []
_name_cols: list[str] = []

_STATE_LOCK = threading.RLock()

# Numeric/string coercion caches. read_datafile loads every column as a string
# (dtype=str), and each figure render re-derives the same numeric coordinate and
# colour columns — and the string id column — many times per call. The dataframe
# is immutable between reloads, so these coercions are stable; cache them keyed by
# column name and clear the caches in reload_state whenever the datafile is swapped.
_num_cache: dict[str, pd.Series] = {}
_id_str_cache: Optional[pd.Series] = None


def numeric_col(col: str) -> pd.Series:
    """Float coercion of a datafile column, computed once per reload and cached."""
    cached = _num_cache.get(col)
    if cached is None:
        s = _ann_df[col] if col in _ann_df.columns else pd.Series(pd.NA, index=_ann_df.index)
        cached = pd.to_numeric(s, errors="coerce")
        _num_cache[col] = cached
    return cached


def id_str() -> pd.Series:
    """The id column as strings, computed once per reload and cached."""
    global _id_str_cache
    if _id_str_cache is None:
        _id_str_cache = _ann_df[_id_col].astype(str)
    return _id_str_cache


def reload_state() -> tuple[str, list]:
    """Reload the datafile into module globals. Returns (status_message, layers).

    The visualizer stores a few derived values as module globals because Dash
    callbacks need fast access to them. Reload builds a complete new state first
    and then swaps all globals while holding one lock, so render callbacks cannot
    observe a half-updated mix of old and new dataframe/axis metadata.

    Layer-ID validation against the freshly loaded datafile happens here too
    (spec §13), not just on the explicit Reload button: every path that calls
    reload_state — startup, manual reload, Boltz completion, RMSD completion —
    keeps logs/layers.json in sync with whatever the datafile now contains.
    layers.json on disk is treated as canonical, per §14.4.
    """
    global _STATE, _ann_df, _types, _col_meta, _id_col, _x_col, _y_col, _query_ids, _name_cols, _id_str_cache
    new_state = load_visualizer_state(ENTRY)
    with _STATE_LOCK:
        _STATE = new_state
        _ann_df = new_state.df
        _types = new_state.types
        _col_meta = new_state.col_meta
        _id_col = new_state.id_col
        _x_col = new_state.x_col
        _y_col = new_state.y_col
        _query_ids = new_state.query_ids
        _name_cols = new_state.name_cols
        _num_cache.clear()
        _REGION_CACHE.clear()
        _id_str_cache = None
        valid_ids = new_state.df[new_state.id_col].astype(str).tolist()

    current_layers = layer_store.read_layers(ENTRY.layers_path)
    cleaned_layers, layer_msg = layer_store.validate_layers(current_layers, valid_ids)
    if cleaned_layers != current_layers:
        layer_store.write_layers(ENTRY.layers_path, cleaned_layers)

    status = " · ".join(x for x in [new_state.warning, layer_msg] if x)
    return status, cleaned_layers


def effective_type(col: str) -> str:
    m = _col_meta.get(col, {})
    return m.get("override") or m.get("type", "skip")


def cols_of_type(*types: str) -> list[str]:
    return [c for c in _col_meta if effective_type(c) in types]


def label_filter_cols() -> list[str]:
    return [c for c in _ann_df.columns if _types.get(c) == TYPE_LABEL]


def slider_config(s: pd.Series):
    numeric = pd.to_numeric(s, errors="coerce").dropna()
    if numeric.empty:
        return 0.0, 1.0, 0.01, {0.0: "0", 1.0: "1"}
    vmin, vmax = float(numeric.min()), float(numeric.max())
    step = (vmax - vmin) / 200.0 if not np.isclose(vmin, vmax) else 1.0
    ticks = np.linspace(vmin, vmax, 5)
    marks = {float(v): f"{v:.2f}".rstrip("0").rstrip(".") for v in ticks}
    return vmin, vmax, step, marks


def apply_filters(pool_df, cont_conditions, bool_conditions, cat_conditions, tag_conditions, id_search_ids=None):
    mask = pd.Series(True, index=pool_df.index)
    for col, lo, hi in cont_conditions:
        if col in pool_df.columns:
            num = numeric_col(col).loc[pool_df.index]
            mask &= num.between(lo, hi, inclusive="both")
    for col, val in bool_conditions:
        if col in pool_df.columns:
            mask &= boolean_mask(pool_df[col]) == val
    for col, vals in cat_conditions:
        if vals and col in pool_df.columns:
            mask &= pool_df[col].astype(str).isin([str(v) for v in vals])
    for col, selected_tags in tag_conditions:
        if selected_tags and col in pool_df.columns:
            selected = set(selected_tags)
            mask &= pool_df[col].apply(lambda cell: bool({t.strip() for t in str(cell).split(",") if t.strip()} & selected) if str(cell) else False)
    if id_search_ids:
        id_match = id_str().loc[pool_df.index].isin([str(x) for x in id_search_ids])
        name_match = pd.Series(False, index=pool_df.index)
        for nc in _name_cols:
            if nc in pool_df.columns:
                name_match |= pool_df[nc].astype(str).isin([str(x) for x in id_search_ids])
        mask &= id_match | name_match
    return mask




def _parse_numeric_value(value, fallback=None):
    """Parse Dash numeric values robustly, including comma-decimal strings.

    Some browser locales display/submit decimal commas in number inputs. The
    sliders always use numeric floats, but the text inputs may arrive as strings
    such as "0,05963". Keep the filter path tolerant so typing in the boxes and
    moving the slider are equivalent.
    """
    if value is None or value == "":
        return fallback
    try:
        return float(str(value).strip().replace(",", "."))
    except (TypeError, ValueError):
        return fallback


def parse_current_conditions(enabled_values, enabled_ids, cont_values, cont_ids, bool_values, bool_ids, cat_values, cat_ids, tag_values, tag_ids, id_enabled=None, id_raw=None):
    enabled_cols = {id_d["col"] for id_d, val in zip(enabled_ids, enabled_values) if val and "on" in val}
    cont_conds = []
    for id_d, v in zip(cont_ids, cont_values):
        if v is None or id_d["col"] not in enabled_cols:
            continue
        lo = _parse_numeric_value(v[0] if len(v) > 0 else None)
        hi = _parse_numeric_value(v[1] if len(v) > 1 else None)
        if lo is None or hi is None:
            continue
        if lo > hi:
            lo, hi = hi, lo
        cont_conds.append((id_d["col"], lo, hi))
    bool_map = {"all": None, "true": True, "false": False}
    bool_conds = [(id_d["col"], bool_map[v]) for id_d, v in zip(bool_ids, bool_values) if id_d["col"] in enabled_cols and v != "all" and bool_map.get(v) is not None]
    cat_conds = [(id_d["col"], vals) for id_d, vals in zip(cat_ids, cat_values) if vals and id_d["col"] in enabled_cols]
    tag_conds = [(id_d["col"], tags) for id_d, tags in zip(tag_ids, tag_values) if tags and id_d["col"] in enabled_cols]
    id_search_active = bool(id_enabled and "on" in id_enabled)
    id_search_ids = [x.strip() for x in (id_raw or "").split(",") if x.strip()] if id_search_active and id_raw else None
    return cont_conds, bool_conds, cat_conds, tag_conds, id_search_ids


def auto_layer_name(cont_conds, bool_conds, cat_conds, tag_conds, color_mode, fixed_color, cont_col, colormap, reversed_, id_search_ids=None):
    parts = []
    if id_search_ids:
        parts.append(f"ID ({len(id_search_ids)} sequence{'s' if len(id_search_ids) != 1 else ''})")
    for col, lo, hi in cont_conds:
        parts.append(f"{col} [{lo:.2f}-{hi:.2f}]")
    for col, val in bool_conds:
        parts.append(f"{col}={'True' if val else 'False'}")
    for col, vals in cat_conds:
        parts.append(f"{col}={','.join(str(v) for v in vals[:3])}" + (f"+{len(vals)-3}" if len(vals) > 3 else ""))
    for col, tags in tag_conds:
        parts.append(f"{col}={','.join(str(t) for t in tags[:3])}" + (f"+{len(tags)-3}" if len(tags) > 3 else ""))
    parts.append(fixed_color if color_mode == "fixed" else f"{cont_col} · {colormap}{' rev' if reversed_ else ''}")
    return " · ".join(parts) if parts else "Layer"


def make_layer(cont_conds, bool_conds, cat_conds, tag_conds, color_mode, fixed_color, cont_col, colormap, reversed_, color_range_mode, name, filter_state, alpha, point_size, marker_symbol, id_search_ids=None):
    mask = apply_filters(_ann_df, cont_conds, bool_conds, cat_conds, tag_conds, id_search_ids)
    sub = _ann_df[mask]
    ids = sub[_id_col].astype(str).tolist()
    return {
        "id": str(uuid.uuid4())[:8],
        "name": name,
        "visible": True,
        "created_utc": layer_store.utc_now(),
        "ids": ids,
        "n_points": len(ids),
        "style": {
            "color_mode": color_mode,
            "fixed_color": fixed_color,
            "cont_col": cont_col,
            "colormap": colormap,
            "reversed": reversed_,
            "color_range": color_range_mode,
            "alpha": alpha,
            "point_size": point_size,
            "marker_symbol": marker_symbol or DEFAULT_SYMBOL,
        },
        "filter_state": filter_state,
    }


def _merge_layer_lists(primary, fallback):
    """Return a de-duplicated layer list, preserving primary order first.

    Dash Store state can occasionally be stale immediately after dynamic panel
    rebuilds. For save-layer operations, combine the browser Store with the
    persisted layers.json so a new save cannot accidentally replace the whole
    saved-layer list with only the new layer.
    """
    out = []
    seen = set()
    for source in (primary or []), (fallback or []):
        for layer in source or []:
            lid = layer.get("id") if isinstance(layer, dict) else None
            if not lid or lid in seen:
                continue
            out.append(layer)
            seen.add(lid)
    return out



def _pattern_click_fired(triggered, clicks, ids):
    """Return True only for a real user click on a pattern-matched button.

    Dash may fire ALL-pattern callbacks when dynamic components are inserted or
    removed. In those cases n_clicks is usually None/0, and writing layers.json
    would persist stale browser state.
    """
    if not isinstance(triggered, dict):
        return False
    lid = triggered.get("lid")
    typ = triggered.get("type")
    for click, id_d in zip(clicks or [], ids or []):
        if not isinstance(id_d, dict):
            continue
        if id_d.get("lid") == lid and id_d.get("type") == typ:
            return bool(click and click > 0)
    return False

def resolve_colormap(colormap, reversed_):
    return f"{colormap}_r" if reversed_ else colormap


def _add_cont_trace(fig, cont_axis_idx, name, ids, x, y, vals, cmin, cmax, colormap, reversed_, point_size, alpha, symbol=DEFAULT_SYMBOL):
    axis_key = f"coloraxis{cont_axis_idx}"
    fig.update_layout(**{axis_key: dict(
        colorscale=resolve_colormap(colormap, reversed_), cmin=cmin, cmax=cmax,
        colorbar=dict(title=dict(text=name, side="right"), len=0.85, y=0.5, x=1.02 + (cont_axis_idx - 1) * 0.08, thickness=15),
    )})
    fig.add_trace(go.Scattergl(
        x=x, y=y, mode="markers",
        marker=dict(size=point_size, color=vals, coloraxis=axis_key, opacity=alpha, symbol=symbol, line=dict(width=0.5, color="black")),
        name=name, customdata=ids, hovertemplate="<b>%{customdata}</b><extra></extra>",
    ))


def _covered(df, x_col, y_col):
    if not x_col or not y_col or x_col not in df.columns or y_col not in df.columns:
        return pd.Series(False, index=df.index)
    x = numeric_col(x_col).loc[df.index]
    y = numeric_col(y_col).loc[df.index]
    return x.notna() & y.notna()


def _add_query_traces(fig, df, id_col, x_col, y_col, marker_size, marker_alpha):
    if COL_QUERY not in df.columns:
        return
    qmask = boolean_mask(df[COL_QUERY]).fillna(False) & _covered(df, x_col, y_col)
    qdf = df[qmask]
    symbols = ["star", "diamond", "cross", "triangle-up", "pentagon", "hexagram"]
    colors = ["#e74c3c", "#f39c12", "#8e44ad", "#16a085", "#2980b9", "#c0392b"]
    for i, (_, row) in enumerate(qdf.iterrows()):
        qid = str(row[id_col])
        fig.add_trace(go.Scattergl(
            x=[row[x_col]], y=[row[y_col]], mode="markers",
            marker=dict(size=marker_size, symbol=symbols[i % len(symbols)], color=colors[i % len(colors)], opacity=marker_alpha, line=dict(width=1.5, color="black")),
            name=f"Query: {qid}", customdata=[qid], hovertemplate="<b>%{customdata}</b><br>Query<extra></extra>",
        ))


def _add_working_filter(fig, plot_df, df, color_mode, fixed_color, cont_col, colormap, reversed_, color_range_mode, alpha, point_size, cont_axis_idx, id_col, x_col, y_col, symbol):
    plot_df = plot_df[_covered(plot_df, x_col, y_col)]
    if plot_df.empty:
        return cont_axis_idx
    if color_mode == "fixed":
        fig.add_trace(go.Scattergl(
            x=numeric_col(x_col).loc[plot_df.index], y=numeric_col(y_col).loc[plot_df.index], mode="markers",
            marker=dict(size=point_size, color=fixed_color, symbol=symbol, opacity=alpha, line=dict(width=0.5, color="black")),
            name="Working filter", customdata=id_str().loc[plot_df.index].tolist(), hovertemplate="<b>%{customdata}</b><extra></extra>",
        ))
    elif color_mode == "continuous" and cont_col and cont_col in plot_df.columns:
        vals = numeric_col(cont_col).loc[plot_df.index]
        has_val = vals.notna()
        no_val = plot_df[~has_val]
        color_df = plot_df[has_val]
        if not no_val.empty:
            fig.add_trace(go.Scattergl(
                x=numeric_col(x_col).loc[no_val.index], y=numeric_col(y_col).loc[no_val.index], mode="markers",
                marker=dict(size=point_size, color="lightgrey", symbol=symbol, opacity=alpha, line=dict(width=0.5, color="#aaa")),
                name="Working filter — no value", customdata=id_str().loc[no_val.index].tolist(), hovertemplate="<b>%{customdata}</b><br>No value<extra></extra>",
            ))
        if color_df.empty:
            return cont_axis_idx
        v = vals[has_val]
        if color_range_mode == "global":
            full = numeric_col(cont_col).dropna()
            cmin, cmax = (float(full.min()), float(full.max())) if not full.empty else (0.0, 1.0)
        else:
            cmin, cmax = (float(v.min()), float(v.max())) if not v.empty else (0.0, 1.0)
        cont_axis_idx += 1
        _add_cont_trace(fig, cont_axis_idx, f"Working filter — {cont_col}", id_str().loc[color_df.index].tolist(), numeric_col(x_col).loc[color_df.index].tolist(), numeric_col(y_col).loc[color_df.index].tolist(), v.tolist(), cmin, cmax, colormap, reversed_, point_size, alpha, symbol)
    return cont_axis_idx


# Qualitative palette for cluster-region fills. Plotly's Dark24: 24 colours all
# at similar saturation, so no cluster reads far brighter/paler than another
# (D3 category20, used previously, pairs each hue with a washed-out tint).
REGION_PALETTE = ["#2E91E5", "#E15F99", "#1CA71C", "#FB0D0D", "#DA16FF",
                  "#B68100", "#750D86", "#EB663B", "#511CFB", "#00A08B",
                  "#FB00D1", "#FC0080", "#B2828D", "#6C7C32", "#778AAE",
                  "#862A16", "#A777F1", "#620042", "#1616A7", "#DA60CA",
                  "#6C4516", "#0D2A63", "#AF0038", "#222A2A"]

# Region rings are expensive (KDE ~0.5-1s); cache them so only a change to the
# column / shape / extent / axes recomputes, not every unrelated slider tick.
# reload_state() clears this when the datafile is reloaded.
_REGION_CACHE = {}


def _hex_to_rgba(hex_color, alpha):
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


# Category labels that are not real clusters (blank cells, HDBSCAN noise).
_NOISE_LABELS = ("", "noise", "nan", "NaN", "<NA>")
NOISE_COLOR = "#b0b0b0"


def _cluster_color_map(region_col):
    """Stable {category: colour} for a cluster column (noise excluded).

    Built from every category in the column (not just the visible/filtered ones)
    so a cluster keeps the same colour across the KDE, hull, and colour-points
    modes and regardless of filtering.
    """
    vals = _ann_df[region_col].astype(str)
    cats = sorted(c for c in vals.unique() if c not in _NOISE_LABELS)
    return {c: REGION_PALETTE[i % len(REGION_PALETTE)] for i, c in enumerate(cats)}


def _region_rings_for(region_col, shape, coverage, tightness, x_col, y_col):
    """(cat, color, rings) per cluster in the current axes, memoized."""
    key = (region_col, shape, round(coverage, 3), round(tightness, 3), x_col, y_col)
    if key in _REGION_CACHE:
        return _REGION_CACHE[key]
    from sse_tools.cluster_regions import cluster_region_rings

    df = _ann_df
    covered = _covered(df, x_col, y_col)
    vals = df[region_col].astype(str)
    cmap = _cluster_color_map(region_col)
    cats = [c for c in sorted(vals[covered].unique()) if c not in _NOISE_LABELS]
    xs_all, ys_all = numeric_col(x_col), numeric_col(y_col)
    out = []
    for cat in cats:
        idx = df.index[covered & (vals == cat)]
        pts = np.column_stack([xs_all.loc[idx].to_numpy(float),
                               ys_all.loc[idx].to_numpy(float)])
        pts = pts[~np.isnan(pts).any(axis=1)]
        if len(pts) < 3:
            continue
        rings = cluster_region_rings(pts, shape, tightness=tightness, coverage=coverage)
        if rings:
            out.append((cat, cmap.get(cat, REGION_PALETTE[0]), rings))
    _REGION_CACHE[key] = out
    return out


def _add_cluster_regions(fig, region_col, shape, opacity, coverage, tightness, x_col, y_col):
    if not region_col or region_col not in _ann_df.columns:
        return
    for cat, color, rings in _region_rings_for(region_col, shape, coverage,
                                               tightness, x_col, y_col):
        fillrgba = _hex_to_rgba(color, opacity)
        # Scattergl (not Scatter) so regions share the WebGL layer with the
        # points; Plotly draws all WebGL above all SVG, so an SVG fill could
        # never sit on top of the gl markers. One trace per ring (grouped so a
        # cluster is a single legend entry) avoids None-separated fill artefacts.
        for j, r in enumerate(rings):
            fig.add_trace(go.Scattergl(
                x=r[:, 0].tolist(), y=r[:, 1].tolist(),
                mode="lines", fill="toself", fillcolor=fillrgba,
                line=dict(color=color, width=1.2), name=f"▣ {cat}",
                legendgroup=f"region-{cat}", showlegend=(j == 0),
                hoverinfo="skip",
            ))


def _add_cluster_points(fig, plot_df, point_size, alpha, symbol, region_col, x_col, y_col):
    """Colour the (filtered) covered points by cluster, replacing the working
    filter. Noise / unlabelled points render grey and beneath the clusters, with
    the same palette the region overlays use so switching modes stays legible.

    Draw order is bottom-up: noise, then clusters largest to smallest, so a small
    cluster is never buried under a big one it overlaps.
    """
    if not region_col or region_col not in plot_df.columns:
        return
    sub = plot_df[_covered(plot_df, x_col, y_col)]
    if sub.empty:
        return
    cmap = _cluster_color_map(region_col)
    vals = sub[region_col].astype(str)
    counts = vals.value_counts()
    xs, ys = numeric_col(x_col), numeric_col(y_col)
    # Plotly stacks traces in the order they are added: noise first (bottom),
    # then clusters largest to smallest, so a small cluster is never buried
    # under a big one it overlaps.
    cats = sorted(vals.unique(), key=lambda c: (c in cmap, -int(counts[c]), c))
    # legendrank keeps the legend readable (clusters by name, noise last)
    # independently of that draw order.
    ranks = {c: 1001 + i for i, c in enumerate(sorted(cmap))}
    for cat in cats:
        idx = sub.index[vals == cat]
        in_map = cat in cmap
        fig.add_trace(go.Scattergl(
            x=xs.loc[idx].tolist(), y=ys.loc[idx].tolist(), mode="markers",
            marker=dict(size=point_size, color=cmap[cat] if in_map else NOISE_COLOR,
                        symbol=symbol, opacity=alpha,
                        line=dict(width=0.3, color="rgba(0,0,0,0.35)")),
            name=f"■ {cat}" if in_map else "■ noise",
            legendgroup=f"cpoint-{cat}", showlegend=True,
            legendrank=ranks.get(cat, 1001 + len(ranks)),
            customdata=id_str().loc[idx].tolist(),
            hovertemplate="<b>%{customdata}</b><extra></extra>",
        ))


def make_figure(cont_conds, bool_conds, cat_conds, tag_conds, color_mode, fixed_color, cont_col, colormap, reversed_, color_range_mode, alpha, point_size, bg_size, marker_size, marker_alpha, marker_mode, wf_position, wf_visible, layers, id_search_ids=None, wf_symbol=DEFAULT_SYMBOL, selection_ids=None, selection_color=DEFAULT_SELECTION_COLOR, region_col=None, region_shape="kde", region_position="below", region_opacity=0.25, region_coverage=0.4, region_tightness=0.25):
    df = _ann_df
    x_col, y_col = _x_col, _y_col
    if df.empty:
        return go.Figure(), 0, 0
    if not x_col or not y_col:
        fig = go.Figure()
        fig.update_layout(template="simple_white", annotations=[dict(text="No coordinates yet — run sse_coordinates.py, then reload.", showarrow=False, xref="paper", yref="paper", x=0.5, y=0.5, font=dict(size=14, color="#aaa"))])
        return fig, 0, 0

    fig = go.Figure()
    has_region_col = bool(region_col) and region_col in df.columns
    draw_regions = has_region_col and region_shape in ("kde", "hull")
    color_points_mode = has_region_col and region_shape == "points"
    if draw_regions and region_position == "below":
        _add_cluster_regions(fig, region_col, region_shape, region_opacity, region_coverage, region_tightness, x_col, y_col)
    coord_mask = _covered(df, x_col, y_col)
    filter_mask = apply_filters(df, cont_conds, bool_conds, cat_conds, tag_conds, id_search_ids)
    plot_df = df[filter_mask]
    n_filtered = int(filter_mask.sum())
    n_covered_filtered = int((filter_mask & coord_mask).sum())

    query_mask = boolean_mask(df[COL_QUERY]).fillna(False) if COL_QUERY in df.columns else pd.Series(False, index=df.index)
    # When the working filter is active, keep only its excluded points in the
    # subdued background trace.  Drawing every point here and relying on the
    # coloured trace to cover matches is unreliable (notably with transparent
    # markers and WebGL rendering) and makes filtered-in points look grey too.
    working_filter_active = bool(cont_conds or bool_conds or cat_conds or tag_conds or id_search_ids)
    bg_mask = coord_mask & ~query_mask
    if working_filter_active and wf_visible:
        bg_mask &= ~filter_mask
    bg_df = df[bg_mask]
    fig.add_trace(go.Scattergl(
        x=numeric_col(x_col).loc[bg_df.index], y=numeric_col(y_col).loc[bg_df.index], mode="markers",
        marker=dict(size=bg_size, color="lightgrey", opacity=0.4), name="All sequences", hoverinfo="skip", showlegend=True,
    ))
    if marker_mode == "bottom":
        _add_query_traces(fig, df, _id_col, x_col, y_col, marker_size, marker_alpha)

    cont_axis_idx = 0
    if wf_visible and wf_position == "bottom":
        if color_points_mode:
            _add_cluster_points(fig, plot_df, point_size, alpha, wf_symbol, region_col, x_col, y_col)
        else:
            cont_axis_idx = _add_working_filter(fig, plot_df, df, color_mode, fixed_color, cont_col, colormap, reversed_, color_range_mode, alpha, point_size, cont_axis_idx, _id_col, x_col, y_col, wf_symbol)

    for layer in reversed([l for l in (layers or []) if l.get("visible", True)]):
        ids = set(str(x) for x in layer.get("ids", []))
        if not ids:
            continue
        sub = df[id_str().isin(ids)]
        sub = sub[_covered(sub, x_col, y_col)]
        if sub.empty:
            continue
        style = layer.get("style", {})
        l_alpha = style.get("alpha", alpha)
        l_size = style.get("point_size", point_size)
        l_symbol = style.get("marker_symbol", DEFAULT_SYMBOL)
        if style.get("color_mode") == "continuous" and style.get("cont_col") in sub.columns:
            ccol = style.get("cont_col")
            vals = numeric_col(ccol).loc[sub.index]
            has_val = vals.notna()
            if has_val.any():
                v = vals[has_val]
                hv_idx = v.index
                cont_axis_idx += 1
                if style.get("color_range") == "global":
                    full = numeric_col(ccol).dropna()
                    cmin, cmax = float(full.min()), float(full.max())
                else:
                    cmin, cmax = float(v.min()), float(v.max())
                _add_cont_trace(fig, cont_axis_idx, layer.get("name", "Layer"), id_str().loc[hv_idx].tolist(), numeric_col(x_col).loc[hv_idx].tolist(), numeric_col(y_col).loc[hv_idx].tolist(), v.tolist(), cmin, cmax, style.get("colormap", DEFAULT_COLORMAP), style.get("reversed", False), l_size, l_alpha, l_symbol)
        else:
            fig.add_trace(go.Scattergl(
                x=numeric_col(x_col).loc[sub.index], y=numeric_col(y_col).loc[sub.index], mode="markers",
                marker=dict(size=l_size, color=style.get("fixed_color", DEFAULT_FIXED_COLOR), symbol=l_symbol, opacity=l_alpha, line=dict(width=0.5, color="black")),
                name=layer.get("name", "Layer"), customdata=id_str().loc[sub.index].tolist(), hovertemplate="<b>%{customdata}</b><extra></extra>",
            ))

    if wf_visible and wf_position == "top":
        if color_points_mode:
            _add_cluster_points(fig, plot_df, point_size, alpha, wf_symbol, region_col, x_col, y_col)
        else:
            cont_axis_idx = _add_working_filter(fig, plot_df, df, color_mode, fixed_color, cont_col, colormap, reversed_, color_range_mode, alpha, point_size, cont_axis_idx, _id_col, x_col, y_col, wf_symbol)
    # "Above" is added after the working-filter trace (the actual data points),
    # not just the grey background, so it sits over every point cloud. Query
    # markers and the selection highlight are drawn next, staying on top.
    if draw_regions and region_position == "above":
        _add_cluster_regions(fig, region_col, region_shape, region_opacity, region_coverage, region_tightness, x_col, y_col)
    if marker_mode == "top":
        _add_query_traces(fig, df, _id_col, x_col, y_col, marker_size, marker_alpha)
    if selection_ids:
        sel_df = df[id_str().isin([str(x) for x in selection_ids])]
        sel_df = sel_df[_covered(sel_df, x_col, y_col)]
        if not sel_df.empty:
            fig.add_trace(go.Scattergl(
                x=numeric_col(x_col).loc[sel_df.index], y=numeric_col(y_col).loc[sel_df.index], mode="markers",
                marker=dict(size=point_size + 6, symbol="circle-open", color=selection_color, opacity=1.0, line=dict(width=2.5, color=selection_color)),
                name=f"Selected ({len(sel_df)})", customdata=id_str().loc[sel_df.index].tolist(), hovertemplate="<b>%{customdata}</b><br>Selected<extra></extra>",
            ))
    fig.update_layout(
        xaxis=dict(title=x_col, showgrid=False),
        yaxis=dict(title=y_col, scaleanchor="x", scaleratio=1, showgrid=False),
        template="simple_white", dragmode="pan", uirevision="sequence-space",
        legend=dict(itemsizing="constant", yanchor="bottom", y=0.01, xanchor="right", x=0.99),
        margin=dict(l=40, r=20, t=40, b=40),
    )
    return fig, n_filtered, n_covered_filtered


def make_details_panel(sequence_id: str):
    if _ann_df.empty:
        return html.Div("No data loaded.")
    row = _ann_df.loc[_ann_df[_id_col].astype(str) == str(sequence_id)]
    if row.empty:
        return html.Div(f"No metadata found for {sequence_id}.")
    row = row.iloc[0]
    groups = {"Identity": [], "Taxonomy": [], "Scores & annotations": [], "Sequence": []}
    for col in _ann_df.columns:
        if col == _id_col or _types.get(col) == TYPE_COORDINATE:
            continue
        val = row.get(col, "")
        val = "" if pd.isna(val) else str(val)
        if not val.strip():
            continue
        key = col.lower()
        if col == COL_SEQ or "sequence" in key:
            group = "Sequence"
        elif any(x in key for x in ("tax", "organism", "species", "genus", "family", "phylum", "kingdom")):
            group = "Taxonomy"
        elif any(x in key for x in ("name", "label", "accession", "description")):
            group = "Identity"
        else:
            group = "Scores & annotations"
        display_val = val
        value_class = "detail-value detail-sequence" if group == "Sequence" else "detail-value"
        groups[group].append(html.Div([
            html.Div(col.replace("_", " "), className="detail-key", title=col),
            html.Div(display_val, className=value_class, title=display_val),
        ], className="detail-row"))

    sections = [
        html.Section([html.H5(name), html.Div(items, className="detail-grid")], className="detail-section")
        for name, items in groups.items() if items
    ]
    return html.Div([
        html.Div([
            html.Div([html.Div("Sequence inspector", className="detail-eyebrow"), html.H3(sequence_id, title=sequence_id)]),
            html.Span("Selected", className="status-pill"),
        ], className="detail-header"),
        *sections,
    ], className="details-panel")


def make_filter_panel():
    if not _col_meta:
        return html.P("No columns loaded.", style={"fontSize": "12px", "color": "var(--text-faint)"})
    df = _ann_df
    children = []
    cont_cols = [c for c in cols_of_type("continuous") if _types.get(c) == TYPE_LABEL]
    bool_cols = [c for c in cols_of_type("boolean") if _types.get(c) == TYPE_LABEL]
    cat_cols = [c for c in cols_of_type("categorical") if _types.get(c) == TYPE_LABEL]
    tag_cols = [c for c in cols_of_type("tag_split") if _types.get(c) == TYPE_LABEL]

    def checkbox(col_id, label):
        return dcc.Checklist(id={"type": "filter-enabled", "col": col_id}, options=[{"label": html.Span(label, style=LABEL_STYLE), "value": "on"}], value=[], inputStyle={"cursor": "pointer"})

    if cont_cols:
        children.append(html.P("Continuous", style=SECTION_STYLE))
        for col in cont_cols:
            vmin, vmax, step, marks = slider_config(df[col])
            inp = {"width": "70px", "fontSize": "11px", "padding": "2px 4px", "border": "1px solid var(--border)", "borderRadius": "3px", "textAlign": "center"}
            children.append(html.Div([checkbox(col, col), html.Div(id={"type": "filter-control", "col": col}, children=[html.Div([dcc.Input(id={"type": "cont-min-input", "col": col}, type="number", value=vmin, step="any", debounce=True, style=inp), html.Span("–", style={"margin": "0 4px", "fontSize": "11px", "color": "var(--text-faint)"}), dcc.Input(id={"type": "cont-max-input", "col": col}, type="number", value=vmax, step="any", debounce=True, style=inp)], style={"display": "flex", "alignItems": "center", "marginBottom": "4px"}), dcc.RangeSlider(id={"type": "cont-slider", "col": col}, min=vmin, max=vmax, step=step, value=[vmin, vmax], marks=marks, allowCross=False, updatemode="mouseup")], style={"display": "none", **CONTROL_WRAPPER})], style={"marginBottom": "6px"}))
    if bool_cols:
        children.append(html.P("Boolean", style=SECTION_STYLE))
        for col in bool_cols:
            children.append(html.Div([checkbox(col, col), html.Div(id={"type": "filter-control", "col": col}, children=[dcc.RadioItems(id={"type": "bool-filter", "col": col}, options=[{"label": " All", "value": "all"}, {"label": " True", "value": "true"}, {"label": " False", "value": "false"}], value="all", labelStyle={"display": "inline-block", "marginRight": "10px", "fontSize": "13px"})], style={"display": "none", **CONTROL_WRAPPER})], style={"marginBottom": "6px"}))
    if cat_cols:
        children.append(html.P("Categorical", style=SECTION_STYLE))
        for col in cat_cols:
            opts = sorted(df[col].replace("", pd.NA).dropna().astype(str).unique().tolist())
            children.append(html.Div([checkbox(col, col), html.Div(id={"type": "filter-control", "col": col}, children=[dcc.Dropdown(id={"type": "cat-filter", "col": col}, options=[{"label": v, "value": v} for v in opts], value=[], multi=True, placeholder="Select values…", style={"fontSize": "12px"})], style={"display": "none", **CONTROL_WRAPPER})], style={"marginBottom": "6px"}))
    if tag_cols:
        children.append(html.P("Tag-split", style=SECTION_STYLE))
        for col in tag_cols:
            tags = _col_meta[col].get("tags", [])
            children.append(html.Div([checkbox(col, col), html.Div(id={"type": "filter-control", "col": col}, children=[dcc.Dropdown(id={"type": "tag-filter", "col": col}, options=[{"label": t, "value": t} for t in tags], value=[], multi=True, placeholder="Select tags…", style={"fontSize": "12px"})], style={"display": "none", **CONTROL_WRAPPER})], style={"marginBottom": "6px"}))
    return children or html.P("No filterable label columns.", style={"fontSize": "12px", "color": "var(--text-faint)"})


def make_col_settings_panel():
    opts = [{"label": "Continuous", "value": "continuous"}, {"label": "Boolean", "value": "boolean"}, {"label": "Categorical", "value": "categorical"}, {"label": "Tag split", "value": "tag_split"}, {"label": "Skip", "value": "skip"}]
    rows = []
    for col, meta in _col_meta.items():
        if _types.get(col) != TYPE_LABEL:
            continue
        rows.append(html.Div([html.Span(col, title=col, style={"fontSize": "11px", "fontWeight": "500", "overflow": "hidden", "textOverflow": "ellipsis", "whiteSpace": "nowrap", "flexGrow": "1"}), dcc.Dropdown(id={"type": "col-override", "col": col}, options=opts, value=effective_type(col), clearable=False, style={"fontSize": "11px", "width": "120px"})], style={"display": "flex", "alignItems": "center", "gap": "6px", "padding": "4px 2px", "borderBottom": "1px solid var(--border-soft)"}))
    return html.Div(rows, style={"maxHeight": "350px", "overflowY": "auto"})


def make_sidebar(layers):
    if not layers:
        return html.Div("No saved layers yet.", style={"fontSize": "12px", "color": "var(--text-faint)", "padding": "8px 0"})
    rows = []
    n = len(layers)
    for i, layer in enumerate(layers):
        lid = layer["id"]
        visible = layer.get("visible", True)
        style = layer.get("style", {})
        swatch_style = {"width": "12px", "height": "12px", "borderRadius": "2px", "border": "1px solid var(--border)", "flexShrink": "0"}
        if style.get("color_mode") == "continuous":
            swatch = html.Div(style={**swatch_style, "background": "linear-gradient(to bottom, #440154, #31688e, #35b779, #fde725)"})
        else:
            swatch = html.Div(style={**swatch_style, "backgroundColor": style.get("fixed_color", DEFAULT_FIXED_COLOR)})
        btn = {"background": "none", "border": "none", "cursor": "pointer", "padding": "1px 3px", "fontSize": "11px", "lineHeight": "1"}
        n_total = len(layer.get("ids", []))
        n_cov = n_total
        if _x_col and _y_col and not _ann_df.empty:
            sub = _ann_df[id_str().isin([str(x) for x in layer.get("ids", [])])]
            n_cov = int(_covered(sub, _x_col, _y_col).sum())
        rows.append(html.Div([swatch, html.Div([html.Div(layer.get("name", "Layer"), title=layer.get("name", "Layer"), style={"fontSize": "11px", "fontWeight": "500", "overflow": "hidden", "textOverflow": "ellipsis", "whiteSpace": "nowrap"}), html.Div(f"{n_cov:,}/{n_total:,} visible here", style={"fontSize": "10px", "color": "var(--text-faint)"})], style={"flexGrow": "1", "minWidth": "0", "margin": "0 4px"}), html.Button(html.Span(className="ic ic-up"), id={"type": "layer-up", "lid": lid}, disabled=i == 0, title="Move up", className="sse-icon-btn", style=btn), html.Button(html.Span(className="ic ic-down"), id={"type": "layer-down", "lid": lid}, disabled=i == n - 1, title="Move down", className="sse-icon-btn", style=btn), html.Button(html.Span(className="ic ic-load"), id={"type": "layer-load", "lid": lid}, title="Load into working filter", className="sse-icon-btn", style={**btn, "color": "var(--accent)"}), html.Button(html.Span(className="ic ic-eye" if visible else "ic ic-eyeoff"), id={"type": "layer-toggle", "lid": lid}, title="Show/hide layer", className="sse-icon-btn", style={**btn, "opacity": "1" if visible else "0.5"}), html.Button(html.Span(className="ic ic-x"), id={"type": "layer-delete", "lid": lid}, title="Delete layer", className="sse-icon-btn", style={**btn, "color": "var(--danger)"})], style={"display": "flex", "alignItems": "center", "padding": "5px 2px", "borderBottom": "1px solid var(--border-soft)", "opacity": "1" if visible else "0.5"}))
    return html.Div(rows, style={"maxHeight": "400px", "overflowY": "auto"})


SSE_INDEX_STRING = """<!DOCTYPE html>
<html lang="en" data-theme="pipeline">
<head>
{%metas%}
<title>{%title%}</title>
{%favicon%}
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Geist:wght@400;500;600;700&family=Geist+Mono:wght@400;500&display=swap" rel="stylesheet">
{%css%}
<style>
:root,
html[data-theme="clean-lab"] {
  --page:        #f5f9f9;
  --surface:     #fffffe;
  --surface-2:   #e3f6f5;
  --border:      #d2e4e4;
  --border-soft: #e8f2f2;
  --text:        #272343;
  --text-muted:  #545b70;
  --text-faint:  #9aa0b0;

  --accent:      #0f7c78;
  --accent-weak: #d3ecec;
  --on-accent:   #fffffe;

  --success:     #2e7d64;
  --warning:     #b4622f;
  --danger:      #c0392b;
  --boltz:       #6b46c1;

  --shadow-card: 0 1px 2px rgba(39,35,67,.05), 0 1px 3px rgba(39,35,67,.04);
  color-scheme: light;
}

html[data-theme="dark-lab"] {
  --page:        #14121d;
  --surface:     #201d30;
  --surface-2:   #2a2743;
  --border:      #38344f;
  --border-soft: #2a2740;
  --text:        #ecebf4;
  /* darkened from the previous grey so it reads on light dcc cells too */
  --text-muted:  #c2c7d6;
  --text-faint:  #8b90a6;

  --accent:      #34c6bd;
  --accent-weak: #143f3d;
  --on-accent:   #15121f;

  --success:     #4bbf8a;
  --warning:     #d98b4a;
  --danger:      #e5695f;
  --boltz:       #a78bfa;

  --shadow-card: 0 1px 2px rgba(0,0,0,.45), 0 1px 3px rgba(0,0,0,.35);
  color-scheme: dark;
}

/* Rose Quartz — light theme built from the "Happy Hues" palette. Raw palette colours
   are noted inline; muted/faint text, weak tints and the shadow are derived from them.
   --on-accent is the dark headline navy (not white) because #a786df is a light purple:
   navy on it is ~6.4:1, white is ~2.8:1. */
html[data-theme="rose-quartz"] {
  --page:        #fec7d7;   /* Background */
  --surface:     #efecf6;   /* pale lavender-grey fields (lightened from Secondary #d9d4e7) */
  --surface-2:   #d9d4e7;   /* Secondary — slightly deeper inset */
  --border:      #c7c0dd;   /* darker lavender so cards read on the tint */
  --border-soft: #ded9ec;
  --text:        #0e172c;   /* Headline / Paragraph */
  --text-muted:  #3d4560;
  --text-faint:  #8a86a0;

  --accent:      #a786df;   /* Tertiary */
  --accent-weak: #e9ddf9;
  --on-accent:   #0e172c;   /* dark navy reads on the light-purple accent */

  --success:     #2f9e73;
  --warning:     #c9722f;
  --danger:      #d6455f;
  --boltz:       #6d3fc4;   /* deeper violet, distinct from --accent */

  --shadow-card: 0 1px 2px rgba(14,23,44,.08), 0 1px 3px rgba(14,23,44,.06);
  color-scheme: light;
}

/* Deep Canopy — dark teal theme from the "Happy Hues" palette. Raw palette colours
   noted inline; faint text, weak tints, surfaces and shadow are derived. --accent is
   the amber button colour, so --on-accent is the near-black teal (dark on light amber). */
html[data-theme="deep-canopy"] {
  --page:        #001e1d;   /* Stroke / Button text (near-black teal) */
  --surface:     #004643;   /* Background (deep-teal cards) */
  --surface-2:   #0a5a54;   /* lighter teal inset */
  --border:      #10695f;
  --border-soft: #0a4a44;
  --text:        #fffffe;   /* Headline */
  --text-muted:  #abd1c6;   /* Paragraph (sage) */
  --text-faint:  #7fa99e;

  --accent:      #f9bc60;   /* Button (amber) */
  --accent-weak: #3f3218;   /* dark amber tint */
  --on-accent:   #001e1d;   /* near-black teal reads on the amber accent */

  --success:     #4fc98a;
  --warning:     #e8964f;
  --danger:      #e16162;   /* Tertiary (coral) */
  --boltz:       #a78bfa;   /* violet, stays distinct on teal */

  --shadow-card: 0 1px 2px rgba(0,0,0,.45), 0 1px 3px rgba(0,0,0,.35);
  color-scheme: dark;
}

/* Pipeline — mirrors the Pipeline Control Center app (deep-teal glass, Geist,
   amber/lime accents). Palette lifted 1:1 from pipeline-ui-poc/app/globals.css
   (--ink/--muted/--teal/--lime/--amber/--red); surfaces are the solid form of
   that app's panel gradient, borders are its teal hairlines. The ambient shell,
   panel gradients and teal-gradient buttons live in the polish block below. */
html[data-theme="pipeline"] {
  --page:        #101615;   /* --canvas (near-black teal behind the shell) */
  --surface:     #06333c;   /* solid of the panel gradient rgba(5,47,56,.83) */
  --surface-2:   #0a4650;   /* lighter teal inset */
  --border:      rgba(96,205,195,.20);   /* --line-strong hairline */
  --border-soft: rgba(96,205,195,.10);   /* --line */
  --text:        #edf8f4;   /* --ink */
  --text-muted:  #9dbeba;   /* --muted, lifted a touch for dcc cells */
  --text-faint:  #6f9691;   /* --faint */

  --accent:      #0eb5a4;   /* --teal */
  --accent-weak: rgba(14,181,164,.16);
  --on-accent:   #ecfffa;

  --lime:        #b8e94f;   /* --lime (the signature Pipeline pop accent) */
  --success:     #b8e94f;   /* positive == lime, as in the Pipeline app */
  --warning:     #e6bf56;   /* --amber */
  --danger:      #ff8a8a;   /* --red */
  --boltz:       #a78bfa;   /* violet, stays distinct on teal */

  /* Near-flat: only a faint inner highlight, so panels read like the Pipeline
     app's cards (flush) rather than floating on a drop shadow. */
  --shadow-card: inset 0 1px rgba(255,255,255,.02);
  color-scheme: dark;
}

/* page + default text follow the theme (covers gutters too) */
html, body { background: var(--page); color: var(--text); }
html, body { margin:0; min-height:100%; }
.sse-root  { color: var(--text);
  font-family: -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; }

/* control panels rendered as cards (top-level <details> in each column) */
.sse-root .sse-col > details {
  background: var(--surface); border: 1px solid var(--border);
  border-radius: 7px; box-shadow: var(--shadow-card);
  padding: 6px 12px 10px; margin-bottom: 12px; }
.sse-root .sse-col > details > summary { margin-bottom: 8px; }

/* labels + option text track the theme (fixes options vanishing in Dark Lab) */
.sse-root label, .sse-root label span { color: var(--text) !important; }

/* Dash input (dash-input-container > dash-input-element + stepper buttons):
   the container is the single white cell; the inner element contributes no box
   (no border AND no focus outline — the outline is what created the inner box);
   steppers hidden. Focus ring goes on the container instead. */
.sse-root [class*="dash-input-container"] {
  background: #ffffff !important; border: 1px solid var(--border) !important;
  border-radius: 4px !important; box-shadow: none !important; }
.sse-root [class*="dash-input-container"]:focus-within {
  outline: 2px solid var(--accent); outline-offset: 0; }
.sse-root [class*="dash-input-element"],
.sse-root input[type=text],
.sse-root input[type=number],
.sse-root input[type=password] {
  background: transparent !important; color: #1a1a1a !important;
  border: none !important; outline: none !important; box-shadow: none !important; }
.sse-root textarea {
  background: #ffffff !important; color: #1a1a1a !important;
  border: 1px solid var(--border) !important; border-radius: 4px !important; }
.sse-root input::placeholder, .sse-root textarea::placeholder { color: #9aa0b0; }
.sse-root [class*="dash-input-stepper"],
.sse-root [class*="dash-stepper"] { display: none !important; }
.sse-root input[type=number] { -moz-appearance: textfield; }
.sse-root input[type=number]::-webkit-inner-spin-button,
.sse-root input[type=number]::-webkit-outer-spin-button { -webkit-appearance: none; margin: 0; }

/* Dash dropdown — closed control: white cell with dark text. Scoped to the
   dash-dropdown-* classes ONLY, so RadioItems/Checklist (which share the generic
   dash-options-list-option class) are left to the theme's label rule. */
[class*="dash-dropdown"] { background: #ffffff !important; }
[class*="dash-dropdown"] * { color: #1a1a1a !important; }
/* Dropdown popup — always a white menu with black text, both themes.
   Each option is a <label> and its text is a descendant node. The theme's generic
   `.sse-root label span { color: var(--text) !important }` (light) also hits the option
   text and, on the white menu, made it faint. Root cause is specificity, not colour: that
   rule is (0,1,2), so any override must exceed it. We anchor on the dropdown-only
   dash-dropdown-options container plus the <label> itself (confirmed present in devtools),
   which reaches (0,2,1) > (0,1,2) and does NOT depend on the option carrying a
   dash-dropdown-option class (builds vary). The container anchor keeps RadioItems /
   Checklist — which live under dash-options-list, not dash-dropdown-options — untouched,
   so they still read as light text on the dark panel. The unscoped variant covers builds
   where the menu is portaled outside .sse-root. */
[class*="dash-dropdown-options"] label,
[class*="dash-dropdown-options"] label *,
.sse-root [class*="dash-dropdown-options"] label,
.sse-root [class*="dash-dropdown-options"] label * {
  color: #1a1a1a !important; }

/* checkbox / radio / range accent */
.sse-root input[type=checkbox],
.sse-root input[type=radio],
.sse-root input[type=range] { accent-color: var(--accent); }
.sse-root :focus-visible { outline: 2px solid var(--accent); outline-offset: 1px; }

/* -------------------------------------------------------------------------
   dcc component internals. !important is used deliberately to beat dcc's
   bundled component CSS. Verify class names in devtools if a control does not
   pick up the theme (rc-slider = .rc-slider-*).
   ------------------------------------------------------------------------- */

/* dcc.Slider / RangeSlider (rc-slider) -> teal accent */
.sse-root .rc-slider-rail   { background: var(--border) !important; }
.sse-root .rc-slider-track  { background: var(--accent) !important; }
.sse-root .rc-slider-handle { border-color: var(--accent) !important; background: var(--surface) !important; }
.sse-root .rc-slider-handle:hover,
.sse-root .rc-slider-handle:active,
.sse-root .rc-slider-handle:focus { border-color: var(--accent) !important; box-shadow: 0 0 0 4px var(--accent-weak) !important; }
.sse-root .rc-slider-dot,
.sse-root .rc-slider-dot-active { border-color: var(--accent) !important; }
.sse-root .rc-slider-mark-text  { color: var(--text-faint); }

/* dash Slider/RangeSlider — Radix build (.dash-slider-* classes) */
.sse-root .dash-slider-track { background: var(--border) !important; }
.sse-root .dash-slider-track > span,
.sse-root .dash-slider-range,
.sse-root [class*="dash-slider"][class*="range"] { background: var(--accent) !important; }
.sse-root .dash-slider-thumb,
.sse-root [class*="dash-slider-thumb"] {
  background: var(--accent) !important; border-color: var(--accent) !important; box-shadow: none !important; }
.sse-root .dash-slider-thumb:focus,
.sse-root [class*="dash-slider-thumb"]:focus { box-shadow: 0 0 0 4px var(--accent-weak) !important; }
.sse-root .dash-slider-mark,
.sse-root .dash-slider-mark-outside-selection { color: var(--text-muted) !important; background: transparent !important; }

/* =========================================================================
   Redesign layer (2026-07): additive polish only — no callback changes.
   Inline mask-image icons (colour follows currentColor), card/summary polish,
   a consistent button hierarchy via :hover (filter/!important beat inline base
   styles without touching them), a plot "hero" card, and collapsible rails.
   ------------------------------------------------------------------------- */
:root {
  --i-chev:  url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='2.4' stroke-linecap='round' stroke-linejoin='round'%3E%3Cpath d='M9 6l6 6-6 6'/%3E%3C/svg%3E");
  --i-up:    url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Cpath d='M12 19V5M6 11l6-6 6 6'/%3E%3C/svg%3E");
  --i-down:  url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Cpath d='M12 5v14M6 13l6 6 6-6'/%3E%3C/svg%3E");
  --i-load:  url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Cpath d='M4 12h12M12 8l4 4-4 4'/%3E%3C/svg%3E");
  --i-eye:   url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Cpath d='M2 12s3.5-7 10-7 10 7 10 7-3.5 7-10 7-10-7-10-7z'/%3E%3Ccircle cx='12' cy='12' r='2.6'/%3E%3C/svg%3E");
  --i-eyeoff:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Cpath d='M4 4l16 16M9.9 5.2A9.6 9.6 0 0 1 12 5c6.5 0 10 7 10 7a15 15 0 0 1-3.4 4M6.1 7.9A15 15 0 0 0 2 12s3.5 7 10 7a9.6 9.6 0 0 0 3.2-.5'/%3E%3C/svg%3E");
  --i-x:     url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 24 24' fill='none' stroke='black' stroke-width='2' stroke-linecap='round' stroke-linejoin='round'%3E%3Cpath d='M6 6l12 12M18 6L6 18'/%3E%3C/svg%3E");
}
.sse-root .ic { display:inline-block; width:14px; height:14px; flex:none;
  background-color: currentColor; vertical-align:-2px;
  -webkit-mask-repeat:no-repeat; mask-repeat:no-repeat;
  -webkit-mask-position:center; mask-position:center;
  -webkit-mask-size:contain; mask-size:contain; }
.sse-root .ic-up{ -webkit-mask-image:var(--i-up); mask-image:var(--i-up); }
.sse-root .ic-down{ -webkit-mask-image:var(--i-down); mask-image:var(--i-down); }
.sse-root .ic-load{ -webkit-mask-image:var(--i-load); mask-image:var(--i-load); }
.sse-root .ic-eye{ -webkit-mask-image:var(--i-eye); mask-image:var(--i-eye); }
.sse-root .ic-eyeoff{ -webkit-mask-image:var(--i-eyeoff); mask-image:var(--i-eyeoff); }
.sse-root .ic-x{ -webkit-mask-image:var(--i-x); mask-image:var(--i-x); }

/* Top-level panels: chevron affordance + hover, marker removed */
.sse-root .sse-col > details > summary { list-style:none; display:flex; align-items:center; gap:7px; }
.sse-root .sse-col > details > summary::-webkit-details-marker { display:none; }
.sse-root .sse-col > details > summary::before {
  content:""; width:13px; height:13px; flex:none; color:var(--text-faint);
  background-color: currentColor;
  -webkit-mask:var(--i-chev) center/contain no-repeat; mask:var(--i-chev) center/contain no-repeat;
  transform: rotate(0deg); transition: transform .18s ease; }
.sse-root .sse-col > details[open] > summary::before { transform: rotate(90deg); }
.sse-root .sse-col > details > summary:hover { color: var(--accent); }
.sse-root .sse-col > details { transition: border-color .15s ease, box-shadow .15s ease; }
.sse-root .sse-col > details:hover { border-color: var(--accent); }

/* Button hierarchy — hover states layered on top of existing inline styles */
.sse-root button { transition: filter .12s ease, background-color .12s ease, border-color .12s ease; }
.sse-root .sse-btn-primary:hover { filter: brightness(1.07); }
.sse-root .sse-btn-primary:active { filter: brightness(.95); }
.sse-root .sse-btn-sec:hover { background: var(--surface-2) !important; border-color: var(--accent) !important; }
.sse-root .sse-icon-btn:hover { background: var(--surface-2) !important; color: var(--accent) !important; }
.sse-root .sse-icon-btn { border-radius:6px; }

/* Coordinate field labels sat flush against each other ("Coordinate systemX axis") */
.sse-root #coord-system-panel > label, .sse-root #coord-free-panel > label {
  display:block; margin:9px 0 3px; }

/* Plot as hero: rounded card around the graph */
.sse-root .sse-plot-card { background: var(--surface); border:1px solid var(--border);
  border-radius:12px; box-shadow: var(--shadow-card); padding:8px 8px 4px; }
.sse-root .sse-plot-card .js-plotly-plot .plot-container { border-radius:8px; }
/* Plotly modebar: follow the theme instead of fixed dark-on-dark icons */
.sse-root .sse-plot-card .modebar { background: transparent !important; }
.sse-root .sse-plot-card .modebar-btn .icon path { fill: var(--text-faint) !important; }
.sse-root .sse-plot-card .modebar-btn.active .icon path,
.sse-root .sse-plot-card .modebar-btn:hover .icon path { fill: var(--accent) !important; }

/* Collapsible rails */
.sse-root .sse-rail { transition: width .2s ease, min-width .2s ease, opacity .18s ease, padding .2s ease; }
.sse-root .rail-toggle { background:none; border:none; cursor:pointer; color:var(--text-faint);
  font-size:16px; line-height:1; padding:4px 7px; border-radius:6px; }
.sse-root .rail-toggle:hover { background:var(--surface-2); color:var(--accent); }

/* Workspace shell and typography */
.sse-root { font-size:14px; }
.sse-root .sse-header { position:sticky; top:0; z-index:20; padding:14px 0 12px;
  background:color-mix(in srgb, var(--page) 92%, transparent); backdrop-filter:blur(12px); }
.sse-root .sse-header h2 { font-size:24px; letter-spacing:-.025em; }
.sse-root .sse-workspace { min-height:calc(100vh - 104px); }
.sse-root .sse-rail { position:sticky; top:88px; max-height:calc(100vh - 104px) !important;
  scrollbar-width:thin; scrollbar-color:var(--border) transparent; }
.sse-root .sse-center { padding-bottom:32px; }
.sse-root .sse-toolbar { min-height:34px; }
.sse-root .selection-toolbar { padding:8px 10px !important; border-radius:9px !important;
  box-shadow:var(--shadow-card); }
.sse-root .plot-status { font-variant-numeric:tabular-nums; }

/* Details become a scannable inspector instead of an undifferentiated table. */
.sse-root .details-card { padding:0 !important; overflow:hidden; }
.sse-root .details-empty { padding:20px; color:var(--text-muted); text-align:center; }
.sse-root .details-panel { padding:18px; }
.sse-root .detail-header { display:flex; align-items:flex-start; justify-content:space-between;
  gap:16px; padding-bottom:14px; border-bottom:1px solid var(--border); }
.sse-root .detail-header h3 { margin:3px 0 0; font-size:18px; line-height:1.25;
  overflow-wrap:anywhere; }
.sse-root .detail-eyebrow { color:var(--text-muted); font-size:11px; font-weight:700;
  letter-spacing:.08em; text-transform:uppercase; }
.sse-root .status-pill { flex:none; padding:4px 9px; border-radius:999px;
  color:var(--accent); background:var(--accent-weak); font-size:11px; font-weight:700; }
.sse-root .detail-section { margin-top:16px; }
.sse-root .detail-section h5 { margin:0 0 8px; color:var(--text-muted); font-size:11px;
  letter-spacing:.07em; text-transform:uppercase; }
.sse-root .detail-grid { display:grid; grid-template-columns:repeat(auto-fit, minmax(210px, 1fr)); gap:1px;
  overflow:hidden; border:1px solid var(--border); border-radius:8px; background:var(--border); }
.sse-root .detail-row { min-width:0; padding:9px 11px; background:var(--surface); }
.sse-root .detail-key { margin-bottom:3px; color:var(--text-muted); font-size:11px;
  font-weight:600; text-transform:capitalize; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
.sse-root .detail-value { line-height:1.4; overflow-wrap:anywhere; }
.sse-root .detail-sequence { max-height:5.6em; overflow:auto; font:11px/1.4 ui-monospace, SFMono-Regular, Menlo, monospace; }

/* Larger targets and clearer control rhythm. */
.sse-root button { min-height:30px; }
.sse-root .sse-icon-btn, .sse-root .rail-toggle { min-width:30px; min-height:30px; }
.sse-root .sse-col > details { border-radius:10px; padding:8px 12px 12px; }
.sse-root .sse-col > details > summary { min-height:30px; font-size:13px; }

@media (max-width: 1100px) {
  .sse-root { padding:14px !important; }
  .sse-root .sse-workspace { display:grid !important; grid-template-columns:minmax(240px, 280px) minmax(0, 1fr); }
  .sse-root #right-rail { position:static; grid-column:1 / -1; width:auto !important; min-width:0 !important;
    max-height:none !important; padding:18px 0 0 !important; display:grid; grid-template-columns:repeat(2, minmax(0, 1fr)); gap:12px; }
  .sse-root #right-rail > * { min-width:0; }
  .sse-root .sse-center { padding-right:0 !important; }
}

@media (max-width: 760px) {
  .sse-root .sse-header { position:static; align-items:flex-start !important; gap:12px; flex-direction:column; }
  .sse-root .sse-header > div:last-child { width:100%; }
  .sse-root .sse-header [class*="dash-dropdown"] { flex:1; }
  .sse-root .sse-workspace { display:flex !important; flex-direction:column; }
  .sse-root #left-rail, .sse-root #right-rail { position:static; width:100% !important; min-width:0 !important;
    max-height:none !important; padding:0 !important; opacity:1 !important; overflow:visible !important; }
  .sse-root #left-rail { order:2; margin-top:14px; }
  .sse-root #right-rail { order:3; display:block; }
  .sse-root .sse-center { order:1; padding:0 !important; }
  .sse-root .rail-toggle { display:none; }
  .sse-root .selection-toolbar { flex-wrap:wrap; gap:6px; }
  .sse-root .selection-toolbar > span { flex-basis:100%; }
  .sse-root .sse-plot-card { padding:4px; }
  .sse-root #latent-graph { height:62vh !important; min-height:460px; }
  .sse-root .detail-grid { grid-template-columns:1fr; }
}

@media (prefers-reduced-motion: reduce) { .sse-root * { transition: none !important; } }

/* =========================================================================
   Pipeline theme polish (data-theme="pipeline" only) — makes the visualizer
   read like the Pipeline Control Center: an ambient rounded "app-shell" with a
   perspective grid + teal glow, glassy panel gradients, teal-gradient primary
   buttons, and the Geist typeface. Scoped to the theme so the other four are
   untouched. Backgrounds use !important to beat the inline var(--surface) fills.
   ------------------------------------------------------------------------- */
html[data-theme="pipeline"] { font-synthesis-weight: none; }
html[data-theme="pipeline"] body {
  background:
    radial-gradient(circle at 16% 6%, rgba(0,118,126,.13), transparent 34%),
    radial-gradient(circle at 88% 26%, rgba(0,89,94,.11), transparent 30%),
    var(--page);
}

/* Rounded, bordered, glowing shell around the whole workspace (~.app-shell). */
html[data-theme="pipeline"] .sse-root {
  position: relative;
  isolation: isolate;
  font-family: "Geist", -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif !important;
  margin: 14px auto !important;
  padding: 20px 26px 40px !important;
  border: 1px solid rgba(91,207,197,.12);
  border-radius: 26px;
  background:
    radial-gradient(circle at 20% 4%, rgba(0,118,126,.15), transparent 33%),
    radial-gradient(circle at 86% 32%, rgba(0,89,94,.13), transparent 29%),
    linear-gradient(145deg, #001e27 0%, #002a32 48%, #001b24 100%);
  box-shadow: 0 18px 65px rgba(0,0,0,.34);
}
/* Ambient perspective grid + glow, fixed to the viewport, behind the content
   (z-index:-1 inside the isolated .sse-root stacking context). */
html[data-theme="pipeline"] .sse-root::before {
  content:""; position: fixed; z-index:-1; left:3%; bottom:-50px;
  width:60%; height:56%; opacity:.5; pointer-events:none;
  background-image:
    linear-gradient(rgba(31,201,190,.05) 1px, transparent 1px),
    linear-gradient(90deg, rgba(31,201,190,.05) 1px, transparent 1px);
  background-size:38px 38px;
  transform: perspective(430px) rotateX(63deg) rotateZ(-9deg);
  transform-origin:50% 100%;
  -webkit-mask-image: linear-gradient(90deg,#000 0%,#000 56%,transparent 94%);
  mask-image: linear-gradient(90deg,#000 0%,#000 56%,transparent 94%);
}
html[data-theme="pipeline"] .sse-root::after {
  content:""; position: fixed; z-index:-1; left:8%; bottom:6%;
  width:30%; height:26%; border-radius:50%; filter: blur(42px);
  background: rgba(0,186,178,.13); pointer-events:none;
}

/* Header: teal-tinted glass strip + glowing title (~.entry-header h1). */
html[data-theme="pipeline"] .sse-root .sse-header {
  background: color-mix(in srgb, #01222b 88%, transparent);
}
html[data-theme="pipeline"] .sse-root .sse-header h2 {
  color: #f3fbf8; letter-spacing: -.03em;
  text-shadow: 0 1px 22px rgba(20,186,174,.12);
}

/* Panels as flat cards (~.panel) — near-uniform teal fill, thin hairline,
   NO drop shadow, so they sit flush like the Pipeline app instead of floating. */
html[data-theme="pipeline"] .sse-root .sse-col > details,
html[data-theme="pipeline"] .sse-root .sse-plot-card,
html[data-theme="pipeline"] .sse-root .details-card,
html[data-theme="pipeline"] .sse-root .selection-toolbar {
  background: linear-gradient(150deg, rgba(7,52,61,.5), rgba(5,44,52,.52)) !important;
  border-color: rgba(96,205,195,.13) !important;
  box-shadow: inset 0 1px rgba(255,255,255,.02) !important;
  -webkit-backdrop-filter: none; backdrop-filter: none;
}
html[data-theme="pipeline"] .sse-root .sse-col > details:hover {
  border-color: rgba(96,205,195,.28) !important;
}
/* Open panel gets a subtle lime spine, echoing the Pipeline app's active nav. */
html[data-theme="pipeline"] .sse-root .sse-col > details[open] {
  box-shadow: inset 2px 0 0 rgba(184,233,79,.55), inset 0 1px rgba(255,255,255,.02) !important;
}

/* Let the white plot fill the card and take its rounded corners, so it reads as
   part of the UI rather than a white square dropped on top. Drop the card padding
   and clip the Plotly surface to the card radius. */
html[data-theme="pipeline"] .sse-root .sse-plot-card {
  padding: 0 !important;
  overflow: hidden;
  border-radius: 14px !important;
}
html[data-theme="pipeline"] .sse-root .sse-plot-card .js-plotly-plot,
html[data-theme="pipeline"] .sse-root .sse-plot-card .plot-container,
html[data-theme="pipeline"] .sse-root .sse-plot-card .svg-container,
html[data-theme="pipeline"] .sse-root .sse-plot-card .main-svg {
  border-radius: 14px !important;
}

/* Lime accents — the pops that make the Pipeline app read as "alive". */
html[data-theme="pipeline"] .sse-root .status-pill {
  color: #dbff8a !important;
  background: rgba(184,233,79,.12) !important;
  box-shadow: inset 0 0 0 1px rgba(184,233,79,.22);
}
html[data-theme="pipeline"] .sse-root .detail-eyebrow,
html[data-theme="pipeline"] .sse-root .detail-section h5 {
  color: #b8e94f !important;
}
/* Checked filters/toggles glow lime (radios stay teal via accent-color below). */
html[data-theme="pipeline"] .sse-root input[type=checkbox] { accent-color: #a9dd44; }

/* Primary buttons -> teal gradient (~.primary-button). */
html[data-theme="pipeline"] .sse-root .sse-btn-primary {
  background: linear-gradient(90deg, #0b958b, #0eb5a4) !important;
  color: #ecfffa !important; border: 0 !important;
  box-shadow: 0 8px 22px rgba(0,169,156,.18), inset 0 1px rgba(255,255,255,.13) !important;
}
html[data-theme="pipeline"] .sse-root .sse-btn-primary:hover {
  background: linear-gradient(90deg, #0eb5a4, #18c2af) !important;
}
/* Secondary buttons -> teal-outlined glass (~.secondary-button). */
html[data-theme="pipeline"] .sse-root .sse-btn-sec {
  background: rgba(5,51,59,.55) !important;
  border-color: rgba(96,205,195,.22) !important;
  color: #bcd2cf !important;
}

/* Dark teal fields to match Pipeline's .field inputs. Overrides the base
   white-cell rules for this theme only. The CLOSED dropdown trigger
   (button.dash-dropdown) and its value go dark; the OPEN options menu keeps
   its white/black treatment (it is a separate dash-dropdown-options node, not
   a descendant of the trigger), so readability is preserved. */
html[data-theme="pipeline"] .sse-root [class*="dash-input-container"],
html[data-theme="pipeline"] .sse-root textarea,
html[data-theme="pipeline"] .sse-root button.dash-dropdown {
  background: rgba(0,32,40,.72) !important;
  border-color: rgba(96,205,195,.24) !important;
}
/* The base [class*="dash-dropdown"]{background:#fff} paints the trigger's inner
   wrapper/value nodes white, hiding the dark button fill above. Make just those
   closed-control nodes transparent (the open menu is dash-dropdown-options and
   is untouched, so it keeps its white/black treatment). */
html[data-theme="pipeline"] .sse-root .dash-dropdown-wrapper,
html[data-theme="pipeline"] .sse-root .dash-dropdown-value,
html[data-theme="pipeline"] .sse-root .dash-dropdown-value-item,
html[data-theme="pipeline"] .sse-root .dash-dropdown-grid-container,
html[data-theme="pipeline"] .sse-root .dash-dropdown-trigger {
  background: transparent !important;
}
html[data-theme="pipeline"] .sse-root input[type=text],
html[data-theme="pipeline"] .sse-root input[type=number],
html[data-theme="pipeline"] .sse-root input[type=password],
html[data-theme="pipeline"] .sse-root [class*="dash-input-element"],
html[data-theme="pipeline"] .sse-root textarea,
html[data-theme="pipeline"] .sse-root button.dash-dropdown,
html[data-theme="pipeline"] .sse-root button.dash-dropdown span {
  color: #dcebe8 !important;
}
html[data-theme="pipeline"] .sse-root [class*="dash-input-container"]:focus-within,
html[data-theme="pipeline"] .sse-root button.dash-dropdown:focus-visible {
  outline: none !important; border-color: var(--accent) !important;
  box-shadow: 0 0 0 3px rgba(14,181,164,.14) !important;
}
</style>
</head>
<body>
{%app_entry%}
<footer>
{%config%}
{%scripts%}
{%renderer%}
</footer>
</body>
</html>
"""


def build_app(entry_arg: str):
    global ENTRY
    ENTRY = resolve_entry(entry_arg)
    warning, loaded_layers = reload_state()

    app = Dash(__name__, suppress_callback_exceptions=True)
    app.index_string = SSE_INDEX_STRING
    app.title = f"SSE — {ENTRY.stem}"

    coord_system_options = [{"label": k, "value": k} for k in _STATE.coord_systems]
    first_system = next(iter(_STATE.coord_systems), None)
    cont_cols = [c for c in cols_of_type("continuous") if _types.get(c) == TYPE_LABEL]

    def coord_axis_options(system):
        cols = _STATE.coord_systems.get(system, []) if system else []
        return [{"label": c, "value": c} for c in cols]

    def all_coord_opts():
        return [{"label": c, "value": c} for c in _STATE.coord_cols]

    app.layout = html.Div([
        dcc.Store(id="data-loaded-store", data=True),
        dcc.Store(id="theme-store", data="pipeline"),
        dcc.Store(id="reload-counter", data=0),
        dcc.Store(id="fixed-color-store", data=DEFAULT_FIXED_COLOR),
        dcc.Store(id="layers-store", data=loaded_layers),
        dcc.Store(id="load-trigger-store", data=None),
        dcc.Store(id="wf-visible-store", data=True),
        dcc.Store(id="coord-cols-store", data={"x": _x_col, "y": _y_col}),
        dcc.Store(id="view-init-store", data=axis_range(_ann_df, _x_col, _y_col)),
        dcc.Store(id="selection-store", data=[]),
        dcc.Store(id="selection-color-store", data=DEFAULT_SELECTION_COLOR),
        dcc.Store(id="filter-pending-store", data=None),
        dcc.Store(id="wheel-zoom-init-store", data=None),
        dcc.Store(id="left-collapsed-store", data=False),
        dcc.Store(id="right-collapsed-store", data=False),
        dcc.Store(id="plot-theme-store", data=None),
        dcc.Download(id="extract-download"),
        dcc.Download(id="figure-download"),

        html.Div([html.Div([html.H2("Sequence Space Explorer", style={"margin": "0", "color": "var(--text)"}), html.Span(id="subtitle-text", children=f"Entry: {ENTRY.stem} · {_ann_df.shape[0]:,} rows · {len(_STATE.coord_cols)} coordinate column(s)", style={"color": "var(--text-muted)", "fontSize": "13px"}), html.Span(id="reload-status", children=warning or "", style={"color": "var(--warning)", "fontSize": "12px", "marginLeft": "10px"})]), html.Div([dcc.Dropdown(id="theme-select", options=[{"label": "Pipeline", "value": "pipeline"}, {"label": "Clean Lab", "value": "clean-lab"}, {"label": "Dark Lab", "value": "dark-lab"}, {"label": "Rose Quartz", "value": "rose-quartz"}, {"label": "Deep Canopy", "value": "deep-canopy"}], value="pipeline", clearable=False, persistence=True, style={"width": "150px", "fontSize": "12px", "marginRight": "10px"}), html.Button("Reload datafile", id="reload-btn", n_clicks=0, className="sse-btn-sec", style={"padding": "7px 12px", "backgroundColor": "var(--surface-2)", "border": "1px solid var(--border)", "borderRadius": "7px", "cursor": "pointer", "fontSize": "12px"})], style={"display": "flex", "alignItems": "center"})], className="sse-header", style={"display": "flex", "alignItems": "center", "justifyContent": "space-between", "marginBottom": "16px", "borderBottom": "1px solid var(--border)"}),

        html.Div([
            html.Div([
                html.Details([html.Summary("Coordinates", style={"fontWeight": "bold", "cursor": "pointer", "marginBottom": "10px"}), dcc.RadioItems(id="coord-mode", options=[{"label": " Coordinate system mode", "value": "system"}, {"label": " Advanced free-axis mode", "value": "free"}], value="system", labelStyle={"display": "block", "fontSize": "12px", "marginBottom": "3px"}), html.Div(id="coord-system-panel", children=[html.Label("Coordinate system", style={"fontSize": "12px"}), dcc.Dropdown(id="coord-system-select", options=coord_system_options, value=first_system, clearable=False, style={"fontSize": "12px", "marginBottom": "4px"}), html.Label("X axis", style={"fontSize": "12px"}), dcc.Dropdown(id="x-axis-system", options=coord_axis_options(first_system), value=_x_col, clearable=False, style={"fontSize": "12px", "marginBottom": "4px"}), html.Label("Y axis", style={"fontSize": "12px"}), dcc.Dropdown(id="y-axis-system", options=coord_axis_options(first_system), value=_y_col, clearable=False, style={"fontSize": "12px"})]), html.Div(id="coord-free-panel", children=[html.Label("X axis", style={"fontSize": "12px"}), dcc.Dropdown(id="x-axis-free", options=all_coord_opts(), value=_x_col, clearable=False, style={"fontSize": "12px", "marginBottom": "4px"}), html.Label("Y axis", style={"fontSize": "12px"}), dcc.Dropdown(id="y-axis-free", options=all_coord_opts(), value=_y_col, clearable=False, style={"fontSize": "12px"})], style={"display": "none"}), html.Div(id="coord-warning", style={"fontSize": "11px", "color": "var(--warning)", "marginTop": "6px"})], open=True, style={"marginBottom": "16px"}),

                html.Details([html.Summary("Appearance", style={"fontWeight": "bold", "cursor": "pointer", "marginBottom": "10px"}), html.P("Filtered points", style={**SECTION_STYLE, "margin": "0 0 4px 0"}), html.Label("Opacity", style={"fontSize": "12px"}), dcc.Slider(id="alpha-slider", min=0.05, max=1.0, step=0.05, value=DEFAULT_ALPHA, marks={0.05: "0.05", 0.5: "0.5", 1.0: "1"}), html.Label("Size", style={"fontSize": "12px"}), dcc.Slider(id="point-size-slider", min=2, max=20, step=1, value=DEFAULT_POINT_SIZE, marks={2: "2", 6: "6", 12: "12", 20: "20"}), html.Label("Shape", style={"fontSize": "12px"}), dcc.Dropdown(id="marker-symbol", options=MARKER_SYMBOL_OPTIONS, value=DEFAULT_SYMBOL, clearable=False, style={"fontSize": "12px", "marginBottom": "6px"}), html.P("Background points", style={**SECTION_STYLE, "margin": "14px 0 4px 0"}), html.Label("Size", style={"fontSize": "12px"}), dcc.Slider(id="bg-size-slider", min=1, max=12, step=1, value=DEFAULT_BG_SIZE, marks={1:"1",4:"4",8:"8",12:"12"}), html.P("Query markers", style={**SECTION_STYLE, "margin": "14px 0 4px 0"}), html.Label("Opacity", style={"fontSize":"12px"}), dcc.Slider(id="marker-alpha-slider", min=0.05, max=1.0, step=0.05, value=DEFAULT_MARKER_ALPHA, marks={0.05:"0.05",0.5:"0.5",1.0:"1"}), html.Label("Size", style={"fontSize":"12px"}), dcc.Slider(id="marker-size-slider", min=6, max=40, step=1, value=DEFAULT_MARKER_SIZE, marks={6:"6",14:"14",28:"28",40:"40"}), html.Label("Position", style={"fontSize":"12px"}), dcc.RadioItems(id="marker-mode", options=[{"label":" On top","value":"top"},{"label":" Below overlays","value":"bottom"},{"label":" Hide","value":"none"}], value=DEFAULT_MARKER_MODE, labelStyle={"display":"block","fontSize":"13px"}), html.P("Working filter position", style={**SECTION_STYLE, "margin": "14px 0 4px 0"}), dcc.RadioItems(id="wf-position", options=[{"label":" On top of saved layers","value":"top"},{"label":" Below saved layers","value":"bottom"}], value=DEFAULT_WF_POSITION, labelStyle={"display":"block","fontSize":"13px"})], open=False, style={"marginBottom": "16px"}),

                html.Details([html.Summary("Colour", style={"fontWeight": "bold", "cursor": "pointer", "marginBottom": "10px"}), dcc.RadioItems(id="color-mode", options=[{"label":" Fixed color","value":"fixed"},{"label":" Continuous","value":"continuous"}], value=DEFAULT_COLOR_MODE, labelStyle={"display":"block","fontSize":"13px"}), html.Div(id="fixed-color-panel", children=[html.Label("Pick a color", style={"fontSize":"12px"}), html.Div([html.Div(id={"type":"color-chip","color":c}, style={"width":"22px","height":"22px","backgroundColor":c,"borderRadius":"3px","cursor":"pointer","border":"2px solid transparent","display":"inline-block","marginRight":"4px","marginBottom":"4px"}) for c in FIXED_COLOR_OPTIONS])]), html.Div(id="continuous-color-panel", children=[html.Label("Color by", style={"fontSize":"12px"}), dcc.Dropdown(id="cont-color-col", options=[{"label":c,"value":c} for c in cont_cols], value=cont_cols[0] if cont_cols else None, clearable=False, style={"fontSize":"13px","marginBottom":"6px"}), html.Label("Colormap", style={"fontSize":"12px"}), dcc.Dropdown(id="colormap-select", options=[{"label":c,"value":c} for c in COLORMAP_OPTIONS], value=DEFAULT_COLORMAP, clearable=False, style={"fontSize":"13px","marginBottom":"6px"}), dcc.Checklist(id="colormap-reversed", options=[{"label":html.Span(" Reverse colormap", style={"fontSize":"12px"}),"value":"reversed"}], value=[]), html.Label("Color range", style={"fontSize":"12px"}), dcc.RadioItems(id="color-range-mode", options=[{"label":" Global","value":"global"},{"label":" Subset","value":"subset"}], value=DEFAULT_COLOR_RANGE, labelStyle={"display":"block","fontSize":"12px"})], style={"display":"none"})], open=False, style={"marginBottom":"16px"}),

                html.Details([html.Summary("Filters", style={"fontWeight":"bold","cursor":"pointer","marginBottom":"10px"}), html.Div(id="filter-panel", children=make_filter_panel())], open=True, style={"marginBottom":"16px"}),
                html.Details([html.Summary("Search by ID", style={"fontWeight":"bold","cursor":"pointer","marginBottom":"10px"}), dcc.Checklist(id="id-search-enabled", options=[{"label":html.Span(" Enable ID/name search", style=LABEL_STYLE),"value":"on"}], value=[]), html.Div(id="id-search-control", children=[dcc.Textarea(id="id-search-input", placeholder="IDs separated by commas", style={"width":"100%","fontSize":"12px","minHeight":"60px","fontFamily":"monospace"}), html.Div(id="id-search-status", style={"fontSize":"11px","marginTop":"4px","color":"var(--text-muted)"})], style={"display":"none", **CONTROL_WRAPPER})], open=False, style={"marginBottom":"16px"}),
                html.Details([html.Summary("Save layer", style={"fontWeight":"bold","cursor":"pointer","marginBottom":"10px"}), dcc.Input(id="layer-name-input", type="text", placeholder="Auto-generated if blank", style={"width":"100%","fontSize":"12px","marginBottom":"8px"}), html.Button("Save layer", id="save-layer-btn", n_clicks=0, className="sse-btn-primary", style={"width":"100%","padding":"8px","backgroundColor":"var(--accent)","color":"var(--on-accent)","border":"none","borderRadius":"7px","cursor":"pointer","fontSize":"13px","fontWeight":"600"}), html.Div(id="save-layer-status", style={"fontSize":"11px","marginTop":"6px","color":"var(--text-muted)"})], open=False, style={"marginBottom":"16px"}),
                html.Details([html.Summary("Cluster regions", style={"fontWeight":"bold","cursor":"pointer","marginBottom":"10px"}),
                    html.P("Show clusters as shaded regions or coloured points. Off until a column is picked.", style={"fontSize":"11px","color":"var(--text-muted)","margin":"0 0 8px 0"}),
                    html.Label("Cluster column", style={"fontSize":"12px"}),
                    dcc.Dropdown(id="cluster-region-col", options=[{"label":c,"value":c} for c in _types if _types.get(c)==TYPE_LABEL and c.endswith("_cluster")], value=None, placeholder="None (off)", clearable=True, style={"fontSize":"12px","marginBottom":"6px"}),
                    html.Label("Style", style={"fontSize":"12px"}),
                    dcc.RadioItems(id="cluster-region-shape", options=[{"label":" Density region (KDE)","value":"kde"},{"label":" Concave hull region","value":"hull"},{"label":" Colour points by cluster","value":"points"}], value="kde", labelStyle={"display":"block","fontSize":"13px"}),
                    html.Label("Position", style={"fontSize":"12px"}),
                    dcc.RadioItems(id="cluster-region-position", options=[{"label":" Below points","value":"below"},{"label":" On top","value":"above"}], value="below", labelStyle={"display":"block","fontSize":"13px"}),
                    html.Label("Fill opacity", style={"fontSize":"12px"}),
                    dcc.Slider(id="cluster-region-opacity", min=0.05, max=0.7, step=0.05, value=0.25, marks={0.05:"0.05",0.35:"0.35",0.7:"0.7"}),
                    html.Label("Density coverage (KDE)", style={"fontSize":"12px"}),
                    dcc.Slider(id="cluster-region-coverage", min=0.1, max=0.9, step=0.05, value=0.4, marks={0.1:"0.1",0.4:"0.4",0.9:"0.9"}),
                    html.Label("Hull tightness", style={"fontSize":"12px"}),
                    dcc.Slider(id="cluster-region-tightness", min=0.0, max=0.8, step=0.05, value=0.25, marks={0:"0",0.25:"0.25",0.8:"0.8"}),
                ], open=False, style={"marginBottom":"16px"}),
                html.Details([html.Summary("Column settings", style={"fontWeight":"bold","cursor":"pointer","marginBottom":"10px"}), html.Button("Rebuild filter panel", id="rebuild-filters-btn", n_clicks=0, style={"width":"100%","padding":"5px","fontSize":"12px","marginBottom":"8px"}), html.Div(id="col-settings-panel", children=make_col_settings_panel())], open=False, style={"marginBottom":"16px"}),
            ], id="left-rail", className="sse-col sse-rail", style={"width": "290px", "minWidth": "270px", "flexShrink": "0", "overflowY": "auto", "maxHeight": "90vh", "paddingRight": "12px"}),

            html.Div([html.Div([html.Button("‹", id="left-collapse-btn", n_clicks=0, className="rail-toggle", title="Collapse controls panel"), html.Span(id="point-count", className="plot-status", style={"fontSize":"12px","color":"var(--text-muted)","marginLeft":"2px"}), html.Button(html.Span(className="ic ic-eye"), id="wf-toggle-btn", n_clicks=0, title="Show/hide working filter", className="sse-icon-btn", style={"background":"none","border":"none","cursor":"pointer","color":"var(--text-muted)","padding":"4px 6px","marginLeft":"6px"}), html.Span(style={"flexGrow":"1"}), html.Button("›", id="right-collapse-btn", n_clicks=0, className="rail-toggle", title="Collapse layers panel")], className="sse-toolbar", style={"display":"flex","alignItems":"center","marginBottom":"6px"}), html.Div([html.Span(id="selection-count", children="No sequences selected", style={"fontSize":"12px","color":"var(--text-muted)","flexGrow":"1"}), html.Div([html.Div(style={"width":"18px","height":"18px","backgroundColor":c,"borderRadius":"50%","cursor":"pointer","display":"inline-block","marginRight":"4px","border":"2px solid var(--border)"}, id={"type":"sel-color-chip","color":c}, title=c) for c in SELECTION_COLOR_OPTIONS], style={"display":"inline-flex","alignItems":"center","marginRight":"6px"}), html.Button("Clear", id="clear-selection-btn", n_clicks=0, className="sse-btn-sec", style={"fontSize":"11px","padding":"3px 8px","marginRight":"4px","backgroundColor":"var(--surface-2)","color":"var(--text)","border":"1px solid var(--border)","borderRadius":"6px","cursor":"pointer"}), html.Button("Use as working filter", id="selection-to-wl-btn", n_clicks=0, className="sse-btn-sec", style={"fontSize":"11px","padding":"3px 8px","marginRight":"4px","backgroundColor":"var(--surface-2)","color":"var(--text)","border":"1px solid var(--border)","borderRadius":"6px","cursor":"pointer"}), html.Button("Export selection for Boltz", id="selection-export-btn", n_clicks=0, title="Save the selected sequences to a cache the pipeline's Boltz-2 module can import", className="sse-btn-sec", style={"fontSize":"11px","padding":"3px 8px","backgroundColor":"var(--surface-2)","color":"var(--text)","border":"1px solid var(--border)","borderRadius":"6px","cursor":"pointer"})], className="selection-toolbar", style={"display":"flex","alignItems":"center","marginBottom":"8px","backgroundColor":"var(--surface)","border":"1px solid var(--border)"}), html.Div(id="selection-export-status", style={"fontSize":"11px","color":"var(--text-muted)","marginBottom":"8px","minHeight":"14px"}), html.Div(dcc.Graph(id="latent-graph", style={"height":"72vh", "minHeight":"540px"}, config={"displaylogo":False,"scrollZoom":False,"responsive":True,"modeBarButtonsToAdd":["lasso2d","select2d"],"toImageButtonOptions":{"format":"png","filename":"sequence-space"}}), className="sse-plot-card"), html.Div(id="click-details", children=html.Div("Select a point to inspect metadata and analysis context.", className="details-empty"), className="details-card", style={"marginTop":"16px","backgroundColor":"var(--surface)","border":"1px solid var(--border)","borderRadius":"10px","fontSize":"13px","minHeight":"60px"}), html.Div(id="load-warning", style={"marginTop":"8px","fontSize":"11px","color":"var(--warning)"})], className="sse-center", style={"flexGrow":"1","minWidth":"0","paddingLeft":"16px","paddingRight":"16px"}),

            html.Div([html.Div([html.Span("Saved layers", style={"fontWeight":"bold","fontSize":"14px","color":"var(--text)"}), html.Button("Clear all", id="clear-layers-btn", n_clicks=0, style={"background":"none","border":"none","color":"var(--danger)","cursor":"pointer","fontSize":"11px","float":"right"})], style={"marginBottom":"8px","borderBottom":"1px solid var(--border)","paddingBottom":"6px"}), html.Div(id="layers-sidebar", children=make_sidebar(loaded_layers)), html.Button("Extract visible layers", id="extract-btn", n_clicks=0, className="sse-btn-primary", style={"width":"100%","padding":"8px","backgroundColor":"var(--accent)","color":"var(--on-accent)","border":"none","borderRadius":"7px","cursor":"pointer","fontSize":"12px","fontWeight":"600","marginBottom":"6px","marginTop":"8px"}), html.Div(id="extract-status", style={"fontSize":"11px","color":"var(--text-muted)","marginBottom":"8px"}), html.Details([
                    html.Summary("Export figure", style={"fontWeight":"600","fontSize":"12px","cursor":"pointer","marginBottom":"6px"}),
                    html.Label("Format", style={"fontSize":"11px"}),
                    dcc.RadioItems(id="export-format", options=EXPORT_FORMAT_OPTIONS, value="png", labelStyle={"display":"inline-block","fontSize":"12px","marginRight":"8px"}, style={"marginBottom":"4px"}),
                    html.Label("Resolution", style={"fontSize":"11px"}),
                    dcc.RadioItems(id="export-dpi", options=EXPORT_DPI_OPTIONS, value=300, labelStyle={"display":"block","fontSize":"12px","marginBottom":"2px"}, style={"marginBottom":"4px"}),
                    dcc.Checklist(id="export-legend", options=[{"label":html.Span(" Include legend", style={"fontSize":"12px"}),"value":"show"}], value=["show"], style={"marginBottom":"4px"}),
                    dcc.Checklist(id="export-transparent", options=[{"label":html.Span(" Transparent background", style={"fontSize":"12px"}),"value":"on"}], value=[], style={"marginBottom":"6px"}),
                    html.Label("Axis colour (hex)", style={"fontSize":"11px"}),
                    dcc.Input(id="export-axis-color", type="text", value="#000000", debounce=True, style={"width":"100%","fontSize":"12px","marginBottom":"4px","boxSizing":"border-box"}),
                    html.Label("Axis label colour (hex)", style={"fontSize":"11px"}),
                    dcc.Input(id="export-label-color", type="text", value="#000000", debounce=True, style={"width":"100%","fontSize":"12px","marginBottom":"4px","boxSizing":"border-box"}),
                    html.Label("Marker edge colour (hex)", style={"fontSize":"11px"}),
                    dcc.Input(id="export-edge-color", type="text", value="#000000", debounce=True, style={"width":"100%","fontSize":"12px","marginBottom":"4px","boxSizing":"border-box"}),
                    html.Label("Background colour (hex, ignored if transparent)", style={"fontSize":"11px"}),
                    dcc.Input(id="export-bg-color", type="text", value="#ffffff", debounce=True, style={"width":"100%","fontSize":"12px","marginBottom":"4px","boxSizing":"border-box"}),
                    html.Div([html.Label("Width × height (px)", style={"fontSize":"11px", "flexGrow":"1"}), html.Button("Reset", id="export-size-reset", n_clicks=0, title="Reset to 1200 × 800 px", style={"fontSize":"10px", "padding":"1px 6px", "border":"1px solid var(--border)", "backgroundColor":"var(--surface-2)", "borderRadius":"3px", "cursor":"pointer"})], style={"display":"flex", "alignItems":"center", "gap":"6px", "marginBottom":"2px"}),
                    html.Div([
                        dcc.Input(id="export-width", type="number", value=EXPORT_DEFAULT_WIDTH, min=300, max=8000, step=50, debounce=True, style={"width":"48%","fontSize":"12px","marginRight":"4%","boxSizing":"border-box"}),
                        dcc.Input(id="export-height", type="number", value=EXPORT_DEFAULT_HEIGHT, min=300, max=8000, step=50, debounce=True, style={"width":"48%","fontSize":"12px","boxSizing":"border-box"}),
                    ], style={"display":"flex","marginBottom":"6px"}),
                    html.Label("Save to", style={"fontSize":"11px"}),
                    dcc.RadioItems(id="export-destination", options=[{"label":html.Span(" Browser download", style={"fontSize":"12px"}),"value":"browser"},{"label":html.Span(" Entry figures/", style={"fontSize":"12px"}),"value":"server"}], value="browser", labelStyle={"display":"block","marginBottom":"2px"}, inputStyle={"cursor":"pointer","marginRight":"4px"}, style={"marginBottom":"6px","marginTop":"2px"}),
                    html.Button("Save figure", id="export-btn", n_clicks=0, className="sse-btn-primary", style={"width":"100%","padding":"8px","backgroundColor":"var(--accent)","color":"var(--on-accent)","border":"none","borderRadius":"7px","cursor":"pointer","fontSize":"12px","fontWeight":"600","marginBottom":"4px"}),
                    html.Div(id="export-status", style={"fontSize":"11px","color":"var(--text-muted)","wordBreak":"break-all"})
                ], open=False, style={"marginBottom":"12px"}),

                html.Details([
                    html.Summary("Structure prediction & RMSD", style={"fontWeight":"600","fontSize":"12px","cursor":"pointer","marginBottom":"8px"}),
                    html.Div([
                        html.P("Boltz-2 structure prediction and RMSD analysis now run in the pipeline app, not here.", style={"fontSize":"11px","color":"var(--text-muted)","margin":"0 0 6px 0"}),
                        html.P("Select the points you want, click “Export selection for Boltz” above, then open the pipeline's “Structure & binding” module to import that selection and run prediction. New pTM / pLDDT and RMSD columns appear here after you reload the datafile.", style={"fontSize":"11px","color":"var(--text-muted)","margin":"0"}),
                    ]),
                ], open=False, style={"marginBottom":"12px"})
            ], id="right-rail", className="sse-col sse-rail", style={"width":"235px","minWidth":"215px","flexShrink":"0","paddingLeft":"12px","overflowY":"auto","maxHeight":"90vh"}),
        ], className="sse-workspace", style={"display":"flex","gap":"0","alignItems":"flex-start"})
    ], id="app-root", className="sse-root", style={"padding":"20px","fontFamily":"sans-serif","maxWidth":"1900px","margin":"0 auto"})

    app.clientside_callback(
        "function(v){ v = v || 'pipeline'; document.documentElement.setAttribute('data-theme', v); return v; }",
        Output("theme-store", "data"),
        Input("theme-select", "value"),
    )

    # Collapsible side rails. Toggling the rail's own style dict from the server
    # risks clobbering it around dynamic panel rebuilds, so we mutate the DOM node's
    # style directly on the client and only round-trip the boolean state + the
    # toggle glyph. The plot column has flexGrow:1, so it reclaims the freed width.
    for side, rail_id, btn_id, store_id, expanded_w, expanded_mw, pad_side, glyph_open, glyph_closed in [
        ("left", "left-rail", "left-collapse-btn", "left-collapsed-store", "290px", "270px", "paddingRight", "‹", "›"),
        ("right", "right-rail", "right-collapse-btn", "right-collapsed-store", "235px", "215px", "paddingLeft", "›", "‹"),
    ]:
        app.clientside_callback(
            f"""
            function(n, collapsed){{
                const now = !collapsed;
                const el = document.getElementById('{rail_id}');
                if (el){{
                    if (now){{
                        el.style.width='0px'; el.style.minWidth='0px';
                        el.style.{pad_side}='0px'; el.style.overflow='hidden'; el.style.opacity='0';
                    }} else {{
                        el.style.width='{expanded_w}'; el.style.minWidth='{expanded_mw}';
                        el.style.{pad_side}='12px'; el.style.overflowY='auto'; el.style.opacity='1';
                    }}
                }}
                return [now, now ? '{glyph_closed}' : '{glyph_open}'];
            }}
            """,
            [Output(store_id, "data"), Output(btn_id, "children")],
            Input(btn_id, "n_clicks"),
            State(store_id, "data"),
            prevent_initial_call=True,
        )

    # Custom throttled scroll-zoom for the latent graph. Plotly's native scrollZoom
    # rubberbands on fast wheel input when a scaleanchor (equal-aspect) constraint is
    # set, because each wheel tick triggers its own constraint-solving relayout and
    # consecutive ticks solve to inconsistent ranges. We disable native scrollZoom and
    # coalesce wheel events into a single relayout per animation frame, scaling both
    # axes by the same factor so the 1:1 aspect ratio stays consistent (no oscillation).
    app.clientside_callback(
        """
        function(fig) {
            const nu = window.dash_clientside.no_update;
            const container = document.getElementById('latent-graph');
            if (!container) return nu;
            if (container._customWheelZoom) return nu;
            container._customWheelZoom = true;

            let pending = 0;      // accumulated zoom exponent (log-factor)
            let scheduled = false;
            let lastEvent = null;

            function applyZoom() {
                scheduled = false;
                const gd = container.classList.contains('js-plotly-plot')
                    ? container : container.querySelector('.js-plotly-plot');
                if (!gd || !gd._fullLayout) { pending = 0; return; }
                const drag = gd.querySelector('.nsewdrag');
                const xa = gd._fullLayout.xaxis, ya = gd._fullLayout.yaxis;
                if (!drag || !xa || !ya || !xa.range || !ya.range) { pending = 0; return; }
                const rect = drag.getBoundingClientRect();
                if (!rect.width || !rect.height) { pending = 0; return; }
                const x0 = xa.range[0], x1 = xa.range[1];
                const y0 = ya.range[0], y1 = ya.range[1];
                const fx = (lastEvent.clientX - rect.left) / rect.width;
                const fy = (lastEvent.clientY - rect.top) / rect.height;
                const cx = x0 + fx * (x1 - x0);          // cursor position in data coords
                const cy = y1 - fy * (y1 - y0);          // (pixel-top corresponds to y-max)
                const f = Math.exp(Math.max(-0.5, Math.min(0.5, pending)));  // clamp per-frame
                pending = 0;
                window.Plotly.relayout(gd, {
                    'xaxis.range': [cx - (cx - x0) * f, cx + (x1 - cx) * f],
                    'yaxis.range': [cy - (cy - y0) * f, cy + (y1 - cy) * f]
                });
            }

            container.addEventListener('wheel', function(e) {
                e.preventDefault();
                e.stopPropagation();
                lastEvent = e;
                let dy = e.deltaY;
                if (e.deltaMode === 1) dy *= 16;               // line units -> approx px
                else if (e.deltaMode === 2) dy *= window.innerHeight;
                pending += dy * 0.0015;   // deltaY>0 (scroll down) -> f>1 -> zoom out
                if (!scheduled) {
                    scheduled = true;
                    window.requestAnimationFrame(applyZoom);
                }
            }, { passive: false, capture: true });

            return nu;
        }
        """,
        Output("wheel-zoom-init-store", "data"),
        Input("latent-graph", "figure"),
    )

    # Theme-aware plot chrome. make_figure renders with template="simple_white"
    # (white paper, dark axes) because it has no idea which theme is active — the
    # theme is a client-only CSS-variable swap. So after every figure render, and on
    # every theme change, restyle the paper/axes/legend from the live CSS variables:
    # this keeps the CSS the single source of truth (exact colour match in all four
    # themes), updates instantly on switch, and avoids re-rendering the ~10k-point
    # figure server-side. Paper/plot backgrounds go transparent so the surface of the
    # .sse-plot-card shows through. Data-trace colours are intentionally left alone.
    app.clientside_callback(
        """
        function(theme, fig) {
            const nu = window.dash_clientside.no_update;
            let tries = 0;
            function apply() {
                const gd = document.getElementById('latent-graph');
                if (!gd) return;
                const plot = gd.classList.contains('js-plotly-plot') ? gd : gd.querySelector('.js-plotly-plot');
                if (!plot || !window.Plotly || !plot._fullLayout) {
                    if (tries++ < 30) window.requestAnimationFrame(apply);
                    return;
                }
                const cs = getComputedStyle(document.documentElement);
                const v = n => cs.getPropertyValue(n).trim();
                const text = v('--text'), muted = v('--text-muted'), border = v('--border'), surface = v('--surface');
                // Pipeline theme: render the plot as a white "figure sheet" so the
                // data stays readable against the dark card. Other themes keep the
                // transparent-paper / CSS-variable treatment.
                const white = theme === 'pipeline';
                const paper = white ? '#ffffff' : 'rgba(0,0,0,0)';
                const axisTxt = white ? '#3b4a55' : muted;
                const titleTxt = white ? '#1f2733' : text;
                const line = white ? '#cbd8d8' : border;
                const legBg = white ? '#ffffff' : surface;
                window.Plotly.relayout(plot, {
                    'paper_bgcolor': paper,
                    'plot_bgcolor': paper,
                    'font.color': axisTxt,
                    'xaxis.color': axisTxt,
                    'yaxis.color': axisTxt,
                    'xaxis.linecolor': line,
                    'yaxis.linecolor': line,
                    'xaxis.zerolinecolor': line,
                    'yaxis.zerolinecolor': line,
                    'xaxis.title.font.color': titleTxt,
                    'yaxis.title.font.color': titleTxt,
                    'legend.bgcolor': legBg,
                    'legend.bordercolor': line,
                    'legend.borderwidth': 1,
                    'legend.font.color': titleTxt
                });
                // The pipeline theme drops the plot-card padding, which resizes the
                // container after Plotly's initial draw; relayout does not remeasure,
                // so force a resize to make the plot fill (and round to) the card.
                if (window.Plotly.Plots && window.Plotly.Plots.resize) window.Plotly.Plots.resize(plot);
            }
            apply();
            return nu;
        }
        """,
        Output("plot-theme-store", "data"),
        Input("theme-store", "data"),
        Input("latent-graph", "figure"),
    )

    # ---------------- callbacks ----------------
    @app.callback(Output("fixed-color-panel", "style"), Output("continuous-color-panel", "style"), Input("color-mode", "value"))
    def toggle_color_panels(mode):
        return ({"display":"block"}, {"display":"none"}) if mode == "fixed" else ({"display":"none"}, {"display":"block"})

    @app.callback(Output("coord-system-panel", "style"), Output("coord-free-panel", "style"), Input("coord-mode", "value"))
    def toggle_coord_mode(mode):
        return ({"display":"block"}, {"display":"none"}) if mode == "system" else ({"display":"none"}, {"display":"block"})

    @app.callback(Output("x-axis-system", "options"), Output("x-axis-system", "value"), Output("y-axis-system", "options"), Output("y-axis-system", "value"), Input("coord-system-select", "value"))
    def update_system_axes(system):
        cols = _STATE.coord_systems.get(system, [])
        opts = [{"label": c, "value": c} for c in cols]
        x = cols[0] if cols else None
        y = cols[1] if len(cols) > 1 else x
        return opts, x, opts, y

    @app.callback(Output("coord-cols-store", "data"), Output("coord-warning", "children"), Output("view-init-store", "data", allow_duplicate=True), Input("coord-mode", "value"), Input("x-axis-system", "value"), Input("y-axis-system", "value"), Input("x-axis-free", "value"), Input("y-axis-free", "value"), prevent_initial_call=True)
    def set_axes(mode, xs, ys, xf, yf):
        global _x_col, _y_col
        x, y = (xf, yf) if mode == "free" else (xs, ys)
        _x_col, _y_col = x, y
        warn = ""
        if mode == "free" and x and y and coordinate_system_key(x) != coordinate_system_key(y):
            warn = "Advanced axis pairing: X and Y come from different coordinate systems."
        return {"x": x, "y": y}, warn, axis_range(_ann_df, x, y)

    @app.callback(Output("filter-panel", "children", allow_duplicate=True), Output("col-settings-panel", "children", allow_duplicate=True), Output("subtitle-text", "children"), Output("reload-status", "children"), Output("reload-counter", "data"), Output("layers-store", "data", allow_duplicate=True), Output("filter-pending-store", "data", allow_duplicate=True), Input("reload-btn", "n_clicks"), State("reload-counter", "data"), State("layers-store", "data"), prevent_initial_call=True)
    def manual_reload(n, counter, current_layers):
        warn, cleaned_layers = reload_state()
        with _STATE_LOCK:
            subtitle = f"Entry: {ENTRY.stem} · {_ann_df.shape[0]:,} rows · {len(_STATE.coord_cols)} coordinate column(s)"
        status = " · ".join(x for x in [warn, "Reloaded datafile."] if x)
        # Clear dynamic/pattern-matched panels first. A second callback rebuilds
        # them after Dash has removed the old ALL-pattern components from the DOM.
        # Updating them in one callback can trigger Dash renderer JS errors
        # ("Cannot read properties of undefined (reading 'props')").
        pending = {"source": "reload", "counter": (counter or 0) + 1}
        return html.Div(), html.Div(), subtitle, status, (counter or 0) + 1, cleaned_layers, pending

    @app.callback(Output("filter-panel", "children", allow_duplicate=True), Output("col-settings-panel", "children", allow_duplicate=True), Output("cont-color-col", "options", allow_duplicate=True), Output("cont-color-col", "value", allow_duplicate=True), Input("filter-pending-store", "data"), prevent_initial_call=True)
    def rebuild_dynamic_panels(pending):
        if not pending:
            return no_update, no_update, no_update, no_update
        with _STATE_LOCK:
            cont = [c for c in cols_of_type("continuous") if _types.get(c) == TYPE_LABEL]
            return make_filter_panel(), make_col_settings_panel(), [{"label": c, "value": c} for c in cont], cont[0] if cont else None

    @app.callback(Output("fixed-color-store", "data"), [Input({"type":"color-chip", "color":c}, "n_clicks") for c in FIXED_COLOR_OPTIONS], prevent_initial_call=True)
    def pick_fixed_color(*_):
        triggered = ctx.triggered_id
        if isinstance(triggered, dict):
            return triggered.get("color", no_update)
        return no_update

    @app.callback(Output({"type":"filter-control", "col":ALL}, "style"), Input({"type":"filter-enabled", "col":ALL}, "value"))
    def toggle_filter_controls(enabled_values):
        return [{**CONTROL_WRAPPER, "display":"block"} if v and "on" in v else {**CONTROL_WRAPPER, "display":"none"} for v in enabled_values]

    @app.callback(
        Output({"type":"cont-min-input", "col":ALL}, "value"),
        Output({"type":"cont-max-input", "col":ALL}, "value"),
        Input({"type":"cont-slider", "col":ALL}, "value"),
        prevent_initial_call=True,
    )
    def sync_cont_slider_to_inputs(values):
        """Keep the numeric range boxes in sync when the slider handles move."""
        mins, maxs = [], []
        for v in values or []:
            if isinstance(v, (list, tuple)) and len(v) >= 2:
                mins.append(v[0])
                maxs.append(v[1])
            else:
                mins.append(no_update)
                maxs.append(no_update)
        return mins, maxs

    @app.callback(
        Output({"type":"cont-slider", "col":ALL}, "value", allow_duplicate=True),
        Input({"type":"cont-min-input", "col":ALL}, "value"),
        Input({"type":"cont-max-input", "col":ALL}, "value"),
        State({"type":"cont-slider", "col":ALL}, "value"),
        State({"type":"cont-slider", "col":ALL}, "min"),
        State({"type":"cont-slider", "col":ALL}, "max"),
        prevent_initial_call=True,
    )
    def sync_cont_inputs_to_slider(min_values, max_values, current_values, slider_mins, slider_maxs):
        """Typing in range boxes updates the real slider value used by filters."""
        out = []
        for mn, mx, cur, smin, smax in zip(min_values or [], max_values or [], current_values or [], slider_mins or [], slider_maxs or []):
            cur = cur if isinstance(cur, (list, tuple)) and len(cur) >= 2 else [smin, smax]
            lo = _parse_numeric_value(mn, cur[0])
            hi = _parse_numeric_value(mx, cur[1])
            lo = max(smin, min(lo, smax)) if lo is not None else cur[0]
            hi = max(smin, min(hi, smax)) if hi is not None else cur[1]
            if lo > hi:
                lo, hi = hi, lo
            out.append([lo, hi])
        return out

    @app.callback(Output("id-search-control", "style"), Input("id-search-enabled", "value"))
    def toggle_id_search(enabled):
        return {"display":"block", **CONTROL_WRAPPER} if enabled and "on" in enabled else {"display":"none", **CONTROL_WRAPPER}

    @app.callback(Output("id-search-status", "children"), Input("id-search-enabled", "value"), Input("id-search-input", "n_blur"), State("id-search-input", "value"))
    def id_search_status(enabled, _n, raw):
        if not enabled or "on" not in enabled or not raw:
            return ""
        queried = [x.strip() for x in raw.split(",") if x.strip()]
        id_match = _ann_df[_id_col].astype(str).isin(queried)
        name_match = pd.Series(False, index=_ann_df.index)
        for nc in _name_cols:
            name_match |= _ann_df[nc].astype(str).isin(queried)
        return f"{int((id_match | name_match).sum())}/{len(queried)} matched across {', '.join([_id_col] + _name_cols)}"

    @app.callback(Output("wf-visible-store", "data"), Output("wf-toggle-btn", "children"), Input("wf-toggle-btn", "n_clicks"), State("wf-visible-store", "data"), prevent_initial_call=True)
    def toggle_wf_visibility(_, currently_visible):
        now = not currently_visible
        return now, html.Span(className="ic ic-eye" if now else "ic ic-eyeoff")

    @app.callback(Output("layers-store", "data"), Output("save-layer-status", "children"), Output("layer-name-input", "value"), Input("save-layer-btn", "n_clicks"), State("layer-name-input", "value"), State("color-mode", "value"), State("fixed-color-store", "data"), State("cont-color-col", "value"), State("colormap-select", "value"), State("colormap-reversed", "value"), State("color-range-mode", "value"), State("alpha-slider", "value"), State("point-size-slider", "value"), State("marker-symbol", "value"), State({"type":"filter-enabled", "col":ALL}, "value"), State({"type":"filter-enabled", "col":ALL}, "id"), State({"type":"cont-slider", "col":ALL}, "value"), State({"type":"cont-slider", "col":ALL}, "id"), State({"type":"bool-filter", "col":ALL}, "value"), State({"type":"bool-filter", "col":ALL}, "id"), State({"type":"cat-filter", "col":ALL}, "value"), State({"type":"cat-filter", "col":ALL}, "id"), State({"type":"tag-filter", "col":ALL}, "value"), State({"type":"tag-filter", "col":ALL}, "id"), State("id-search-enabled", "value"), State("id-search-input", "value"), State("layers-store", "data"), prevent_initial_call=True)
    def save_layer_cb(n_clicks, custom_name, color_mode, fixed_color, cont_col, colormap, colormap_reversed, color_range_mode, alpha, point_size, marker_symbol, enabled_values, enabled_ids, cont_values, cont_ids, bool_values, bool_ids, cat_values, cat_ids, tag_values, tag_ids, id_enabled, id_raw, existing_layers):
        if not n_clicks or _ann_df.empty:
            return no_update, no_update, no_update
        cont_conds, bool_conds, cat_conds, tag_conds, id_search_ids = parse_current_conditions(enabled_values, enabled_ids, cont_values, cont_ids, bool_values, bool_ids, cat_values, cat_ids, tag_values, tag_ids, id_enabled, id_raw)
        reversed_ = bool(colormap_reversed and "reversed" in colormap_reversed)
        name = (custom_name or "").strip() or auto_layer_name(cont_conds, bool_conds, cat_conds, tag_conds, color_mode, fixed_color, cont_col, colormap, reversed_, id_search_ids)
        fs = {"cont": [[c, lo, hi] for c, lo, hi in cont_conds], "bool": [[c, v] for c, v in bool_conds], "cat": [[c, vs] for c, vs in cat_conds], "tag": [[c, ts] for c, ts in tag_conds], "id_search_ids": id_search_ids or [], "color_mode": color_mode, "fixed_color": fixed_color, "cont_col": cont_col, "colormap": colormap, "reversed": reversed_, "color_range": color_range_mode, "alpha": alpha, "point_size": point_size, "marker_symbol": marker_symbol or DEFAULT_SYMBOL}
        layer = make_layer(cont_conds, bool_conds, cat_conds, tag_conds, color_mode, fixed_color, cont_col, colormap, reversed_, color_range_mode, name, fs, alpha, point_size, marker_symbol, id_search_ids)

        # Use the persisted file as a fallback when the browser Store state is
        # stale/empty. Without this, a new save can overwrite layers.json with
        # only the newly created layer, making it look as if previous layers
        # were replaced.
        persisted_layers = layer_store.read_layers(ENTRY.layers_path)
        # File state is canonical here. Browser Store can be stale around dynamic
        # panel rebuilds, so prefer layers.json and use Store only as fallback.
        base_layers = _merge_layer_lists(persisted_layers, existing_layers or [])
        updated = [layer] + base_layers
        layer_store.write_layers(ENTRY.layers_path, updated)
        return updated, f"✓ Saved '{name}' ({layer['n_points']:,} sequences)", ""

    @app.callback(Output("layers-store", "data", allow_duplicate=True), Input({"type":"layer-up", "lid":ALL}, "n_clicks"), Input({"type":"layer-down", "lid":ALL}, "n_clicks"), State({"type":"layer-up", "lid":ALL}, "id"), State({"type":"layer-down", "lid":ALL}, "id"), State("layers-store", "data"), prevent_initial_call=True)
    def reorder_layers(up_clicks, down_clicks, up_ids, down_ids, layers):
        triggered = ctx.triggered_id
        clicks = up_clicks if isinstance(triggered, dict) and triggered.get("type") == "layer-up" else down_clicks
        ids = up_ids if isinstance(triggered, dict) and triggered.get("type") == "layer-up" else down_ids
        if not _pattern_click_fired(triggered, clicks, ids):
            return no_update
        base = _merge_layer_lists(layers or [], layer_store.read_layers(ENTRY.layers_path))
        if not base:
            return no_update
        lid, typ = triggered["lid"], triggered["type"]
        out = list(base)
        idx = next((i for i, l in enumerate(out) if l.get("id") == lid), None)
        if idx is None:
            return no_update
        if typ == "layer-up" and idx > 0:
            out[idx], out[idx - 1] = out[idx - 1], out[idx]
        elif typ == "layer-down" and idx < len(out) - 1:
            out[idx], out[idx + 1] = out[idx + 1], out[idx]
        else:
            return no_update
        layer_store.write_layers(ENTRY.layers_path, out)
        return out

    @app.callback(Output("layers-store", "data", allow_duplicate=True), Input({"type":"layer-toggle", "lid":ALL}, "n_clicks"), State({"type":"layer-toggle", "lid":ALL}, "id"), State("layers-store", "data"), prevent_initial_call=True)
    def toggle_layer(clicks, ids, layers):
        triggered = ctx.triggered_id
        if not _pattern_click_fired(triggered, clicks, ids):
            return no_update
        base = _merge_layer_lists(layers or [], layer_store.read_layers(ENTRY.layers_path))
        if not base:
            return no_update
        out = [{**l, "visible": not l.get("visible", True)} if l.get("id") == triggered["lid"] else l for l in base]
        layer_store.write_layers(ENTRY.layers_path, out)
        return out

    @app.callback(Output("layers-store", "data", allow_duplicate=True), Input({"type":"layer-delete", "lid":ALL}, "n_clicks"), State({"type":"layer-delete", "lid":ALL}, "id"), State("layers-store", "data"), prevent_initial_call=True)
    def delete_layer(clicks, ids, layers):
        triggered = ctx.triggered_id
        if not _pattern_click_fired(triggered, clicks, ids):
            return no_update
        base = _merge_layer_lists(layers or [], layer_store.read_layers(ENTRY.layers_path))
        if not base:
            return no_update
        out = [l for l in base if l.get("id") != triggered["lid"]]
        layer_store.write_layers(ENTRY.layers_path, out)
        return out

    @app.callback(Output("layers-store", "data", allow_duplicate=True), Input("clear-layers-btn", "n_clicks"), prevent_initial_call=True)
    def clear_layers(_):
        layer_store.write_layers(ENTRY.layers_path, [])
        return []

    @app.callback(Output("layers-sidebar", "children"), Input("layers-store", "data"), Input("coord-cols-store", "data"))
    def render_sidebar(layers, _coord):
        return make_sidebar(layers or [])

    @app.callback(Output("load-trigger-store", "data"), Input({"type":"layer-load", "lid":ALL}, "n_clicks"), State({"type":"layer-load", "lid":ALL}, "id"), State("layers-store", "data"), prevent_initial_call=True)
    def layer_load(clicks, ids, layers):
        triggered = ctx.triggered_id
        # Only a real click on the explicit "load into working filter" button
        # may restore a layer. Reordering/toggling/deleting re-renders the
        # sidebar, which inserts a fresh set of ALL-pattern layer-load buttons;
        # Dash may fire this callback during that insertion with n_clicks=None/0.
        # Without this guard, the newly topmost layer can be mistaken for a load
        # action and its filters get applied to the working layer.
        if not _pattern_click_fired(triggered, clicks, ids):
            return no_update
        if not layers:
            return no_update
        layer = next((l for l in layers if l.get("id") == triggered["lid"]), None)
        return layer.get("filter_state") if layer else no_update

    @app.callback(Output("color-mode", "value"), Output("fixed-color-store", "data", allow_duplicate=True), Output("cont-color-col", "value", allow_duplicate=True), Output("colormap-select", "value"), Output("colormap-reversed", "value"), Output("color-range-mode", "value"), Output("alpha-slider", "value"), Output("point-size-slider", "value"), Output("marker-symbol", "value", allow_duplicate=True), Output({"type":"filter-enabled", "col":ALL}, "value"), Output({"type":"cont-slider", "col":ALL}, "value"), Output({"type":"bool-filter", "col":ALL}, "value"), Output({"type":"cat-filter", "col":ALL}, "value"), Output({"type":"tag-filter", "col":ALL}, "value"), Output("id-search-enabled", "value"), Output("id-search-input", "value"), Input("load-trigger-store", "data"), State({"type":"filter-enabled", "col":ALL}, "id"), State({"type":"cont-slider", "col":ALL}, "id"), State({"type":"bool-filter", "col":ALL}, "id"), State({"type":"cat-filter", "col":ALL}, "id"), State({"type":"tag-filter", "col":ALL}, "id"), prevent_initial_call=True)
    def restore_layer(fs, enabled_ids, cont_ids, bool_ids, cat_ids, tag_ids):
        if not fs:
            return (no_update,)*9 + ([no_update]*len(enabled_ids), [no_update]*len(cont_ids), [no_update]*len(bool_ids), [no_update]*len(cat_ids), [no_update]*len(tag_ids), no_update, no_update)
        cont_lookup = {c: [lo, hi] for c, lo, hi in fs.get("cont", [])}
        bool_lookup = {c: ("true" if v else "false") for c, v in fs.get("bool", [])}
        cat_lookup = {c: vs for c, vs in fs.get("cat", [])}
        tag_lookup = {c: ts for c, ts in fs.get("tag", [])}
        enabled = [["on"] if d["col"] in cont_lookup or d["col"] in bool_lookup or d["col"] in cat_lookup or d["col"] in tag_lookup else [] for d in enabled_ids]

        def _cont_default(col):
            vmin, vmax, _step, _marks = slider_config(numeric_col(col))
            return [vmin, vmax]
        cont = [cont_lookup[d["col"]] if d["col"] in cont_lookup else _cont_default(d["col"]) for d in cont_ids]
        bools = [bool_lookup.get(d["col"], "all") for d in bool_ids]
        cats = [cat_lookup.get(d["col"], []) for d in cat_ids]
        tags = [tag_lookup.get(d["col"], []) for d in tag_ids]
        saved_ids = fs.get("id_search_ids", [])
        return fs.get("color_mode", DEFAULT_COLOR_MODE), fs.get("fixed_color", DEFAULT_FIXED_COLOR), fs.get("cont_col"), fs.get("colormap", DEFAULT_COLORMAP), ["reversed"] if fs.get("reversed") else [], fs.get("color_range", DEFAULT_COLOR_RANGE), fs.get("alpha", DEFAULT_ALPHA), fs.get("point_size", DEFAULT_POINT_SIZE), fs.get("marker_symbol", DEFAULT_SYMBOL), enabled, cont, bools, cats, tags, ["on"] if saved_ids else [], ", ".join(saved_ids)

    @app.callback(Output({"type":"cat-filter", "col":ALL}, "options"), Input({"type":"cat-filter", "col":ALL}, "value"), Input({"type":"cat-filter", "col":ALL}, "id"), Input({"type":"filter-enabled", "col":ALL}, "value"), Input({"type":"filter-enabled", "col":ALL}, "id"), Input({"type":"cont-slider", "col":ALL}, "value"), Input({"type":"cont-slider", "col":ALL}, "id"), Input({"type":"bool-filter", "col":ALL}, "value"), Input({"type":"bool-filter", "col":ALL}, "id"), Input({"type":"tag-filter", "col":ALL}, "value"), Input({"type":"tag-filter", "col":ALL}, "id"), prevent_initial_call=True)
    def narrow_cat_options(cat_values, cat_ids, enabled_values, enabled_ids, cont_values, cont_ids, bool_values, bool_ids, tag_values, tag_ids):
        if _ann_df.empty or not cat_ids:
            return [no_update] * len(cat_ids)
        enabled_cols = {d["col"] for d, v in zip(enabled_ids, enabled_values) if v and "on" in v}
        cont_conds = [(d["col"], v[0], v[1]) for d, v in zip(cont_ids, cont_values) if v is not None and d["col"] in enabled_cols]
        bool_map = {"all": None, "true": True, "false": False}
        bool_conds = [(d["col"], bool_map[v]) for d, v in zip(bool_ids, bool_values) if d["col"] in enabled_cols and v != "all" and bool_map.get(v) is not None]
        tag_conds = [(d["col"], t) for d, t in zip(tag_ids, tag_values) if t and d["col"] in enabled_cols]
        cat_all = [(d["col"], v) for d, v in zip(cat_ids, cat_values) if v and d["col"] in enabled_cols]
        out = []
        for d, val in zip(cat_ids, cat_values):
            col = d["col"]
            cat_loo = [(c, v) for c, v in cat_all if c != col]
            mask = apply_filters(_ann_df, cont_conds, bool_conds, cat_loo, tag_conds)
            counts = _ann_df.loc[mask, col].astype(str).value_counts()
            selected = {str(v) for v in (val or [])}
            live = sorted(counts[counts > 0].index.tolist())
            opts = [{"label": f"{v}  ({int(counts[v]):,})", "value": v} for v in live]
            opts += [{"label": f"{v}  —  0", "value": v} for v in sorted(s for s in selected if int(counts.get(s, 0)) == 0)]
            out.append(opts)
        return out

    @app.callback(Output("latent-graph", "figure"), Output("point-count", "children"), Input("data-loaded-store", "data"), Input("reload-counter", "data"), Input("coord-cols-store", "data"), Input("alpha-slider", "value"), Input("point-size-slider", "value"), Input("bg-size-slider", "value"), Input("marker-size-slider", "value"), Input("marker-alpha-slider", "value"), Input("marker-mode", "value"), Input("wf-position", "value"), Input("color-mode", "value"), Input("fixed-color-store", "data"), Input("cont-color-col", "value"), Input("colormap-select", "value"), Input("colormap-reversed", "value"), Input("color-range-mode", "value"), Input({"type":"filter-enabled", "col":ALL}, "value"), Input({"type":"filter-enabled", "col":ALL}, "id"), Input({"type":"cont-slider", "col":ALL}, "value"), Input({"type":"cont-slider", "col":ALL}, "id"), Input({"type":"bool-filter", "col":ALL}, "value"), Input({"type":"bool-filter", "col":ALL}, "id"), Input({"type":"cat-filter", "col":ALL}, "value"), Input({"type":"cat-filter", "col":ALL}, "id"), Input({"type":"tag-filter", "col":ALL}, "value"), Input({"type":"tag-filter", "col":ALL}, "id"), Input("layers-store", "data"), Input("wf-visible-store", "data"), Input("id-search-enabled", "value"), Input("id-search-input", "n_blur"), Input("marker-symbol", "value"), Input("selection-store", "data"), Input("selection-color-store", "data"), Input("cluster-region-col", "value"), Input("cluster-region-shape", "value"), Input("cluster-region-position", "value"), Input("cluster-region-opacity", "value"), Input("cluster-region-coverage", "value"), Input("cluster-region-tightness", "value"), State("id-search-input", "value"))
    def update_figure(_loaded, _reload, _coord, alpha, point_size, bg_size, marker_size, marker_alpha, marker_mode, wf_position, color_mode, fixed_color, cont_col, colormap, colormap_reversed, color_range_mode, enabled_values, enabled_ids, cont_values, cont_ids, bool_values, bool_ids, cat_values, cat_ids, tag_values, tag_ids, layers, wf_visible, id_enabled, _blur, wf_symbol, selection_ids, selection_color, region_col, region_shape, region_position, region_opacity, region_coverage, region_tightness, id_raw):
        cont_conds, bool_conds, cat_conds, tag_conds, id_search_ids = parse_current_conditions(enabled_values, enabled_ids, cont_values, cont_ids, bool_values, bool_ids, cat_values, cat_ids, tag_values, tag_ids, id_enabled, id_raw)
        with _STATE_LOCK:
            fig, n_filtered, n_cov = make_figure(cont_conds, bool_conds, cat_conds, tag_conds, color_mode, fixed_color, cont_col, colormap or DEFAULT_COLORMAP, bool(colormap_reversed and "reversed" in colormap_reversed), color_range_mode or DEFAULT_COLOR_RANGE, alpha or DEFAULT_ALPHA, point_size or DEFAULT_POINT_SIZE, bg_size or DEFAULT_BG_SIZE, marker_size or DEFAULT_MARKER_SIZE, marker_alpha or DEFAULT_MARKER_ALPHA, marker_mode or DEFAULT_MARKER_MODE, wf_position or DEFAULT_WF_POSITION, True if wf_visible is None else wf_visible, layers or [], id_search_ids, wf_symbol or DEFAULT_SYMBOL, selection_ids or [], selection_color or DEFAULT_SELECTION_COLOR, region_col, region_shape or "kde", region_position or "below", region_opacity if region_opacity is not None else 0.25, region_coverage if region_coverage is not None else 0.4, region_tightness if region_tightness is not None else 0.25)
            total = len(_ann_df)
        n_vis = len([l for l in (layers or []) if l.get("visible", True)])
        count = f"Showing {n_filtered:,} / {total:,} filtered sequences · {n_cov:,} have coordinates here"
        if n_vis:
            count += f" · {n_vis} saved layer(s) visible"
        return fig, count

    @app.callback(Output("click-details", "children"), Output("selection-store", "data", allow_duplicate=True), Input("latent-graph", "clickData"), State("selection-store", "data"), prevent_initial_call=True)
    def handle_click(click_data, current_selection):
        if not click_data or not click_data.get("points"):
            return "Click a point to show details here.", no_update
        cd = click_data["points"][0].get("customdata")
        if cd is None:
            return "No data for this point.", no_update
        seq_id = cd if isinstance(cd, str) else cd[0]
        sel = list(current_selection or [])
        if seq_id in sel:
            sel.remove(seq_id)
        else:
            sel.append(seq_id)
        with _STATE_LOCK:
            details = make_details_panel(seq_id)
        return details, sel

    @app.callback(Output("selection-count", "children"), Input("selection-store", "data"))
    def selection_count(sel):
        n = len(sel or [])
        return "No sequences selected" if n == 0 else f"{n} sequence{'s' if n != 1 else ''} selected"

    @app.callback(Output("selection-color-store", "data"), [Input({"type":"sel-color-chip", "color":c}, "n_clicks") for c in SELECTION_COLOR_OPTIONS], prevent_initial_call=True)
    def pick_selection_color(*_):
        triggered = ctx.triggered_id
        return triggered.get("color") if isinstance(triggered, dict) else no_update

    @app.callback(Output("selection-store", "data", allow_duplicate=True), Input("clear-selection-btn", "n_clicks"), prevent_initial_call=True)
    def clear_selection(_):
        return []

    @app.callback(Output("id-search-enabled", "value", allow_duplicate=True), Output("id-search-input", "value", allow_duplicate=True), Input("selection-to-wl-btn", "n_clicks"), State("selection-store", "data"), prevent_initial_call=True)
    def selection_to_wl(n, sel):
        return (["on"], ", ".join(sel)) if n and sel else (no_update, no_update)

    @app.callback(Output("selection-export-status", "children"), Input("selection-export-btn", "n_clicks"), State("selection-store", "data"), prevent_initial_call=True)
    def export_selection(n, sel):
        # Write the selected sequences to a timestamped selection cache on disk.
        # The pipeline's Boltz-2 module (scripts/sse_boltz.py) imports these.
        if not n:
            return no_update
        ids = [str(x) for x in (sel or [])]
        if not ids:
            return "No sequences selected to export."
        if COL_SEQ not in _ann_df.columns:
            return f"Cannot export: no '{COL_SEQ}' column in this entry."
        with _STATE_LOCK:
            rows = _ann_df[id_str().isin(ids)]
            seq_by_id = {str(r[_id_col]): str(r[COL_SEQ]).strip() for _, r in rows.iterrows()}
        sequences = [{"id": sid, "sequence": seq_by_id.get(sid, "")} for sid in ids]
        missing = [s["id"] for s in sequences if not s["sequence"]]
        try:
            path = selection_cache.write_selection(ENTRY.selections_dir, ENTRY.stem, sequences)
        except Exception as exc:
            return f"Export failed: {exc}"
        payload = selection_cache.read_selection(path)
        msg = f"Exported {payload.get('count', 0)} sequence(s) → selections/{path.name}"
        if missing:
            msg += f" · {len(missing)} skipped (no sequence)"
        return msg

    @app.callback(Output("extract-download", "data"), Output("extract-status", "children"), Input("extract-btn", "n_clicks"), State("layers-store", "data"), prevent_initial_call=True)
    def extract_layers(n, layers):
        visible = [l for l in (layers or []) if l.get("visible", True)]
        if not n or not visible:
            return no_update, "No visible layers to extract."
        all_ids = []
        for layer in visible:
            all_ids.extend(layer.get("ids", []))
        unique_ids = list(dict.fromkeys([str(x) for x in all_ids]))
        if not unique_ids:
            return no_update, "No sequences found in visible layers."
        out = _ann_df[_ann_df[_id_col].astype(str).isin(unique_ids)].copy()
        order = {sid: i for i, sid in enumerate(unique_ids)}
        out = out.sort_values(by=_id_col, key=lambda s: s.astype(str).map(order)).reset_index(drop=True)
        out = out.drop(columns=[c for c in out.columns if _types.get(c) == TYPE_COORDINATE])
        run_id = str(uuid.uuid4())[:8]
        filename = f"sse_extraction_{run_id}.csv"
        path = ENTRY.logs_dir / f"extraction_{run_id}.csv"
        out.to_csv(path, index=False)
        return dcc.send_data_frame(out.to_csv, filename, index=False), f"{len(out):,} unique sequences · saved log copy to {path.name}"

    def _valid_hex(value, default):
        val = (value or default).strip()
        if not val.startswith("#"):
            val = "#" + val
        if re.fullmatch(r"#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6})", val):
            return val
        return default

    def _render_figure_bytes(fig, fmt, width, height, scale):
        import plotly.io as pio
        kwargs = {} if fmt == "svg" else {"width": width, "height": height, "scale": scale}
        return pio.to_image(fig, format=fmt, **kwargs)

    def _export_with_timeout(fig, fmt, width, height, scale, timeout=45):
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(_render_figure_bytes, fig, fmt, width, height, scale)
            try:
                return future.result(timeout=timeout)
            except concurrent.futures.TimeoutError:
                raise TimeoutError(
                    f"kaleido timed out after {timeout}s. Ensure kaleido>=1.0.0 is installed."
                )

    @app.callback(
        Output("export-width", "value"),
        Output("export-height", "value"),
        Input("export-size-reset", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_export_size(n):
        if not n:
            return no_update, no_update
        return EXPORT_DEFAULT_WIDTH, EXPORT_DEFAULT_HEIGHT

    @app.callback(
        Output("export-status", "children"),
        Output("figure-download", "data"),
        Input("export-btn", "n_clicks"),
        State("latent-graph", "figure"),
        State("export-format", "value"),
        State("export-dpi", "value"),
        State("export-legend", "value"),
        State("export-transparent", "value"),
        State("export-axis-color", "value"),
        State("export-label-color", "value"),
        State("export-edge-color", "value"),
        State("export-bg-color", "value"),
        State("export-width", "value"),
        State("export-height", "value"),
        State("export-destination", "value"),
        prevent_initial_call=True,
    )
    def export_figure(n, figure, fmt, dpi, legend_val, transparent_val,
                      axis_color, label_color, edge_color, bg_color,
                      width, height, dest):
        if not n or figure is None:
            return no_update, no_update
        try:
            import datetime
            fig = go.Figure(figure)
            fmt = fmt or "png"
            dpi = dpi or 300
            width = int(width or EXPORT_DEFAULT_WIDTH)
            height = int(height or EXPORT_DEFAULT_HEIGHT)
            axis_col = _valid_hex(axis_color, "#000000")
            label_col = _valid_hex(label_color, axis_col)
            edge_col = _valid_hex(edge_color, "#000000")
            bg_col = _valid_hex(bg_color, "#ffffff")
            transparent = bool(transparent_val and "on" in transparent_val)
            bg = "rgba(0,0,0,0)" if transparent else bg_col

            fig.update_layout(
                showlegend=bool(legend_val and "show" in legend_val),
                paper_bgcolor=bg,
                plot_bgcolor=bg,
                font=dict(color=label_col),
                legend=dict(font=dict(color=label_col)),
                xaxis=dict(
                    color=axis_col,
                    linecolor=axis_col,
                    tickcolor=axis_col,
                    tickfont=dict(color=label_col),
                    title=dict(font=dict(color=label_col)),
                    showgrid=False,
                    zeroline=False,
                ),
                yaxis=dict(
                    color=axis_col,
                    linecolor=axis_col,
                    tickcolor=axis_col,
                    tickfont=dict(color=label_col),
                    title=dict(font=dict(color=label_col)),
                    showgrid=False,
                    zeroline=False,
                ),
            )

            for trace in fig.data:
                marker = getattr(trace, "marker", None)
                if marker is not None:
                    if getattr(marker, "line", None) is not None:
                        trace.marker.line.color = edge_col
                    else:
                        trace.marker.line = dict(color=edge_col, width=0.5)

            scale = dpi / 96
            img_bytes = _export_with_timeout(fig, fmt, width, height, scale)
            filename = f"sse_figure_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.{fmt}"
            if dest == "server":
                ENTRY.figures_dir.mkdir(parents=True, exist_ok=True)
                out = ENTRY.figures_dir / filename
                out.write_bytes(img_bytes)
                return f"Saved to {out}", no_update
            return f"Downloaded {filename}", dcc.send_bytes(img_bytes, filename)
        except ImportError:
            return 'Export failed: kaleido is not installed. Run: pip install "kaleido>=1.0.0"', no_update
        except TimeoutError as exc:
            return f"Export failed: {exc}", no_update
        except Exception as exc:
            return f"Export failed: {exc}", no_update

    @app.callback(Output("filter-panel", "children", allow_duplicate=True), Output("col-settings-panel", "children", allow_duplicate=True), Output("filter-pending-store", "data", allow_duplicate=True), Input("rebuild-filters-btn", "n_clicks"), State({"type":"col-override", "col":ALL}, "value"), State({"type":"col-override", "col":ALL}, "id"), prevent_initial_call=True)
    def rebuild_filters(n, override_values, override_ids):
        if not n:
            return no_update, no_update, no_update
        for id_d, val in zip(override_ids, override_values):
            col = id_d["col"]
            if col in _col_meta:
                _col_meta[col]["override"] = val
                if val == "tag_split" and not _col_meta[col].get("tags"):
                    tags = set()
                    for v in _ann_df[col].dropna().astype(str):
                        tags.update(t.strip() for t in v.split(",") if t.strip())
                    _col_meta[col]["tags"] = sorted(tags)
        return html.Div(), html.Div(), {"source": "rebuild", "counter": n}


    return app


def main(argv=None):
    ap = argparse.ArgumentParser(description="Open the SSE visualizer for one entry.")
    ap.add_argument("entry", help="Entry stem, entry directory, or .sse.tsv datafile path")
    ap.add_argument("--port", type=int, default=PORT)
    args = ap.parse_args(argv)
    try:
        app = build_app(args.entry)
    except SSEError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"Starting Sequence Space Explorer for {ENTRY.stem} on http://127.0.0.1:{args.port}")
    # use_reloader=False: the reloader forks a second process that actually holds
    # the port, which the pipeline runner can't see or stop, so it would survive
    # shutdown and keep the port bound. One process shuts down cleanly.
    app.run(debug=True, port=args.port, use_reloader=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
