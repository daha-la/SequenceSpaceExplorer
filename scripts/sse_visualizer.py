"""SSE visualizer: Dash viewer for one entry-centered SSE datafile.

Run:
    python scripts/sse_visualizer.py <entry-stem|entry-dir|datafile.sse.tsv>

The app reads exactly one .sse.tsv datafile, using the Type row as the contract:
id / label / coordinate. Structure prediction and RMSD are exposed in the UI,
but the heavy logic lives in sse_tools.boltz and sse_tools.rmsd.
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
from sse_tools import jobs as job_store
from sse_tools import boltz as boltz_backend
from sse_tools import rmsd as rmsd_backend

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
SECTION_STYLE = {"fontWeight": "600", "fontSize": "12px", "margin": "10px 0 6px 0", "color": "#555"}
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

_boltz_executor = concurrent.futures.ThreadPoolExecutor(max_workers=2)
_boltz_futures: dict[str, concurrent.futures.Future] = {}
_STATE_LOCK = threading.RLock()


def reload_state() -> str:
    """Reload the datafile into module globals. Returns a warning/status string.

    The visualizer stores a few derived values as module globals because Dash
    callbacks need fast access to them. Reload builds a complete new state first
    and then swaps all globals while holding one lock, so render callbacks cannot
    observe a half-updated mix of old and new dataframe/axis metadata.
    """
    global _STATE, _ann_df, _types, _col_meta, _id_col, _x_col, _y_col, _query_ids, _name_cols
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
    return new_state.warning


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
            num = pd.to_numeric(pool_df[col], errors="coerce")
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
        id_match = pool_df[_id_col].astype(str).isin([str(x) for x in id_search_ids])
        name_match = pd.Series(False, index=pool_df.index)
        for nc in _name_cols:
            if nc in pool_df.columns:
                name_match |= pool_df[nc].astype(str).isin([str(x) for x in id_search_ids])
        mask &= id_match | name_match
    return mask




def _coerce_float(value, default: float) -> float:
    """Parse Dash numeric input robustly, including comma decimal locales."""
    if value is None or value == "":
        return default
    if isinstance(value, str):
        value = value.strip().replace(",", ".")
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _coerce_int(value, default: int) -> int:
    try:
        return int(float(str(value).strip().replace(",", ".")))
    except (TypeError, ValueError):
        return default

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
    x = pd.to_numeric(df[x_col], errors="coerce")
    y = pd.to_numeric(df[y_col], errors="coerce")
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
            x=pd.to_numeric(plot_df[x_col], errors="coerce"), y=pd.to_numeric(plot_df[y_col], errors="coerce"), mode="markers",
            marker=dict(size=point_size, color=fixed_color, symbol=symbol, opacity=alpha, line=dict(width=0.5, color="black")),
            name="Working filter", customdata=plot_df[id_col].astype(str).tolist(), hovertemplate="<b>%{customdata}</b><extra></extra>",
        ))
    elif color_mode == "continuous" and cont_col and cont_col in plot_df.columns:
        vals = pd.to_numeric(plot_df[cont_col], errors="coerce")
        has_val = vals.notna()
        no_val = plot_df[~has_val]
        color_df = plot_df[has_val]
        if not no_val.empty:
            fig.add_trace(go.Scattergl(
                x=pd.to_numeric(no_val[x_col], errors="coerce"), y=pd.to_numeric(no_val[y_col], errors="coerce"), mode="markers",
                marker=dict(size=point_size, color="lightgrey", symbol=symbol, opacity=alpha, line=dict(width=0.5, color="#aaa")),
                name="Working filter — no value", customdata=no_val[id_col].astype(str).tolist(), hovertemplate="<b>%{customdata}</b><br>No value<extra></extra>",
            ))
        if color_df.empty:
            return cont_axis_idx
        v = vals[has_val]
        if color_range_mode == "global":
            full = pd.to_numeric(df[cont_col], errors="coerce").dropna()
            cmin, cmax = (float(full.min()), float(full.max())) if not full.empty else (0.0, 1.0)
        else:
            cmin, cmax = (float(v.min()), float(v.max())) if not v.empty else (0.0, 1.0)
        cont_axis_idx += 1
        _add_cont_trace(fig, cont_axis_idx, f"Working filter — {cont_col}", color_df[id_col].astype(str).tolist(), pd.to_numeric(color_df[x_col], errors="coerce").tolist(), pd.to_numeric(color_df[y_col], errors="coerce").tolist(), v.tolist(), cmin, cmax, colormap, reversed_, point_size, alpha, symbol)
    return cont_axis_idx


def make_figure(cont_conds, bool_conds, cat_conds, tag_conds, color_mode, fixed_color, cont_col, colormap, reversed_, color_range_mode, alpha, point_size, bg_size, marker_size, marker_alpha, marker_mode, wf_position, wf_visible, layers, id_search_ids=None, wf_symbol=DEFAULT_SYMBOL, selection_ids=None, selection_color=DEFAULT_SELECTION_COLOR):
    df = _ann_df
    x_col, y_col = _x_col, _y_col
    if df.empty:
        return go.Figure(), 0, 0
    if not x_col or not y_col:
        fig = go.Figure()
        fig.update_layout(template="simple_white", annotations=[dict(text="No coordinates yet — run sse_coordinates.py, then reload.", showarrow=False, xref="paper", yref="paper", x=0.5, y=0.5, font=dict(size=14, color="#aaa"))])
        return fig, 0, 0

    fig = go.Figure()
    coord_mask = _covered(df, x_col, y_col)
    filter_mask = apply_filters(df, cont_conds, bool_conds, cat_conds, tag_conds, id_search_ids)
    plot_df = df[filter_mask]
    n_filtered = int(filter_mask.sum())
    n_covered_filtered = int((filter_mask & coord_mask).sum())

    query_mask = boolean_mask(df[COL_QUERY]).fillna(False) if COL_QUERY in df.columns else pd.Series(False, index=df.index)
    bg_df = df[coord_mask & ~query_mask]
    fig.add_trace(go.Scattergl(
        x=pd.to_numeric(bg_df[x_col], errors="coerce"), y=pd.to_numeric(bg_df[y_col], errors="coerce"), mode="markers",
        marker=dict(size=bg_size, color="lightgrey", opacity=0.4), name="All sequences", hoverinfo="skip", showlegend=True,
    ))
    if marker_mode == "bottom":
        _add_query_traces(fig, df, _id_col, x_col, y_col, marker_size, marker_alpha)

    cont_axis_idx = 0
    if wf_visible and wf_position == "bottom":
        cont_axis_idx = _add_working_filter(fig, plot_df, df, color_mode, fixed_color, cont_col, colormap, reversed_, color_range_mode, alpha, point_size, cont_axis_idx, _id_col, x_col, y_col, wf_symbol)

    for layer in reversed([l for l in (layers or []) if l.get("visible", True)]):
        ids = set(str(x) for x in layer.get("ids", []))
        if not ids:
            continue
        sub = df[df[_id_col].astype(str).isin(ids)]
        sub = sub[_covered(sub, x_col, y_col)]
        if sub.empty:
            continue
        style = layer.get("style", {})
        l_alpha = style.get("alpha", alpha)
        l_size = style.get("point_size", point_size)
        l_symbol = style.get("marker_symbol", DEFAULT_SYMBOL)
        if style.get("color_mode") == "continuous" and style.get("cont_col") in sub.columns:
            ccol = style.get("cont_col")
            vals = pd.to_numeric(sub[ccol], errors="coerce")
            has_val = vals.notna()
            if has_val.any():
                v = vals[has_val]
                cont_axis_idx += 1
                cmin, cmax = (float(v.min()), float(v.max())) if style.get("color_range") != "global" else (
                    float(pd.to_numeric(df[ccol], errors="coerce").dropna().min()),
                    float(pd.to_numeric(df[ccol], errors="coerce").dropna().max()),
                )
                _add_cont_trace(fig, cont_axis_idx, layer.get("name", "Layer"), sub.loc[has_val, _id_col].astype(str).tolist(), pd.to_numeric(sub.loc[has_val, x_col], errors="coerce").tolist(), pd.to_numeric(sub.loc[has_val, y_col], errors="coerce").tolist(), v.tolist(), cmin, cmax, style.get("colormap", DEFAULT_COLORMAP), style.get("reversed", False), l_size, l_alpha, l_symbol)
        else:
            fig.add_trace(go.Scattergl(
                x=pd.to_numeric(sub[x_col], errors="coerce"), y=pd.to_numeric(sub[y_col], errors="coerce"), mode="markers",
                marker=dict(size=l_size, color=style.get("fixed_color", DEFAULT_FIXED_COLOR), symbol=l_symbol, opacity=l_alpha, line=dict(width=0.5, color="black")),
                name=layer.get("name", "Layer"), customdata=sub[_id_col].astype(str).tolist(), hovertemplate="<b>%{customdata}</b><extra></extra>",
            ))

    if wf_visible and wf_position == "top":
        cont_axis_idx = _add_working_filter(fig, plot_df, df, color_mode, fixed_color, cont_col, colormap, reversed_, color_range_mode, alpha, point_size, cont_axis_idx, _id_col, x_col, y_col, wf_symbol)
    if marker_mode == "top":
        _add_query_traces(fig, df, _id_col, x_col, y_col, marker_size, marker_alpha)
    if selection_ids:
        sel_df = df[df[_id_col].astype(str).isin([str(x) for x in selection_ids])]
        sel_df = sel_df[_covered(sel_df, x_col, y_col)]
        if not sel_df.empty:
            fig.add_trace(go.Scattergl(
                x=pd.to_numeric(sel_df[x_col], errors="coerce"), y=pd.to_numeric(sel_df[y_col], errors="coerce"), mode="markers",
                marker=dict(size=point_size + 6, symbol="circle-open", color=selection_color, opacity=1.0, line=dict(width=2.5, color=selection_color)),
                name=f"Selected ({len(sel_df)})", customdata=sel_df[_id_col].astype(str).tolist(), hovertemplate="<b>%{customdata}</b><br>Selected<extra></extra>",
            ))
    fig.update_layout(
        xaxis=dict(title=x_col, constrain="domain"),
        yaxis=dict(title=y_col, scaleanchor="x", scaleratio=1),
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
    rows = [html.Tr([html.Td("ID", style={"fontWeight": "bold", "padding": "4px 10px"}), html.Td(sequence_id, style={"padding": "4px 10px"})])]
    for col in _ann_df.columns:
        if col == _id_col or _types.get(col) == TYPE_COORDINATE:
            continue
        val = row.get(col, "")
        val = "" if pd.isna(val) else str(val)
        rows.append(html.Tr([html.Td(col, style={"fontWeight": "bold", "padding": "4px 10px", "verticalAlign": "top"}), html.Td(val, style={"padding": "4px 10px", "wordBreak": "break-all"})]))
    return html.Div([html.H4(f"Selected: {sequence_id}", style={"marginBottom": "10px"}), html.Table(rows, style={"borderCollapse": "collapse", "width": "100%", "fontSize": "13px"})])


def make_filter_panel():
    if not _col_meta:
        return html.P("No columns loaded.", style={"fontSize": "12px", "color": "#aaa"})
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
            inp = {"width": "70px", "fontSize": "11px", "padding": "2px 4px", "border": "1px solid #ddd", "borderRadius": "3px", "textAlign": "center"}
            children.append(html.Div([checkbox(col, col), html.Div(id={"type": "filter-control", "col": col}, children=[html.Div([dcc.Input(id={"type": "cont-min-input", "col": col}, type="number", value=vmin, step="any", debounce=True, style=inp), html.Span("–", style={"margin": "0 4px", "fontSize": "11px", "color": "#777"}), dcc.Input(id={"type": "cont-max-input", "col": col}, type="number", value=vmax, step="any", debounce=True, style=inp)], style={"display": "flex", "alignItems": "center", "marginBottom": "4px"}), dcc.RangeSlider(id={"type": "cont-slider", "col": col}, min=vmin, max=vmax, step=step, value=[vmin, vmax], marks=marks, allowCross=False, updatemode="mouseup")], style={"display": "none", **CONTROL_WRAPPER})], style={"marginBottom": "6px"}))
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
    return children or html.P("No filterable label columns.", style={"fontSize": "12px", "color": "#aaa"})


def make_col_settings_panel():
    opts = [{"label": "Continuous", "value": "continuous"}, {"label": "Boolean", "value": "boolean"}, {"label": "Categorical", "value": "categorical"}, {"label": "Tag split", "value": "tag_split"}, {"label": "Skip", "value": "skip"}]
    rows = []
    for col, meta in _col_meta.items():
        if _types.get(col) != TYPE_LABEL:
            continue
        rows.append(html.Div([html.Span(col, title=col, style={"fontSize": "11px", "fontWeight": "500", "overflow": "hidden", "textOverflow": "ellipsis", "whiteSpace": "nowrap", "flexGrow": "1"}), dcc.Dropdown(id={"type": "col-override", "col": col}, options=opts, value=effective_type(col), clearable=False, style={"fontSize": "11px", "width": "120px"})], style={"display": "flex", "alignItems": "center", "gap": "6px", "padding": "4px 2px", "borderBottom": "1px solid #f5f5f5"}))
    return html.Div(rows, style={"maxHeight": "350px", "overflowY": "auto"})


def make_sidebar(layers):
    if not layers:
        return html.Div("No saved layers yet.", style={"fontSize": "12px", "color": "#aaa", "padding": "8px 0"})
    rows = []
    n = len(layers)
    for i, layer in enumerate(layers):
        lid = layer["id"]
        visible = layer.get("visible", True)
        style = layer.get("style", {})
        swatch_style = {"width": "12px", "height": "12px", "borderRadius": "2px", "border": "1px solid #ccc", "flexShrink": "0"}
        if style.get("color_mode") == "continuous":
            swatch = html.Div(style={**swatch_style, "background": "linear-gradient(to bottom, #440154, #31688e, #35b779, #fde725)"})
        else:
            swatch = html.Div(style={**swatch_style, "backgroundColor": style.get("fixed_color", DEFAULT_FIXED_COLOR)})
        btn = {"background": "none", "border": "none", "cursor": "pointer", "padding": "1px 3px", "fontSize": "11px", "lineHeight": "1"}
        n_total = len(layer.get("ids", []))
        n_cov = n_total
        if _x_col and _y_col and not _ann_df.empty:
            sub = _ann_df[_ann_df[_id_col].astype(str).isin([str(x) for x in layer.get("ids", [])])]
            n_cov = int(_covered(sub, _x_col, _y_col).sum())
        rows.append(html.Div([swatch, html.Div([html.Div(layer.get("name", "Layer"), title=layer.get("name", "Layer"), style={"fontSize": "11px", "fontWeight": "500", "overflow": "hidden", "textOverflow": "ellipsis", "whiteSpace": "nowrap"}), html.Div(f"{n_cov:,}/{n_total:,} visible here", style={"fontSize": "10px", "color": "#aaa"})], style={"flexGrow": "1", "minWidth": "0", "margin": "0 4px"}), html.Button("↑", id={"type": "layer-up", "lid": lid}, disabled=i == 0, style=btn), html.Button("↓", id={"type": "layer-down", "lid": lid}, disabled=i == n - 1, style=btn), html.Button("⤴", id={"type": "layer-load", "lid": lid}, title="Load into working filter", style={**btn, "color": "#2980b9"}), html.Button("👁" if visible else "🚫", id={"type": "layer-toggle", "lid": lid}, style={**btn, "opacity": "1" if visible else "0.4"}), html.Button("✕", id={"type": "layer-delete", "lid": lid}, style={**btn, "color": "#e74c3c"})], style={"display": "flex", "alignItems": "center", "padding": "5px 2px", "borderBottom": "1px solid #f0f0f0", "opacity": "1" if visible else "0.5"}))
    return html.Div(rows, style={"maxHeight": "400px", "overflowY": "auto"})


def current_job_records():
    data = job_store.read_jobs(ENTRY.jobs_path, mark_stale=False)
    return data.get("boltz", {}), data.get("rmsd", {})


def make_job_table(jobs: dict):
    if not jobs:
        return html.Div("No Boltz jobs yet.", style={"fontSize": "11px", "color": "#aaa"})
    rows = []
    colors = {"queued":"#7f8c8d", "msa":"#2980b9", "predicting":"#e67e22", "done":"#27ae60", "cached":"#27ae60", "error":"#e74c3c", "interrupted":"#e67e22"}
    for job in sorted(jobs.values(), key=lambda j: j.get("updated_utc", ""), reverse=True):
        status = job.get("status", "")
        sid = job.get("sequence_id", "")
        kind = job.get("kind", "apo")
        typ = "apo" if kind == "apo" else f"holo · {job.get('smiles_label') or job.get('smiles_hash', '')}"
        rows.append(html.Tr([html.Td(sid[:16] + ("…" if len(sid) > 16 else ""), title=sid, style={"fontFamily": "monospace", "fontSize": "10px", "padding": "3px 5px"}), html.Td(typ, title=typ, style={"fontSize": "9px", "padding": "3px 5px"}), html.Td(status, style={"fontSize": "10px", "fontWeight": "600", "color": colors.get(status, "#555"), "padding": "3px 5px"}), html.Td(f"{job.get('ptm'):.3f}" if isinstance(job.get("ptm"), (int, float)) else "—", style={"fontFamily": "monospace", "fontSize": "10px", "textAlign": "right", "padding": "3px 5px"}), html.Td(f"{job.get('plddt'):.1f}" if isinstance(job.get("plddt"), (int, float)) else "—", style={"fontFamily": "monospace", "fontSize": "10px", "textAlign": "right", "padding": "3px 5px"})], style={"borderBottom": "1px solid #f5f5f5"}))
    th = {"padding": "3px 5px", "fontSize": "10px", "color": "#7f8c8d", "fontWeight": "600", "textAlign": "left"}
    return html.Table([html.Thead(html.Tr([html.Th("ID", style=th), html.Th("Type", style=th), html.Th("Status", style=th), html.Th("pTM", style={**th, "textAlign": "right"}), html.Th("pLDDT", style={**th, "textAlign": "right"})])), html.Tbody(rows)], style={"width": "100%", "borderCollapse": "collapse"})


def boltz_summary(jobs: dict) -> str:
    if not jobs:
        return ""
    counts = {}
    for j in jobs.values():
        counts[j.get("status", "unknown")] = counts.get(j.get("status", "unknown"), 0) + 1
    order = ["done", "cached", "predicting", "msa", "queued", "interrupted", "error"]
    return " · ".join(f"{counts[k]} {k}" for k in order if k in counts)


def make_rmsd_results_table(results: list[dict]):
    if not results:
        return html.Div()
    results = sorted(results, key=lambda r: (r.get("method", ""), np.isnan(r.get("rmsd", np.nan)) if isinstance(r.get("rmsd"), float) else False, r.get("rmsd") if isinstance(r.get("rmsd"), (int, float)) and not np.isnan(r.get("rmsd")) else 9999))
    rows = []
    for r in results:
        rmsd = r.get("rmsd")
        rmsd_s = f"{rmsd:.3f}" if isinstance(rmsd, (int, float)) and not np.isnan(rmsd) else "—"
        rows.append(html.Tr([html.Td(r.get("query_id", ""), style={"fontFamily": "monospace", "fontSize": "10px", "padding": "3px 5px"}), html.Td(str(r.get("query_rank", 0)), style={"fontSize": "10px", "textAlign": "right", "padding": "3px 5px"}), html.Td(str(r.get("n_aligned", "")), style={"fontSize": "10px", "textAlign": "right", "padding": "3px 5px"}), html.Td(rmsd_s, style={"fontFamily": "monospace", "fontWeight": "600", "fontSize": "10px", "textAlign": "right", "padding": "3px 5px"}), html.Td(r.get("method", ""), style={"fontSize": "9px", "padding": "3px 5px"}), html.Td("✓ cached" if r.get("cached") else "", style={"fontSize": "9px", "color": "#27ae60", "padding": "3px 5px"})], style={"borderBottom": "1px solid #f5f5f5"}))
    th = {"padding": "3px 5px", "fontSize": "10px", "color": "#7f8c8d", "fontWeight": "600", "textAlign": "left"}
    return html.Table([html.Thead(html.Tr([html.Th("Query", style=th), html.Th("Rank", style={**th, "textAlign": "right"}), html.Th("Aligned", style={**th, "textAlign": "right"}), html.Th("RMSD Å", style={**th, "textAlign": "right"}), html.Th("Method", style=th), html.Th("", style=th)])), html.Tbody(rows)], style={"width": "100%", "borderCollapse": "collapse"})


def build_app(entry_arg: str):
    global ENTRY
    ENTRY = resolve_entry(entry_arg)
    warning = reload_state()
    loaded_layers, layer_msg = layer_store.validate_layers(layer_store.read_layers(ENTRY.layers_path), _ann_df[_id_col].astype(str).tolist())
    if layer_msg:
        layer_store.write_layers(ENTRY.layers_path, loaded_layers)
    job_store.read_jobs(ENTRY.jobs_path, mark_stale=True)

    app = Dash(__name__, suppress_callback_exceptions=True)
    app.title = f"SSE — {ENTRY.stem}"

    coord_system_options = [{"label": k, "value": k} for k in _STATE.coord_systems]
    first_system = next(iter(_STATE.coord_systems), None)
    cont_cols = [c for c in cols_of_type("continuous") if _types.get(c) == TYPE_LABEL]

    def coord_axis_options(system):
        cols = _STATE.coord_systems.get(system, []) if system else []
        return [{"label": c, "value": c} for c in cols]

    def all_coord_opts():
        return [{"label": c, "value": c} for c in _STATE.coord_cols]

    _boltz_jobs, _rmsd_jobs = current_job_records()

    app.layout = html.Div([
        dcc.Store(id="data-loaded-store", data=True),
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
        dcc.Store(id="boltz-key-valid-store", data=False),
        dcc.Store(id="boltz-clicked-id-store", data=None),
        dcc.Interval(id="boltz-interval", interval=3000, n_intervals=0, disabled=True),
        dcc.Download(id="extract-download"),
        dcc.Download(id="figure-download"),

        html.Div([html.Div([html.H2("Sequence Space Explorer", style={"margin": "0", "color": "#2c3e50"}), html.Span(id="subtitle-text", children=f"Entry: {ENTRY.stem} · {_ann_df.shape[0]:,} rows · {len(_STATE.coord_cols)} coordinate column(s)", style={"color": "#7f8c8d", "fontSize": "13px"})]), html.Button("Reload datafile", id="reload-btn", n_clicks=0, style={"padding": "6px 10px", "backgroundColor": "#ecf0f1", "border": "1px solid #bdc3c7", "borderRadius": "4px", "cursor": "pointer", "fontSize": "12px"})], style={"display": "flex", "alignItems": "center", "justifyContent": "space-between", "marginBottom": "16px", "borderBottom": "2px solid #ecf0f1", "paddingBottom": "10px"}),

        html.Div([
            html.Div([
                html.Details([html.Summary("Entry", style={"fontWeight": "bold", "cursor": "pointer", "marginBottom": "10px"}), html.Div([html.Div(f"Entry: {ENTRY.stem}", style={"fontSize": "12px"}), html.Div(str(ENTRY.datafile_path), style={"fontSize": "10px", "color": "#888", "wordBreak": "break-all"}), html.Div(id="reload-status", children=warning or layer_msg or "", style={"fontSize": "11px", "color": "#e67e22", "marginTop": "6px"})])], open=True, style={"marginBottom": "16px"}),

                html.Details([html.Summary("Coordinates", style={"fontWeight": "bold", "cursor": "pointer", "marginBottom": "10px"}), dcc.RadioItems(id="coord-mode", options=[{"label": " Coordinate system mode", "value": "system"}, {"label": " Advanced free-axis mode", "value": "free"}], value="system", labelStyle={"display": "block", "fontSize": "12px", "marginBottom": "3px"}), html.Div(id="coord-system-panel", children=[html.Label("Coordinate system", style={"fontSize": "12px"}), dcc.Dropdown(id="coord-system-select", options=coord_system_options, value=first_system, clearable=False, style={"fontSize": "12px", "marginBottom": "4px"}), html.Label("X axis", style={"fontSize": "12px"}), dcc.Dropdown(id="x-axis-system", options=coord_axis_options(first_system), value=_x_col, clearable=False, style={"fontSize": "12px", "marginBottom": "4px"}), html.Label("Y axis", style={"fontSize": "12px"}), dcc.Dropdown(id="y-axis-system", options=coord_axis_options(first_system), value=_y_col, clearable=False, style={"fontSize": "12px"})]), html.Div(id="coord-free-panel", children=[html.Label("X axis", style={"fontSize": "12px"}), dcc.Dropdown(id="x-axis-free", options=all_coord_opts(), value=_x_col, clearable=False, style={"fontSize": "12px", "marginBottom": "4px"}), html.Label("Y axis", style={"fontSize": "12px"}), dcc.Dropdown(id="y-axis-free", options=all_coord_opts(), value=_y_col, clearable=False, style={"fontSize": "12px"})], style={"display": "none"}), html.Div(id="coord-warning", style={"fontSize": "11px", "color": "#e67e22", "marginTop": "6px"})], open=True, style={"marginBottom": "16px"}),

                html.Details([html.Summary("Appearance", style={"fontWeight": "bold", "cursor": "pointer", "marginBottom": "10px"}), html.P("Filtered points", style={**SECTION_STYLE, "margin": "0 0 4px 0"}), html.Label("Opacity", style={"fontSize": "12px"}), dcc.Slider(id="alpha-slider", min=0.05, max=1.0, step=0.05, value=DEFAULT_ALPHA, marks={0.05: "0.05", 0.5: "0.5", 1.0: "1"}), html.Label("Size", style={"fontSize": "12px"}), dcc.Slider(id="point-size-slider", min=2, max=20, step=1, value=DEFAULT_POINT_SIZE, marks={2: "2", 6: "6", 12: "12", 20: "20"}), html.Label("Shape", style={"fontSize": "12px"}), dcc.Dropdown(id="marker-symbol", options=MARKER_SYMBOL_OPTIONS, value=DEFAULT_SYMBOL, clearable=False, style={"fontSize": "12px", "marginBottom": "6px"}), html.P("Background points", style={**SECTION_STYLE, "margin": "14px 0 4px 0"}), html.Label("Size", style={"fontSize": "12px"}), dcc.Slider(id="bg-size-slider", min=1, max=12, step=1, value=DEFAULT_BG_SIZE, marks={1:"1",4:"4",8:"8",12:"12"}), html.P("Query markers", style={**SECTION_STYLE, "margin": "14px 0 4px 0"}), html.Label("Opacity", style={"fontSize":"12px"}), dcc.Slider(id="marker-alpha-slider", min=0.05, max=1.0, step=0.05, value=DEFAULT_MARKER_ALPHA, marks={0.05:"0.05",0.5:"0.5",1.0:"1"}), html.Label("Size", style={"fontSize":"12px"}), dcc.Slider(id="marker-size-slider", min=6, max=40, step=1, value=DEFAULT_MARKER_SIZE, marks={6:"6",14:"14",28:"28",40:"40"}), html.Label("Position", style={"fontSize":"12px"}), dcc.RadioItems(id="marker-mode", options=[{"label":" On top","value":"top"},{"label":" Below overlays","value":"bottom"},{"label":" Hide","value":"none"}], value=DEFAULT_MARKER_MODE, labelStyle={"display":"block","fontSize":"13px"}), html.P("Working filter position", style={**SECTION_STYLE, "margin": "14px 0 4px 0"}), dcc.RadioItems(id="wf-position", options=[{"label":" On top of saved layers","value":"top"},{"label":" Below saved layers","value":"bottom"}], value=DEFAULT_WF_POSITION, labelStyle={"display":"block","fontSize":"13px"})], open=False, style={"marginBottom": "16px"}),

                html.Details([html.Summary("Colour", style={"fontWeight": "bold", "cursor": "pointer", "marginBottom": "10px"}), dcc.RadioItems(id="color-mode", options=[{"label":" Fixed color","value":"fixed"},{"label":" Continuous","value":"continuous"}], value=DEFAULT_COLOR_MODE, labelStyle={"display":"block","fontSize":"13px"}), html.Div(id="fixed-color-panel", children=[html.Label("Pick a color", style={"fontSize":"12px"}), html.Div([html.Div(id={"type":"color-chip","color":c}, style={"width":"22px","height":"22px","backgroundColor":c,"borderRadius":"3px","cursor":"pointer","border":"2px solid transparent","display":"inline-block","marginRight":"4px","marginBottom":"4px"}) for c in FIXED_COLOR_OPTIONS])]), html.Div(id="continuous-color-panel", children=[html.Label("Color by", style={"fontSize":"12px"}), dcc.Dropdown(id="cont-color-col", options=[{"label":c,"value":c} for c in cont_cols], value=cont_cols[0] if cont_cols else None, clearable=False, style={"fontSize":"13px","marginBottom":"6px"}), html.Label("Colormap", style={"fontSize":"12px"}), dcc.Dropdown(id="colormap-select", options=[{"label":c,"value":c} for c in COLORMAP_OPTIONS], value=DEFAULT_COLORMAP, clearable=False, style={"fontSize":"13px","marginBottom":"6px"}), dcc.Checklist(id="colormap-reversed", options=[{"label":html.Span(" Reverse colormap", style={"fontSize":"12px"}),"value":"reversed"}], value=[]), html.Label("Color range", style={"fontSize":"12px"}), dcc.RadioItems(id="color-range-mode", options=[{"label":" Global","value":"global"},{"label":" Subset","value":"subset"}], value=DEFAULT_COLOR_RANGE, labelStyle={"display":"block","fontSize":"12px"})], style={"display":"none"})], open=True, style={"marginBottom":"16px"}),

                html.Details([html.Summary("Filters", style={"fontWeight":"bold","cursor":"pointer","marginBottom":"10px"}), html.Div(id="filter-panel", children=make_filter_panel())], open=True, style={"marginBottom":"16px"}),
                html.Details([html.Summary("Search by ID", style={"fontWeight":"bold","cursor":"pointer","marginBottom":"10px"}), dcc.Checklist(id="id-search-enabled", options=[{"label":html.Span(" Enable ID/name search", style=LABEL_STYLE),"value":"on"}], value=[]), html.Div(id="id-search-control", children=[dcc.Textarea(id="id-search-input", placeholder="IDs separated by commas", style={"width":"100%","fontSize":"12px","minHeight":"60px","fontFamily":"monospace"}), html.Div(id="id-search-status", style={"fontSize":"11px","marginTop":"4px","color":"#555"})], style={"display":"none", **CONTROL_WRAPPER})], open=False, style={"marginBottom":"16px"}),
                html.Details([html.Summary("Save layer", style={"fontWeight":"bold","cursor":"pointer","marginBottom":"10px"}), dcc.Input(id="layer-name-input", type="text", placeholder="Auto-generated if blank", style={"width":"100%","fontSize":"12px","marginBottom":"8px"}), html.Button("Save layer", id="save-layer-btn", n_clicks=0, style={"width":"100%","padding":"6px","backgroundColor":"#2c3e50","color":"white","border":"none","borderRadius":"4px","cursor":"pointer","fontSize":"13px"}), html.Div(id="save-layer-status", style={"fontSize":"11px","marginTop":"6px","color":"#555"})], open=True, style={"marginBottom":"16px"}),
                html.Details([html.Summary("Column settings", style={"fontWeight":"bold","cursor":"pointer","marginBottom":"10px"}), html.Button("Rebuild filter panel", id="rebuild-filters-btn", n_clicks=0, style={"width":"100%","padding":"5px","fontSize":"12px","marginBottom":"8px"}), html.Div(id="col-settings-panel", children=make_col_settings_panel())], open=False, style={"marginBottom":"16px"}),
            ], style={"width": "290px", "minWidth": "270px", "flexShrink": "0", "overflowY": "auto", "maxHeight": "90vh", "paddingRight": "12px", "borderRight": "1px solid #ecf0f1"}),

            html.Div([html.Div([html.Span(id="point-count", style={"fontSize":"12px","color":"#7f8c8d"}), html.Button("👁", id="wf-toggle-btn", n_clicks=0, title="Show/hide working filter", style={"background":"none","border":"none","cursor":"pointer","fontSize":"14px","marginLeft":"8px"})], style={"display":"flex","alignItems":"center","marginBottom":"6px"}), html.Div([html.Span(id="selection-count", children="No sequences selected", style={"fontSize":"12px","color":"#7f8c8d","flexGrow":"1"}), html.Div([html.Div(style={"width":"16px","height":"16px","backgroundColor":c,"borderRadius":"3px","cursor":"pointer","display":"inline-block","marginRight":"3px","border":"2px solid #ccc"}, id={"type":"sel-color-chip","color":c}, title=c) for c in SELECTION_COLOR_OPTIONS], style={"display":"inline-flex","alignItems":"center","marginRight":"6px"}), html.Button("Clear", id="clear-selection-btn", n_clicks=0, style={"fontSize":"11px","padding":"2px 6px","marginRight":"4px"}), html.Button("→ Working layer", id="selection-to-wl-btn", n_clicks=0, style={"fontSize":"11px","padding":"2px 6px"})], style={"display":"flex","alignItems":"center","marginBottom":"6px","padding":"4px 6px","backgroundColor":"#f8f9fa","borderRadius":"4px","border":"1px solid #ecf0f1"}), dcc.Graph(id="latent-graph", style={"height":"75vh"}, config={"displaylogo":False,"scrollZoom":True,"modeBarButtonsToAdd":["lasso2d","select2d"]}), html.Div(id="click-details", children="Click a point to show details here.", style={"marginTop":"16px","padding":"14px","border":"1px solid #ecf0f1","borderRadius":"8px","fontSize":"13px","minHeight":"60px"}), html.Div(id="load-warning", style={"marginTop":"8px","fontSize":"11px","color":"#e67e22"})], style={"flexGrow":"1","minWidth":"0","paddingLeft":"16px","paddingRight":"16px"}),

            html.Div([html.Div([html.Span("Saved layers", style={"fontWeight":"bold","fontSize":"14px","color":"#2c3e50"}), html.Button("Clear all", id="clear-layers-btn", n_clicks=0, style={"background":"none","border":"none","color":"#e74c3c","cursor":"pointer","fontSize":"11px","float":"right"})], style={"marginBottom":"8px","borderBottom":"1px solid #ecf0f1","paddingBottom":"6px"}), html.Div(id="layers-sidebar", children=make_sidebar(loaded_layers)), html.Button("⬇ Extract visible layers", id="extract-btn", n_clicks=0, style={"width":"100%","padding":"5px","backgroundColor":"#2980b9","color":"white","border":"none","borderRadius":"4px","cursor":"pointer","fontSize":"12px","marginBottom":"6px","marginTop":"8px"}), html.Div(id="extract-status", style={"fontSize":"11px","color":"#555","marginBottom":"8px"}), html.Details([
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
                    html.Label("Gridline colour (hex)", style={"fontSize":"11px"}),
                    dcc.Input(id="export-grid-color", type="text", value="#000000", debounce=True, style={"width":"100%","fontSize":"12px","marginBottom":"4px","boxSizing":"border-box"}),
                    html.Label("Marker edge colour (hex)", style={"fontSize":"11px"}),
                    dcc.Input(id="export-edge-color", type="text", value="#000000", debounce=True, style={"width":"100%","fontSize":"12px","marginBottom":"4px","boxSizing":"border-box"}),
                    html.Label("Background colour (hex, ignored if transparent)", style={"fontSize":"11px"}),
                    dcc.Input(id="export-bg-color", type="text", value="#ffffff", debounce=True, style={"width":"100%","fontSize":"12px","marginBottom":"4px","boxSizing":"border-box"}),
                    html.Div([html.Label("Width × height (px)", style={"fontSize":"11px", "flexGrow":"1"}), html.Button("Reset", id="export-size-reset", n_clicks=0, title="Reset to 1200 × 800 px", style={"fontSize":"10px", "padding":"1px 6px", "border":"1px solid #bdc3c7", "backgroundColor":"#ecf0f1", "borderRadius":"3px", "cursor":"pointer"})], style={"display":"flex", "alignItems":"center", "gap":"6px", "marginBottom":"2px"}),
                    html.Div([
                        dcc.Input(id="export-width", type="number", value=EXPORT_DEFAULT_WIDTH, min=300, max=8000, step=50, debounce=True, style={"width":"48%","fontSize":"12px","marginRight":"4%","boxSizing":"border-box"}),
                        dcc.Input(id="export-height", type="number", value=EXPORT_DEFAULT_HEIGHT, min=300, max=8000, step=50, debounce=True, style={"width":"48%","fontSize":"12px","boxSizing":"border-box"}),
                    ], style={"display":"flex","marginBottom":"6px"}),
                    html.Label("Save to", style={"fontSize":"11px"}),
                    dcc.RadioItems(id="export-destination", options=[{"label":html.Span(" Browser download", style={"fontSize":"12px"}),"value":"browser"},{"label":html.Span(" Entry figures/", style={"fontSize":"12px"}),"value":"server"}], value="browser", labelStyle={"display":"block","marginBottom":"2px"}, inputStyle={"cursor":"pointer","marginRight":"4px"}, style={"marginBottom":"6px","marginTop":"2px"}),
                    html.Button("📷 Save figure", id="export-btn", n_clicks=0, style={"width":"100%","padding":"5px","backgroundColor":"#27ae60","color":"white","border":"none","borderRadius":"4px","cursor":"pointer","fontSize":"12px","marginBottom":"4px"}),
                    html.Div(id="export-status", style={"fontSize":"11px","color":"#555","wordBreak":"break-all"})
                ], open=False, style={"marginBottom":"10px","borderTop":"1px solid #ecf0f1","paddingTop":"8px"}),

                html.Details([html.Summary("Boltz-2 structure prediction", style={"fontWeight":"600","fontSize":"12px","cursor":"pointer","marginBottom":"8px"}), html.Label("NVIDIA API key", style={"fontSize":"11px"}), dcc.Input(id="boltz-api-key", type="password", placeholder="nvapi-…", style={"width":"100%","fontSize":"11px","fontFamily":"monospace","marginBottom":"4px"}), html.Button("Check API key", id="boltz-check-key-btn", n_clicks=0, style={"width":"100%","padding":"4px","fontSize":"11px","marginBottom":"4px"}), html.Div(id="boltz-key-status", style={"fontSize":"11px","marginBottom":"8px","minHeight":"14px"}), dcc.Checklist(id="boltz-use-msa", options=[{"label":html.Span(" Use MSA (recommended)", style={"fontSize":"12px"}),"value":"on"}], value=["on"], style={"marginBottom":"8px"}), html.Label("Substrate SMILES (optional)", style={"fontSize":"11px"}), dcc.Textarea(id="boltz-smiles", placeholder="One SMILES per line", style={"width":"100%","fontSize":"11px","fontFamily":"monospace","resize":"vertical","minHeight":"58px"}), html.Label("Ligand label (optional)", style={"fontSize":"11px"}), dcc.Input(id="boltz-smiles-label", type="text", placeholder="e.g. UDP-Glc", style={"width":"100%","fontSize":"11px","marginBottom":"6px"}), html.Details([
                    html.Summary("Prediction parameters", style={"fontSize":"11px","cursor":"pointer","color":"#555","marginBottom":"6px"}),
                    html.Div([
                        html.Div([html.Label("Recycling steps", style={"fontSize":"11px","lineHeight":"22px"}), dcc.Input(id="boltz-recycling-steps", type="number", value=3, min=1, max=10, step=1, debounce=True, style={"width":"72px","fontSize":"11px","padding":"2px 4px","boxSizing":"border-box","textAlign":"right"})], style={"display":"grid","gridTemplateColumns":"1fr 76px","alignItems":"center","gap":"6px","marginBottom":"4px"}),
                        html.Div([html.Label("Sampling steps", style={"fontSize":"11px","lineHeight":"22px"}), dcc.Input(id="boltz-sampling-steps", type="number", value=200, min=10, max=500, step=10, debounce=True, style={"width":"72px","fontSize":"11px","padding":"2px 4px","boxSizing":"border-box","textAlign":"right"})], style={"display":"grid","gridTemplateColumns":"1fr 76px","alignItems":"center","gap":"6px","marginBottom":"4px"}),
                        html.Div([html.Label("Diffusion samples", style={"fontSize":"11px","lineHeight":"22px"}), dcc.Input(id="boltz-diffusion-samples", type="number", value=5, min=1, max=10, step=1, debounce=True, style={"width":"72px","fontSize":"11px","padding":"2px 4px","boxSizing":"border-box","textAlign":"right"})], style={"display":"grid","gridTemplateColumns":"1fr 76px","alignItems":"center","gap":"6px","marginBottom":"4px"}),
                        html.Div([html.Label("Step scale", style={"fontSize":"11px","lineHeight":"22px"}), dcc.Input(id="boltz-step-scale", type="number", value=1.638, min=0.1, max=5, step=0.001, debounce=True, style={"width":"72px","fontSize":"11px","padding":"2px 4px","boxSizing":"border-box","textAlign":"right"})], style={"display":"grid","gridTemplateColumns":"1fr 76px","alignItems":"center","gap":"6px","marginBottom":"2px"}),
                    ], style={"width":"100%"}),
                ], open=False, style={"marginBottom":"8px"}), dcc.Checklist(id="boltz-force-rerun", options=[{"label":html.Span(" Force re-run (ignore cache)", style={"fontSize":"12px"}),"value":"on"}], value=[]), html.Button("⚗ Send to Boltz-2", id="boltz-submit-btn", n_clicks=0, disabled=True, style={"width":"100%","padding":"6px","backgroundColor":"#8e44ad","color":"white","border":"none","borderRadius":"4px","cursor":"not-allowed","fontSize":"12px","marginBottom":"4px","opacity":"0.4"}), html.Div(id="boltz-submit-status", style={"fontSize":"11px","color":"#555","marginBottom":"6px"}), html.Div(id="boltz-summary", children=boltz_summary(_boltz_jobs), style={"fontSize":"11px","color":"#7f8c8d","marginBottom":"6px","fontStyle":"italic"}), html.Div(id="boltz-job-table", children=make_job_table(_boltz_jobs), style={"fontSize":"11px","overflowX":"auto"})], open=False, style={"marginTop":"12px","borderTop":"1px solid #ecf0f1","paddingTop":"8px"}),

                html.Details([html.Summary("RMSD analysis", style={"fontWeight":"600","fontSize":"12px","cursor":"pointer","marginBottom":"8px"}), html.Div(id="rmsd-struct-list", style={"fontSize":"11px","color":"#555","marginBottom":"4px"}), html.Label("Reference structure", style={"fontSize":"11px"}), dcc.Dropdown(id="rmsd-reference-select", placeholder="Select reference…", options=[], clearable=False, style={"fontSize":"11px","marginBottom":"4px"}), html.Div([html.Label("Reference rank", style={"fontSize":"11px","marginRight":"6px"}), dcc.Input(id="rmsd-ref-rank", type="number", value=0, min=0, step=1, debounce=True, style={"width":"55px","fontSize":"11px"})], style={"display":"flex","alignItems":"center","marginBottom":"8px"}), html.Label("Scope", style={"fontSize":"11px"}), dcc.RadioItems(id="rmsd-scope", options=[{"label":html.Span(" All completed apo structures", style={"fontSize":"12px"}),"value":"all"},{"label":html.Span(" Selected sequences only", style={"fontSize":"12px"}),"value":"selected"}], value="all", labelStyle={"display":"block","fontSize":"12px"}), html.Details([html.Summary("Advanced options — per-sequence rank", style={"fontSize":"11px","cursor":"pointer","color":"#555"}), html.Div(id="rmsd-rank-overrides", style={"fontSize":"11px"}), dcc.Store(id="rmsd-rank-store", data={}), dcc.Store(id="rmsd-seq-store", data=[])], open=False, style={"marginBottom":"8px"}), html.Label("Alignment method", style={"fontSize":"11px"}), dcc.RadioItems(id="rmsd-method", options=[{"label":html.Span(" Sequence-guided (super)", style={"fontSize":"12px"}),"value":"seq"},{"label":html.Span(" Structure-based (CE)", style={"fontSize":"12px"}),"value":"ce"},{"label":html.Span(" Both", style={"fontSize":"12px"}),"value":"both"}], value="seq", labelStyle={"display":"block","fontSize":"12px"}), html.Button("Calculate RMSDs", id="rmsd-calc-btn", n_clicks=0, style={"width":"100%","padding":"6px","backgroundColor":"#2980b9","color":"white","border":"none","borderRadius":"4px","cursor":"pointer","fontSize":"12px","marginBottom":"4px"}), html.Div(id="rmsd-status", style={"fontSize":"11px","color":"#555","marginBottom":"6px","minHeight":"14px"}), html.Div(id="rmsd-results-table", style={"fontSize":"11px","overflowX":"auto"})], open=False, style={"marginTop":"12px","borderTop":"1px solid #ecf0f1","paddingTop":"8px"})
            ], style={"width":"235px","minWidth":"215px","flexShrink":"0","borderLeft":"1px solid #ecf0f1","paddingLeft":"12px","overflowY":"auto","maxHeight":"90vh"}),
        ], style={"display":"flex","gap":"0","alignItems":"flex-start"})
    ], style={"padding":"20px","fontFamily":"sans-serif","maxWidth":"1900px","margin":"0 auto"})

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
        warn = reload_state()
        with _STATE_LOCK:
            valid_ids = _ann_df[_id_col].astype(str).tolist()
            subtitle = f"Entry: {ENTRY.stem} · {_ann_df.shape[0]:,} rows · {len(_STATE.coord_cols)} coordinate column(s)"
        cleaned_layers, layer_msg = layer_store.validate_layers(current_layers or [], valid_ids)
        if cleaned_layers != (current_layers or []):
            layer_store.write_layers(ENTRY.layers_path, cleaned_layers)
        status = " · ".join(x for x in [warn, layer_msg, "Reloaded datafile."] if x)
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
        return now, "👁" if now else "🚫"

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
        cont = [cont_lookup.get(d["col"], [slider_config(_ann_df[d["col"]])[0], slider_config(_ann_df[d["col"]])[1]]) for d in cont_ids]
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

    @app.callback(Output("latent-graph", "figure"), Output("point-count", "children"), Input("data-loaded-store", "data"), Input("reload-counter", "data"), Input("coord-cols-store", "data"), Input("alpha-slider", "value"), Input("point-size-slider", "value"), Input("bg-size-slider", "value"), Input("marker-size-slider", "value"), Input("marker-alpha-slider", "value"), Input("marker-mode", "value"), Input("wf-position", "value"), Input("color-mode", "value"), Input("fixed-color-store", "data"), Input("cont-color-col", "value"), Input("colormap-select", "value"), Input("colormap-reversed", "value"), Input("color-range-mode", "value"), Input({"type":"filter-enabled", "col":ALL}, "value"), Input({"type":"filter-enabled", "col":ALL}, "id"), Input({"type":"cont-slider", "col":ALL}, "value"), Input({"type":"cont-slider", "col":ALL}, "id"), Input({"type":"bool-filter", "col":ALL}, "value"), Input({"type":"bool-filter", "col":ALL}, "id"), Input({"type":"cat-filter", "col":ALL}, "value"), Input({"type":"cat-filter", "col":ALL}, "id"), Input({"type":"tag-filter", "col":ALL}, "value"), Input({"type":"tag-filter", "col":ALL}, "id"), Input("layers-store", "data"), Input("wf-visible-store", "data"), Input("id-search-enabled", "value"), Input("id-search-input", "n_blur"), Input("marker-symbol", "value"), Input("selection-store", "data"), Input("selection-color-store", "data"), State("id-search-input", "value"))
    def update_figure(_loaded, _reload, _coord, alpha, point_size, bg_size, marker_size, marker_alpha, marker_mode, wf_position, color_mode, fixed_color, cont_col, colormap, colormap_reversed, color_range_mode, enabled_values, enabled_ids, cont_values, cont_ids, bool_values, bool_ids, cat_values, cat_ids, tag_values, tag_ids, layers, wf_visible, id_enabled, _blur, wf_symbol, selection_ids, selection_color, id_raw):
        cont_conds, bool_conds, cat_conds, tag_conds, id_search_ids = parse_current_conditions(enabled_values, enabled_ids, cont_values, cont_ids, bool_values, bool_ids, cat_values, cat_ids, tag_values, tag_ids, id_enabled, id_raw)
        with _STATE_LOCK:
            fig, n_filtered, n_cov = make_figure(cont_conds, bool_conds, cat_conds, tag_conds, color_mode, fixed_color, cont_col, colormap or DEFAULT_COLORMAP, bool(colormap_reversed and "reversed" in colormap_reversed), color_range_mode or DEFAULT_COLOR_RANGE, alpha or DEFAULT_ALPHA, point_size or DEFAULT_POINT_SIZE, bg_size or DEFAULT_BG_SIZE, marker_size or DEFAULT_MARKER_SIZE, marker_alpha or DEFAULT_MARKER_ALPHA, marker_mode or DEFAULT_MARKER_MODE, wf_position or DEFAULT_WF_POSITION, True if wf_visible is None else wf_visible, layers or [], id_search_ids, wf_symbol or DEFAULT_SYMBOL, selection_ids or [], selection_color or DEFAULT_SELECTION_COLOR)
            total = len(_ann_df)
        n_vis = len([l for l in (layers or []) if l.get("visible", True)])
        count = f"Showing {n_filtered:,} / {total:,} filtered sequences · {n_cov:,} have coordinates here"
        if n_vis:
            count += f" · {n_vis} saved layer(s) visible"
        return fig, count

    @app.callback(Output("click-details", "children"), Output("selection-store", "data", allow_duplicate=True), Output("boltz-clicked-id-store", "data"), Input("latent-graph", "clickData"), State("selection-store", "data"), prevent_initial_call=True)
    def handle_click(click_data, current_selection):
        if not click_data or not click_data.get("points"):
            return "Click a point to show details here.", no_update, no_update
        cd = click_data["points"][0].get("customdata")
        if cd is None:
            return "No data for this point.", no_update, no_update
        seq_id = cd if isinstance(cd, str) else cd[0]
        sel = list(current_selection or [])
        if seq_id in sel:
            sel.remove(seq_id)
        else:
            sel.append(seq_id)
        with _STATE_LOCK:
            details = make_details_panel(seq_id)
        return details, sel, seq_id

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
        State("export-grid-color", "value"),
        State("export-edge-color", "value"),
        State("export-bg-color", "value"),
        State("export-width", "value"),
        State("export-height", "value"),
        State("export-destination", "value"),
        prevent_initial_call=True,
    )
    def export_figure(n, figure, fmt, dpi, legend_val, transparent_val,
                      axis_color, label_color, grid_color, edge_color, bg_color,
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
            grid_col = _valid_hex(grid_color, "#000000")
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
                    gridcolor=grid_col,
                    zerolinecolor=grid_col,
                ),
                yaxis=dict(
                    color=axis_col,
                    linecolor=axis_col,
                    tickcolor=axis_col,
                    tickfont=dict(color=label_col),
                    title=dict(font=dict(color=label_col)),
                    gridcolor=grid_col,
                    zerolinecolor=grid_col,
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

    @app.callback(Output("boltz-key-status", "children"), Output("boltz-key-status", "style"), Output("boltz-key-valid-store", "data"), Input("boltz-check-key-btn", "n_clicks"), State("boltz-api-key", "value"), prevent_initial_call=True)
    def check_key(n, key):
        ok, msg = boltz_backend.validate_api_key(key or "")
        color = "#27ae60" if ok else "#e74c3c"
        return ("✓ " if ok else "✗ ") + msg, {"fontSize":"11px","marginBottom":"8px","minHeight":"14px","color":color}, ok

    @app.callback(Output("boltz-key-valid-store", "data", allow_duplicate=True), Output("boltz-key-status", "children", allow_duplicate=True), Input("boltz-api-key", "value"), prevent_initial_call=True)
    def invalidate_key(_):
        return False, ""

    @app.callback(Output("boltz-submit-btn", "disabled"), Output("boltz-submit-btn", "style"), Input("boltz-key-valid-store", "data"), Input("boltz-clicked-id-store", "data"))
    def boltz_button_state(valid, clicked):
        base = {"width":"100%","padding":"6px","backgroundColor":"#8e44ad","color":"white","border":"none","borderRadius":"4px","fontSize":"12px","marginBottom":"4px"}
        enabled = bool(valid) and bool(clicked)
        return (not enabled, {**base, "cursor":"pointer" if enabled else "not-allowed", "opacity":"1" if enabled else "0.4"})

    @app.callback(Output("boltz-submit-status", "children"), Output("boltz-interval", "disabled"), Input("boltz-submit-btn", "n_clicks"), State("boltz-clicked-id-store", "data"), State("boltz-api-key", "value"), State("boltz-key-valid-store", "data"), State("boltz-use-msa", "value"), State("boltz-force-rerun", "value"), State("boltz-smiles", "value"), State("boltz-smiles-label", "value"), State("boltz-recycling-steps", "value"), State("boltz-sampling-steps", "value"), State("boltz-diffusion-samples", "value"), State("boltz-step-scale", "value"), prevent_initial_call=True)
    def submit_boltz(n, clicked_id, key, valid, use_msa, force_val, smiles, smiles_label, recycling, sampling, diffusion, scale):
        if not n:
            return no_update, no_update
        if not valid:
            return "API key not validated.", no_update
        if not clicked_id:
            return "Click a sequence first.", no_update
        rows = _ann_df[_ann_df[_id_col].astype(str) == str(clicked_id)]
        if rows.empty or COL_SEQ not in _ann_df.columns:
            return "Selected sequence not found or Sequence column missing.", no_update
        sequence = str(rows.iloc[0][COL_SEQ]).strip()
        params = boltz_backend.BoltzParams(_coerce_int(recycling, 3), _coerce_int(sampling, 200), _coerce_int(diffusion, 5), _coerce_float(scale, 1.638))
        try:
            job, should_run, msg = boltz_backend.submit_or_cache(ENTRY, str(clicked_id), sequence, api_key=key or "", smiles=smiles or "", smiles_label=smiles_label or "", use_msa=bool(use_msa and "on" in use_msa), params=params, force=bool(force_val and "on" in force_val))
            if should_run:
                fut = _boltz_executor.submit(boltz_backend.run_prediction, ENTRY, job["job_key"])
                _boltz_futures[job["job_key"]] = fut
            return f"{msg} {clicked_id}", False
        except Exception as exc:
            return f"Boltz submit failed: {exc}", no_update

    @app.callback(
        Output("boltz-job-table", "children"),
        Output("boltz-summary", "children"),
        Output("boltz-interval", "disabled", allow_duplicate=True),
        Output("reload-counter", "data", allow_duplicate=True),
        Output("filter-panel", "children", allow_duplicate=True),
        Output("col-settings-panel", "children", allow_duplicate=True),
        Output("filter-pending-store", "data", allow_duplicate=True),
        Input("boltz-interval", "n_intervals"),
        State("reload-counter", "data"),
        prevent_initial_call=True,
    )
    def poll_boltz(_n, counter):
        # Surface worker exceptions as error jobs. When a worker finishes, reload
        # the datafile and rebuild dynamic filter/colour controls so newly written
        # Boltz pTM/pLDDT columns immediately become available.
        completed_any = False
        for key, fut in list(_boltz_futures.items()):
            if fut.done():
                completed_any = True
                try:
                    fut.result()
                except Exception as exc:
                    job_store.update_job(ENTRY.jobs_path, "boltz", key, status="error", error=str(exc))
                _boltz_futures.pop(key, None)
        with _STATE_LOCK:
            before_cols = set(_ann_df.columns)
        warn = reload_state()
        with _STATE_LOCK:
            after_cols = set(_ann_df.columns)
        completed_any = completed_any or (before_cols != after_cols)
        jobs, _ = current_job_records()
        active = any(j.get("status") in {"queued", "msa", "predicting"} for j in jobs.values())
        pending = {"source": "boltz", "counter": (counter or 0) + 1} if completed_any else no_update
        cleared = html.Div() if completed_any else no_update
        return make_job_table(jobs), boltz_summary(jobs), not active, (counter or 0) + 1, cleared, cleared, pending

    @app.callback(Output("boltz-job-table", "children", allow_duplicate=True), Output("boltz-summary", "children", allow_duplicate=True), Input("boltz-submit-status", "children"), prevent_initial_call=True)
    def refresh_boltz_table(_):
        jobs, _r = current_job_records()
        return make_job_table(jobs), boltz_summary(jobs)

    @app.callback(
        Output("rmsd-reference-select", "options"),
        Output("rmsd-reference-select", "value"),
        Output("rmsd-rank-overrides", "children"),
        Output("rmsd-struct-list", "children"),
        Output("rmsd-rank-store", "data"),
        Output("rmsd-seq-store", "data"),
        Input("boltz-interval", "n_intervals"),
        Input("reload-counter", "data"),
        Input("boltz-submit-status", "children"),
        State("rmsd-rank-store", "data"),
        State("rmsd-reference-select", "value"),
    )
    def refresh_rmsd_structs(_i, _r, _s, rank_store, current_ref):
        seqs = rmsd_backend.list_apo_structures(ENTRY)
        if not seqs:
            empty = html.Div("No predicted apo structures yet.", style={"fontSize":"11px","color":"#aaa"})
            return [], None, empty, empty, {}, []

        rank_store = {str(k): int(v or 0) for k, v in (rank_store or {}).items()}
        seq_ids = [str(s["id"]) for s in seqs]
        max_by_id = {str(s["id"]): int(s.get("max_rank", 0) or 0) for s in seqs}

        # Preserve existing user-entered ranks when the structure list refreshes;
        # only initialize newly discovered structures to rank 0. Clamp vanished or
        # out-of-range values rather than resetting the whole store.
        clean_store = {}
        for sid in seq_ids:
            clean_store[sid] = max(0, min(int(rank_store.get(sid, 0) or 0), max_by_id[sid]))

        opts = [{"label": sid, "value": sid} for sid in seq_ids]
        selected_ref = current_ref if current_ref in seq_ids else seq_ids[0]

        overrides = []
        for s in seqs:
            sid = str(s["id"])
            max_rank = max_by_id[sid]
            overrides.append(html.Div([
                html.Span(sid, style={"fontFamily":"monospace","fontSize":"10px","flexGrow":"1","overflow":"hidden","textOverflow":"ellipsis","whiteSpace":"nowrap","maxWidth":"120px"}),
                html.Span(f"0-{max_rank}", style={"fontSize":"9px","color":"#aaa","marginRight":"4px"}),
                dcc.Input(
                    id={"type":"rmsd-rank", "sid":sid},
                    type="number",
                    value=clean_store[sid],
                    min=0, max=max_rank, step=1, debounce=True,
                    style={"width":"45px","fontSize":"10px","padding":"2px 4px","border":"1px solid #ddd","borderRadius":"3px"},
                ),
            ], style={"display":"grid","gridTemplateColumns":"minmax(0, 1fr) 28px 52px","alignItems":"center","gap":"4px","marginBottom":"3px"}))

        return opts, selected_ref, overrides, html.Div(f"{len(seqs)} apo structure(s) available.", style={"fontSize":"11px","color":"#555"}), clean_store, seq_ids

    @app.callback(
        Output("rmsd-rank-store", "data", allow_duplicate=True),
        Input({"type":"rmsd-rank", "sid":ALL}, "value"),
        State({"type":"rmsd-rank", "sid":ALL}, "id"),
        State("rmsd-rank-store", "data"),
        prevent_initial_call=True,
    )
    def update_rmsd_rank_store(values, ids, current_store):
        if values is None or ids is None:
            return no_update
        store = dict(current_store or {})
        changed = False
        for value, id_d in zip(values, ids):
            sid = str(id_d.get("sid", ""))
            if not sid:
                continue
            try:
                rank = max(0, int(value or 0))
            except (TypeError, ValueError):
                rank = 0
            if store.get(sid) != rank:
                store[sid] = rank
                changed = True
        return store if changed else no_update

    @app.callback(
        Output("rmsd-status", "children"),
        Output("rmsd-results-table", "children"),
        Output("filter-panel", "children", allow_duplicate=True),
        Output("reload-counter", "data", allow_duplicate=True),
        Output("filter-pending-store", "data", allow_duplicate=True),
        Input("rmsd-calc-btn", "n_clicks"),
        State("rmsd-reference-select", "value"),
        State("rmsd-ref-rank", "value"),
        State("rmsd-rank-store", "data"),
        State({"type":"rmsd-rank", "sid":ALL}, "value"),
        State({"type":"rmsd-rank", "sid":ALL}, "id"),
        State("rmsd-method", "value"),
        State("rmsd-scope", "value"),
        State("selection-store", "data"),
        State("reload-counter", "data"),
        prevent_initial_call=True,
    )
    def calc_rmsd(n, ref_id, ref_rank, rank_store, rank_values, rank_ids, method, scope, selection, counter):
        if not n:
            return no_update, no_update, no_update, no_update, no_update
        if not ref_id:
            return "Select a reference first.", no_update, no_update, no_update, no_update
        methods = ["seq", "ce"] if method == "both" else [method or "seq"]
        rank_store = dict(rank_store or {})
        # Use the live input values from the DOM at click time as the source of
        # truth, so clicking Calculate immediately after editing a rank does not
        # depend on a blur/Enter event updating rmsd-rank-store first.
        for value, id_d in zip(rank_values or [], rank_ids or []):
            sid = str((id_d or {}).get("sid", ""))
            if not sid:
                continue
            try:
                rank_store[sid] = max(0, int(value or 0))
            except (TypeError, ValueError):
                rank_store[sid] = 0
        query_ids = [str(x) for x in (selection or [])] if scope == "selected" else None
        if scope == "selected" and not query_ids:
            return "No selected sequences to compare.", no_update, no_update, no_update, no_update
        try:
            res = rmsd_backend.calculate_rmsds(ENTRY, ref_id, int(ref_rank or 0), query_ids=query_ids, query_rank_map=rank_store or {}, methods=methods)
            reload_state()
            status = f"Done — {res['n_new']} calculated, {res['n_cached']} from cache. Columns: {', '.join(res['columns'])}"
            pending = {"source": "rmsd", "counter": (counter or 0) + 1}
            return status, make_rmsd_results_table(res["results"]), html.Div(), (counter or 0) + 1, pending
        except Exception as exc:
            return f"RMSD failed: {exc}", no_update, no_update, no_update, no_update

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
    app.run(debug=True, port=args.port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
