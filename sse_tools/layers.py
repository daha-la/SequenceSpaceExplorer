"""Persistent saved-layer helpers for the SSE visualizer."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_layers(path) -> list:
    p = Path(path)
    if not p.exists():
        return []
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return []
    return data if isinstance(data, list) else []


def write_layers(path, layers: list) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(p) + ".tmp")
    tmp.write_text(json.dumps(layers or [], indent=2), encoding="utf-8")
    os.replace(tmp, p)


def validate_layers(layers: list, valid_ids: Iterable[str]) -> tuple[list, str]:
    """Keep layers loadable after datafile edits; remove IDs that no longer exist."""
    valid = set(str(x) for x in valid_ids)
    changed = False
    missing_total = 0
    out = []
    for layer in layers or []:
        layer = dict(layer)
        ids = [str(x) for x in layer.get("ids", [])]
        kept = [x for x in ids if x in valid]
        missing_total += len(ids) - len(kept)
        if kept != ids:
            layer["ids"] = kept
            changed = True
        out.append(layer)
    msg = f"Removed {missing_total} missing sequence ID(s) from saved layers." if missing_total else ""
    return out, msg
