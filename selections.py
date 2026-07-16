"""Selection-cache helpers shared by the visualizer and the Boltz-2 pipeline.

The visualizer lets a user pick points in the plot and export the underlying
sequences to a cache on disk. The Boltz-2 pipeline module (scripts/sse_boltz.py)
later imports that cache and runs structure prediction + RMSD. This module owns
the on-disk schema so the two sides never drift.

Each cache is a timestamped JSON file under entries/<entry>/selections/:

    {
      "version": 1,
      "entry": "<stem>",
      "created_utc": "<iso8601>",
      "count": <n>,
      "sequences": [{"id": "<seq id>", "sequence": "<AA>"}, ...]
    }
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Optional

SELECTION_VERSION = 1
SELECTION_PREFIX = "selection_"
SELECTION_SUFFIX = ".json"


def _now_stamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def write_selection(selections_dir: Path, entry_stem: str,
                    sequences: Iterable[dict]) -> Path:
    """Write a timestamped selection cache and return its path.

    `sequences` is an iterable of {"id", "sequence"} dicts. Entries missing an
    id or a sequence are dropped, and duplicate ids keep their first occurrence.
    """
    selections_dir = Path(selections_dir)
    selections_dir.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    cleaned: list[dict] = []
    for item in sequences:
        sid = str(item.get("id", "")).strip()
        seq = str(item.get("sequence", "")).strip()
        if not sid or not seq or sid in seen:
            continue
        seen.add(sid)
        cleaned.append({"id": sid, "sequence": seq})

    payload = {
        "version": SELECTION_VERSION,
        "entry": entry_stem,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "count": len(cleaned),
        "sequences": cleaned,
    }
    path = selections_dir / f"{SELECTION_PREFIX}{_now_stamp()}{SELECTION_SUFFIX}"
    # Avoid clobbering a same-second export by appending a disambiguating suffix.
    n = 1
    while path.exists():
        path = selections_dir / f"{SELECTION_PREFIX}{_now_stamp()}_{n}{SELECTION_SUFFIX}"
        n += 1
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def list_selections(selections_dir: Path) -> list[Path]:
    """Return selection cache paths, newest first (by modified time)."""
    selections_dir = Path(selections_dir)
    if not selections_dir.exists():
        return []
    files = [p for p in selections_dir.glob(f"{SELECTION_PREFIX}*{SELECTION_SUFFIX}")
             if p.is_file()]
    return sorted(files, key=lambda p: p.stat().st_mtime, reverse=True)


def resolve_selection(selections_dir: Path, name: Optional[str] = None) -> Optional[Path]:
    """Resolve a selection reference to a concrete path.

    `name` may be an absolute/relative path, a bare filename inside
    selections_dir, or None (returns the most recent cache). Returns None when
    nothing matches.
    """
    selections_dir = Path(selections_dir)
    if name:
        candidate = Path(name)
        if candidate.is_file():
            return candidate
        candidate = selections_dir / name
        if candidate.is_file():
            return candidate
        if not name.endswith(SELECTION_SUFFIX):
            candidate = selections_dir / f"{name}{SELECTION_SUFFIX}"
            if candidate.is_file():
                return candidate
        return None
    picks = list_selections(selections_dir)
    return picks[0] if picks else None


def read_selection(path: Path) -> dict:
    """Load a selection cache, returning the parsed payload dict."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Malformed selection cache: {path}")
    data.setdefault("sequences", [])
    return data
