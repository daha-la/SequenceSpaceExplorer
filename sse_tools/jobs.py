"""Persistent job-state helpers for visualizer-launched Boltz/RMSD work."""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

TERMINAL_STATUSES = {"done", "cached", "error", "interrupted"}
NONTERMINAL_STATUSES = {"queued", "msa", "predicting", "running"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _empty() -> dict:
    return {"version": 1, "boltz": {}, "rmsd": {}}


def read_jobs(path, *, mark_stale: bool = True) -> dict:
    p = Path(path)
    if not p.exists():
        return _empty()
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            data = _empty()
    except Exception:
        data = _empty()
    data.setdefault("version", 1)
    data.setdefault("boltz", {})
    data.setdefault("rmsd", {})
    if mark_stale:
        changed = False
        for section in ("boltz", "rmsd"):
            for key, rec in list(data.get(section, {}).items()):
                if isinstance(rec, dict) and rec.get("status") in NONTERMINAL_STATUSES:
                    rec["status"] = "interrupted"
                    rec["error"] = rec.get("error") or "Job was non-terminal when the app last stopped."
                    rec["updated_utc"] = utc_now()
                    changed = True
        if changed:
            write_jobs(path, data)
    return data


def write_jobs(path, jobs: dict) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    jobs = jobs or _empty()
    jobs.setdefault("version", 1)
    jobs.setdefault("boltz", {})
    jobs.setdefault("rmsd", {})
    tmp = Path(str(p) + ".tmp")
    tmp.write_text(json.dumps(jobs, indent=2), encoding="utf-8")
    os.replace(tmp, p)


def upsert_job(path, section: str, key: str, record: dict) -> dict:
    jobs = read_jobs(path, mark_stale=False)
    jobs.setdefault(section, {})[key] = {**record, "updated_utc": utc_now()}
    write_jobs(path, jobs)
    return jobs


def update_job(path, section: str, key: str, **fields) -> dict:
    jobs = read_jobs(path, mark_stale=False)
    rec = jobs.setdefault(section, {}).setdefault(key, {"job_key": key})
    rec.update(fields)
    rec["updated_utc"] = utc_now()
    write_jobs(path, jobs)
    return rec
