"""Boltz-2 backend for the SSE visualizer.

This module contains the network/MSA/file/datafile work. Dash code should only
call these functions and render the returned/persisted job records.
"""

from __future__ import annotations

import csv
import hashlib
import io
import os
import re
import tarfile
import threading
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import pandas as pd
import requests

from .common import COL_ID, COL_SEQ, TYPE_LABEL, merge_columns, read_datafile
from .jobs import update_job, upsert_job, utc_now
from .visualizer_state import EntryContext

BOLTZ2_URL = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
MSA_API_BASE = "https://api.colabfold.com"
MSA_SUBMIT_INTERVAL = 3.0
MSA_TIMEOUT = 600
MSA_MAX_SEQS = 1000

BOLTZ_LOG_FIELDS = [
    "run_id", "timestamp", "sequence_id", "sequence",
    "kind", "smiles", "smiles_label", "smiles_hash",
    "status", "ptm", "plddt", "msa_used", "cif_paths", "error",
    "recycling_steps", "sampling_steps", "diffusion_samples", "step_scale",
]

_APO_SCALAR_COLS = [
    "boltz_apo_ptm",
    "boltz_apo_plddt",
]


def _safe_col_token(value: str, fallback: str = "ligand") -> str:
    """Return a compact, column-safe token for ligand-specific Boltz columns."""
    token = re.sub(r"[^0-9A-Za-z]+", "_", str(value or "").strip()).strip("_").lower()
    return token[:40] or fallback


def scalar_columns_for_job(job: dict) -> list[str]:
    """Scalar datafile columns written for a Boltz job.

    Apo predictions have one canonical column set per sequence. Holo predictions
    can be repeated with different ligands, so their columns include a ligand
    label when supplied, otherwise the SMILES hash.
    """
    if job.get("kind") == "holo":
        tag = _safe_col_token(job.get("smiles_label") or job.get("smiles_hash") or "ligand")
        prefix = f"boltz_holo_{tag}"
        return [f"{prefix}_ptm", f"{prefix}_plddt"]
    return list(_APO_SCALAR_COLS)

_msa_submit_lock = threading.Lock()
_msa_last_submit = 0.0


@dataclass
class BoltzParams:
    recycling_steps: int = 3
    sampling_steps: int = 200
    diffusion_samples: int = 5
    step_scale: float = 1.638

    def as_payload(self) -> dict:
        return {
            "recycling_steps": int(self.recycling_steps),
            "sampling_steps": int(self.sampling_steps),
            "diffusion_samples": int(self.diffusion_samples),
            "step_scale": float(self.step_scale),
        }

    def as_record(self) -> dict:
        return self.as_payload()


def smiles_hash(smiles: str) -> str:
    lines = [s.strip() for s in (smiles or "").splitlines() if s.strip()]
    if not lines:
        return "apo"
    return hashlib.md5("\n".join(lines).encode("utf-8")).hexdigest()[:10]


def sanitize_label(label: str, max_len: int = 40) -> str:
    import re
    cleaned = re.sub(r"[\s/\\'\"<>:|?*]+", "-", (label or "").strip())
    cleaned = cleaned.strip("-. ")
    return cleaned[:max_len]


def holo_folder_name(seq_id: str, smiles_label: str, h: str) -> str:
    lbl = sanitize_label(smiles_label)
    return f"{seq_id}__{lbl}__{h}" if lbl else f"{seq_id}__{h}"


def boltz_log_path(entry: EntryContext) -> Path:
    return entry.logs_dir / "boltz_log.csv"


def ensure_boltz_log(entry: EntryContext) -> None:
    p = boltz_log_path(entry)
    p.parent.mkdir(parents=True, exist_ok=True)
    if not p.exists():
        with open(p, "w", newline="", encoding="utf-8") as fh:
            csv.DictWriter(fh, fieldnames=BOLTZ_LOG_FIELDS).writeheader()


def append_boltz_log(entry: EntryContext, job: dict) -> None:
    ensure_boltz_log(entry)
    params = job.get("params", {})
    row = {
        "run_id": job.get("run_id", ""),
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "sequence_id": job.get("sequence_id", job.get("id", "")),
        "sequence": job.get("sequence", ""),
        "kind": job.get("kind", "apo"),
        "smiles": job.get("smiles", ""),
        "smiles_label": job.get("smiles_label", ""),
        "smiles_hash": job.get("smiles_hash", "apo"),
        "status": job.get("status", ""),
        "ptm": job.get("ptm", ""),
        "plddt": job.get("plddt", ""),
        "msa_used": job.get("msa_used", ""),
        "cif_paths": "; ".join(str(p) for p in job.get("cif_paths", [])),
        "error": job.get("error", ""),
        "recycling_steps": params.get("recycling_steps", ""),
        "sampling_steps": params.get("sampling_steps", ""),
        "diffusion_samples": params.get("diffusion_samples", ""),
        "step_scale": params.get("step_scale", ""),
    }
    with open(boltz_log_path(entry), "a", newline="", encoding="utf-8") as fh:
        csv.DictWriter(fh, fieldnames=BOLTZ_LOG_FIELDS).writerow(row)


def validate_api_key(api_key: str) -> tuple[bool, str]:
    if not api_key or not api_key.strip():
        return False, "Enter a key first."
    payload = {
        "polymers": [{"id": "A", "molecule_type": "protein", "sequence": "MAST"}],
        "recycling_steps": 1,
        "sampling_steps": 1,
        "diffusion_samples": 1,
        "step_scale": 1.638,
    }
    headers = {"Authorization": f"Bearer {api_key.strip()}", "Content-Type": "application/json"}
    try:
        r = requests.post(BOLTZ2_URL, headers=headers, json=payload, timeout=30)
        if r.status_code in (200, 422):
            return True, "Valid"
        if r.status_code in (401, 403):
            return False, "Invalid key."
        return False, f"HTTP {r.status_code}."
    except requests.exceptions.Timeout:
        return False, "Timed out."
    except Exception as exc:
        return False, str(exc)


def seq_hash(seq: str) -> str:
    return hashlib.md5(seq.strip().upper().encode()).hexdigest()[:12]


def trim_msa(msa: str) -> str:
    out, count = [], 0
    for line in msa.splitlines():
        if line.startswith(">"):
            if count >= MSA_MAX_SEQS:
                break
            count += 1
        out.append(line)
    return "\n".join(out)


def generate_msa(entry: EntryContext, sequence: str) -> Optional[str]:
    global _msa_last_submit
    cache_file = entry.msa_cache_dir / f"{seq_hash(sequence)}.a3m"
    if cache_file.exists():
        return cache_file.read_text(encoding="utf-8")
    with _msa_submit_lock:
        wait = MSA_SUBMIT_INTERVAL - (time.time() - _msa_last_submit)
        if wait > 0:
            time.sleep(wait)
        _msa_last_submit = time.time()
    hdrs = {"User-Agent": "sse-boltz2/1.0"}
    try:
        r = requests.post(
            f"{MSA_API_BASE}/ticket/msa",
            data={"q": f">query\n{sequence}\n", "mode": "env"},
            headers=hdrs,
            timeout=30,
        )
        if r.status_code == 429:
            time.sleep(int(r.headers.get("Retry-After", 90)))
            r = requests.post(
                f"{MSA_API_BASE}/ticket/msa",
                data={"q": f">query\n{sequence}\n", "mode": "env"},
                headers=hdrs,
                timeout=30,
            )
        r.raise_for_status()
        ticket = r.json()["id"]
        deadline = time.time() + MSA_TIMEOUT
        while time.time() < deadline:
            time.sleep(10)
            rp = requests.get(f"{MSA_API_BASE}/ticket/{ticket}", headers=hdrs, timeout=15)
            rp.raise_for_status()
            status = rp.json().get("status", "")
            if status == "COMPLETE":
                break
            if status in ("ERROR", "UNKNOWN"):
                return None
        else:
            return None
        rd = requests.get(f"{MSA_API_BASE}/result/download/{ticket}", headers=hdrs, timeout=120)
        rd.raise_for_status()
        msa_text = None
        with tarfile.open(fileobj=io.BytesIO(rd.content), mode="r:gz") as tar:
            for member in tar.getmembers():
                if member.name.endswith(".a3m"):
                    f = tar.extractfile(member)
                    if f:
                        msa_text = f.read().decode("utf-8")
                    break
        if msa_text:
            cache_file.parent.mkdir(parents=True, exist_ok=True)
            cache_file.write_text(msa_text, encoding="utf-8")
            return trim_msa(msa_text)
        return None
    except Exception:
        return None


def build_payload(sequence: str, msa: Optional[str], params: BoltzParams, smiles: str = "") -> dict:
    polymer = {"id": "A", "molecule_type": "protein", "sequence": sequence}
    if msa:
        polymer["msa"] = {"colabfold": {"a3m": {"format": "a3m", "alignment": msa}}}
    payload = {"polymers": [polymer], **params.as_payload()}
    smiles_list = [s.strip() for s in (smiles or "").splitlines() if s.strip()]
    if smiles_list:
        chain_ids = [chr(ord("B") + i) for i in range(len(smiles_list))]
        payload["ligands"] = [
            {"id": cid, "molecule_type": "smiles", "smiles": smi}
            for cid, smi in zip(chain_ids, smiles_list)
        ]
    return payload


def structure_dir(entry: EntryContext, job: dict) -> Path:
    seq_id = job["sequence_id"]
    if job.get("kind") == "holo":
        folder = holo_folder_name(seq_id, job.get("smiles_label", ""), job.get("smiles_hash", ""))
        return entry.structures_dir / "holo" / folder
    return entry.structures_dir / "apo" / seq_id


def check_cache(entry: EntryContext, seq_id: str, sequence: str, smiles: str = "") -> Optional[dict]:
    p = boltz_log_path(entry)
    if not p.exists():
        return None
    want_hash = smiles_hash(smiles)
    try:
        with open(p, newline="", encoding="utf-8") as fh:
            match = None
            for row in csv.DictReader(fh):
                row_hash = row.get("smiles_hash") or smiles_hash(row.get("smiles", ""))
                if (row.get("sequence_id") == seq_id and row.get("sequence") == sequence
                        and row_hash == want_hash and row.get("status") in ("done", "cached")):
                    match = row
            return match
    except Exception:
        return None


def make_job_record(sequence_id: str, sequence: str, *, api_key: str, smiles: str = "",
                    smiles_label: str = "", use_msa: bool = True,
                    params: Optional[BoltzParams] = None) -> dict:
    params = params or BoltzParams()
    smiles = (smiles or "").strip()
    kind = "holo" if smiles else "apo"
    h = smiles_hash(smiles)
    job_key = sequence_id if kind == "apo" else f"{sequence_id}::{h}"
    return {
        "job_key": job_key,
        "sequence_id": sequence_id,
        "sequence": sequence,
        "kind": kind,
        "smiles": smiles,
        "smiles_label": (smiles_label or "").strip() if kind == "holo" else "",
        "smiles_hash": h,
        "status": "queued",
        "ptm": None,
        "plddt": None,
        "msa_used": None,
        "cif_paths": [],
        "error": "",
        "params": params.as_record(),
        "api_key": api_key.strip(),
        "use_msa": bool(use_msa),
        "submitted_utc": utc_now(),
    }


def _datafile_scalar(value, *, digits: int | None = None) -> str:
    """Return a TSV-safe string for a scalar value written into the SSE datafile.

    read_datafile() returns string/object columns. Assigning Python floats into
    those columns can raise on newer pandas versions (or when warnings are
    treated strictly). Keep all Boltz writeback values as strings; the
    visualizer will still classify numeric-looking strings as continuous after
    reload.
    """
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, (int, float)) and digits is not None:
        return f"{float(value):.{digits}f}"
    return str(value)


def _obsolete_boltz_datafile_columns(columns) -> list[str]:
    """Status/MSA bookkeeping columns that should no longer live in the SSE datafile.

    Detailed job state remains in jobs.json and boltz_log.csv. The datafile only
    keeps prediction scores that are useful for colouring/filtering.
    """
    obsolete = []
    for col in columns:
        name = str(col)
        if name in {"boltz_apo_status", "boltz_apo_msa_used"}:
            obsolete.append(name)
        elif re.match(r"^boltz_holo_.+_(status|msa_used)$", name):
            obsolete.append(name)
    return obsolete


def write_boltz_scalars(entry: EntryContext, job: dict) -> None:
    """Append/replace Boltz score columns for one sequence in the SSE datafile.

    The SSE datafile only receives pTM and pLDDT columns. Status/MSA bookkeeping
    stays in jobs.json and boltz_log.csv, and obsolete status/MSA columns from
    older app versions are removed during writeback.
    """
    cols = scalar_columns_for_job(job)
    df, _types = read_datafile(entry.datafile_path)
    base = pd.DataFrame({COL_ID: df[COL_ID].astype(str)})
    for col in cols:
        base[col] = df[col].astype(str) if col in df.columns else ""

    m = base[COL_ID] == str(job["sequence_id"])
    base.loc[m, cols[0]] = _datafile_scalar(job.get("ptm"), digits=6)
    base.loc[m, cols[1]] = _datafile_scalar(job.get("plddt"), digits=3)

    new_types = {col: TYPE_LABEL for col in cols}
    merge_columns(
        entry.datafile_path,
        base[[COL_ID] + cols],
        new_types,
        manifest_path=entry.manifest_path if entry.manifest_path.exists() else None,
        provenance_source="sse",
        tool="sse_boltz",
        notes="Boltz-2 pTM/pLDDT score columns written by the SSE visualizer.",
        force=True,
        drop_columns=_obsolete_boltz_datafile_columns(df.columns),
    )


def write_apo_scalars(entry: EntryContext, job: dict) -> None:
    """Backward-compatible alias for older callers."""
    write_boltz_scalars(entry, job)


def run_prediction(entry: EntryContext, job_key: str) -> dict:
    """Run one persisted job to terminal state. Intended for a background thread."""
    from .jobs import read_jobs

    jobs = read_jobs(entry.jobs_path, mark_stale=False)
    job = dict(jobs.get("boltz", {}).get(job_key, {}))
    if not job:
        return {"job_key": job_key, "status": "error", "error": "Job not found."}

    update_job(entry.jobs_path, "boltz", job_key, status="msa" if job.get("use_msa") else "predicting")
    msa = None
    if job.get("use_msa"):
        msa = generate_msa(entry, job["sequence"])
        if msa is None:
            time.sleep(90)
            msa = generate_msa(entry, job["sequence"])
    update_job(entry.jobs_path, "boltz", job_key, status="predicting", msa_used=msa is not None)

    params = BoltzParams(**{k: job.get("params", {}).get(k, getattr(BoltzParams(), k)) for k in BoltzParams().__dataclass_fields__})
    payload = build_payload(job["sequence"], msa, params, job.get("smiles", ""))
    headers = {"Authorization": f"Bearer {job.get('api_key', '')}", "Content-Type": "application/json"}
    out_dir = structure_dir(entry, job)
    out_dir.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now().strftime("%Y%m%d%H%M%S")

    final = {**job, "run_id": run_id, "msa_used": msa is not None}
    for attempt in range(3):
        try:
            r = requests.post(BOLTZ2_URL, headers=headers, json=payload, timeout=1800)
            if r.status_code == 200:
                data = r.json()
                cif_paths = []
                for i, structure in enumerate(data.get("structures", [])):
                    cif_path = out_dir / f"{job['sequence_id']}_Rank_{i}.cif"
                    cif_path.write_text(structure["structure"], encoding="utf-8")
                    cif_paths.append(str(cif_path.resolve()))
                final.update({
                    "status": "done",
                    "ptm": data.get("ptm_scores", [None])[0],
                    "plddt": data.get("complex_plddt_scores", [None])[0],
                    "cif_paths": cif_paths,
                    "error": "",
                })
                break
            if r.status_code == 422 and msa is not None:
                msa = None
                final["msa_used"] = False
                payload = build_payload(job["sequence"], None, params, job.get("smiles", ""))
                continue
            if r.status_code == 429:
                time.sleep(60 * (attempt + 1))
                continue
            time.sleep(10)
            final.update({"status": "error", "error": f"HTTP {r.status_code}: {r.text[:300]}"})
            break
        except requests.exceptions.Timeout:
            time.sleep(10)
        except Exception as exc:
            final.update({"status": "error", "error": str(exc)})
            break
    else:
        final.update({"status": "error", "error": "All retries exhausted"})

    update_job(entry.jobs_path, "boltz", job_key, **final)
    append_boltz_log(entry, final)
    if final.get("status") in {"done", "error"}:
        try:
            write_boltz_scalars(entry, final)
        except Exception as exc:
            update_job(entry.jobs_path, "boltz", job_key, error=f"Datafile writeback failed: {exc}")
    return final


def submit_or_cache(entry: EntryContext, sequence_id: str, sequence: str, *, api_key: str,
                    smiles: str = "", smiles_label: str = "", use_msa: bool = True,
                    params: Optional[BoltzParams] = None, force: bool = False) -> tuple[dict, bool, str]:
    """Create/update a job record. Returns (job, should_run, message)."""
    job = make_job_record(sequence_id, sequence, api_key=api_key, smiles=smiles,
                          smiles_label=smiles_label, use_msa=use_msa, params=params)
    if not force:
        cached = check_cache(entry, sequence_id, sequence, smiles)
        if cached:
            job.update({
                "status": "cached",
                "ptm": float(cached["ptm"]) if cached.get("ptm") else None,
                "plddt": float(cached["plddt"]) if cached.get("plddt") else None,
                "msa_used": cached.get("msa_used") == "True",
                "cif_paths": [p for p in cached.get("cif_paths", "").split("; ") if p],
                "error": "",
            })
            upsert_job(entry.jobs_path, "boltz", job["job_key"], job)
            write_boltz_scalars(entry, job)
            return job, False, "Loaded from cache."
    upsert_job(entry.jobs_path, "boltz", job["job_key"], job)
    return job, True, "Queued."
