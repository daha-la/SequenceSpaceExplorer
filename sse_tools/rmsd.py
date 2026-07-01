"""RMSD backend for the SSE visualizer."""

from __future__ import annotations

import csv
import math
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .common import COL_ID, TYPE_LABEL, merge_columns, read_datafile
from .jobs import upsert_job, utc_now
from .visualizer_state import EntryContext

RMSD_LOG_FIELDS = [
    "reference_id", "query_id", "reference_rank", "query_rank",
    "n_aligned_residues", "rmsd", "method", "timestamp",
]


def rmsd_log_path(entry: EntryContext) -> Path:
    return entry.logs_dir / "rmsd_log.csv"


def ensure_rmsd_log(entry: EntryContext) -> None:
    p = rmsd_log_path(entry)
    p.parent.mkdir(parents=True, exist_ok=True)
    if not p.exists():
        with open(p, "w", newline="", encoding="utf-8") as fh:
            csv.DictWriter(fh, fieldnames=RMSD_LOG_FIELDS).writeheader()


def append_rmsd_rows(entry: EntryContext, rows: list[dict]) -> None:
    ensure_rmsd_log(entry)
    with open(rmsd_log_path(entry), "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=RMSD_LOG_FIELDS)
        for row in rows:
            w.writerow(row)


def read_rmsd_log(entry: EntryContext) -> pd.DataFrame:
    p = rmsd_log_path(entry)
    if not p.exists():
        return pd.DataFrame(columns=RMSD_LOG_FIELDS)
    try:
        return pd.read_csv(p, dtype=str)
    except Exception:
        return pd.DataFrame(columns=RMSD_LOG_FIELDS)


def list_apo_structures(entry: EntryContext) -> list[dict]:
    root = entry.structures_dir / "apo"
    if not root.exists():
        return []
    out = []
    for seq_dir in sorted(root.iterdir()):
        if not seq_dir.is_dir():
            continue
        cifs = sorted(seq_dir.glob("*_Rank_*.cif"))
        if cifs:
            out.append({"id": seq_dir.name, "cif_paths": [str(c) for c in cifs], "max_rank": len(cifs) - 1})
    return out


def cif_path_for(entry: EntryContext, seq_id: str, rank: int) -> Optional[Path]:
    p = entry.structures_dir / "apo" / seq_id / f"{seq_id}_Rank_{rank}.cif"
    return p if p.exists() else None


def rmsd_column(reference_id: str, reference_rank: int, method: str) -> str:
    safe = str(reference_id).replace(" ", "_")
    return f"RMSD_vs_{safe}_r{int(reference_rank)}_{method}"


def cached(entry: EntryContext, ref_id, qry_id, ref_rank, qry_rank, method) -> Optional[tuple[float, int]]:
    df = read_rmsd_log(entry)
    if df.empty:
        return None
    if "method" not in df.columns:
        df["method"] = "seq"
    df["method"] = df["method"].fillna("seq").replace("", "seq")
    mask = (
        (df["reference_id"] == str(ref_id))
        & (df["query_id"] == str(qry_id))
        & (df["reference_rank"] == str(ref_rank))
        & (df["query_rank"] == str(qry_rank))
        & (df["method"] == method)
    )
    hits = df[mask]
    if hits.empty:
        return None
    row = hits.iloc[-1]
    try:
        rmsd = float(row["rmsd"])
        n_al = int(row["n_aligned_residues"])
        return (rmsd, n_al) if n_al > 0 and not np.isnan(rmsd) else None
    except Exception:
        return None


def parse_ca_coords(cif_path: str) -> tuple[list[str], np.ndarray]:
    three_to_one = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
        "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
        "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    }
    seq, coords, seen = [], [], {}
    idx = {"group":None, "atom_id":None, "comp_id":None, "seq_id":None, "chain":None, "x":None, "y":None, "z":None}
    in_atom_loop = False
    col_order = []
    try:
        with open(cif_path, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip()
                if line.strip() == "loop_":
                    col_order = []
                    in_atom_loop = False
                    continue
                if line.startswith("_atom_site."):
                    tag = line.strip().split(".")[1]
                    col_order.append(tag)
                    in_atom_loop = True
                    continue
                if not in_atom_loop:
                    continue
                if line.startswith("_") or line.startswith("#") or not line.strip():
                    in_atom_loop = False
                    col_order = []
                    continue
                if idx["x"] is None and col_order:
                    tag_map = {t: i for i, t in enumerate(col_order)}
                    idx["group"] = tag_map.get("group_PDB")
                    idx["atom_id"] = tag_map.get("label_atom_id")
                    idx["comp_id"] = tag_map.get("label_comp_id")
                    idx["seq_id"] = tag_map.get("label_seq_id")
                    idx["chain"] = tag_map.get("label_asym_id")
                    idx["x"] = tag_map.get("Cartn_x")
                    idx["y"] = tag_map.get("Cartn_y")
                    idx["z"] = tag_map.get("Cartn_z")
                parts = line.split()
                if not parts:
                    continue
                if idx["group"] is not None and idx["group"] < len(parts) and parts[idx["group"]] != "ATOM":
                    continue
                if idx["atom_id"] is not None and idx["atom_id"] < len(parts):
                    if parts[idx["atom_id"]] != "CA":
                        continue
                else:
                    continue
                chain = parts[idx["chain"]] if idx["chain"] is not None else "A"
                seq_id = parts[idx["seq_id"]] if idx["seq_id"] is not None else "0"
                key = (chain, seq_id)
                if key in seen:
                    continue
                seen[key] = True
                try:
                    x = float(parts[idx["x"]]); y = float(parts[idx["y"]]); z = float(parts[idx["z"]])
                except (TypeError, ValueError, IndexError):
                    continue
                res = parts[idx["comp_id"]] if idx["comp_id"] is not None else "UNK"
                seq.append(three_to_one.get(res, "X"))
                coords.append([x, y, z])
    except Exception:
        pass
    return seq, np.array(coords, dtype=float) if coords else np.empty((0, 3))


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    P = P - P.mean(axis=0)
    Q = Q - Q.mean(axis=0)
    H = Q.T @ P
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt(((P - (Q @ R.T)) ** 2).sum(axis=1).mean()))


def compute_rmsd_seq(ref_cif: str, qry_cif: str) -> tuple[float, int]:
    from Bio.Align import PairwiseAligner
    ref_seq, ref_coords = parse_ca_coords(ref_cif)
    qry_seq, qry_coords = parse_ca_coords(qry_cif)
    if len(ref_seq) == 0 or len(qry_seq) == 0:
        return float("nan"), 0
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    try:
        best = next(iter(aligner.align("".join(ref_seq), "".join(qry_seq))))
    except StopIteration:
        return float("nan"), 0
    ref_idx, qry_idx = [], []
    for (r0, r1), (q0, q1) in zip(best.aligned[0], best.aligned[1]):
        for k in range(min(r1 - r0, q1 - q0)):
            ref_idx.append(r0 + k)
            qry_idx.append(q0 + k)
    if len(ref_idx) < 3:
        return float("nan"), 0
    P = ref_coords[ref_idx]
    Q = qry_coords[qry_idx]
    mask = np.ones(len(P), dtype=bool)
    for _ in range(10):
        Pm, Qm = P[mask], Q[mask]
        Pc, Qc = Pm - Pm.mean(0), Qm - Qm.mean(0)
        H = Qc.T @ Pc
        U, S, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        Pf = P - Pm.mean(0)
        Qf = Q - Qm.mean(0)
        dists = np.sqrt(((Pf - Qf @ R.T) ** 2).sum(axis=1))
        cutoff = dists[mask].mean() + 2.0 * dists[mask].std()
        new_mask = dists <= cutoff
        if new_mask.sum() == mask.sum():
            break
        mask = new_mask
    if int(mask.sum()) < 3:
        return float("nan"), 0
    return kabsch_rmsd(P[mask], Q[mask]), int(mask.sum())


def compute_rmsd_ce(ref_cif: str, qry_cif: str) -> tuple[float, int]:
    from Bio.PDB.cealign import CEAligner
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser(QUIET=True)
    try:
        ref_struct = parser.get_structure("ref", ref_cif)
        qry_struct = parser.get_structure("qry", qry_cif)
        aligner = CEAligner()
        aligner.set_reference(ref_struct)
        aligner.align(qry_struct)
        n_aligned = len(aligner._superimposer.reference_coords)
        return float(aligner.rms), n_aligned
    except Exception:
        return float("nan"), 0


def _tsv_safe_scalar(value) -> str:
    """Return a string-safe scalar for SSE datafile writeback.

    SSE datafiles are read and written as text. Assigning Python floats into
    existing string columns can fail on newer pandas versions, so RMSD values
    are normalized to strings before merge_columns(). Numeric-looking strings
    are still classified as continuous by the visualizer after reload.
    """
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    if isinstance(value, (float, np.floating)):
        if not np.isfinite(float(value)):
            return ""
        return f"{float(value):.6g}"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return str(value)


def write_rmsd_columns(entry: EntryContext, reference_id: str, reference_rank: int, methods: list[str]) -> None:
    df, types = read_datafile(entry.datafile_path)
    out = pd.DataFrame({COL_ID: df[COL_ID].astype(str)})
    log_df = read_rmsd_log(entry)
    if log_df.empty:
        return
    if "method" not in log_df.columns:
        log_df["method"] = "seq"
    log_df["method"] = log_df["method"].fillna("seq").replace("", "seq")
    log_df["rmsd"] = pd.to_numeric(log_df["rmsd"], errors="coerce")
    new_cols = []
    for method in methods:
        col = rmsd_column(reference_id, reference_rank, method)
        if col in df.columns:
            out[col] = df[col].astype(str).where(df[col].notna(), "")
        else:
            out[col] = ""
        m = (
            (log_df["reference_id"] == str(reference_id))
            & (log_df["reference_rank"] == str(reference_rank))
            & (log_df["method"] == method)
        )
        col_data = (
            log_df[m]
            .sort_values("timestamp")
            .drop_duplicates("query_id", keep="last")
            .set_index("query_id")["rmsd"]
        )
        mapped = out[COL_ID].map(col_data)
        good = mapped.notna()
        if good.any():
            out.loc[good, col] = mapped[good].map(_tsv_safe_scalar).astype(str)
        out[col] = out[col].map(_tsv_safe_scalar).astype(str)
        new_cols.append(col)
    merge_columns(
        entry.datafile_path,
        out[[COL_ID] + new_cols],
        {c: TYPE_LABEL for c in new_cols},
        manifest_path=entry.manifest_path if entry.manifest_path.exists() else None,
        provenance_source="sse",
        tool="sse_rmsd",
        notes=f"RMSD values against {reference_id} rank {reference_rank}.",
        force=True,
    )


def calculate_rmsds(entry: EntryContext, reference_id: str, reference_rank: int = 0,
                    query_ids: list[str] | None = None,
                    query_rank_map: dict[str, int] | None = None,
                    methods: list[str] | None = None,
                    force: bool = False) -> dict:
    methods = methods or ["seq"]
    query_rank_map = query_rank_map or {}
    ref_cif = cif_path_for(entry, reference_id, int(reference_rank))
    if ref_cif is None:
        raise FileNotFoundError(f"Rank {reference_rank} CIF not found for {reference_id}.")
    all_structs = list_apo_structures(entry)
    if query_ids is not None:
        wanted = set(str(x) for x in query_ids)
        queries = [s for s in all_structs if s["id"] in wanted and s["id"] != reference_id]
    else:
        queries = [s for s in all_structs if s["id"] != reference_id]
    results = []
    new_rows = []
    for s in queries:
        qry_id = s["id"]
        qry_rank = int(query_rank_map.get(qry_id, 0))
        qry_cif = cif_path_for(entry, qry_id, qry_rank)
        if qry_cif is None:
            continue
        for method in methods:
            was_cached = False
            hit = None if force else cached(entry, reference_id, qry_id, reference_rank, qry_rank, method)
            if hit is not None:
                rmsd, n_al = hit
                was_cached = True
            else:
                if method == "ce":
                    rmsd, n_al = compute_rmsd_ce(str(ref_cif), str(qry_cif))
                else:
                    rmsd, n_al = compute_rmsd_seq(str(ref_cif), str(qry_cif))
                new_rows.append({
                    "reference_id": str(reference_id),
                    "query_id": str(qry_id),
                    "reference_rank": str(reference_rank),
                    "query_rank": str(qry_rank),
                    "n_aligned_residues": str(n_al),
                    "rmsd": str(round(rmsd, 4)) if not np.isnan(rmsd) else "",
                    "method": method,
                    "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                })
            results.append({"query_id": qry_id, "query_rank": qry_rank, "method": method,
                            "n_aligned": n_al, "rmsd": rmsd, "cached": was_cached})
    if new_rows:
        append_rmsd_rows(entry, new_rows)
    if results:
        write_rmsd_columns(entry, reference_id, reference_rank, methods)
    key = f"ref={reference_id}|rank={reference_rank}|methods={','.join(methods)}|scope={'selected' if query_ids else 'all'}"
    upsert_job(entry.jobs_path, "rmsd", key, {
        "job_key": key,
        "status": "done",
        "reference_id": reference_id,
        "reference_rank": int(reference_rank),
        "methods": methods,
        "scope": "selected" if query_ids else "all",
        "n_results": len(results),
        "n_new": len(new_rows),
        "updated_utc": utc_now(),
    })
    return {"results": results, "n_new": len(new_rows), "n_cached": sum(1 for r in results if r["cached"]),
            "columns": [rmsd_column(reference_id, reference_rank, m) for m in methods]}
