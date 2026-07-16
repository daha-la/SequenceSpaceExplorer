"""Structure input for the Foldseek-only embedders (ProstT5, SaProt).

Shared sub-pipeline: read Cα coordinates (tCa) from the original Foldseek JSON,
reconstruct the backbone (ca_reconstruct), convert to a 3Di sequence (mini3di),
and cache the result per Foldseek search so switching models or resuming does
not repeat it. ESM-C does not use any of this.

tCa is read from the JSON in initial_files/, not from the datafile — coordinates
would bloat the datafile, and creation deliberately leaves them out. The lookup
is keyed on (id, sequence) = (target, tSeq).
"""

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from ..common import abort, COL_ID, COL_SEQ
from .ca_reconstruct import reconstruct_backbone


# ---------------------------------------------------------------- tCa lookup

def _parse_tca(tca: str):
    if not tca:
        return None
    try:
        vals = np.fromstring(tca, sep=",", dtype=np.float64)
    except Exception:
        return None
    if vals.size == 0 or vals.size % 3 != 0 or not np.isfinite(vals).all():
        return None
    return vals.reshape(-1, 3)


def load_tca_lookup(json_path) -> dict:
    """{(target, tSeq): (N,3) ndarray} from the Foldseek JSON."""
    with open(json_path, encoding="utf-8") as fh:
        data = json.load(fh)
    lookup = {}
    for res in data[0].get("results", []):
        for hits in (res.get("alignments") or {}).values():
            for h in hits:
                key = (h.get("target", ""), h.get("tSeq", "") or "")
                if key in lookup:
                    continue
                coords = _parse_tca(h.get("tCa", ""))
                if coords is not None:
                    lookup[key] = coords
    return lookup


# ------------------------------------------------------- backbone -> 3Di

def _init_encoder():
    try:
        import mini3di
    except ImportError:
        abort("mini3di is required for structure embedders "
              "(pip install mini3di).")
    return mini3di.Encoder()


def _chain_from_backbone(rec):
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Atom import Atom
    chain = Chain("A")
    L = rec["CA"].shape[0]
    for i in range(L):
        residue = Residue((" ", i + 1, " "), "GLY", "")
        for name in ("N", "CA", "C", "CB"):
            coord = rec[name][i].astype(np.float64)
            residue.add(Atom(name, coord, 0.0, 1.0, " ", name, i + 1,
                             element=name[0]))
        chain.add(residue)
    return chain


def _coords_to_3di(encoder, ca):
    if ca.shape[0] < 3:
        return None
    try:
        rec = reconstruct_backbone(ca)
        states = encoder.encode_chain(_chain_from_backbone(rec))
        return encoder.build_sequence(states)
    except Exception:
        return None


# ----------------------------------------------------------------- 3Di cache

def _cache_path(entry_dir: Path, json_path) -> Path:
    cache_dir = Path(entry_dir) / "embeddings" / "3di_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    stem = Path(json_path).stem
    return cache_dir / f"{stem}_3di.tsv"


def _load_cache(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except Exception:
        print(f"  warning: unreadable 3Di cache {path}; recomputing.")
        return {}
    return {(r["id"], r["Sequence"]): r["di_seq"]
            for _, r in df.iterrows() if r.get("di_seq")}


def _save_cache(path: Path, cache: dict):
    rows = [{"id": tid, "Sequence": seq, "di_seq": di}
            for (tid, seq), di in cache.items()]
    pd.DataFrame(rows, columns=["id", "Sequence", "di_seq"]).to_csv(
        path, sep="\t", index=False)


# ----------------------------------------------------------- public entry point

def compute_3di(df: pd.DataFrame, ctx, max_residues: int):
    """Return ({id: di_seq} for rows that succeeded, skip_report).

    skip_report buckets: no_tca, too_long, threedi_fail. Rows are keyed on
    (id, Sequence) against the tCa lookup and the cache.
    """
    if ctx.foldseek_json is None or not Path(ctx.foldseek_json).exists():
        abort("structure embedder needs the source Foldseek JSON; none found. "
              "Pass --foldseek-json or keep it in initial_files/.")

    print(f"  indexing tCa from {ctx.foldseek_json} ...")
    tca = load_tca_lookup(ctx.foldseek_json)
    cache_path = _cache_path(ctx.entry_dir, ctx.foldseek_json)
    cache = _load_cache(cache_path)
    n_cache_start = len(cache)

    encoder = None
    di_by_id, skip = {}, {"no_tca": [], "too_long": [], "threedi_fail": []}

    for _, r in tqdm(df.iterrows(), total=len(df), desc="3Di conversion"):
        rid, seq = r[COL_ID], r[COL_SEQ]
        key = (rid, seq)
        di = cache.get(key)
        if di is None:
            coords = tca.get(key)
            if coords is None:
                skip["no_tca"].append(rid)
                continue
            if coords.shape[0] > max_residues:
                skip["too_long"].append(rid)
                continue
            if encoder is None:
                encoder = _init_encoder()
            di = _coords_to_3di(encoder, coords)
            if di is None:
                skip["threedi_fail"].append(rid)
                continue
            cache[key] = di
        di_by_id[rid] = di

    if len(cache) > n_cache_start:
        _save_cache(cache_path, cache)
    return di_by_id, skip
