"""Embedder contract + shared embedding machinery.

An embedder turns an entry's rows into an embedding matrix (one row per
sequence, columns = embedding dimensions). Embedders differ in their model and
their *input requirement* (ESM-C needs the AA sequence; structure embedders need
a 3Di sequence derived from Foldseek coordinates). The streaming-to-disk loop
with ID-based resume, pooling, and device selection are shared here so the
per-model code is just the forward pass.

Heavy imports (torch, model libraries) are lazy — inside the methods that need
them — so the CLI, the reducers, and the merge path all work without a GPU stack
installed; only actually running an embedder pulls torch in.

CONTRACT
--------
Subclass Embedder and set `name` (the --embedder key) and `requires_structure`.
Implement:
  tag(args)             -> str : coordinate-system tag encoding everything that
                                 changes the embedding identity (model variant,
                                 pooling). Becomes the column-name prefix.
  prepare(df, ctx, args)-> (entries, skip_report)
                                 entries: list of dicts, each at least {"id": ...}
                                 plus whatever encode_batch needs.
                                 skip_report: {reason: [ids]} for rows that cannot
                                 be embedded (empty unless the source is partial).
  load_model(args)      -> model_ctx : load model/tokenizer once.
  encode_batch(model_ctx, batch, args) -> list[np.ndarray] : pooled vector per
                                 entry in `batch`, aligned to it.
The default embed() runs the shared streaming loop over these.
"""

import math
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..common import abort


@dataclass
class EmbedContext:
    """What an embedder may need beyond the table itself."""
    entry_dir: Path
    stem: str
    source_type: str
    source_file: str
    initial_files_dir: Path
    foldseek_json: Optional[Path] = None   # resolved for structure embedders


# --------------------------------------------------------------- torch helpers

# --------------------------------------------------------------- torch helpers

def resolve_device(requested: str) -> str:
    """Resolve --device to a concrete backend. 'auto' picks cuda > mps > cpu;
    an explicit choice that is unavailable aborts (rather than silently
    downgrading). 'mps' is Apple Silicon's Metal backend."""
    import torch
    req = (requested or "auto").lower()
    has_cuda = torch.cuda.is_available()
    has_mps = (getattr(torch.backends, "mps", None) is not None
               and torch.backends.mps.is_available())
    if req == "auto":
        return "cuda" if has_cuda else ("mps" if has_mps else "cpu")
    if req == "cuda" and not has_cuda:
        abort("--device cuda requested but no CUDA device is available.")
    if req == "mps" and not has_mps:
        abort("--device mps requested but no MPS (Apple Silicon) GPU is available.")
    if req not in ("cuda", "mps", "cpu"):
        abort(f"unknown --device {requested!r} (choose auto/cuda/mps/cpu).")
    return req


def use_fp16(device: str) -> bool:
    """Half precision only on CUDA. fp16 matmuls are unsupported or pathologically
    slow on CPU, and flaky on MPS, so those stay fp32."""
    return device == "cuda"


def warn_if_slow(device: str, n: int):
    if device == "cpu" and n > 1000:
        print(f"  WARNING: embedding {n} sequences on CPU can take hours. "
              f"A CUDA machine is far faster for datasets this size.", flush=True)


def to_numpy_fp32(t):
    import torch
    return t.detach().to(device="cpu", dtype=torch.float32).numpy()


def pool(emb, how: str, mask=None):
    """Pool a (L, D) tensor to (D,). If `mask` is given, pool over True positions
    (falling back to all positions if the mask is empty)."""
    valid = emb
    if mask is not None:
        sel = emb[mask]
        valid = sel if sel.numel() > 0 else emb
    if how == "mean":
        return valid.mean(dim=0)
    if how == "max":
        return valid.max(dim=0).values
    if how == "min":
        return valid.min(dim=0).values
    abort(f"unsupported pooling {how!r}")


# --------------------------------------------------------------- the contract

class Embedder:
    name: str = ""
    requires_structure: bool = False

    def tag(self, args) -> str:
        raise NotImplementedError

    def prepare(self, df: pd.DataFrame, ctx: EmbedContext, args):
        raise NotImplementedError

    def load_model(self, args):
        raise NotImplementedError

    def encode_batch(self, model_ctx, batch, args):
        raise NotImplementedError

    def embed(self, entries, out_path, args) -> pd.DataFrame:
        """Shared streaming loop with ID-based resume. Returns the embedding
        matrix (column 'ID' + integer-named dimension columns)."""
        return _run_streaming(self, entries, out_path, args)


# --------------------------------------------------------- streaming + resume

def _read_done_ids(path):
    existing = pd.read_csv(path, sep="\t", usecols=["ID"])
    return set(existing["ID"].astype(str))


def _resume(tmp_path, out_path, entries, *, force_embedding=False):
    """Seed resume state from a completed matrix first, then a .part file.

    A completed <tag>.emb.tsv is a reusable cache across reducer changes. A
    .part file is only an interrupted-run resume source. If new rows have been
    added since the completed matrix was written, copy the completed matrix to
    .part and append only the missing IDs.
    """
    if force_embedding:
        for p in (tmp_path, out_path):
            if os.path.exists(p):
                os.remove(p)
        return set(), False

    wanted = {str(e["id"]) for e in entries}

    if os.path.exists(out_path):
        done = _read_done_ids(out_path)
        covered = wanted & done
        missing = wanted - done
        if not missing:
            print(f"  using cached embeddings: {len(covered)} already embedded.")
            return done, True
        shutil.copyfile(out_path, tmp_path)
        print(f"  extending cached embeddings: {len(covered)} done, {len(missing)} missing.")
        return done, True

    if os.path.exists(tmp_path):
        done = _read_done_ids(tmp_path)
        print(f"  resuming interrupted run: {len(wanted & done)} already embedded.")
        return done, True

    return set(), False

def _run_streaming(embedder: Embedder, entries, out_path, args) -> pd.DataFrame:
    from tqdm import tqdm

    tmp = str(out_path) + ".part"
    force_embedding = bool(getattr(args, "force_embedding", False))
    done, wrote_header = _resume(tmp, out_path, entries,
                                  force_embedding=force_embedding)
    remaining = [e for e in entries if str(e["id"]) not in done]
    if not remaining:
        # Complete-cache path: _resume may have reused out_path directly, in
        # which case no .part file exists. Interrupted-run path: the .part file
        # may itself already cover every requested ID, so promote it to the
        # completed cache.
        if os.path.exists(tmp):
            os.replace(tmp, out_path)
        elif not os.path.exists(out_path):
            abort("embedding cache state is inconsistent: all IDs are marked "
                  "done but neither the completed embedding file nor the .part "
                  "file exists.")
        return pd.read_csv(out_path, sep="\t")

    import torch
    model_ctx = embedder.load_model(args)
    warn_if_slow(model_ctx.get("device", "cpu"), len(remaining))
    buffer, processed = [], 0
    state = {"wrote_header": wrote_header}

    def flush():
        if not buffer:
            return
        pd.DataFrame(buffer).to_csv(tmp, sep="\t", index=False, mode="a",
                                    header=not state["wrote_header"])
        state["wrote_header"] = True
        buffer.clear()

    bs = args.batch_size
    with torch.inference_mode():
        for start in tqdm(range(0, len(remaining), bs),
                          total=math.ceil(len(remaining) / bs),
                          desc=f"Embedding ({embedder.name})"):
            batch = remaining[start:start + bs]
            vecs = embedder.encode_batch(model_ctx, batch, args)
            for e, vec in zip(batch, vecs):
                row = {"ID": str(e["id"])}
                row.update({str(i): float(v) for i, v in enumerate(vec)})
                buffer.append(row)
                processed += 1
                if processed % args.write_every == 0:
                    flush()
    flush()
    os.replace(tmp, out_path)
    print(f"  embeddings written to {out_path}")
    return pd.read_csv(out_path, sep="\t")
