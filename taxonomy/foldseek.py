"""Foldseek taxonomy strategy.

Resolves each row's NCBI taxId from what Foldseek/creation already carried
into the datafile (the `Databases` and `taxId` source labels, spec §5.2) --
no protein lookup needed, unlike the em strategy.

Three resolution paths, dispatched per row by `Databases` membership:

  A. taxId already in the file (afdb*, pdb100, BFVD, bfmd, cath50)
     -> used directly.
  B. gmgcl_id -> the GMGC unigene API returns a taxonomy list whose 'id' is
     an NCBI taxId. The row's id carries a Foldseek '_trun_<n>' suffix that
     must be stripped first to recover the GMGC unigene id.
  C. mgnify_esm30 -> MGnify protein (MGYP) accessions have no precomputed
     per-protein NCBI taxId; flagged via the NO_TAXONOMY sentinel.

NOTE on `Databases` as a tag_split column: SSE's own Foldseek reader dedups
multiple hit rows for the same target into one row, aggregating every
database the target appeared in in a comma-joined `Databases` label (spec
§5.2) -- so a single row can carry more than one of the tags above (e.g. a
target present in both `afdb50` and `gmgcl_id`). This module treats
`Databases` as a set of tags, not a scalar, and resolves in the priority
order A > B > C per row: an existing taxId wins if present, GMGC is tried
next, and a row is only flagged `no_taxonomy` if none of its tags offer any
other path. This priority is a judgment call, not something the original
one-row-per-hit script had to make -- confirm it matches your intent.
"""
import json
import re
import sys
import time

from .base import NO_TAXONOMY, http_request

NAME = "foldseek"

GMGC_BATCH_URL = "https://gmgc.embl.de/api/v1.0/unigenes/unigene"

# Foldseek 'Databases' tags that need the GMGC API path
GMGC_DBS = {"gmgcl_id"}
# tags with no taxonomy path at all
NO_TAX_DBS = {"mgnify_esm30"}

# ranks used only to pick the deepest GMGC taxonomy assignment
_GMGC_RANK_ORDER = ["superkingdom", "phylum", "class", "order",
                    "family", "genus", "species"]

# strips Foldseek's truncation marker, e.g. "..._trun_0" -> "..."
_TRUN_RE = re.compile(r"_trun_\d+$")


def detect(df, types) -> bool:
    """Positive signal: the datafile carries the Foldseek-only source labels
    laid down at creation (spec §5.2)."""
    return "Databases" in df.columns and "taxId" in df.columns


def _db_tags(value) -> set:
    return {t.strip() for t in str(value or "").split(",") if t.strip()}


def _strip_trun(target_id: str) -> str:
    return _TRUN_RE.sub("", target_id)


def _fetch_gmgc_taxids(unigene_ids, gmgc_batch=50, delay=0.34):
    """Resolve GMGC unigene IDs -> NCBI taxId via the GMGC batch API.

    Returns {unigene_id: taxid}, picking the deepest-rank assignment the
    API's taxonomy list offers for each unigene.
    """
    out = {}
    unigene_ids = sorted(unigene_ids)
    n_batches = (len(unigene_ids) + gmgc_batch - 1) // gmgc_batch
    for bi, i in enumerate(range(0, len(unigene_ids), gmgc_batch), start=1):
        sub = unigene_ids[i:i + gmgc_batch]
        print(f"  gmgc batch {bi}/{n_batches} ({len(sub)} unigenes)...",
              file=sys.stderr, flush=True)
        body = json.dumps({"names": sub}).encode("utf-8")
        text = http_request(GMGC_BATCH_URL, data=body,
                            headers={"Content-Type": "application/json"})
        if not text:
            continue
        try:
            data = json.loads(text)
        except json.JSONDecodeError:
            continue
        records = data.values() if isinstance(data, dict) else data
        for rec in records:
            if not isinstance(rec, dict):
                continue
            name = rec.get("name") or rec.get("query")
            tax = rec.get("taxonomy") or []
            if not name or not tax:
                continue
            best = None
            for t in tax:
                rank = (t.get("rank") or "").lower()
                tid = t.get("id")
                if not tid:
                    continue
                if rank in _GMGC_RANK_ORDER:
                    if best is None or _GMGC_RANK_ORDER.index(rank) >= _GMGC_RANK_ORDER.index(best[0]):
                        best = (rank, tid)
                elif best is None:
                    best = ("", tid)
            if best:
                out[name] = str(best[1])
        time.sleep(delay)
    return out


def resolve_taxids(ids, df, id_col, key_params, args, batch=100):
    """{id: taxid | NO_TAXONOMY} using the datafile's own Databases/taxId
    columns. `key_params`/`batch` are accepted for interface parity with the
    em strategy but unused here -- nothing in this path calls NCBI efetch.
    """
    ids = set(str(i) for i in ids)
    sub = df[df[id_col].astype(str).isin(ids)]

    out = {}
    gmgc_unigene_by_id = {}
    for _, row in sub.iterrows():
        rid = str(row[id_col])
        tags = _db_tags(row.get("Databases", ""))

        tid = str(row.get("taxId", "") or "").strip()
        if tid and tid not in ("0", "NA", "N/A"):
            out[rid] = tid
            continue
        if tags & GMGC_DBS:
            gmgc_unigene_by_id[rid] = _strip_trun(rid)
            continue
        if tags & NO_TAX_DBS:
            out[rid] = NO_TAXONOMY
        # else: no recognized tag and no taxId -> left unresolved

    if gmgc_unigene_by_id:
        gmgc_batch = getattr(args, "gmgc_batch", 50)
        gmgc_map = _fetch_gmgc_taxids(set(gmgc_unigene_by_id.values()), gmgc_batch=gmgc_batch)
        for rid, unigene in gmgc_unigene_by_id.items():
            if unigene in gmgc_map:
                out[rid] = gmgc_map[unigene]

    return out
