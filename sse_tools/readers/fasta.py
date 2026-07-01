"""FASTA reader (.fasta/.faa/.fa/.txt)."""

from collections import Counter

import pandas as pd

from .base import ReaderResult
from ..common import abort, COL_ID, COL_SEQ

# UniProt/UniRef-style db tags recognised in headers (db|ACCESSION|NAME).
# Extend if your FASTA sources use other tags.
KNOWN_DB_TAGS = {"sp", "tr", "up", "ur", "ur100", "ur90", "ur50"}


def _parse_fasta(path):
    """Minimal FASTA parser. Returns list of (header_without_'>', sequence)."""
    records, header, chunks = [], None, []
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            elif line.strip():
                chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def read_fasta(path, args) -> ReaderResult:
    """FASTA.

    ID resolution per spec §5.3: first whitespace token of the header; if it
    splits on '|' into a known db tag + accession, take the accession, else the
    first pipe-field. Collisions are resolved in a second pass by appending the
    next unused pipe-field; unresolved collisions abort (caught downstream by the
    uniqueness check). The post-whitespace remainder of the header is kept as a
    `Description` label. `--query` matches the FULL header.
    """
    records = _parse_fasta(path)
    if not records:
        abort("--source fasta: no records found (no sequence information).")

    # Pass 1: candidate ID + fallback (next unused pipe-field).
    cand_ids, fallbacks, headers, descs, seqs = [], [], [], [], []
    for header, seq in records:
        first_ws = header.split(None, 1)
        token = first_ws[0] if first_ws else header
        desc = first_ws[1] if len(first_ws) > 1 else ""
        fields = token.split("|")
        if fields[0] in KNOWN_DB_TAGS and len(fields) >= 2:
            cand = fields[1]
            fallback = fields[2] if len(fields) >= 3 and fields[2] else None
        else:
            cand = fields[0]
            fallback = fields[1] if len(fields) >= 2 and fields[1] else None
        cand_ids.append("_".join(cand.split()))  # whitespace guard
        fallbacks.append("_".join(fallback.split()) if fallback else None)
        headers.append(header)
        descs.append(desc)
        seqs.append(seq)

    # Pass 2: append fallback for colliding candidates.
    counts = Counter(cand_ids)
    final_ids = []
    for cand, fb in zip(cand_ids, fallbacks):
        final_ids.append(f"{cand}_{fb}" if counts[cand] > 1 and fb else cand)

    df = pd.DataFrame({
        COL_ID: final_ids,
        COL_SEQ: seqs,
        "Description": descs,
        "_full_header": headers,
    })
    return ReaderResult(table=df, id_col=COL_ID, seq_col=COL_SEQ, source="fasta",
                        match_col="_full_header", match_label="full header",
                        auto_query=None, notes={"records": len(records)})
