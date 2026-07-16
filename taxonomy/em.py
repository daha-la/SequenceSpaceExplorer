"""Generic / EnzymeMiner taxonomy strategy.

Resolves each row's NCBI taxId via efetch protein on the datafile's id
column. This is the universal fallback: it needs nothing beyond id + network
access, so it works for any datafile whose id column holds NCBI protein
accessions -- not only EnzymeMiner entries.
"""
import sys
import time
import xml.etree.ElementTree as ET

from .base import eutils

NAME = "em"


def detect(df, types) -> bool:
    """EM/generic has no positive signal of its own -- it is only chosen
    when no other strategy's detect() fires (see taxonomy/__init__.py's
    detect_strategy). This always returns False so a strategy with actual
    evidence (e.g. foldseek's Databases/taxId columns) is preferred in auto
    mode; --strategy em still works regardless of what this returns.
    """
    return False


def resolve_taxids(ids, df, id_col, key_params, args, batch=100):
    """efetch protein (GBSeq XML) -> {id: taxid}.

    Ids are matched both with and without a version suffix, since some
    sources carry versioned accessions (P12345.2) and some don't.
    """
    ids = sorted(set(str(i) for i in ids))
    delay = 0.12 if key_params.get("api_key") else 0.4  # respect NCBI rate limits
    n_batches = (len(ids) + batch - 1) // batch
    by_acc = {}
    for bi, i in enumerate(range(0, len(ids), batch), start=1):
        sub = ids[i:i + batch]
        print(f"  em batch {bi}/{n_batches} ({len(sub)} ids)...",
              file=sys.stderr, flush=True)
        params = dict(key_params, db="protein", id=",".join(sub),
                      rettype="gb", retmode="xml")
        xml = eutils("efetch", params)
        if xml:
            try:
                root = ET.fromstring(xml)
            except ET.ParseError:
                root = None
            if root is not None:
                for seq in root.iter("GBSeq"):
                    primary = seq.findtext("GBSeq_primary-accession")
                    accver = seq.findtext("GBSeq_accession-version")
                    acc = accver or primary
                    taxid = None
                    for qual in seq.iter("GBQualifier"):
                        if qual.findtext("GBQualifier_name") == "db_xref":
                            v = qual.findtext("GBQualifier_value") or ""
                            if v.startswith("taxon:"):
                                taxid = v.split(":", 1)[1]
                    if acc and taxid:
                        by_acc[acc] = taxid
                        by_acc[acc.split(".")[0]] = taxid  # also index without version
        time.sleep(delay)

    out = {}
    for i in ids:
        if i in by_acc:
            out[i] = by_acc[i]
        elif i.split(".")[0] in by_acc:
            out[i] = by_acc[i.split(".")[0]]
    return out
