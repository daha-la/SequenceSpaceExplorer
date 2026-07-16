"""Shared NCBI taxonomy lookup primitives.

Strategy-agnostic: the E-utilities HTTP primitive, lineage expansion, rank
aliases, and the merged column set. Each strategy module (em.py, foldseek.py)
only supplies how to get from a datafile row to an NCBI taxId; everything
past "I have a taxId" is shared here. Mirrors the readers/ split (creation,
spec §5.0) and the embedders//reducers/ split (coordinates, spec §11A): one
shared foundation, one file per strategy, an explicit registry in __init__.py.
"""
import http.client
import socket
import sys
import time
import xml.etree.ElementTree as ET
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# NCBI is renaming the top rank "superkingdom" -> "domain"; accept either
# and normalize both onto the "superkingdom" output column.
RANKS = ["superkingdom", "phylum", "class", "order",
         "family", "genus", "species"]
RANK_ALIASES = {r: r for r in RANKS}
RANK_ALIASES["domain"] = "superkingdom"

# label columns merged into the datafile by every strategy
TAX_COLS = ["tax_status", "taxid", "tax_organism"] + RANKS

# sentinel returned by a strategy's resolve_taxids() for a row that has no
# taxonomy path at all (e.g. an MGnify hit under the foldseek strategy).
# Distinct from "absent from the returned dict", which means "unresolved,
# try again later" rather than "known to have nothing."
NO_TAXONOMY = "__no_taxonomy__"


def http_request(url, data=None, headers=None, retries=5):
    """POST (if data) or GET a URL, return decoded text. Retries w/ backoff
    on HTTP, URL, connection-reset, and incomplete-read errors.
    """
    for attempt in range(retries):
        try:
            req = Request(url, data=data, headers=headers or {})
            with urlopen(req, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except (HTTPError, URLError, http.client.IncompleteRead,
                http.client.HTTPException, socket.error) as e:
            wait = 2 ** attempt
            sys.stderr.write(f"  request error ({e}); retry in {wait}s\n")
            time.sleep(wait)
    sys.stderr.write(f"  giving up after {retries} attempts: {url}\n")
    return None


def eutils(endpoint, params):
    """POST to an NCBI E-utilities endpoint (large id lists go in the body)."""
    return http_request(f"{EUTILS}/{endpoint}.fcgi",
                         data=urlencode(params).encode("utf-8"))


def fetch_lineages(taxids, key_params, batch=100, delay=0.34):
    """efetch taxonomy -> {taxid: {tax_organism, rank: name, ...}}.

    Shared by every strategy: once a strategy has produced a taxId, expanding
    it to a full lineage no longer depends on where the taxId came from.
    """
    out = {}
    taxids = sorted(set(taxids))
    n_batches = (len(taxids) + batch - 1) // batch
    for bi, i in enumerate(range(0, len(taxids), batch), start=1):
        sub = taxids[i:i + batch]
        print(f"  lineage batch {bi}/{n_batches} ({len(sub)} taxids)...",
              file=sys.stderr, flush=True)
        params = dict(key_params, db="taxonomy", id=",".join(sub), retmode="xml")
        xml = eutils("efetch", params)
        if not xml:
            continue
        try:
            root = ET.fromstring(xml)
        except ET.ParseError:
            continue
        for tax in root.findall("Taxon"):  # top-level records only
            tid = tax.findtext("TaxId")
            sciname = tax.findtext("ScientificName")
            if not tid or not sciname:
                continue
            rec = {"tax_organism": sciname}
            own = RANK_ALIASES.get(tax.findtext("Rank"))
            if own:
                rec[own] = sciname
            lineage_ex = tax.find("LineageEx")
            if lineage_ex is not None:
                for t in lineage_ex.findall("Taxon"):
                    r = RANK_ALIASES.get(t.findtext("Rank"))
                    n = t.findtext("ScientificName")
                    if r and n:
                        rec[r] = n
            out[tid] = rec
        time.sleep(delay)
    return out
