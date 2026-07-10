#!/usr/bin/env python3
"""Merge NCBI taxonomy columns into an SSE entry's datafile.

Resolves each row's NCBI taxId using a strategy chosen by --strategy (or
auto-detected from the datafile's columns), then expands every resolved
taxId into a full lineage (superkingdom -> species) via NCBI Entrez efetch,
and merges the result into the entry's .sse.tsv as label columns through
common.merge_columns (the only write path onto a datafile).

Strategies (sse_tools/taxonomy/):
  em        efetch protein on the id column. Universal; works for any
            datafile whose id holds an NCBI protein accession.
  foldseek  uses the Databases/taxId source labels already carried into the
            datafile by the Foldseek reader (spec §5.2) -- existing taxId,
            GMGC unigene lookup, or a no-taxonomy flag for MGnify hits.
  auto (default)  detect_strategy() picks foldseek if the datafile carries
            Databases+taxId columns, else falls back to em.

Batches are cached to a per-entry tmp folder as they complete, so an
interrupted run can be resumed by re-running the same command. Once a run
finishes, the tmp cache is deleted -- but a plain rerun still resumes
correctly, because already-resolved rows are re-seeded from the datafile's
own taxonomy columns (if present) rather than refetched. --force wipes the
cache, ignores whatever's already on the datafile, refetches every row from
scratch, and overwrites existing taxonomy columns. --retry-failed additionally
re-attempts rows previously marked taxid_unresolved/lineage_unresolved
(transient failures) without touching rows that already succeeded or that are
genuinely no_taxonomy.

Usage:
    python fetch_taxonomy.py ENTRY --email you@dtu.dk
    python fetch_taxonomy.py ENTRY --email you@dtu.dk --strategy foldseek
    python fetch_taxonomy.py ENTRY --email you@dtu.dk --retry-failed
"""
import argparse
import csv
import os
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sse_tools import common
from sse_tools import taxonomy
from sse_tools.taxonomy.base import NO_TAXONOMY, RANKS, TAX_COLS, fetch_lineages

# tax_status values that represent a definitive, non-retryable outcome.
# Anything else (empty, taxid_unresolved, lineage_unresolved) is retryable.
_TERMINAL_STATUSES = {"ok", "no_taxonomy"}


def load_cache(cache_path, out_cols):
    if not cache_path.exists():
        return []
    with open(cache_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    for r in rows:
        for c in out_cols:
            r.setdefault(c, "")
    return rows


def save_cache(cache_path, rows, out_cols):
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _merge_and_cleanup(datafile_path, manifest_path, tmp_dir, df, id_col,
                       out_cols, row_by_id, strategy_name, args, merge_force):
    new_df = pd.DataFrame([{c: row_by_id[i][c] for c in out_cols} for i in df[id_col]])
    new_types = {id_col: common.TYPE_ID, **{c: common.TYPE_LABEL for c in TAX_COLS}}

    common.merge_columns(
        datafile_path, new_df, new_types,
        manifest_path=manifest_path if manifest_path.exists() else None,
        provenance_source="sse", tool="fetch_taxonomy", version="2.1",
        params=f"strategy={strategy_name},batch={args.batch}",
        id_col=id_col, force=merge_force,
    )

    shutil.rmtree(tmp_dir, ignore_errors=True)

    counts = {}
    for r in row_by_id.values():
        counts[r["tax_status"]] = counts.get(r["tax_status"], 0) + 1
    print(f"done; merged taxonomy columns into {datafile_path}")
    for st in sorted(counts):
        print(f"  {st or '(blank)'}: {counts[st]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("entry", help="entry stem (in entries/) or path to a .sse.tsv")
    ap.add_argument("--email", required=True, help="NCBI requires an email")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"),
                    help="optional NCBI API key (raises rate limit to 10/s)")
    ap.add_argument("--strategy", choices=["auto"] + list(taxonomy.REGISTRY),
                    default="auto",
                    help="taxId-resolution strategy, or auto-detect from the datafile")
    ap.add_argument("--batch", type=int, default=100,
                    help="ids / taxIds per NCBI request batch")
    ap.add_argument("--gmgc-batch", type=int, default=50,
                    help="unigene ids per GMGC API request (foldseek strategy)")
    ap.add_argument("--force", action="store_true",
                    help="wipe any cached progress, ignore the datafile's existing "
                         "taxonomy columns, and refetch every row from scratch")
    ap.add_argument("--retry-failed", action="store_true",
                    help="also re-attempt rows previously marked taxid_unresolved "
                         "or lineage_unresolved, without refetching rows that already "
                         "succeeded (ok) or are genuinely no_taxonomy")
    ap.add_argument("--entries-dir", default=None)
    args = ap.parse_args()

    datafile_path = common.resolve_entry_path(args.entry, args.entries_dir)
    stem = datafile_path.stem.removesuffix(".sse")
    manifest_path = datafile_path.parent / "logs" / f"{stem}.sse.manifest.json"
    tmp_dir = datafile_path.parent / "tmp_taxonomy"
    cache_path = tmp_dir / "taxonomy_cache.tsv"

    if args.force and tmp_dir.exists():
        shutil.rmtree(tmp_dir)

    df, types = common.read_datafile(datafile_path)
    id_col = common.id_column(types)
    out_cols = [id_col] + TAX_COLS
    existing_tax_cols = [c for c in TAX_COLS if c in df.columns]

    strategy_name = args.strategy
    if strategy_name == "auto":
        strategy_name = taxonomy.detect_strategy(df, types)
    strategy = taxonomy.REGISTRY[strategy_name]
    print(f"taxonomy strategy: {strategy_name}"
          f"{' (auto-detected)' if args.strategy == 'auto' else ''}")

    # Seed already-resolved rows from the datafile's own taxonomy columns, if
    # any -- this is what makes a plain rerun after a *successful* prior run
    # resume instead of refetching everything (the tmp cache is deleted on
    # success, so it alone can't carry that state across runs). Skipped
    # entirely under --force: every row starts from scratch.
    datafile_seed = {}
    if existing_tax_cols and not args.force:
        for _, row in df[[id_col] + existing_tax_cols].iterrows():
            rid = str(row[id_col])
            rec = {c: "" for c in out_cols}
            rec[id_col] = rid
            for c in existing_tax_cols:
                v = row[c]
                rec[c] = "" if pd.isna(v) else str(v)
            datafile_seed[rid] = rec

    rows = load_cache(cache_path, out_cols)
    cached_ids = {r[id_col] for r in rows}
    for i in df[id_col].astype(str):
        if i in cached_ids:
            continue
        rows.append(datafile_seed.get(i) or {**{c: "" for c in out_cols}, id_col: i})
    row_by_id = {r[id_col]: r for r in rows}

    key_params = {"email": args.email, "tool": "sse_taxonomy"}
    if args.api_key:
        key_params["api_key"] = args.api_key

    def _needs_resolution(rid):
        status = row_by_id[rid]["tax_status"]
        if not status:
            return True
        if args.retry_failed and status not in _TERMINAL_STATUSES:
            return True
        return False

    todo_ids = [i for i in df[id_col].astype(str) if _needs_resolution(i)]
    print(f"{len(df)} rows total, {len(df) - len(todo_ids)} already resolved, "
          f"{len(todo_ids)} to resolve")

    # Reconstructing over already-existing taxonomy columns is always safe
    # here: every previously-resolved row was seeded above and carried
    # through untouched, so the merge writes back the same values plus
    # whatever's newly resolved -- never a destructive overwrite. --force
    # forces it too, for the case where existing_tax_cols is empty because
    # this is a genuinely first run.
    merge_force = args.force or bool(existing_tax_cols)

    if not todo_ids:
        _merge_and_cleanup(datafile_path, manifest_path, tmp_dir, df, id_col,
                            out_cols, row_by_id, strategy_name, args, merge_force)
        return

    # Rows already carrying a taxid from a prior lineage_unresolved attempt
    # only need lineage retried, not a redundant taxId re-resolution.
    todo_lineage_only = [i for i in todo_ids
                         if row_by_id[i]["tax_status"] == "lineage_unresolved"
                         and row_by_id[i].get("taxid")]
    todo_resolve = [i for i in todo_ids if i not in set(todo_lineage_only)]

    # Clear stale status before re-resolving this pass -- otherwise a retried
    # row that succeeds still carries its old taxid_unresolved/
    # lineage_unresolved status and gets silently excluded from `pending`
    # below (which requires a blank tax_status to mean "not yet handled
    # this pass").
    for i in todo_ids:
        row_by_id[i]["tax_status"] = ""

    # --- step 1: resolve a taxId (or NO_TAXONOMY) for every pending row ---
    if todo_resolve:
        print(f"resolving taxIds for {len(todo_resolve)} row(s) via '{strategy_name}'...")
        resolved = strategy.resolve_taxids(todo_resolve, df, id_col, key_params, args,
                                           batch=args.batch)
        for i in todo_resolve:
            val = resolved.get(i)
            if val is None:
                row_by_id[i]["tax_status"] = "taxid_unresolved"
            elif val == NO_TAXONOMY:
                row_by_id[i]["tax_status"] = "no_taxonomy"
            else:
                row_by_id[i]["taxid"] = val
    save_cache(cache_path, rows, out_cols)  # persist taxId resolution before lineage step

    # --- step 2: expand all resolved taxIds into full lineages ------------
    pending = [i for i in todo_ids if row_by_id[i].get("taxid")
               and not row_by_id[i].get("tax_status")]
    all_taxids = {row_by_id[i]["taxid"] for i in pending}
    print(f"expanding {len(all_taxids)} unique taxId(s) into lineages...")
    lineages = fetch_lineages(all_taxids, key_params, batch=args.batch)

    # --- step 3: write lineage columns, saving after every batch ----------
    for start in range(0, len(pending), args.batch):
        for i in pending[start:start + args.batch]:
            tid = row_by_id[i]["taxid"]
            lin = lineages.get(tid)
            if not lin:
                row_by_id[i]["tax_status"] = "lineage_unresolved"
                continue
            row_by_id[i]["tax_status"] = "ok"
            row_by_id[i]["tax_organism"] = lin.get("tax_organism", "")
            for rk in RANKS:
                row_by_id[i][rk] = lin.get(rk, "")
        save_cache(cache_path, rows, out_cols)

    _merge_and_cleanup(datafile_path, manifest_path, tmp_dir, df, id_col,
                        out_cols, row_by_id, strategy_name, args, merge_force)


if __name__ == "__main__":
    main()
