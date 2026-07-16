#!/usr/bin/env python3
"""sse_boltz.py — Boltz-2 structure prediction + RMSD for an entry.

Imports a selection cache exported from the visualizer (a timestamped JSON file
under entries/<entry>/selections/), predicts a structure for each selected
sequence with the NVIDIA-hosted Boltz-2 model, optionally with one or more
substrate SMILES (holo prediction), and optionally computes pairwise RMSDs
between the predicted apo structures.

This is the pipeline-side driver for the workflow that used to live in the Dash
visualizer. The heavy lifting is unchanged and reused from sse_tools.boltz and
sse_tools.rmsd: predictions write pTM/pLDDT columns into the entry .sse.tsv, save
.cif structures under entries/<entry>/structures/, and log to
logs/boltz_log.csv; RMSDs write RMSD_vs_* columns and logs/rmsd_log.csv.

The NVIDIA API key is read from the BOLTZ_API_KEY environment variable (the
pipeline UI injects it as a per-run secret).

examples:
  sse_boltz.py akr --list-selections
  sse_boltz.py akr                                   # newest selection, apo
  sse_boltz.py akr --selection selection_20260715_101500.json
  sse_boltz.py akr --smiles "OC[C@H]1O..." --smiles-label UDP-Glc
  sse_boltz.py akr --rmsd --rmsd-reference OleD_S1 --rmsd-method both
"""

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sse_tools.common import SSEError, abort  # noqa: E402
from sse_tools import selections as selection_cache  # noqa: E402

TOOL_NAME = "sse_boltz"
TOOL_VERSION = "1.0.0"


def _check_dependencies() -> None:
    """Fail early with an install hint if a runtime dependency is missing."""
    missing = []
    try:
        import requests  # noqa: F401
    except ImportError:
        missing.append("requests")
    try:
        import Bio  # noqa: F401
    except ImportError:
        missing.append("biopython")
    if missing:
        abort(
            "Missing required package(s): "
            + ", ".join(missing)
            + f".\n  Install them into this environment with: pip install {' '.join(missing)}"
        )


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Boltz-2 structure prediction + RMSD for a selection of sequences.")
    ap.add_argument("entry", help="Entry stem, entry directory, or .sse.tsv datafile path")
    ap.add_argument("--selection", default=None,
                    help="Selection cache filename (in selections/) or path. "
                         "Default: the most recent selection cache.")
    ap.add_argument("--list-selections", action="store_true",
                    help="List available selection caches for the entry and exit.")

    # Boltz-2 parameters (mirror sse_tools.boltz.BoltzParams).
    ap.add_argument("--smiles", default="",
                    help="Substrate SMILES for holo prediction; one per line "
                         "(use \\n) or a single string. Empty = apo prediction.")
    ap.add_argument("--smiles-label", default="",
                    help="Human-readable ligand label used in holo output columns/folders.")
    ap.add_argument("--no-msa", action="store_true",
                    help="Skip ColabFold MSA generation (faster, lower accuracy).")
    ap.add_argument("--recycling-steps", type=int, default=3)
    ap.add_argument("--sampling-steps", type=int, default=200)
    ap.add_argument("--diffusion-samples", type=int, default=5)
    ap.add_argument("--step-scale", type=float, default=1.638)
    ap.add_argument("--force", action="store_true",
                    help="Re-run predictions even if a cached result exists.")

    # RMSD options.
    ap.add_argument("--rmsd", action="store_true",
                    help="After prediction, compute RMSDs between predicted apo structures.")
    ap.add_argument("--rmsd-reference", default=None,
                    help="Reference sequence id for RMSD (required with --rmsd).")
    ap.add_argument("--rmsd-ref-rank", type=int, default=0,
                    help="Reference structure rank (default 0).")
    ap.add_argument("--rmsd-method", choices=["seq", "ce", "both"], default="seq",
                    help="Alignment method: sequence-guided, structure-based (CE), or both.")
    ap.add_argument("--rmsd-scope", choices=["all", "selected"], default="all",
                    help="Compare against all apo structures or only the selected ids.")

    ap.add_argument("--entries-dir", default=None,
                    help="Override the entries directory (defaults to <repo>/entries).")
    return ap


def _print_selection_list(entry) -> None:
    picks = selection_cache.list_selections(entry.selections_dir)
    if not picks:
        print(f"No selection caches found in {entry.selections_dir}")
        print("Export a selection from the visualizer first "
              "(select points, then 'Export selection for Boltz').")
        return
    print(f"{len(picks)} selection cache(s) in {entry.selections_dir} (newest first):")
    for p in picks:
        try:
            payload = selection_cache.read_selection(p)
            n = payload.get("count", len(payload.get("sequences", [])))
            created = payload.get("created_utc", "")
            print(f"  {p.name}  ({n} sequence(s), created {created})")
        except Exception:
            print(f"  {p.name}  (unreadable)")


def main(argv=None):
    args = build_parser().parse_args(argv)

    # resolve_entry is safe to import unconditionally (only pandas/numpy). The
    # boltz backend imports `requests` at module scope, so it is imported only
    # after the dependency check, and --list-selections skips it entirely.
    from sse_tools.visualizer_state import resolve_entry

    entries_dir = Path(args.entries_dir) if args.entries_dir else None

    try:
        entry = resolve_entry(args.entry, entries_dir=entries_dir)

        if args.list_selections:
            _print_selection_list(entry)
            return 0

        _check_dependencies()
        from sse_tools import boltz as boltz_backend
        from sse_tools import rmsd as rmsd_backend

        # Resolve which selection cache to run.
        path = selection_cache.resolve_selection(entry.selections_dir, args.selection)
        if path is None:
            if args.selection:
                abort(f"Selection cache not found: {args.selection} "
                      f"(looked in {entry.selections_dir}).")
            abort(f"No selection caches in {entry.selections_dir}. Export one from "
                  "the visualizer first ('Export selection for Boltz').")
        payload = selection_cache.read_selection(path)
        sequences = payload.get("sequences", [])
        if not sequences:
            abort(f"Selection cache {path.name} contains no sequences.")

        # API key.
        api_key = os.environ.get("BOLTZ_API_KEY", "").strip()
        if not api_key:
            abort("BOLTZ_API_KEY is not set. Enter the NVIDIA API key in the "
                  "pipeline module, or export BOLTZ_API_KEY in your shell.")
        ok, msg = boltz_backend.validate_api_key(api_key)
        if not ok:
            abort(f"Boltz-2 API key check failed: {msg}")
        print("Boltz-2 API key: valid")

        params = boltz_backend.BoltzParams(
            recycling_steps=args.recycling_steps,
            sampling_steps=args.sampling_steps,
            diffusion_samples=args.diffusion_samples,
            step_scale=args.step_scale,
        )
        smiles = (args.smiles or "").replace("\\n", "\n")
        kind = "holo" if smiles.strip() else "apo"
        use_msa = not args.no_msa

        print(f"selection: {path.name} — {len(sequences)} sequence(s), kind={kind}, "
              f"msa={'on' if use_msa else 'off'}, force={args.force}")

        selected_ids = [str(item.get("id", "")).strip() for item in sequences]
        selected_ids = [sid for sid in selected_ids if sid]

        n_done = n_cached = n_error = 0
        for i, item in enumerate(sequences, start=1):
            seq_id = str(item.get("id", "")).strip()
            sequence = str(item.get("sequence", "")).strip()
            if not seq_id or not sequence:
                print(f"[{i}/{len(sequences)}] skipped (missing id or sequence)")
                continue
            print(f"[{i}/{len(sequences)}] {seq_id} — submitting…", flush=True)
            try:
                job, should_run, submit_msg = boltz_backend.submit_or_cache(
                    entry, seq_id, sequence, api_key=api_key,
                    smiles=smiles, smiles_label=args.smiles_label,
                    use_msa=use_msa, params=params, force=args.force)
                if should_run:
                    final = boltz_backend.run_prediction(entry, job["job_key"])
                else:
                    final = job
                status = final.get("status", "")
                if status in ("done",):
                    n_done += 1
                elif status in ("cached",):
                    n_cached += 1
                else:
                    n_error += 1
                ptm = final.get("ptm")
                plddt = final.get("plddt")
                score = ""
                if ptm is not None or plddt is not None:
                    score = f" pTM={ptm} pLDDT={plddt}"
                err = f" — {final.get('error')}" if final.get("error") else ""
                print(f"[{i}/{len(sequences)}] {seq_id} — {status}{score}{err}")
            except Exception as exc:  # keep going on a per-sequence failure
                n_error += 1
                print(f"[{i}/{len(sequences)}] {seq_id} — error: {exc}", file=sys.stderr)

        print(f"predictions: {n_done} done, {n_cached} cached, {n_error} error(s)")

        # RMSD.
        if args.rmsd:
            if not args.rmsd_reference:
                abort("--rmsd requires --rmsd-reference <sequence id>.")
            methods = ["seq", "ce"] if args.rmsd_method == "both" else [args.rmsd_method]
            query_ids = selected_ids if args.rmsd_scope == "selected" else None
            print(f"RMSD: reference={args.rmsd_reference} rank={args.rmsd_ref_rank} "
                  f"methods={','.join(methods)} scope={args.rmsd_scope}")
            res = rmsd_backend.calculate_rmsds(
                entry, args.rmsd_reference, int(args.rmsd_ref_rank),
                query_ids=query_ids, methods=methods)
            print(f"RMSD: {res['n_new']} calculated, {res['n_cached']} from cache")
            if res.get("columns"):
                print("RMSD columns: " + ", ".join(res["columns"]))

    except (SSEError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"structures -> {entry.structures_dir}")
    print(f"datafile   -> {entry.datafile_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
