#!/usr/bin/env python3
"""
sse_coordinates.py — add coordinate columns to an SSE datafile.

Embeds the entry's sequences (or structures) and reduces them to coordinates,
then merges the result into the datafile as `coordinate`-typed columns. The
embedder and reducer are independent plug-ins (see sse_tools/embedders/ and
sse_tools/reducers/), selected with --embedder / --reducer.

Coordinate columns are named `<tag>_<LABEL><n>` where `tag` encodes everything
that defines the embedding (model variant, pooling) and `LABEL` is the reducer's
component label, e.g. `esmc600m_mean_PC1`, `saprot_mean_UMAP1`. This keeps
multiple coordinate systems distinct in one datafile and groupable in the
visualizer. `--label` overrides the tag.

Rows the embedder cannot process (e.g. a Foldseek hit with no Cα, or the query)
cause a loud abort by default; `--include_empty` proceeds and leaves those rows
with empty coordinate cells (unplottable in that system).
"""

import argparse
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sse_tools.common import (SSEError, abort, COL_ID, TYPE_COORDINATE,
                              read_datafile, read_manifest, merge_columns,
                              entry_dir, ENTRIES_DIR, INITIAL_FILES_DIR)
from sse_tools.embedders import REGISTRY as EMBEDDERS, EmbedContext
from sse_tools.reducers import REGISTRY as REDUCERS

TOOL_NAME = "sse_coordinates"
TOOL_VERSION = "1.0.0"
DATAFILE_SUFFIX = ".sse.tsv"


def resolve_entry(entry: str, entries_dir: Path):
    """Accept an entry stem or a path to the datafile. Returns (datafile, dir, stem)."""
    p = Path(entry)
    if p.is_file():
        edir = p.parent
        stem = p.name[:-len(DATAFILE_SUFFIX)] if p.name.endswith(DATAFILE_SUFFIX) else p.stem
        return p, edir, stem
    edir = entry_dir(entry, entries_dir)
    datafile = edir / f"{entry}{DATAFILE_SUFFIX}"
    if not datafile.exists():
        abort(f"datafile not found: {datafile}. Pass an entry stem (looked up in "
              f"{entries_dir}) or a path to the .sse.tsv.")
    return datafile, edir, entry


def resolve_foldseek_json(manifest, edir, initial_dir, override):
    if override:
        p = Path(override)
        if not p.exists():
            abort(f"--foldseek-json not found: {p}")
        return p
    src = (manifest or {}).get("source_file")
    if src:
        cand = Path(initial_dir) / src
        if cand.exists():
            return cand
    return None


def report_skips(skip, total):
    n = sum(len(v) for v in skip.values())
    if n == 0:
        return None
    lines = [f"{n} of {total} rows cannot be embedded:"]
    for reason, ids in skip.items():
        if ids:
            sample = ", ".join(map(str, ids[:5]))
            lines.append(f"  {reason}: {len(ids)}  (e.g. {sample}"
                         f"{' ...' if len(ids) > 5 else ''})")
    return "\n".join(lines)




def coordinate_system_columns(types: dict, tag: str, reducer_label: str) -> list:
    """Existing coordinate columns belonging to this tag+reducer system."""
    prefix = f"{tag}_{reducer_label}"
    pattern = re.compile(rf"^{re.escape(prefix)}\d+$")
    return [c for c, t in types.items()
            if t == TYPE_COORDINATE and pattern.match(str(c))]



def build_parser():
    p = argparse.ArgumentParser(
        prog="sse_coordinates.py",
        description="Embed an SSE datafile and merge coordinate columns into it.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="examples:\n"
               "  sse_coordinates.py akr\n"
               "  sse_coordinates.py oleD --embedder saprot --include_empty\n"
               "  sse_coordinates.py akr --embedder esmc --esmc-model esmc_300m "
               "--pooling max --n-components 20\n",
    )
    p.add_argument("entry", help="entry stem (looked up in entries/) or a path to the .sse.tsv.")
    p.add_argument("--embedder", default="esmc", choices=sorted(EMBEDDERS),
                   help="embedding model (default: esmc).")
    p.add_argument("--reducer", default="pca", choices=sorted(REDUCERS),
                   help="dimensionality reduction (default: pca).")
    p.add_argument("--pooling", default="mean", choices=["mean", "max", "min"],
                   help="pooling over residues (default: mean).")
    p.add_argument("--device", default="auto",
                   choices=["auto", "cuda", "mps", "cpu"],
                   help="compute device (default: auto -> cuda > mps > cpu). "
                        "'mps' is Apple Silicon. An explicit unavailable device aborts.")
    p.add_argument("--esmc-model", default="esmc_600m",
                   choices=["esmc_300m", "esmc_600m"], help="ESM-C variant.")
    p.add_argument("--prostt5-checkpoint", default="Rostlab/ProstT5")
    p.add_argument("--saprot-checkpoint", default="westlake-repl/SaProt_650M_AF2")
    p.add_argument("--n-components", type=int, default=10,
                   help="number of coordinate components (default: 10; "
                        "use --n-components 2 for a UMAP landscape).")
    p.add_argument("--umap-neighbors", type=int, default=15,
                   help="UMAP n_neighbors (default: 15).")
    p.add_argument("--umap-min-dist", type=float, default=0.1,
                   help="UMAP min_dist (default: 0.1).")
    p.add_argument("--umap-metric", default="euclidean",
                   help="UMAP metric (default: euclidean).")
    p.add_argument("--tsne-perplexity", type=float, default=30.0,
                   help="t-SNE perplexity (default: 30; auto-capped below n_samples).")
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--write-every", type=int, default=1000)
    p.add_argument("--max-residues", type=int, default=1500,
                   help="structure embedders: skip targets longer than this.")
    p.add_argument("--foldseek-json", help="structure embedders: source Foldseek JSON "
                                           "(default: from the manifest, in initial_files/).")
    p.add_argument("--label", help="override the coordinate-system tag (column prefix).")
    p.add_argument("--include_empty", action="store_true",
                   help="proceed even if some rows cannot be embedded (they get "
                        "empty coordinate cells instead of aborting).")
    p.add_argument("--rereduce", action="store_true",
                   help="reuse the cached embedding matrix, rerun the reducer, and overwrite "
                        "an existing coordinate system with the same tag/reducer prefix.")
    p.add_argument("--reembed", action="store_true",
                   help="discard cached embeddings for this tag, re-embed from scratch, rerun "
                        "the reducer, and overwrite colliding coordinate columns.")
    p.add_argument("--force", action="store_true",
                   help="overwrite an existing coordinate system by reusing cached embeddings "
                        "and rerunning only the reducer; equivalent to --rereduce. Use "
                        "--reembed to recompute embeddings.")
    p.add_argument("--entries-dir")
    p.add_argument("--initial-files-dir")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    entries_dir = Path(args.entries_dir) if args.entries_dir else ENTRIES_DIR
    initial_dir = Path(args.initial_files_dir) if args.initial_files_dir else INITIAL_FILES_DIR

    try:
        datafile, edir, stem = resolve_entry(args.entry, entries_dir)
        df, types = read_datafile(datafile)
        manifest_path = edir / "logs" / f"{stem}{DATAFILE_SUFFIX[:-4]}.manifest.json"
        manifest = read_manifest(manifest_path) if manifest_path.exists() else {}

        embedder = EMBEDDERS[args.embedder]
        reducer = REDUCERS[args.reducer]
        tag = args.label or embedder.tag(args)
        overwrite_coordinates = bool(args.force or args.reembed or args.rereduce)
        force_embedding = bool(args.reembed)

        existing_system_cols = coordinate_system_columns(types, tag, reducer.label)
        if existing_system_cols and not overwrite_coordinates:
            abort(f"coordinate system already exists for prefix "
                  f"'{tag}_{reducer.label}' ({len(existing_system_cols)} column(s): "
                  f"{existing_system_cols[:5]}"
                  f"{' ...' if len(existing_system_cols) > 5 else ''}). "
                  "Use --rereduce or --force to reuse cached embeddings and replace it, "
                  "or --reembed to recompute embeddings as well.")

        # Structure embedders are Foldseek-only and need the source JSON.
        foldseek_json = None
        if embedder.requires_structure:
            if manifest.get("source_type") != "fs":
                abort(f"--embedder {args.embedder} is structure-based and needs a "
                      f"Foldseek entry (source_type 'fs'); this entry is "
                      f"'{manifest.get('source_type')}'. Use --embedder esmc.")
            foldseek_json = resolve_foldseek_json(manifest, edir, initial_dir,
                                                  args.foldseek_json)

        ctx = EmbedContext(entry_dir=edir, stem=stem,
                           source_type=manifest.get("source_type", ""),
                           source_file=manifest.get("source_file", ""),
                           initial_files_dir=initial_dir, foldseek_json=foldseek_json)

        print(f"preparing {args.embedder} input for {len(df)} rows ...")
        entries, skip = embedder.prepare(df, ctx, args)

        msg = report_skips(skip, len(df))
        if msg and not args.include_empty:
            abort(msg + "\nRe-run with --include_empty to embed the rest and leave "
                  "these rows without coordinates in this system.")
        if msg:
            print(msg)
        if not entries:
            abort("no rows could be embedded; nothing to do.")

        emb_dir = edir / "embeddings"
        emb_dir.mkdir(parents=True, exist_ok=True)
        emb_path = emb_dir / f"{tag}.emb.tsv"

        args.force_embedding = force_embedding
        if force_embedding:
            print(f"embedding {len(entries)} rows -> {tag} (cache ignored) ...")
        else:
            print(f"embedding {len(entries)} rows -> {tag} (reuse cache when possible) ...")
        matrix = embedder.embed(entries, emb_path, args)

        ids = matrix["ID"].astype(str).values
        X = matrix.drop(columns=["ID"]).to_numpy()
        coords, meta = reducer.reduce(X, args.n_components, args)
        k = meta.get("k", coords.shape[1])

        fig_path = edir / "figures" / f"{tag}_{reducer.name}.png"
        fig_path.parent.mkdir(parents=True, exist_ok=True)
        reducer.figure(meta, fig_path)

        coord_cols = [f"{tag}_{reducer.label}{i+1}" for i in range(k)]
        coord_df = pd.DataFrame(coords[:, :k], columns=coord_cols)
        coord_df.insert(0, COL_ID, ids)
        new_types = {c: TYPE_COORDINATE for c in coord_cols}

        params = (f"embedder={args.embedder};model="
                  f"{args.esmc_model if args.embedder == 'esmc' else args.embedder};"
                  f"pooling={args.pooling};reducer={args.reducer};n={k}")
        if existing_system_cols and overwrite_coordinates:
            print(f"replacing existing coordinate system '{tag}_{reducer.label}' "
                  f"({len(existing_system_cols)} old column(s)).")
        merge_columns(datafile, coord_df, new_types, manifest_path=manifest_path,
                      provenance_source="sse", tool=TOOL_NAME, version=TOOL_VERSION,
                      params=params, force=overwrite_coordinates,
                      drop_columns=existing_system_cols if overwrite_coordinates else None)

        _write_log(edir / "logs" / f"{tag}_{reducer.name}.log", args, tag,
                   len(df), len(entries), skip, k, datafile)

    except SSEError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    print(f"merged {k} coordinate column(s) '{tag}_{reducer.label}1..{k}' into {datafile}")
    return 0


def _write_log(path, args, tag, n_rows, n_embedded, skip, k, datafile):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"tool        : {TOOL_NAME} {TOOL_VERSION}",
        f"timestamp   : {datetime.now(timezone.utc).isoformat()}",
        f"datafile    : {datafile}",
        f"tag         : {tag}",
        f"embedder    : {args.embedder}  pooling={args.pooling}",
        f"reducer     : {args.reducer}  components={k}",
        f"rows        : {n_rows}  embedded={n_embedded}  "
        f"skipped={sum(len(v) for v in skip.values())}",
    ]
    for reason, ids in skip.items():
        if ids:
            lines.append(f"  {reason}: {len(ids)}")
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
