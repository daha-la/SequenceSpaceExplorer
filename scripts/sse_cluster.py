#!/usr/bin/env python3
"""sse_cluster.py — add cluster-membership columns to an SSE datafile.

Clusters the entry's cached ESM-C embeddings and merges three `label` columns
per sequence:

  <tag>_<method>_cluster          cluster id (or 'noise' for HDBSCAN outliers)
  <tag>_<method>_representative   True for the medoid (most central) of a cluster
  <tag>_<method>_dist_to_center   distance to the cluster centroid (empty = noise)

The cluster column colours points in the visualizer; the representative flag
marks the sequence nearest each cluster's centre - the "most average" member,
a natural candidate to characterise the group.

Geometry follows the embedding-cache contract (sse_tools/common.resolve_embedding_cache):
by default it reads the L2-normalized cache when one exists, else the raw cache;
--raw forces raw. Clustering runs on a PCA reduction of those vectors by default
(--pca-dims / --pca-variance) because it denoises k-means and is essential for
HDBSCAN in high dimensions; --space full clusters on the embeddings directly.

examples:
  sse_cluster.py akr --clusterer kmeans                 # auto-k by silhouette
  sse_cluster.py akr --clusterer kmeans --k 12
  sse_cluster.py akr --clusterer hdbscan --min-cluster-size 50
  sse_cluster.py akr --clusterer kmeans --space full --raw --force
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sse_tools.common import (SSEError, abort, COL_ID, TYPE_LABEL, ENTRIES_DIR,
                              read_datafile, resolve_entry_path,
                              resolve_embedding_cache, merge_columns)
from sse_tools.esmc_distance import read_embedding_matrix
from sse_tools.clusterers import REGISTRY as CLUSTERERS
from sse_tools import cluster_analysis

TOOL_NAME = "sse_cluster"
TOOL_VERSION = "1.0.0"
DATAFILE_SUFFIX = ".sse.tsv"


def reduce_for_clustering(X, args):
    """Return the matrix to cluster on, plus a description for provenance."""
    if args.space == "full":
        return X, {"space": "full", "dims": X.shape[1], "variance": None}

    from sklearn.decomposition import PCA
    if args.pca_variance is not None:
        if not 0 < args.pca_variance <= 1:
            abort("--pca-variance must be in (0, 1].")
        pca = PCA(n_components=args.pca_variance, random_state=0)
    else:
        k = min(args.pca_dims, X.shape[0], X.shape[1])
        pca = PCA(n_components=k, random_state=0)
    Xr = pca.fit_transform(X.astype(np.float64))
    variance = float(pca.explained_variance_ratio_.sum())
    print(f"  PCA: {X.shape[1]} -> {Xr.shape[1]} dims "
          f"({variance * 100:.1f}% variance retained)")
    return Xr, {"space": "pca", "dims": Xr.shape[1], "variance": variance}


def cluster_geometry(X, labels):
    """Per-point distance to its cluster centroid, and the medoid flag per cluster.

    The medoid (nearest real point to the centroid) is the representative. Noise
    points (label -1) get no centre: distance NaN, representative False.
    """
    labels = np.asarray(labels)
    dist = np.full(len(labels), np.nan, dtype=float)
    representative = np.zeros(len(labels), dtype=bool)
    for c in sorted(set(labels)):
        if c == -1:
            continue
        idx = np.where(labels == c)[0]
        centroid = X[idx].mean(axis=0)
        d = np.linalg.norm(X[idx] - centroid, axis=1)
        dist[idx] = d
        representative[idx[int(np.argmin(d))]] = True
    return representative, dist


def label_strings(labels):
    return ["noise" if int(c) == -1 else str(int(c)) for c in labels]


def cluster_size_summary(labels):
    labels = np.asarray(labels)
    lines = []
    for c in sorted(set(labels)):
        name = "noise" if c == -1 else f"cluster {c}"
        lines.append(f"    {name}: {int((labels == c).sum())}")
    return "\n".join(lines)


def build_parser():
    p = argparse.ArgumentParser(
        prog="sse_cluster.py",
        description="Cluster an entry's cached ESM-C embeddings and merge "
                    "cluster columns into its datafile.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__[__doc__.index("examples:"):],
    )
    p.add_argument("entry", help="entry stem (looked up in entries/) or a path to the .sse.tsv.")
    p.add_argument("--clusterer", default="kmeans", choices=sorted(CLUSTERERS),
                   help="clustering method (default: kmeans).")
    p.add_argument("--embedding", help="embedding tag or path (auto if exactly one exists).")
    p.add_argument("--raw", action="store_true",
                   help="cluster on the raw cache even if an L2-normalized one exists.")

    p.add_argument("--space", default="pca", choices=["pca", "full"],
                   help="cluster on a PCA reduction (default) or the full embeddings.")
    p.add_argument("--pca-dims", type=int, default=50,
                   help="PCA components when --space pca (default: 50).")
    p.add_argument("--pca-variance", type=float,
                   help="instead of --pca-dims, keep enough PCs for this variance "
                        "fraction, e.g. 0.95.")

    # k-means
    p.add_argument("--k", type=int, help="k-means: fixed number of clusters "
                   "(default: auto-select by silhouette).")
    p.add_argument("--k-min", type=int, default=2, help="k-means auto: smallest k (default: 2).")
    p.add_argument("--k-max", type=int, default=20, help="k-means auto: largest k (default: 20).")

    # hdbscan
    p.add_argument("--min-cluster-size", type=int, default=50,
                   help="HDBSCAN: smallest group counted as a cluster (default: 50).")
    p.add_argument("--min-samples", type=int,
                   help="HDBSCAN: conservativeness; higher = more noise "
                        "(default: same as min-cluster-size).")

    p.add_argument("--label", help="override the tag (column prefix); default: the embedding tag.")
    p.add_argument("--force", action="store_true",
                   help="overwrite existing cluster columns for this tag+method.")

    # Tier-2 analysis (runs after the columns merge; outputs regenerate each run)
    p.add_argument("--no-analysis", action="store_true",
                   help="skip the Tier-2 cluster analysis (summary/enrichment/"
                        "shortlist/representatives/diagnostics files).")
    p.add_argument("--analysis-top-n", type=int, default=cluster_analysis.DEFAULT_TOP_N,
                   help="representatives nearest each cluster centre to export "
                        f"(default: {cluster_analysis.DEFAULT_TOP_N}).")
    p.add_argument("--fdr", type=float, default=cluster_analysis.DEFAULT_FDR,
                   help=f"enrichment FDR threshold (default: {cluster_analysis.DEFAULT_FDR}).")
    p.add_argument("--entries-dir")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    entries_dir = Path(args.entries_dir) if args.entries_dir else ENTRIES_DIR

    try:
        datafile = resolve_entry_path(args.entry, entries_dir)
        df, types = read_datafile(datafile)

        emb_path, normalized = resolve_embedding_cache(
            datafile.parent, args.embedding, prefer_normalized=not args.raw)
        print(f"using {'L2-normalized' if normalized else 'raw'} embeddings: {emb_path}")
        matrix, dimensions = read_embedding_matrix(emb_path)

        tag = args.label or emb_path.name.removesuffix(".emb.tsv")
        clusterer = CLUSTERERS[args.clusterer]
        prefix = f"{tag}_{clusterer.label}_"
        out_cols = [f"{prefix}cluster", f"{prefix}representative",
                    f"{prefix}dist_to_center"]

        existing = [c for c in out_cols if c in types]
        if existing and not args.force:
            abort(f"cluster columns already exist for '{prefix}*': {existing}. "
                  "Use --force to overwrite.")

        ids = matrix["ID"].astype(str).to_numpy()
        X = matrix[dimensions].to_numpy(dtype=np.float64)
        Xc, space_meta = reduce_for_clustering(X, args)

        labels, meta = clusterer.cluster(Xc, args)
        representative, dist = cluster_geometry(Xc, labels)

        result = pd.DataFrame({
            COL_ID: ids,
            out_cols[0]: label_strings(labels),
            out_cols[1]: representative,
            out_cols[2]: dist,
        })
        new_types = {c: TYPE_LABEL for c in out_cols}

        n_clusters = meta.get("k", meta.get("n_clusters", 0))
        method_params = (f"k={meta['k']};selected={meta['selected']}"
                         if args.clusterer == "kmeans"
                         else f"min_cluster_size={meta['min_cluster_size']};"
                              f"min_samples={meta['min_samples']}")
        sil = meta.get("silhouette", float("nan"))
        params = (f"embedding={tag};normalize={'l2' if normalized else 'none'};"
                  f"clusterer={args.clusterer};space={space_meta['space']};"
                  f"dims={space_meta['dims']};{method_params};"
                  f"n_clusters={n_clusters};silhouette={sil:.3f}")

        stem = datafile.name.removesuffix(DATAFILE_SUFFIX)
        manifest = datafile.parent / "logs" / f"{stem}.sse.manifest.json"
        merge_columns(
            datafile, result, new_types,
            manifest_path=manifest if manifest.exists() else None,
            provenance_source="sse", tool=TOOL_NAME, version=TOOL_VERSION,
            params=params,
            notes=(f"{args.clusterer} clusters in "
                   f"{'L2-normalized ' if normalized else ''}ESM-C "
                   f"{space_meta['space']} space; representative = cluster medoid"),
            force=args.force, drop_columns=existing or None,
        )

        analysis_written = {}
        if not args.no_analysis:
            # Tier-2: regenerate the quantitative analysis alongside the columns so
            # the two never drift. Namespaced by tag+method so kmeans/hdbscan coexist.
            out_dir = datafile.parent / "cluster_analysis" / f"{tag}_{clusterer.label}"
            analysis_written = cluster_analysis.analyze_clusters(
                df, ids, Xc, labels, dist, representative,
                out_dir=out_dir, top_n=args.analysis_top_n, fdr=args.fdr)

    except (SSEError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"merged {len(out_cols)} cluster column(s) into {datafile}")
    print(f"  {n_clusters} clusters over {len(result)}/{len(df)} embedded rows"
          + (f", {meta['n_noise']} noise" if args.clusterer == "hdbscan" else ""))
    print(cluster_size_summary(labels))
    if analysis_written:
        print(f"Tier-2 analysis -> {out_dir}")
        for name in sorted(analysis_written):
            print(f"    {name}")
    elif not args.no_analysis:
        print("Tier-2 analysis: no real clusters to analyse (skipped).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
