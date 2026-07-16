"""Tier-2 quantitative analysis of an SSE clustering.

Where sse_cluster.py adds the three per-sequence `label` columns (cluster id,
medoid flag, distance-to-centre), this module turns those clusters into the
downstream artefacts a user actually reasons about: what each cluster *is*, what
is statistically special about it, and a small diverse panel of real sequences
per cluster to seed an MSA or structure run.

It runs as the final step of sse_cluster.py (so the outputs are regenerated every
time a clustering is (re)computed and can never drift out of sync with the
columns in the datafile), but everything here is pure library code driven through
`analyze_clusters` — no argparse, no datafile writes. Outputs are written to
    entries/<entry>/cluster_analysis/<tag>_<method>/
namespaced by clustering so kmeans and hdbscan results coexist; each run
overwrites its own directory.

Three files:
  cluster_profiles.tsv   one row per cluster (+ a 'background' row for the whole
                         embedded set, and a 'noise' row if any): the diagnostic
                         numbers (size, fraction, silhouette, mean dist-to-centre)
                         AND the feature profile (median of each curated numeric
                         feature, dominant category of each curated categorical).
  enrichment.tsv         what is *significantly* special about each cluster vs the
                         background: numeric via Mann-Whitney U, categorical via
                         one-sided Fisher exact, Benjamini-Hochberg FDR across all
                         tests, only rows passing the FDR threshold are written.
  representatives.tsv    the 5 sequences nearest each cluster centre (a diverse,
  representatives.fasta  non-redundant panel for Boltz / wet-lab), with sequences.
                         Rank 1 (is_medoid=True) is the cluster medoid.

Geometry note: distances (cohesion, representative ranking) are computed in the
SAME reduced space the clustering ran in (the PCA/full matrix passed in), so they
are consistent with the medoid/dist_to_center columns sse_cluster.py wrote.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from .common import COL_ID, COL_SEQ

# --------------------------------------------------------------------- curated features
# Enrichment and profiling run only over these columns (chosen over auto-detection
# so free-text / near-unique columns like Annotation, Organism, GI, Swiss-Prot,
# Essential residues never generate meaningless tests). Edit these lists to change
# what the analysis reports. Columns absent from a given datafile are skipped
# silently; a curated categorical column that turns out to be near-unique in the
# data is also skipped (see MAX_CATEGORIES).

CURATED_NUMERIC = [
    "Sequence length", "length", "Optimum pH", "Melting temperature",
    "Optimum temp.", "Relative aggregation propensity",
    "MW", "pI", "GRAVY", "aromaticity", "instability_index", "net_charge_pH7",
    "acidic_count", "basic_count", "acidic_ratio", "basic_ratio",
    "ED_RK_ratio", "ED_IK_ratio",
    "boltz_apo_ptm", "boltz_apo_plddt",
    "Identity closest query", "Identity closest known", "Identity closest all",
]

CURATED_CATEGORICAL = [
    "Kingdom", "Solubility", "Salinity", "Temp. range", "Biotic relationship",
    "Disease", "Transmembrane", "tax_status",
    "superkingdom", "phylum", "class", "order", "family", "genus",
    "PSPG 1", "PSPG 2",
]

# Tunables (kept here rather than plumbed through the CLI to keep sse_cluster.py's
# surface small; they are stable analysis choices, not per-run knobs).
NA_TOKENS = {"", "na", "nan", "none", "null", "-", "n/a"}
MAX_CATEGORIES = 100      # skip a categorical feature with more distinct values
MIN_CLUSTER_FOR_TEST = 5  # do not run enrichment tests on tinier clusters
MIN_CATEGORY_COUNT = 3    # smallest in-cluster category count worth a Fisher test
SILHOUETTE_CAP = 4000     # cap points fed to silhouette_samples (O(N^2) memory)
DEFAULT_TOP_N = 5         # representatives nearest centre per cluster
DEFAULT_FDR = 0.05


# ----------------------------------------------------------------- numeric coercion

def _to_numeric(series: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    """Coerce a raw string column to float; return (values, present_mask).

    Empty strings and common NA tokens become NaN; anything unparseable is NaN
    too. `present_mask` marks the finite entries, so callers never test on NaN.
    """
    cleaned = series.astype(str).str.strip()
    cleaned = cleaned.where(~cleaned.str.lower().isin(NA_TOKENS), other=np.nan)
    vals = pd.to_numeric(cleaned, errors="coerce").to_numpy(dtype=float)
    return vals, np.isfinite(vals)


def _categories(series: pd.Series) -> pd.Series:
    """Normalized categorical values with NA tokens dropped (empty = absent)."""
    cleaned = series.astype(str).str.strip()
    return cleaned.where(~cleaned.str.lower().isin(NA_TOKENS))


# ------------------------------------------------------------------ analysis frame

def build_frame(df, ids, labels, dist, representative):
    """Align datafile metadata to the embedded/clustered rows.

    `ids`/`labels`/`dist`/`representative` are in embedding-matrix order (the
    order sse_cluster.py clustered). Returns a DataFrame in that order carrying
    the metadata plus `_cluster` (int, -1 = noise), `_dist`, `_representative`.
    """
    meta = df.drop_duplicates(subset=COL_ID).set_index(COL_ID)
    # reset_index restores COL_ID as a column (index name is COL_ID after set_index).
    frame = meta.reindex([str(i) for i in ids]).reset_index()
    frame["_cluster"] = np.asarray(labels, dtype=int)
    frame["_dist"] = np.asarray(dist, dtype=float)
    frame["_representative"] = np.asarray(representative, dtype=bool)
    return frame


def present_features(frame):
    """Curated feature columns actually usable in this datafile.

    Returns (numeric, categorical). A categorical column with more than
    MAX_CATEGORIES distinct values is dropped as too fine-grained to enrich.
    """
    numeric = [c for c in CURATED_NUMERIC if c in frame.columns]
    categorical = []
    for c in CURATED_CATEGORICAL:
        if c not in frame.columns:
            continue
        if _categories(frame[c]).nunique(dropna=True) <= MAX_CATEGORIES:
            categorical.append(c)
    return numeric, categorical


# ------------------------------------------------------------------- silhouette

def per_cluster_silhouette(X, labels):
    """Mean silhouette per cluster (noise excluded), sampled for large N.

    silhouette_samples is O(N^2) in memory, so above SILHOUETTE_CAP points we
    score a cluster-stratified sample. Returns {cluster_id: mean_silhouette};
    clusters absent from the sample or with <2 surviving clusters get NaN.
    """
    from sklearn.metrics import silhouette_samples

    labels = np.asarray(labels)
    mask = labels != -1
    idx = np.where(mask)[0]
    if idx.size < 2 or len(set(labels[mask])) < 2:
        return {int(c): float("nan") for c in set(labels) if c != -1}

    if idx.size > SILHOUETTE_CAP:
        rng = np.random.default_rng(0)
        # Stratified: keep small clusters represented rather than swamped.
        per = max(2, SILHOUETTE_CAP // len(set(labels[mask])))
        keep = []
        for c in set(labels[mask]):
            members = idx[labels[idx] == c]
            take = members if members.size <= per else rng.choice(members, per, replace=False)
            keep.append(take)
        idx = np.concatenate(keep)

    sub_labels = labels[idx]
    if len(set(sub_labels)) < 2:
        return {int(c): float("nan") for c in set(labels) if c != -1}
    sample_sil = silhouette_samples(X[idx], sub_labels)
    out = {}
    for c in set(labels):
        if c == -1:
            continue
        vals = sample_sil[sub_labels == c]
        out[int(c)] = float(vals.mean()) if vals.size else float("nan")
    return out


# ----------------------------------------------------------------- cluster profiles

def _cluster_ids(labels):
    """Real cluster ids (noise excluded), ascending."""
    return sorted(int(c) for c in set(np.asarray(labels)) if c != -1)


def profile_table(frame, X, numeric, categorical):
    """One row per cluster: diagnostics + feature profile in a single table.

    Rows are the 'background' (whole embedded set), then each real cluster, then a
    'noise' row if the clusterer produced any. Columns:
      - diagnostics: size, fraction, silhouette (mean, per real cluster only),
        mean_dist_to_center (per real cluster only) — the numbers used to judge
        whether k / min_cluster_size is sensible;
      - profile: median of each curated numeric feature and the dominant category
        (with its fraction) of each curated categorical feature.
    silhouette / mean_dist_to_center are blank for background and noise (a
    within-cluster mean has no meaning there); their feature profile is still
    computed, so noise vs background vs clusters can be compared at a glance.
    """
    labels = frame["_cluster"].to_numpy()
    sil = per_cluster_silhouette(X, labels)
    total = len(frame)

    groups = [("background", frame, None)]
    for c in _cluster_ids(labels):
        groups.append((c, frame[frame["_cluster"] == c], c))
    noise = frame[frame["_cluster"] == -1]
    if len(noise):
        groups.append(("noise", noise, None))

    rows = []
    for name, sub, cid in groups:
        row = {"cluster": name, "size": len(sub),
               "fraction": round(len(sub) / total, 4)}
        if cid is None:                                   # background / noise
            row["silhouette"] = ""
            row["mean_dist_to_center"] = ""
        else:
            row["silhouette"] = round(sil.get(cid, float("nan")), 4)
            row["mean_dist_to_center"] = round(float(np.nanmean(sub["_dist"])), 4)
        for col in numeric:
            vals, mask = _to_numeric(sub[col])
            row[f"{col} (median)"] = (round(float(np.median(vals[mask])), 3)
                                      if mask.any() else "")
        for col in categorical:
            cats = _categories(sub[col]).dropna()
            if len(cats):
                top, n = cats.value_counts().index[0], cats.value_counts().iloc[0]
                row[f"{col} (top)"] = f"{top} ({n / len(cats):.0%})"
            else:
                row[f"{col} (top)"] = ""
        rows.append(row)
    return pd.DataFrame(rows)


# -------------------------------------------------------------------- enrichment

def _bh_fdr(pvals):
    """Benjamini-Hochberg q-values (scipy ships this since 1.11)."""
    from scipy.stats import false_discovery_control
    if not len(pvals):
        return np.array([])
    return false_discovery_control(np.asarray(pvals, dtype=float), method="bh")


def enrichment_table(frame, numeric, categorical, fdr=DEFAULT_FDR):
    """Significant per-cluster enrichment vs background (item 2).

    Numeric feature per cluster -> Mann-Whitney U (two-sided) of in-cluster vs
    the rest. Categorical feature per (cluster, category with >= MIN_CATEGORY_COUNT
    members) -> one-sided Fisher exact for over-representation. All raw p-values
    share one Benjamini-Hochberg correction; only rows with q < `fdr` are
    returned, most significant first.
    """
    from scipy.stats import mannwhitneyu, fisher_exact

    labels = frame["_cluster"].to_numpy()
    clusters = [c for c in _cluster_ids(labels)
                if (labels == c).sum() >= MIN_CLUSTER_FOR_TEST]
    tests = []  # each: dict of output fields, plus 'p' used for FDR

    for col in numeric:
        vals, present = _to_numeric(frame[col])
        for c in clusters:
            in_mask = (labels == c) & present
            out_mask = (labels != c) & (labels != -1) & present
            x, y = vals[in_mask], vals[out_mask]
            if len(x) < MIN_CLUSTER_FOR_TEST or len(y) < MIN_CLUSTER_FOR_TEST:
                continue
            try:
                mw = mannwhitneyu(x, y, alternative="two-sided")
            except ValueError:                            # all-identical values
                continue
            p = float(mw.pvalue)
            med_in, med_out = float(np.median(x)), float(np.median(y))
            # Rank-biserial effect size from U: sign shows direction, |r| the size.
            rbc = 2 * float(mw.statistic) / (len(x) * len(y)) - 1
            tests.append({
                "cluster": c, "feature": col, "test": "mann-whitney",
                "category": "", "direction": "higher" if med_in > med_out else "lower",
                "in_cluster": round(med_in, 3), "background": round(med_out, 3),
                "effect_size": round(rbc, 3), "n_in_cluster": int(len(x)),
                "p": p,
            })

    for col in categorical:
        cats = _categories(frame[col])
        known = cats.notna().to_numpy()
        for c in clusters:
            in_mask = (labels == c) & known
            out_mask = (labels != c) & (labels != -1) & known
            n_in, n_out = int(in_mask.sum()), int(out_mask.sum())
            if n_in < MIN_CLUSTER_FOR_TEST or n_out < MIN_CLUSTER_FOR_TEST:
                continue
            in_vals, out_vals = cats[in_mask], cats[out_mask]
            for value, a in in_vals.value_counts().items():
                if a < MIN_CATEGORY_COUNT:
                    continue
                b = n_in - a
                cc = int((out_vals == value).sum())
                d = n_out - cc
                # One-sided: is `value` over-represented inside this cluster?
                odds, p = fisher_exact([[a, b], [cc, d]], alternative="greater")
                frac_in, frac_out = a / n_in, cc / n_out if n_out else 0.0
                if frac_in <= frac_out:                   # only report enrichment
                    continue
                tests.append({
                    "cluster": c, "feature": col, "test": "fisher",
                    "category": value, "direction": "enriched",
                    "in_cluster": round(frac_in, 3), "background": round(frac_out, 3),
                    "effect_size": (round(float(odds), 3) if np.isfinite(odds) else "inf"),
                    "n_in_cluster": int(a), "p": float(p),
                })

    if not tests:
        return pd.DataFrame(columns=[
            "cluster", "feature", "test", "category", "direction", "in_cluster",
            "background", "effect_size", "n_in_cluster", "p_value", "q_value"])

    q = _bh_fdr([t["p"] for t in tests])
    out = pd.DataFrame(tests)
    out["p_value"] = out["p"]
    out["q_value"] = q
    out = out.drop(columns="p")
    out = out[out["q_value"] < fdr].copy()
    out["p_value"] = out["p_value"].map(lambda v: f"{v:.2e}")
    out["q_value"] = out["q_value"].map(lambda v: f"{v:.2e}")
    order = ["cluster", "feature", "test", "category", "direction", "in_cluster",
             "background", "effect_size", "n_in_cluster", "p_value", "q_value"]
    return out.sort_values("q_value").reset_index(drop=True)[order]


# ------------------------------------------------------------------ representatives

def representatives_table(frame, top_n=DEFAULT_TOP_N):
    """The `top_n` sequences nearest each cluster centre.

    A small diverse-but-central panel per cluster (the medoid plus its nearest
    neighbours) — a non-redundant set to seed Boltz / wet-lab. Returns
    (table, fasta_string); sequences are included when the datafile has them.
    """
    seq_available = COL_SEQ in frame.columns
    rows, fasta = [], []
    for c in _cluster_ids(frame["_cluster"].to_numpy()):
        sub = frame[frame["_cluster"] == c].sort_values("_dist").head(top_n)
        for rank, (_, m) in enumerate(sub.iterrows(), start=1):
            seq = str(m[COL_SEQ]) if seq_available else ""
            rows.append({
                "cluster": c, "rank": rank, "id": str(m[COL_ID]),
                "dist_to_center": round(float(m["_dist"]), 4),
                "is_medoid": bool(m["_representative"]),
                "sequence": seq,
            })
            if seq:
                fasta.append(f">{m[COL_ID]} cluster={c} rank={rank} "
                             f"dist={float(m['_dist']):.4f}\n{seq}")
    return pd.DataFrame(rows), "\n".join(fasta) + ("\n" if fasta else "")


# --------------------------------------------------------------------- orchestrator

def analyze_clusters(df, ids, X, labels, dist, representative, *, out_dir,
                     top_n=DEFAULT_TOP_N, fdr=DEFAULT_FDR):
    """Run every Tier-2 analysis and write the three output files to `out_dir`.

    Inputs mirror what sse_cluster.py already holds after clustering:
      df               the datafile rows (metadata + sequences), any order.
      ids              embedding-matrix IDs, in clustered order.
      X                the reduced matrix that was clustered (for distances).
      labels/dist/representative  the per-row cluster outputs (matrix order).

    Returns a dict of {name: written Path} for the caller to print.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    frame = build_frame(df, ids, labels, dist, representative)
    numeric, categorical = present_features(frame)

    if not _cluster_ids(labels):
        # No real clusters (e.g. HDBSCAN found only noise): nothing to analyse.
        return {}

    written = {}

    def _write(name, table):
        path = out_dir / name
        table.to_csv(path, sep="\t", index=False)
        written[name] = path

    _write("cluster_profiles.tsv", profile_table(frame, X, numeric, categorical))
    _write("enrichment.tsv", enrichment_table(frame, numeric, categorical, fdr))

    rep_table, rep_fasta = representatives_table(frame, top_n)
    _write("representatives.tsv", rep_table)
    if rep_fasta:
        (out_dir / "representatives.fasta").write_text(rep_fasta, encoding="utf-8")
        written["representatives.fasta"] = out_dir / "representatives.fasta"

    return written
