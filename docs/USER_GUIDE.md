# Sequence Space Explorer: practical user guide

This guide takes you from a sequence dataset to an annotated, clustered, interactive sequence-space landscape. It focuses on a reliable default workflow first, then shows where to make deliberate choices.

For every command below, run it from the repository root:

```bash
cd /path/to/SequenceSpaceExplorer-main
```

Use `python` instead of `python3` if that is how Python is exposed in your environment.

## What the pipeline produces

The pipeline maintains one central entry at `entries/<name>/`. Its main file is `<name>.sse.tsv`, a typed table with one row per sequence. Tools add columns to this table without removing or reordering its rows.

A typical workflow is:

1. Create an entry from EnzymeMiner/generic TSV, FASTA, or Foldseek JSON.
2. Add taxonomy and external measurements if available.
3. Embed the sequences and calculate PCA, UMAP, or t-SNE coordinates.
4. Optionally calculate embedding distances and clusters.
5. Explore and select sequences in the visualizer, exporting a selection for structure work.
6. Optionally predict structures with Boltz-2 and calculate RMSDs via the pipeline module.

## 1. Prepare the Python environment

Create a Python 3.10 or newer environment and install the provided dependency set:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The dependency file includes all three embedding backends. Choose a different PyTorch build when required for your CUDA, Apple Silicon, or CPU environment. Static figure export also requires compatible Plotly and Kaleido versions.

You can confirm that the scripts start correctly with:

```bash
python scripts/sse_initialization.py --help
python scripts/sse_coordinates.py --help
python scripts/sse_visualizer.py --help
```

## 2. Create an entry

Put a source file in `initial_files/`, or pass any accessible file path directly.

### From EnzymeMiner or another tab-separated table

The default ID and sequence columns are `Accession` and `Sequence`:

```bash
python scripts/sse_initialization.py my_enzymes.tsv \
  --source em \
  --name my_enzymes
```

If the column names differ:

```bash
python scripts/sse_initialization.py my_table.tsv \
  --source em \
  --id_col ProteinID \
  --seq_col Seq \
  --name my_enzymes
```

If a `Closest query` column is present, referenced sequences are marked automatically as queries. Override that behavior with one or more explicit IDs:

```bash
python scripts/sse_initialization.py my_table.tsv \
  --source em \
  --query P12345 Q67890 \
  --name my_enzymes
```

### From FASTA

```bash
python scripts/sse_initialization.py proteins.fasta \
  --source fasta \
  --name my_enzymes
```

FASTA has no automatic query detection. To mark a query, supply its complete header without the leading `>`:

```bash
python scripts/sse_initialization.py proteins.fasta \
  --source fasta \
  --query "sp|P12345|OLED_STRAN OleD reference" \
  --name my_enzymes
```

### From Foldseek web-server JSON

```bash
python scripts/sse_initialization.py foldseek_search.json \
  --source fs \
  --name my_foldseek_entry
```

The reader deduplicates targets found in multiple databases and records those database names as tags. A 100% self-hit matching the query sequence becomes the query row. If none is present, a synthetic query row is created from the search query.

### Check the result

After successful creation, inspect:

```text
entries/my_enzymes/
├── my_enzymes.sse.tsv
├── external_data/
├── figures/
├── logs/
├── msa_cache/
└── structures/
```

The initializer also computes sequence length, charge/composition measures, molecular weight, pI, aromaticity, instability index, and GRAVY. Rows with unusable non-query sequences are dropped and recorded in the manifest. An unusable query or duplicate resolved ID aborts creation without writing a partial entry.

Do not rerun initialization to add data. Use the enrichment tools below. `--force` on initialization deletes and rebuilds the entire existing entry.

## 3. Add annotations before clustering

This step is optional, but doing it before clustering lets the Tier-2 cluster analysis test those annotations for enrichment.

### Add NCBI taxonomy

```bash
python scripts/fetch_taxonomy.py my_enzymes \
  --email you@example.org
```

The default `--strategy auto` uses Foldseek-specific information when `Databases` and `taxId` columns are present; otherwise it looks up the datafile IDs as NCBI protein accessions.

Useful rerun modes are:

```bash
# Retry only unresolved network/API results
python scripts/fetch_taxonomy.py my_enzymes \
  --email you@example.org \
  --retry-failed

# Ignore previous results and fetch every row again
python scripts/fetch_taxonomy.py my_enzymes \
  --email you@example.org \
  --force
```

An NCBI API key can be passed with `--api-key` or the `NCBI_API_KEY` environment variable.

### Merge experimental or external metadata

Place the file under `entries/my_enzymes/external_data/`, then run:

```bash
python scripts/merge_external.py my_enzymes measurements.csv
```

The first external column is treated as its ID by default. Choose specific columns when needed:

```bash
python scripts/merge_external.py my_enzymes measurements.csv \
  --id-col Accession \
  --columns pI,Melting_temperature
```

If the external IDs use a different naming scheme, provide a two-column translation table. Its first column contains SSE datafile IDs and its second contains external IDs:

```bash
python scripts/merge_external.py my_enzymes measurements.csv \
  --translator id_translator.tsv
```

External rows without a matching SSE ID are dropped with a warning. SSE rows without external data remain in place with empty new cells.

Most merged data should use the default `--type label`. Use `--type coordinate` only for columns that should be selectable as plot axes.

## 4. Generate embeddings and coordinates

### Recommended first run: ESM-C plus PCA

```bash
python scripts/sse_coordinates.py my_enzymes \
  --embedder esmc \
  --esmc-model esmc_600m \
  --pooling mean \
  --reducer pca \
  --n-components 10
```

This creates:

- a raw embedding cache under `embeddings/`;
- an L2-normalized cache under `embeddings/normalized/`;
- coordinate columns such as `esmc600m_mean_PC1` in the SSE datafile;
- a diagnostic figure under `figures/`;
- a run log under `logs/`.

Normalization is enabled by default and is the recommended setting. Downstream distances and clustering preferentially use the same normalized cache.

### Add a two-dimensional UMAP using the same cache

```bash
python scripts/sse_coordinates.py my_enzymes \
  --embedder esmc \
  --reducer umap \
  --n-components 2 \
  --umap-neighbors 15 \
  --umap-min-dist 0.1
```

The matching ESM-C cache is reused; the model is not run again. PCA and UMAP coordinate systems can coexist in the datafile.

### Add t-SNE

```bash
python scripts/sse_coordinates.py my_enzymes \
  --embedder esmc \
  --reducer tsne \
  --n-components 2 \
  --tsne-perplexity 30 \
  --tsne-pca 50
```

The default PCA pre-reduction denoises and accelerates t-SNE. Set `--tsne-pca 0` to disable it.

### Choose a compute device

`--device auto` tries CUDA, then Apple MPS, then CPU. You can request one explicitly:

```bash
python scripts/sse_coordinates.py my_enzymes --device cuda
```

An explicitly requested unavailable device aborts. Large embedding jobs can be resumed because completed IDs are written to a `.part` cache during the run.

### Structure-aware coordinates

ProstT5 and SaProt are available only for entries created from Foldseek JSON because they require Foldseek target C-alpha coordinates:

```bash
python scripts/sse_coordinates.py my_foldseek_entry \
  --embedder prostt5 \
  --reducer umap \
  --n-components 2 \
  --include_empty
```

```bash
python scripts/sse_coordinates.py my_foldseek_entry \
  --embedder saprot \
  --reducer umap \
  --n-components 2 \
  --include_empty
```

The structure path is:

```text
Foldseek tCa → reconstructed backbone → 3Di sequence → model embedding
```

The 3Di representation is cached. By default, any unprocessable row aborts the coordinate run. `--include_empty` proceeds with usable structures and leaves missing coordinate cells for the rest.

### Safely rerun coordinates

- A different reducer produces a separate coordinate system and normally reuses the embedding cache.
- `--rereduce` or `--force` reuses the embedding cache and replaces an existing coordinate system with the same prefix.
- `--reembed` discards that embedding cache, recomputes it, and replaces colliding coordinates.
- `--label NAME` gives the embedding/coordinate system a custom prefix.
- `scripts/sse_coordinates_25.py` runs the normal pipeline on the first 25 rows with an isolated `_first25` suffix. It is for development checks, not final analysis.

## 5. Add distances to reference sequences

After an embedding cache exists:

```bash
python scripts/sse_esmc_distance.py my_enzymes
```

By default, one Euclidean-distance column is added for every row marked `query=True`. Choose references explicitly with:

```bash
python scripts/sse_esmc_distance.py my_enzymes \
  --query-id OleD_S1 AgepGT_S3
```

If more than one raw embedding cache exists, select one:

```bash
python scripts/sse_esmc_distance.py my_enzymes \
  --embedding esmc600m_mean
```

The normalized cache is preferred. Pass `--raw` only when raw vector geometry is intentional.

## 6. Cluster the embedding space

Clustering uses the embedding cache, not the displayed two-dimensional coordinates. By default it first reduces embeddings to 50 PCA components.

### K-means with automatic cluster count

```bash
python scripts/sse_cluster.py my_enzymes \
  --clusterer kmeans
```

The tool tests `k=2..20` and keeps the best silhouette score. Set a known value directly with:

```bash
python scripts/sse_cluster.py my_enzymes \
  --clusterer kmeans \
  --k 12
```

K-means assigns every embedded sequence to a cluster and is most appropriate when compact, roughly spherical groups are a useful model.

### HDBSCAN

```bash
python scripts/sse_cluster.py my_enzymes \
  --clusterer hdbscan \
  --min-cluster-size 50
```

HDBSCAN discovers the cluster count and can mark sequences as `noise`. Increase `--min-samples` for more conservative clusters and more noise; reduce `--min-cluster-size` if no clusters are found.

### Clustering geometry choices

```bash
# Retain enough PCs for 95% of variance
python scripts/sse_cluster.py my_enzymes \
  --clusterer kmeans \
  --pca-variance 0.95

# Cluster the complete embedding directly
python scripts/sse_cluster.py my_enzymes \
  --clusterer kmeans \
  --space full
```

The default PCA space is usually preferable, especially for HDBSCAN.

Each run adds three columns:

- `<tag>_<method>_cluster`;
- `<tag>_<method>_representative`, marking the medoid;
- `<tag>_<method>_dist_to_center`.

It also creates `cluster_analysis/<tag>_<method>/` containing cluster profiles, significant enrichment results, and central representative sequences in TSV and FASTA formats. Use `--no-analysis` to skip those files.

If taxonomy, external measurements, or Boltz scores are added later and should participate in enrichment, rerun the same clustering with `--force`.

## 7. Explore the entry interactively

Start the visualizer:

```bash
python scripts/sse_visualizer.py my_enzymes
```

Open the local URL printed by Dash, normally using port 8051. Choose another port if necessary:

```bash
python scripts/sse_visualizer.py my_enzymes --port 8060
```

### Suggested exploration workflow

1. In **Coordinates**, choose a coordinate system and its X and Y axes. Use free-axis mode only when intentionally mixing systems.
2. In **Filters**, enable relevant numeric, boolean, categorical, or tag filters.
3. Use **Colour** to apply a fixed color or a continuous data column.
4. If cluster columns exist, use **Cluster regions** to show KDE regions, concave hulls, or cluster-colored points.
5. Lasso or box-select interesting sequences on the plot.
6. Convert the selection into the working filter or save it as a persistent layer.
7. Click an individual point to inspect all its metadata.
8. Extract visible layers as CSV or export the plot as PNG, SVG, or PDF.

If a column is classified incorrectly, change its role under **Column settings** and rebuild the filter panel.

## 8. Predict structures with Boltz-2

Boltz-2 structure prediction and RMSD run as a pipeline module (`scripts/sse_boltz.py`), driven by a selection you export from the visualizer. This is a two-step workflow: select in the visualizer, then run in the pipeline.

**Step 1 — export a selection from the visualizer.** Select the points you want (lasso/box select, or click individual points), then click **Export selection for Boltz** in the selection toolbar above the plot. This writes a timestamped cache to `entries/<entry>/selections/selection_<timestamp>.json` containing the selected ids and sequences.

**Step 2 — run the Boltz-2 module.** In the pipeline control center, open the **Structure & binding** module (or run the CLI directly), which requires a valid NVIDIA API key:

1. Choose the exported selection (defaults to the most recent).
2. Enter the NVIDIA API key.
3. Keep **Generate MSA** enabled for the recommended path, or disable it for single-sequence predictions.
4. Leave the SMILES field empty for apo predictions, or enter one or more substrate SMILES (one per line) for holo complexes; optionally label the ligand.
5. Optionally adjust recycling, sampling, diffusion, and step-scale parameters.
6. Optionally enable **Compute RMSDs after prediction** (see below).
7. Run the module.

From the command line the same run is:

```bash
BOLTZ_API_KEY=nvapi-… python scripts/sse_boltz.py my_enzymes \
  --smiles "OC[C@H]1O…" --smiles-label UDP-Glc
python scripts/sse_boltz.py my_enzymes --list-selections   # inspect available caches
```

The module predicts a structure for each selected sequence, saving ranked CIF files under `structures/apo/` or `structures/holo/` and writing pTM and pLDDT score columns back into the SSE datafile. Previously completed matching predictions are reused unless **Force re-run** (`--force`) is selected. Reload the datafile in the visualizer to see the new score columns. When the job finishes, use **Open structures folder** on the job panel to browse the saved `.cif` files and logs.

## 9. Calculate RMSDs

RMSD analysis becomes useful after at least two apo structures exist, and runs as part of the same Boltz-2 module (enable **Measure structural RMSD**, or pass `--rmsd`):

1. Choose the **reference structure** — the one every other structure is aligned to. In the UI this is a dropdown of the analyzed sequences (the ones in the chosen selection); on the CLI pass `--rmsd-reference <id>`. Set its rank with `--rmsd-ref-rank`.
2. Choose what to compare against: all predicted apo structures (`--rmsd-scope all`) or only this selection's sequences (`--rmsd-scope selected`).
3. Choose sequence-guided superposition, structure-based CE alignment, or both (`--rmsd-method`).

RMSD is computed by Kabsch superposition of aligned Cα atoms. The resulting `RMSD_vs_<reference>_r<rank>_<method>` columns are appended to the entry's `.sse.tsv` (and logged to `logs/rmsd_log.csv`), and become available for filtering and coloring after you reload the datafile.

## 10. A recommended end-to-end command sequence

For a generic EnzymeMiner-style dataset:

```bash
# 1. Create the entry
python scripts/sse_initialization.py my_enzymes.tsv \
  --source em \
  --name my_enzymes

# 2. Add annotations used by later cluster enrichment
python scripts/fetch_taxonomy.py my_enzymes \
  --email you@example.org
python scripts/merge_external.py my_enzymes measurements.csv

# 3. Create a general coordinate basis and a 2D landscape
python scripts/sse_coordinates.py my_enzymes \
  --embedder esmc \
  --reducer pca \
  --n-components 10
python scripts/sse_coordinates.py my_enzymes \
  --embedder esmc \
  --reducer umap \
  --n-components 2

# 4. Add reference distances and compare clustering assumptions
python scripts/sse_esmc_distance.py my_enzymes
python scripts/sse_cluster.py my_enzymes --clusterer kmeans
python scripts/sse_cluster.py my_enzymes \
  --clusterer hdbscan \
  --min-cluster-size 50

# 5. Explore and select; use "Export selection for Boltz" to stage sequences
python scripts/sse_visualizer.py my_enzymes

# 6. Predict structures and calculate RMSDs for the exported selection
BOLTZ_API_KEY=nvapi-… python scripts/sse_boltz.py my_enzymes --rmsd --rmsd-reference OleD_S1
```

## Safety and reproducibility rules

- Keep a query/reference sequence usable; initialization will not silently drop it.
- Prefer default L2 normalization so coordinates, distances, and clustering share one geometry.
- Avoid running two independent datafile-writing commands on the same entry simultaneously.
- Treat initialization `--force` as destructive: it recreates the whole entry.
- Other `--force` flags generally replace only the named output columns.
- Do not edit the Type row casually. It controls which columns are identifiers, labels, and coordinates.
- Use the manifest and tool logs to recover the model, normalization, reduction, and clustering choices used for each output.

For exact contracts, defaults, output names, and troubleshooting, see [PIPELINE_REFERENCE.md](PIPELINE_REFERENCE.md).
