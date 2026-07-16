# Sequence Space Explorer: pipeline reference

This document describes the architecture, data contracts, script behavior, choices, outputs, caches, and failure modes of the Sequence Space Explorer pipeline. For a shorter task-oriented walkthrough, see [USER_GUIDE.md](USER_GUIDE.md).

## 1. Architectural overview

The pipeline is entry-centered. An entry is a directory containing one typed `.sse.tsv` datafile plus caches, figures, structures, logs, and derived analyses.

The `.sse.tsv` file is the source of truth for sequence-level data. Initialization creates its rows. Every subsequent column-producing tool uses the shared additive merge contract in `sse_tools/common.py`.

The main stages are:

| Stage | Entry point | Primary result |
|---|---|---|
| Entry creation | `scripts/sse_initialization.py` | Typed SSE table and manifest |
| Annotation | `scripts/fetch_taxonomy.py` | Taxonomy label columns |
| External merge | `scripts/merge_external.py` | User-provided label or coordinate columns |
| Embedding and reduction | `scripts/sse_coordinates.py` | Embedding cache and coordinate columns |
| Query distance | `scripts/sse_esmc_distance.py` | Distance-to-query label columns |
| Clustering | `scripts/sse_cluster.py` | Cluster, representative, and center-distance columns |
| Structure prediction and comparison | `scripts/sse_boltz.py` (`sse_tools/boltz.py`, `sse_tools/rmsd.py`) | CIF structures, pTM/pLDDT columns, RMSD columns |
| Exploration | `scripts/sse_visualizer.py` | Interactive plot, layers, selections, selection export, and exports |

The tools are intentionally modular:

- readers live in `sse_tools/readers/`;
- embedders live in `sse_tools/embedders/`;
- reducers live in `sse_tools/reducers/`;
- clusterers live in `sse_tools/clusterers/`;
- taxonomy strategies live in `sse_tools/taxonomy/`.

Each family has an explicit registry in its `__init__.py`. Adding a plugin generally means implementing its base contract and registering one name; the orchestration script does not need source-specific branching.

## 2. Entry and datafile contract

### 2.1 Entry resolution

Most commands accept `ENTRY` in any of these forms:

- a direct path to a `.sse.tsv` file;
- an entry directory containing exactly one `.sse.tsv` file;
- a bare stem resolved as `entries/<stem>/<stem>.sse.tsv`.

Commands with `--entries-dir` can use another entry root.

### 2.2 Entry directory layout

The complete layout grows as tools are used:

```text
entries/<stem>/
├── <stem>.sse.tsv
├── embeddings/
│   ├── <embedding-tag>.emb.tsv
│   ├── <embedding-tag>.emb.tsv.part
│   ├── normalized/
│   │   └── <embedding-tag>.emb.tsv
│   └── 3di_cache/
│       └── <foldseek-source>_3di.tsv
├── external_data/
├── figures/
│   ├── <tag>_<reducer>.png
│   └── sse_figure_<timestamp>.<format>
├── cluster_analysis/
│   └── <tag>_<clusterer>/
│       ├── cluster_profiles.tsv
│       ├── enrichment.tsv
│       ├── representatives.tsv
│       └── representatives.fasta
├── structures/
│   ├── apo/<sequence-id>/<sequence-id>_Rank_<n>.cif
│   └── holo/<sequence-id>__<ligand>__<hash>/...
├── msa_cache/
│   └── <sequence-hash>.a3m
├── selections/
│   └── selection_<timestamp>.json
└── logs/
    ├── <stem>.sse.manifest.json
    ├── <tag>_<reducer>.log
    ├── boltz_log.csv
    ├── rmsd_log.csv
    ├── jobs.json
    └── layers.json
```

The `.part` embedding file exists only while an interrupted or active embedding run is resumable.

### 2.3 SSE TSV structure

The physical file contains:

1. a normal TSV header row;
2. a Type row;
3. data rows.

Allowed Type tokens are:

| Type | Meaning |
|---|---|
| `id` | The unique row identifier. Exactly one is required and it must be the first column. |
| `label` | Metadata, numeric values, categories, booleans, tags, sequences, and analysis results. |
| `coordinate` | Numeric columns that can be used as plot axes. |

Reserved columns are:

- `id` — resolved unique identifier;
- `Sequence` — amino-acid sequence;
- `query` — query/reference membership.

All data are read as strings. Numeric-looking label values are classified as continuous by the visualizer.

### 2.4 Additive merge invariants

Every column-adding tool uses `merge_columns()` with these guarantees:

- base rows are never reordered, removed, or multiplied;
- both base and incoming IDs must be unique;
- the join is one-to-one and left-sided;
- missing incoming values become empty cells;
- reserved columns cannot be overwritten or dropped;
- a column-name collision aborts unless the tool authorizes replacement through `force=True`;
- replacing a complete coordinate system can remove obsolete component columns atomically;
- the datafile is written to a temporary path and moved into place;
- manifest provenance and column coverage are refreshed when a manifest exists.

The write lock protects concurrent threads inside one process, such as completing visualizer jobs. It does not coordinate separate CLI processes. Do not run independent writers against the same entry simultaneously.

### 2.5 Manifest

`logs/<stem>.sse.manifest.json` contains entry-level facts and one provenance record per column:

- source file and source type;
- creation tool and version;
- creation timestamp;
- ID-resolution notes;
- kept, dropped, and query row counts;
- name, type, source, tool, version, parameters, notes, and populated coverage for each column;
- the most recent modifying tool and timestamp.

The visualizer uses the Type row, not the manifest, as its runtime contract.

## 3. Installation and runtime dependencies

The code requires Python 3.10 or newer. Dependency groups are:

| Capability | Packages used by the code |
|---|---|
| Core table handling | `numpy`, `pandas`, `biopython` |
| Reduction and clustering | `scikit-learn`, `scipy`, `matplotlib` |
| UMAP | `umap-learn` |
| Visualizer | `dash`, `plotly` |
| Static visualizer export | `kaleido`, compatible with the installed Plotly version |
| HTTP and progress | `requests`, `tqdm` |
| ESM-C | `torch`, an `esm` build providing `esm.models.esmc` |
| ProstT5/SaProt | `transformers`; ProstT5 tokenization commonly also needs `sentencepiece` |
| Foldseek structure conversion | `mini3di`, Biopython structure modules |

Model weights are obtained through the model libraries on first use. Taxonomy, ColabFold MSA generation, and Boltz-2 require network access.

## 4. Initialization reference

### 4.1 Command

```text
python scripts/sse_initialization.py INPUT --source {em,fasta,fs} [options]
```

The initializer is bootstrap-only. If the target entry already exists, it aborts unless `--force` is supplied. Initialization `--force` deletes and rebuilds the entire entry.

### 4.2 Source readers

#### `--source em`

Input is a tab-separated EnzymeMiner or generic table.

Options:

- `--id_col`, default `Accession`;
- `--seq_col`, default `Sequence`.

All non-internal source columns are retained as labels. If `Closest query` exists, every ID appearing as a value in that column is automatically marked as a query.

An explicit `--query` matches values in the raw ID column and replaces automatic query detection entirely.

After unusable sequences are dropped, every remaining `Closest query` reference must still identify a surviving row. A dangling reference aborts creation.

#### `--source fasta`

Supported content is ordinary FASTA; the extension itself is not used for parsing.

ID resolution uses the first whitespace-separated token of the full header:

- for recognized database tags such as `sp|ACCESSION|NAME`, the accession is used;
- otherwise the first pipe-separated field is used;
- if candidate IDs collide, the next available pipe field is appended;
- unresolved collisions are caught by the downstream uniqueness check.

The remaining header text becomes `Description`. `--query` must match the complete original header exactly, excluding `>`.

FASTA has no source-derived automatic query.

#### `--source fs`

Input is Foldseek web-server JSON with a top-level result block containing `queries` and `results`.

The reader:

- creates one row per unique target;
- deduplicates targets across databases;
- aggregates database membership into the comma-separated `Databases` label;
- carries `taxName`, `taxId`, and `Description` when present;
- detects a self-hit when `seqId == 100` and the target sequence matches the query sequence;
- otherwise creates a synthetic query from the Foldseek query record.

`--query` matches raw Foldseek target IDs and overrides self-hit/synthetic query detection.

### 4.3 Query behavior

`--query VALUE [VALUE ...]` is exact and source-specific:

| Reader | Match field |
|---|---|
| `em` | Raw configured ID column |
| `fasta` | Complete FASTA header |
| `fs` | Raw Foldseek target |

Any unmatched query value aborts. When one value matches multiple rows, all matches are flagged, though ID uniqueness normally prevents this after resolution.

### 4.4 Sequence validation

A usable sequence is non-empty and contains only:

- the 20 standard amino acids;
- ambiguous `U`, `B`, `Z`, or `X`.

Rules are:

- an unusable query aborts;
- an unusable non-query row is dropped and reported in the manifest;
- duplicate resolved IDs abort;
- absence of all sequence information aborts.

### 4.5 Computed sequence features

Count-based features are available for every usable sequence:

- `length`;
- `acidic_count` and `basic_count`;
- `acidic_ratio` and `basic_ratio`;
- `ED_RK_ratio`;
- `ED_IK_ratio`.

ProtParam-derived features are computed only when the sequence contains the 20 standard amino acids:

- `net_charge_pH7`;
- `MW`;
- `pI`;
- `aromaticity`;
- `instability_index`;
- `GRAVY`.

They remain empty for sequences containing ambiguous residues.

### 4.6 Other initialization options

| Option | Behavior |
|---|---|
| `--name STEM` | Overrides the entry name; default is the source filename stem. |
| `--entries-dir DIR` | Overrides the output root. |
| `--initial-files-dir DIR` | Overrides the lookup directory for a bare input filename. |
| `--force` | Deletes and recreates an existing entry. |

Creation is atomic: the entire entry is assembled in a temporary directory and moved into place only after validation and writing succeed.

## 5. Taxonomy reference

### 5.1 Command

```text
python scripts/fetch_taxonomy.py ENTRY --email ADDRESS [options]
```

Output label columns are:

- `tax_status`;
- `taxid`;
- `tax_organism`;
- `superkingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`.

NCBI's newer `domain` rank is normalized into `superkingdom`.

### 5.2 Strategy selection

`--strategy auto` chooses `foldseek` when the datafile contains both `Databases` and `taxId`; otherwise it chooses `em`.

#### `em` strategy

Each datafile ID is sent to NCBI protein `efetch`. Both versioned and unversioned accession forms are matched. This strategy works for any entry whose IDs are NCBI protein accessions, not just EnzymeMiner input.

#### `foldseek` strategy

Per row, resolution priority is:

1. use an existing non-empty `taxId`;
2. for `gmgcl_id`, remove a Foldseek `_trun_<n>` suffix and query the GMGC unigene API;
3. for `mgnify_esm30`, record the definitive `no_taxonomy` status;
4. otherwise leave the row unresolved.

When a row belongs to several Foldseek databases, existing taxId wins over GMGC, which wins over the no-taxonomy path.

After taxId resolution, all strategies use NCBI taxonomy `efetch` to expand the lineage.

### 5.3 Status values

| Status | Meaning |
|---|---|
| `ok` | TaxId and lineage resolved. |
| `no_taxonomy` | Source is known not to have an NCBI taxonomy path. |
| `taxid_unresolved` | A taxId could not be obtained; potentially retryable. |
| `lineage_unresolved` | A taxId exists but lineage retrieval failed; potentially retryable. |

### 5.4 Resume and overwrite behavior

- Completed batches are cached under `tmp_taxonomy/` during an active run.
- After success the temporary cache is deleted.
- A normal later run seeds completed rows from existing taxonomy columns.
- `--retry-failed` retries unresolved statuses without refetching successful or definitive no-taxonomy rows.
- `--force` ignores existing results, clears temporary state, and refetches every row.

`--batch` controls NCBI batch size. `--gmgc-batch` controls GMGC batch size. `--api-key` or `NCBI_API_KEY` raises the NCBI rate limit.

## 6. External-data merge reference

### 6.1 Command

```text
python scripts/merge_external.py ENTRY FILE [options]
```

A bare external filename is looked up in `entries/<entry>/external_data/`. `.csv` implies comma separation and `.tsv` implies tab separation. Other extensions require `--delimiter`.

### 6.2 ID matching

`--id-col` selects the incoming ID column; the first column is the default.

Matching uses exact strings. There is no whitespace trimming, version removal, or case normalization.

`--translator FILE` supplies a table whose first two columns are interpreted positionally:

1. SSE datafile ID;
2. external ID.

External IDs must be unique in the translator. Blank IDs abort. External rows absent from the translator or SSE table are warned about and dropped.

### 6.3 Column choices

- `--columns a,b,c` selects a subset; otherwise all non-ID columns are merged.
- `--type label` is the default.
- `--type coordinate` marks every selected incoming column as a plot coordinate.
- `--force` replaces colliding non-reserved columns.

Incoming duplicate IDs abort because they would multiply base rows.

## 7. Embedding and coordinate reference

### 7.1 Command

```text
python scripts/sse_coordinates.py ENTRY [options]
```

The embedder and reducer are independent plugin choices. Coordinate columns are named:

```text
<embedding-tag>_<reducer-label><component-number>
```

Examples:

- `esmc600m_mean_PC1`;
- `esmc600m_mean_UMAP1`;
- `prostt5_mean_TSNE1`.

`--label` overrides the embedding tag.

### 7.2 Embedders

| Embedder | Input | Default tag | Restrictions |
|---|---|---|---|
| `esmc` | `Sequence` | `esmc600m_mean` with defaults | Works on every valid entry. |
| `prostt5` | Foldseek-derived 3Di | `prostt5_mean` | Foldseek entry and source JSON required. |
| `saprot` | Paired amino acid + 3Di tokens | `saprot_mean` | Foldseek entry and source JSON required; AA and 3Di lengths must match. |

Common choices:

- `--pooling mean|max|min`, default `mean`;
- `--device auto|cuda|mps|cpu`, default `auto`;
- `--batch-size`, default 32;
- `--write-every`, default 1000.

ESM-C choices:

- `--esmc-model esmc_300m|esmc_600m`, default `esmc_600m`.

Structure-model checkpoint defaults:

- ProstT5: `Rostlab/ProstT5`;
- SaProt: `westlake-repl/SaProt_650M_AF2`.

They can be overridden with `--prostt5-checkpoint` and `--saprot-checkpoint`.

### 7.3 Structure-input subpipeline

For ProstT5 and SaProt, the source manifest must say `source_type: fs`. The original Foldseek JSON is resolved from:

1. explicit `--foldseek-json`; or
2. the manifest's source filename under `initial_files/`.

Target `tCa` strings are parsed into C-alpha coordinates keyed by `(target ID, target sequence)`. For each row:

1. C-alpha coordinates are looked up;
2. targets longer than `--max-residues` are skipped, default 1500;
3. N, C, and C-beta positions are reconstructed from the C-alpha trace;
4. mini3di converts the backbone into a 3Di sequence;
5. the 3Di result is cached by ID and sequence.

Skip categories are `no_tca`, `too_long`, `threedi_fail`, and for SaProt `len_mismatch`.

The default response to any skipped row is to abort. `--include_empty` embeds successful rows and merges empty coordinates for skipped rows.

### 7.4 Streaming embedding cache

Raw caches use `embeddings/<tag>.emb.tsv` with columns `ID`, `0`, `1`, and so on.

The shared streaming loop:

- reuses a complete cache when it covers all requested IDs;
- extends a complete cache when new IDs are missing;
- resumes an interrupted `.part` file;
- appends buffered rows every `--write-every` sequences;
- promotes the `.part` file atomically when complete;
- deletes/recomputes the cache only for `--reembed`.

### 7.5 Normalization

`--normalize` is on by default. Each vector is divided by its L2 norm. The same normalized matrix is:

- passed to the reducer;
- written to `embeddings/normalized/<tag>.emb.tsv`;
- preferred by distance and clustering tools.

`--no-normalize` reduces raw vectors and does not write a normalized cache. If an old normalized sibling already exists, downstream tools will still prefer it unless it is removed or `--raw` is used. The coordinate script prints a warning for this mixed-geometry situation.

### 7.6 Reducers

#### PCA

`--reducer pca` produces up to `--n-components` components, bounded by sample count and embedding width. It reports explained variance and writes a bar/cumulative diagnostic plot.

#### UMAP

`--reducer umap` choices are:

- `--n-components`, default 10;
- `--umap-neighbors`, default 15 and capped below the sample count;
- `--umap-min-dist`, default 0.1;
- `--umap-metric`, default `euclidean`.

A two-dimensional landscape normally uses `--n-components 2`.

#### t-SNE

`--reducer tsne` choices are:

- `--n-components`, default inherited as 10, though 2 or 3 is normally appropriate;
- `--tsne-perplexity`, default 30 and automatically capped below sample count;
- `--tsne-pca`, default 50; zero disables PCA pre-reduction.

Barnes-Hut is used for at most three output components; higher dimensions use the exact method.

### 7.7 Rerun semantics

If columns already match `<tag>_<reducer-label><n>`, the default is to abort.

| Option | Embedding work | Coordinate work |
|---|---|---|
| New reducer | Reuse matching cache | Add a distinct coordinate system |
| `--rereduce` | Reuse cache | Replace matching coordinate system |
| `--force` | Reuse cache | Same effect as `--rereduce` |
| `--reembed` | Delete and recompute tag cache | Replace matching coordinate system |

When replacement changes the number of components, old surplus columns are removed as part of the same merge.

### 7.8 Development wrapper

`scripts/sse_coordinates_25.py` calls the normal coordinate pipeline with the hidden `--limit 25` option. It appends `_first25` to the tag so a partial cache cannot be mistaken for a full embedding.

## 8. Embedding-distance reference

### 8.1 Command

```text
python scripts/sse_esmc_distance.py ENTRY [options]
```

Despite the script name, the implementation reads the common embedding-cache contract and can operate on a selected compatible cache.

### 8.2 Cache resolution

`--embedding` accepts:

- an embedding tag;
- an embedding-cache filename;
- a direct path.

If omitted, exactly one raw cache must exist. Multiple caches require explicit selection. When a normalized sibling exists, it is used unless `--raw` is given.

### 8.3 Query resolution and output

`--query-id` accepts one or more explicit SSE IDs. Without it, every truthy `query` row is used. Missing IDs abort.

For each query, Euclidean distance is calculated in the full chosen embedding space and written as:

```text
<embedding-tag>_distance_to_<safe-query-id>
```

The query's own value is zero. `--force` replaces existing distance columns.

## 9. Clustering reference

### 9.1 Command

```text
python scripts/sse_cluster.py ENTRY [options]
```

Cache selection and normalized/raw behavior match the distance tool.

### 9.2 Clustering space

`--space pca` is the default. It denoises k-means and is especially important for HDBSCAN in wide embeddings.

PCA size is chosen by either:

- `--pca-dims`, default 50; or
- `--pca-variance`, a fraction in `(0, 1]`, for example 0.95.

`--space full` skips this PCA and clusters the embedding vectors directly.

### 9.3 K-means

With `--clusterer kmeans`:

- `--k N` fixes the cluster count and requires at least two;
- without `--k`, the tool sweeps `--k-min` through `--k-max`, defaults 2 through 20;
- each candidate is scored by silhouette;
- silhouette computation samples deterministically above 2,000 non-noise points;
- every sequence is assigned to a cluster.

K-means uses 10 initializations and `random_state=0`.

### 9.4 HDBSCAN

With `--clusterer hdbscan`:

- `--min-cluster-size`, default 50, is the smallest accepted group;
- `--min-samples` controls conservativeness and defaults to the scikit-learn behavior associated with the configured cluster size;
- unassigned points receive label `-1`, written as `noise`;
- the number of clusters is discovered rather than specified.

If no clusters are found, reduce `min-cluster-size` or `min-samples`.

### 9.5 Per-sequence outputs

The output prefix is the embedding tag unless `--label` overrides it. Three label columns are written:

| Suffix | Meaning |
|---|---|
| `_cluster` | Integer cluster label or `noise`. |
| `_representative` | `True` for the real point nearest the cluster centroid. |
| `_dist_to_center` | Euclidean distance to the centroid in the actual clustering space; empty for noise. |

The representative is therefore a medoid-like central sequence, not a synthetic centroid.

K-means and HDBSCAN use different prefixes and can coexist. `--force` replaces matching columns for the same tag and method.

### 9.6 Tier-2 cluster analysis

Unless `--no-analysis` is supplied, analysis is regenerated after clustering so it cannot drift from the cluster assignments.

Outputs are namespaced under `cluster_analysis/<tag>_<method>/`:

#### `cluster_profiles.tsv`

Rows include the complete embedded background, each real cluster, and noise when present. It contains:

- size and fraction;
- mean per-cluster silhouette;
- mean distance to center;
- medians for curated numeric features;
- dominant values for curated categorical features.

#### `enrichment.tsv`

Numeric features use two-sided Mann-Whitney U tests. Categorical values use one-sided Fisher exact tests for over-representation. Benjamini-Hochberg false-discovery correction is applied across tests, and only results passing `--fdr` are written.

Very small clusters and rare categories are excluded from testing. High-cardinality/free-text metadata are not tested automatically.

#### `representatives.tsv` and `.fasta`

The closest `--analysis-top-n` sequences to each center are exported; default is five. Rank one is the medoid.

Curated analysis columns are defined in `sse_tools/cluster_analysis.py`. They include built-in sequence properties, selected EnzymeMiner fields, taxonomy ranks, and Boltz apo scores when present.

## 10. Visualizer reference

### 10.1 Command

```text
python scripts/sse_visualizer.py ENTRY [--port 8051]
```

The visualizer reads one entry at a time. It creates missing `logs`, `figures`, `structures`, `msa_cache`, and `selections` directories. Reloading validates saved-layer IDs against the current datafile.

Selected points can be exported for structure work: the **Export selection for Boltz** button in the selection toolbar writes a timestamped cache to `selections/selection_<timestamp>.json` (containing the selected ids and sequences), which the Boltz-2 module (section 11) then imports. Structure prediction and RMSD themselves no longer run inside the visualizer; only the selection is chosen here.

### 10.2 Coordinate systems

Coordinate columns are grouped by their shared prefix and numbered component suffix.

Modes are:

- **Coordinate system mode** — X and Y must come from one detected coordinate system;
- **Advanced free-axis mode** — any two coordinate columns can be mixed.

Rows lacking either selected coordinate are not plotted in that system.

### 10.3 Label classification

Label columns are classified for UI controls as:

- continuous;
- boolean;
- categorical;
- comma-separated tag set;
- skipped.

Numeric parsing requires at least 80% of non-empty values to parse. Very high-cardinality, date-like, numeric-list, or long near-unique text columns are normally skipped. Some high-cardinality categorical columns can be rescued when they cleanly nest with another category.

**Column settings** can override the classification. **Rebuild filter panel** applies the override.

### 10.4 Filters and search

Filters are enabled individually and combined with logical AND:

- continuous ranges;
- boolean true/false;
- selected categorical values;
- intersection with selected tags;
- optional comma-separated ID/name search.

Name search also checks label columns whose names suggest names, labels, aliases, accessions, or genes.

### 10.5 Appearance and coloring

The working filter supports:

- fixed color; or
- continuous coloring by a numeric label;
- selectable Plotly color maps and reversal;
- global or current-subset color limits;
- point size, opacity, and symbol.

Background point size, query marker size/opacity/position, working-filter draw order, and theme are independently configurable.

### 10.6 Selections and layers

Plotly lasso and box selection populate a persistent selection set. A selection can:

- be cleared;
- receive a display color;
- be converted into the working ID filter.

A saved layer records:

- its matched IDs;
- filter state;
- name;
- fixed or continuous color settings;
- opacity, size, and symbol;
- creation time and visibility.

Layers are stored in `logs/layers.json`. They can be shown, hidden, loaded into the working filter, or deleted.

**Extract visible layers** produces a de-duplicated CSV of sequences appearing in visible layers and saves a log copy in the entry.

### 10.7 Cluster overlays

Any label column ending in `_cluster` is offered as a cluster source. Noise and blank categories are not treated as real regions.

Display choices are:

- KDE highest-density regions with configurable coverage;
- concave hull regions with configurable tightness;
- direct point coloring by cluster.

Regions can appear above or below points and have configurable opacity. Geometry is computed in the currently selected X/Y coordinate space and cached for the current settings.

### 10.8 Figure export

Formats are PNG, SVG, and PDF. Choices include:

- 150, 300, or 600 DPI scaling;
- legend inclusion;
- transparent or colored background;
- axis, label, and marker-edge colors;
- width and height;
- browser download or saving under `figures/`.

Export uses Plotly's static image engine and therefore requires a working Plotly/Kaleido combination.

## 11. Boltz-2 reference

Boltz prediction runs as a pipeline module, `scripts/sse_boltz.py`, over a selection cache exported from the visualizer. A validated NVIDIA API key is required, read from the `BOLTZ_API_KEY` environment variable (the pipeline UI injects it as a per-run secret). The heavy logic is unchanged and lives in `sse_tools/boltz.py`.

### 11.0 Command

```text
BOLTZ_API_KEY=nvapi-… python scripts/sse_boltz.py ENTRY \
  [--selection NAME] [--smiles S] [--smiles-label L] [--no-msa] \
  [--recycling-steps N] [--sampling-steps N] [--diffusion-samples N] [--step-scale F] [--force] \
  [--rmsd --rmsd-reference ID --rmsd-ref-rank N --rmsd-method {seq,ce,both} --rmsd-scope {all,selected}]
python scripts/sse_boltz.py ENTRY --list-selections
```

`--selection` takes a filename in `selections/` or a path; it defaults to the most recent exported cache. The module iterates the selection's sequences, predicting each in turn and printing per-sequence status. It validates that `requests` and Biopython are importable before starting.

### 11.1 Job choices

| Choice | Values |
|---|---|
| MSA | ColabFold MSA enabled by default, or no MSA |
| Complex | Apo when SMILES is empty; holo when one or more SMILES lines are supplied |
| Ligand naming | Optional display/column label; otherwise a SMILES hash is used |
| Recycling steps | Default 3 |
| Sampling steps | Default 200 |
| Diffusion samples | Default 5 |
| Step scale | Default 1.638 |
| Cache | Reuse matching completed job unless Force re-run is selected |

MSAs are cached by sequence under `msa_cache/` and capped at 1,000 sequences. If the Boltz API rejects an MSA payload with HTTP 422, the backend retries without the MSA. HTTP 429 responses back off and retry.

### 11.2 Job identity and persistence

- Apo jobs are keyed by sequence ID.
- Holo jobs are keyed by sequence ID plus a hash of normalized SMILES lines.
- Active state is stored in `logs/jobs.json`.
- Completed/error records are appended to `logs/boltz_log.csv`.
- Matching cached records can restore results without a new API request.

### 11.3 Outputs

Ranked CIF files are stored under the appropriate apo/holo directory.

The datafile receives only analysis-relevant scalar columns:

- apo: `boltz_apo_ptm`, `boltz_apo_plddt`;
- holo: `boltz_holo_<ligand-token>_ptm`, `boltz_holo_<ligand-token>_plddt`.

Detailed status and MSA bookkeeping remain in job/log files. Reload the datafile in the visualizer to surface the new score columns for filtering and coloring.

## 12. RMSD reference

RMSD runs as part of the same Boltz-2 module (`--rmsd`, or the module's "Compute RMSDs after prediction" toggle) and operates only on completed apo structures under `structures/apo/`.

### 12.1 Choices

- reference sequence (`--rmsd-reference`);
- reference prediction rank, default 0 (`--rmsd-ref-rank`);
- scope: all completed apo structures or only the exported selection (`--rmsd-scope`);
- optional prediction-rank override for each query;
- method: sequence-guided, CE structure alignment, or both.

The reference is excluded from the query list.

### 12.2 Methods

#### Sequence-guided (`seq`)

The backend parses C-alpha coordinates and residue identities from CIF files, globally aligns the sequences, iteratively rejects large structural outliers, and applies Kabsch superposition. Fewer than three aligned residues yields no RMSD.

#### Structure-based (`ce`)

Biopython's `CEAligner` performs a structure-derived alignment and reports RMSD over its aligned reference coordinates.

### 12.3 Cache and outputs

Each `(reference, query, reference rank, query rank, method)` result is cached in `logs/rmsd_log.csv`.

Columns are named:

```text
RMSD_vs_<reference-id>_r<reference-rank>_<method>
```

Results are merged with replacement enabled so later calculations can fill additional rows in the same column. The visualizer reloads after calculation.

## 13. Recommended decision rules

### 13.1 Which source reader?

- Use `em` for EnzymeMiner or any TSV with named ID and sequence columns.
- Use `fasta` when sequence plus header-derived identity is sufficient.
- Use `fs` for Foldseek web JSON, especially when structure-aware embeddings or Foldseek-specific taxonomy are desired.

### 13.2 Which embedder?

- Start with ESM-C for a general sequence-similarity landscape.
- Use ProstT5 when structural state encoded as 3Di should dominate.
- Use SaProt when both amino-acid and structural-token context should contribute.
- Compare coordinate systems rather than assuming one representation is universally superior.

### 13.3 Which reducer?

- PCA is deterministic, global, and useful for variance inspection and downstream denoising.
- UMAP is usually the most convenient nonlinear 2D neighborhood landscape.
- t-SNE emphasizes local neighborhoods; global inter-cluster distances should not be over-interpreted.

### 13.4 Which clusterer?

- K-means partitions every point and gives stable centroid-based representatives when roughly compact groups are acceptable.
- HDBSCAN is useful when clusters may be irregular and some points should remain unassigned.
- Run both when cluster stability is an analytical question, not merely a plotting choice.

### 13.5 When should annotations be added?

Add taxonomy and external metadata before clustering if they should appear in profiles or enrichment tests. If Boltz scores are generated later, rerun clustering with `--force` to regenerate analysis against the expanded datafile.

## 14. Reproducibility and provenance

Record or preserve:

- source file and explicit query choices;
- model/checkpoint, pooling, device, and normalization;
- embedding tag and cache used;
- reducer parameters and component count;
- clustering cache geometry, PCA/full space, and method parameters;
- taxonomy strategy and retry/force mode;
- Boltz MSA, ligand, and sampling parameters;
- RMSD reference, ranks, scope, and alignment method.

The manifest records column-producing tool parameters. Coordinate logs record run-level embedding/reduction facts. Cluster analysis is namespaced by embedding tag and clusterer. Boltz/RMSD logs retain job-specific details.

For an externally reproducible analysis, also pin the Python package versions and model revisions; the repository currently does not ship a lock file.

## 15. Troubleshooting

### Entry already exists

Initialization is not an update mechanism. Use taxonomy, external merge, coordinate, distance, or clustering tools to add columns. Use initialization `--force` only when full deletion and recreation are intended.

### Query not found

Query matching is exact and occurs against the source-specific raw field. For FASTA, supply the entire header rather than the resolved accession.

### Duplicate IDs

Resolve duplicates at the source. The pipeline deliberately refuses to invent identity during a merge because duplicates could multiply rows.

### No coordinate columns in the visualizer

Run `sse_coordinates.py`, or merge genuinely numeric coordinates with `merge_external.py --type coordinate`. Merely numeric `label` columns can color/filter but do not become axes.

### Coordinate system already exists

Use:

- `--rereduce` or coordinate `--force` to reuse embeddings and recalculate the reducer;
- `--reembed` when the embeddings themselves must change;
- another reducer or label to keep both systems.

### Structure embedder rejects the entry

ProstT5 and SaProt require an entry initialized with `--source fs` and access to the original Foldseek JSON containing `tCa`. Pass `--foldseek-json` if the manifest source cannot be resolved.

### Some structure rows cannot be embedded

Review the printed skip categories. Increase `--max-residues` if hardware permits, restore missing Foldseek coordinates, or use `--include_empty` to create partial coordinates.

### Multiple embedding caches found

Distance and clustering only auto-select when exactly one raw cache exists. Pass `--embedding <tag>` to make the intended representation explicit.

### Coordinates and clustering disagree

Check normalization. Coordinates may have been generated with `--no-normalize` while an older normalized sibling remained available to downstream tools. Also remember clustering normally uses PCA-reduced full embeddings, not the displayed 2D reducer coordinates.

### HDBSCAN finds only noise

Use PCA clustering space, reduce `--min-cluster-size`, or reduce `--min-samples`. Confirm the selected embedding cache represents the intended geometry.

### Taxonomy rows remain unresolved

Use `--retry-failed` for transient taxId/lineage failures. Confirm that generic-entry IDs are valid NCBI protein accessions or choose the correct strategy explicitly.

### Static figure export fails

Verify that Kaleido imports and is compatible with the installed Plotly version. Browser-based plot interaction can work even when static export is misconfigured.

### Boltz module reports no selection

Export a selection first: in the visualizer, select points and click **Export selection for Boltz**, then choose that cache in the module (or let it default to the most recent). `python scripts/sse_boltz.py ENTRY --list-selections` shows what is available.

### Boltz module aborts on the API key

The key is read from `BOLTZ_API_KEY`. In the pipeline UI, enter it in the module's API-key field; from the shell, export it before running. The module validates the key against the Boltz endpoint before starting and aborts with a clear message if it is missing or rejected.

### RMSD has no candidates

RMSD only scans completed apo structures. Holo structures are not included. Predict at least two apo sequences and confirm the selected ranks exist.

### Unexpected empty values after a merge

ID matching is exact. Check version suffixes, whitespace, and case. Use an explicit translator for external data rather than editing SSE IDs after entry creation.

## 16. Extension contracts

### Reader

Implement `reader(path, args) -> ReaderResult` with a resolved table, ID/sequence columns, source tag, query match field, optional automatic query mask, and manifest notes. Register it in `sse_tools/readers/__init__.py`.

### Embedder

Subclass `Embedder`, define a tag, prepare entries, load a model, and encode batches. Set `requires_structure` when Foldseek structure preparation is mandatory. Register it in `sse_tools/embedders/__init__.py`.

### Reducer

Subclass `Reducer`, define `name` and component `label`, implement `reduce()`, and optionally implement a diagnostic `figure()`. Register it in `sse_tools/reducers/__init__.py`.

### Clusterer

Subclass `Clusterer`, define `name` and `label`, and return integer labels plus metadata from `cluster()`. Use `-1` for noise. Centroids, representatives, and distances are handled by the orchestration script. Register it in `sse_tools/clusterers/__init__.py`.

### Taxonomy strategy

Implement positive-shape detection when possible and `resolve_taxids()` returning taxId, the no-taxonomy sentinel, or no result for unresolved rows. Register it in `sse_tools/taxonomy/__init__.py`. Lineage expansion is shared.
