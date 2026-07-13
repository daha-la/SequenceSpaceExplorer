# scripts/

The CLI entry points for SSE — everything here is meant to be run directly
(`python scripts/<name>.py ...`), from the repo root, inside the [conda
environment](../env/requirement.yaml). Each script owns one step of the
pipeline described in the [top-level README](../README.md); the logic behind
each lives in [`sse_tools/`](../sse_tools/), which these scripts import from.

Every script that takes an `ENTRY` argument accepts it three ways: a bare
entry stem (looked up in `entries/`), a path to an entry directory, or a
direct path to its `.sse.tsv` file. Run any script with `--help` for its
full, authoritative flag reference — this document explains what each script
does and why, not a copy of its `--help` output.

| Script | Role |
|---|---|
| [`sse_initialization.py`](#sse_initializationpy) | Create a new entry from a raw source file. |
| [`sse_coordinates.py`](#sse_coordinatespy) | Embed an entry's sequences/structures and add coordinate columns. |
| [`fetch_taxonomy.py`](#fetch_taxonomypy) | Merge NCBI taxonomy lineages into an entry. |
| [`merge_external.py`](#merge_externalpy) | Merge your own external data into an entry. |
| [`sse_visualizer.py`](#sse_visualizerpy) | Launch the interactive Dash viewer for one entry. |

## `sse_initialization.py`

Turns one raw source file into a new entry: `entries/<stem>/<stem>.sse.tsv`
plus its subfolders and a provenance manifest. This is the only script that
*creates* an entry — every other script only adds to one that already
exists. It's also **bootstrap-only**: you run it once per source file, and
never re-run it to "refresh" an entry. If you get new or corrected data for
an existing entry, that's a job for `merge_external.py` (or, for a full
do-over, `--force`).

### What it does

1. Reads the source file with the reader selected by `--source` (see
   [`initial_files/README.md`](../initial_files/README.md) for what each
   reader expects) and resolves it to a flat table of IDs and sequences.
2. Validates the result: IDs must be unique, every row must have a usable
   sequence (only the 20 standard amino acids plus a handful of ambiguity
   codes the embedders tolerate), and query rows specifically may never be
   dropped for an unusable sequence — that aborts instead, since a reference
   protein silently disappearing would be worse than stopping.
3. Drops any other row with an unusable sequence, and reports how many
   (recorded in the manifest, not just printed).
4. Computes a fixed set of sequence-derived feature columns for every
   surviving row — length, charge/composition ratios, molecular weight, pI,
   aromaticity, instability index, GRAVY — via
   `sse_tools.compute_seq_features`.
5. Writes the datafile (`id`, `Sequence`, `query`, then the source's own
   columns, then the computed features) and a manifest recording where every
   column came from.
6. Creates the entry's subfolders (`external_data/`, `structures/`,
   `figures/`, `msa_cache/`, `logs/`) alongside it.

The whole entry is built in a temporary directory first and only moved into
`entries/<stem>/` once everything succeeded, so a failed or interrupted run
never leaves a half-written entry behind.

### Usage

```bash
python scripts/sse_initialization.py <input_file> --source {em,fasta,fs} [options]
```

`input_file` is a bare filename (looked up in `initial_files/`) or a path.
`--source` is mandatory — there's no format auto-detection, because guessing
wrong silently would be worse than requiring one flag.

| Flag | Meaning |
|---|---|
| `--source` | Required. Which reader to use: `em` (tabular/EnzymeMiner-style), `fasta`, or `fs` (Foldseek webserver JSON). |
| `--id_col`, `--seq_col` | `em` source only. Override the default `Accession`/`Sequence` column names. |
| `--query VALUE [VALUE ...]` | Mark specific rows as queries (reference sequences), overriding whatever the reader auto-detects. Matched against the reader's raw identifying field (raw `--id_col` value for `em`, full header for `fasta`, raw target for `fs`); an unmatched value aborts the run rather than silently doing nothing. |
| `--name` | Override the entry's stem (default: the source filename's stem). |
| `--force` | Delete and fully rebuild an entry that already exists at that name. |
| `--entries-dir`, `--initial-files-dir` | Override the default `entries/`/`initial_files/` locations (mainly useful for testing). |

### Examples

```bash
# EnzymeMiner-style tabular export, default column names
python scripts/sse_initialization.py my_enzymes.tsv --source em

# Foldseek search result, entry named explicitly rather than after the file
python scripts/sse_initialization.py search.json --source fs --name oleD_hits

# FASTA, marking one specific record as the query
python scripts/sse_initialization.py seqs.fasta --source fasta \
    --query 'sp|P12345|OLED_STRAN OleD'

# Tabular source with non-default column names, rebuilding an existing entry
python scripts/sse_initialization.py /path/to/data.tsv --source em \
    --id_col ProteinID --seq_col Seq --force
```

## `sse_coordinates.py`

Embeds an entry's sequences (or structures) with a model and reduces the
result to a handful of plottable coordinates, then merges those columns into
the datafile. Run it once per **coordinate system** you want available in the
viewer — you can run it repeatedly with different `--embedder`/`--reducer`
combinations against the same entry, and switch between the results live in
the visualizer.

### What it does

1. Picks the embedder (`--embedder`) and reducer (`--reducer`), and derives a
   **tag** from the embedder's settings (e.g. `esmc600m_mean` — model size +
   pooling), or uses `--label` if you gave one explicitly. The tag becomes the
   column-name prefix, so `esmc600m_mean` + PCA produces
   `esmc600m_mean_PC1`, `esmc600m_mean_PC2`, .... This is what keeps multiple
   coordinate systems on one entry distinct and independently switchable.
2. Refuses to run if that exact tag+reducer combination already exists as
   columns on the datafile, unless you pass `--rereduce`/`--reembed`/`--force`
   — coordinates are never silently overwritten.
3. Asks the embedder to `prepare()` its input from the datafile's rows.
   Structure-based embedders (`prostt5`, `saprot`) only work on **Foldseek
   entries** (they need the Cα coordinates from the original search JSON to
   derive a 3Di sequence) and abort immediately on a non-Foldseek entry,
   pointing you at `--embedder esmc` instead.
4. Embeds every row in batches, streaming results to
   `entries/<stem>/embeddings/<tag>.emb.tsv` as it goes, so an interrupted run
   resumes instead of restarting, and a later run with a different reducer
   reuses the same cached matrix instead of re-embedding.
5. Reduces the embedding matrix to `--n-components` dimensions, writes a
   diagnostic figure to `entries/<stem>/figures/<tag>_<reducer>.png` (a scree
   plot for PCA; a 2D projection scatter for UMAP/t-SNE, which have no
   explained-variance equivalent), and merges the coordinate columns into the
   datafile.
6. Writes a run log to `entries/<stem>/logs/<tag>_<reducer>.log`.

Rows the embedder can't process (e.g. a Foldseek hit with no usable Cα trace,
or a sequence longer than `--max-residues`) cause a loud abort by default,
since a silently incomplete coordinate system is worse than stopping;
`--include_empty` proceeds instead and leaves those rows with blank
coordinate cells (unplottable in that system only — they're unaffected in
every other column and coordinate system).

### Usage

```bash
python scripts/sse_coordinates.py <entry> [options]
```

| Flag | Meaning |
|---|---|
| `--embedder {esmc,prostt5,saprot}` | Embedding model. `esmc` works on any entry (uses the `Sequence` column). `prostt5`/`saprot` are Foldseek-only and structure-based. Default `esmc`. |
| `--reducer {pca,umap,tsne}` | Dimensionality reduction. Default `pca`. |
| `--pooling {mean,max,min}` | How per-residue embeddings are pooled into one vector per sequence. Default `mean`. |
| `--n-components N` | Coordinate dimensions to keep (default 10; use `2` for a plain UMAP/t-SNE landscape — PCA is usually left higher so you can pick different axis pairs in the viewer). |
| `--device {auto,cuda,mps,cpu}` | Compute device. `auto` (default) picks CUDA, then Apple Silicon MPS, then CPU; an explicit choice that isn't available aborts rather than silently falling back. |
| `--esmc-model {esmc_300m,esmc_600m}` | ESM-C variant. Default `esmc_600m`. |
| `--prostt5-checkpoint`, `--saprot-checkpoint` | Override the HuggingFace checkpoint used by those embedders. |
| `--umap-neighbors`, `--umap-min-dist`, `--umap-metric` | UMAP-specific parameters. |
| `--tsne-perplexity` | t-SNE perplexity (auto-capped below the sample count). |
| `--max-residues N` | Structure embedders only: skip targets longer than this rather than reconstructing an unreasonably large backbone. Default 1500. |
| `--foldseek-json PATH` | Structure embedders only: override the source JSON (default: resolved from the manifest, looked up in `initial_files/`). |
| `--label` | Override the auto-derived coordinate-system tag. |
| `--include_empty` | Proceed even if some rows can't be embedded, leaving them without coordinates in this system. |
| `--rereduce` | Reuse the cached embedding matrix, rerun only the reducer, and replace an existing coordinate system with the same tag+reducer. |
| `--reembed` | Discard cached embeddings for this tag, re-embed from scratch, rerun the reducer, and replace an existing coordinate system. |
| `--force` | Same as `--rereduce` (reuse cached embeddings, replace the coordinate system). |
| `--batch-size`, `--write-every` | Embedding batch size and how often partial results are flushed to disk. |
| `--entries-dir`, `--initial-files-dir` | Override the default `entries/`/`initial_files/` locations. |

### Examples

```bash
# Default: ESM-C + PCA, 10 components
python scripts/sse_coordinates.py akr

# SaProt structure embedding, proceeding past rows with no usable structure
python scripts/sse_coordinates.py oleD --embedder saprot --include_empty

# ESM-C 300M, max pooling, 20 PCs
python scripts/sse_coordinates.py akr --embedder esmc --esmc-model esmc_300m \
    --pooling max --n-components 20

# Replace an existing UMAP layout, reusing the cached ESM-C embeddings
python scripts/sse_coordinates.py akr --embedder esmc --reducer umap --rereduce
```

## `fetch_taxonomy.py`

Resolves an NCBI taxId for every row and expands it into a full lineage
(superkingdom → species), merging the result into the datafile as label
columns: `tax_status`, `taxid`, `tax_organism`, and one column per rank
(`superkingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`).
Needs network access and an email address (an NCBI E-utilities requirement,
not an SSE one — NCBI asks for it to be able to contact you if a script is
hammering their servers by mistake).

### What it does

1. Picks a **strategy** — how to get from a row to a taxId — either the one
   you name with `--strategy` or auto-detected from the datafile's columns
   (`foldseek` if it carries the `Databases`/`taxId` labels the Foldseek
   reader leaves behind, otherwise `em`, which resolves a taxId by looking up
   the id column via NCBI `efetch` and works on any datafile whose id is an
   NCBI protein accession).
2. For a `foldseek` entry, resolves per-row by what's already on the
   datafile — no network lookup needed if a `taxId` label column is already
   present; a `gmgcl_id` hit instead queries the GMGC unigene API; an
   `mgnify_esm30`-only hit has no taxonomy path at all and is flagged
   `no_taxonomy` rather than retried forever.
3. Once every row has a taxId (or a `no_taxonomy` verdict), expands the
   unique set of taxIds into full lineages via one shared NCBI `efetch`
   call per batch — lineage expansion doesn't depend on which strategy
   produced the taxId.
4. Merges the result in via `merge_columns`, the same additive write path
   every other tool uses.

### Resuming and re-running

This script is built to survive being interrupted and to make "did I already
do this?" cheap to answer, since a full run over a large entry can take a
while:

- Progress is cached to `entries/<stem>/tmp_taxonomy/` as each batch
  completes. Killing the script mid-run and re-running the same command
  picks up where it left off.
- Once a run finishes successfully, that tmp cache is deleted — but a later
  plain rerun still resumes correctly, because already-resolved rows are
  re-seeded from the taxonomy columns already sitting on the datafile,
  rather than refetched. In other words, a finished entry never gets
  refetched by accident; you have to ask for that explicitly.
- `--retry-failed` additionally re-attempts rows previously marked
  `taxid_unresolved`/`lineage_unresolved` (i.e. transient lookup failures)
  without touching rows that already succeeded (`ok`) or are genuinely
  without taxonomy (`no_taxonomy`).
- `--force` wipes the tmp cache, ignores whatever's already on the datafile,
  and refetches every row from scratch, overwriting existing taxonomy
  columns. Use this if you suspect the existing taxonomy data is stale or
  wrong, not for routine re-runs.

### Usage

```bash
python scripts/fetch_taxonomy.py <entry> --email you@example.com [options]
```

| Flag | Meaning |
|---|---|
| `--email` | Required by NCBI E-utilities. |
| `--api-key` | Optional NCBI API key (raises the rate limit from ~3/s to 10/s); also read from the `NCBI_API_KEY` environment variable. |
| `--strategy {auto,em,foldseek}` | taxId-resolution strategy. Default `auto` (detects from the datafile's columns). |
| `--batch N` | IDs/taxIds per NCBI request batch. Default 100. |
| `--gmgc-batch N` | Unigene IDs per GMGC API request batch (`foldseek` strategy only). Default 50. |
| `--retry-failed` | Also re-attempt rows previously marked as a transient failure, without refetching already-succeeded or genuinely-no-taxonomy rows. |
| `--force` | Wipe cached progress and existing taxonomy columns; refetch everything from scratch. |
| `--entries-dir` | Override the default `entries/` location. |

### Examples

```bash
# Auto-detect the strategy
python scripts/fetch_taxonomy.py my_hits --email you@dtu.dk

# Force the foldseek strategy explicitly
python scripts/fetch_taxonomy.py my_hits --email you@dtu.dk --strategy foldseek

# Re-attempt only the rows that failed transiently last time
python scripts/fetch_taxonomy.py my_hits --email you@dtu.dk --retry-failed
```

## `merge_external.py`

Merges your own external data — experimental measurements, lab notes,
anything not produced by SSE itself — into an entry's datafile as new label
(or coordinate) columns. This is the general-purpose escape hatch: whatever
data you have that doesn't fit one of the other tools' specific jobs goes in
through here.

### What it does

1. Reads `external_file` (a bare filename is looked up in the entry's own
   `external_data/` folder, the per-entry equivalent of the top-level
   `initial_files/`) and picks its id column — the first column by default,
   or `--id-col`.
2. Optionally narrows to specific columns with `--columns`; otherwise every
   non-id column is merged in.
3. Optionally remaps ids first via `--translator`: a two-column table
   (positionally read — column 1 is the **datafile** id, column 2 is the
   **external** id — the header names themselves don't matter) for when your
   external data uses a different identifier scheme than the entry (e.g. lab
   clone codes instead of accessions). Without a translator, external ids are
   matched against the datafile's id column by exact string equality.
4. Drops external rows whose (possibly translated) id has no match on the
   datafile, with a warning listing which ones — this is expected and not an
   error, since your external dataset will often cover more or different
   sequences than one particular entry. Datafile rows with no matching
   external row are kept, just with the new columns left blank; the merge
   never removes or reorders datafile rows.
5. Merges via `merge_columns`, the same additive write path every other tool
   uses — a name collision with an existing column aborts unless `--force`.

### Usage

```bash
python scripts/merge_external.py <entry> <external_file> [options]
```

| Flag | Meaning |
|---|---|
| `--id-col` | Column in `external_file` holding its id. Default: the first column. |
| `--columns` | Comma-separated list of external columns to merge (default: all columns except the id column). |
| `--translator` | Path (or bare filename, looked up in `external_data/`) to a datafile-id ↔ external-id mapping table, for when the two sources don't share an id scheme. |
| `--type {label,coordinate}` | Type token applied to every merged column. Default `label` — only use `coordinate` if the external data genuinely is plottable coordinates. |
| `--delimiter` | Override delimiter auto-detection (`.tsv` → tab, `.csv` → comma); applies to both `external_file` and `--translator`. |
| `--force` | Overwrite a column that already exists on the datafile instead of aborting. |
| `--entries-dir` | Override the default `entries/` location. |

### Examples

```bash
# Merge every column from a CSV, matched directly by id
python scripts/merge_external.py my_hits dummy_biophysical.csv

# Merge specific columns, remapping lab clone codes onto datafile ids first
python scripts/merge_external.py my_hits dummy_lab_codes.tsv \
    --id-col clone_id --columns activity_units,notes \
    --translator id_translator.tsv

# Re-merge after correcting the source file
python scripts/merge_external.py my_hits dummy_biophysical.csv --force
```

## `sse_visualizer.py`

Launches the interactive Dash app: a local web server showing one entry's
sequence space as a browsable, filterable plot, with structure prediction
(Boltz-2) and structural alignment (RMSD) built in. This is the payoff of
every other script in this folder — everything they add to a datafile
becomes something you can filter, colour, or plot by here.

This script itself has very little to document — it's a thin CLI wrapper
that resolves the entry and starts the server — because nearly everything
interesting about it is the app's UI, not its command line. That UI (every
panel, filter, and tool in the app) has its own document:
**[`docs/visualizer_guide.md`](../docs/visualizer_guide.md)**
(also available as a Word document,
[`docs/SSE_Visualizer_User_Guide.docx`](../docs/SSE_Visualizer_User_Guide.docx)
— see [`docs/README.md`](../docs/README.md)).

### Usage

```bash
python scripts/sse_visualizer.py <entry> [--port 8051]
```

| Flag | Meaning |
|---|---|
| `--port` | Port to serve on. Default `8051`. Change this if you're running the visualizer for more than one entry at once. |

Once running, open the printed `http://127.0.0.1:<port>` address in a
browser. If the entry has no coordinate columns yet, the app still opens —
it shows a message telling you to run `sse_coordinates.py` first, rather
than erroring out.

### Example

```bash
# Default port
python scripts/sse_visualizer.py my_hits

# Running two entries side by side
python scripts/sse_visualizer.py my_hits --port 8051
python scripts/sse_visualizer.py other_hits --port 8052
```
