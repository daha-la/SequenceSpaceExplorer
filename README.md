# Sequence Space Explorer (SSE)

SSE turns a table of protein sequences (a homolog search, an EnzymeMiner
selection table, a Foldseek hit list, a FASTA file, ...) into an interactive,
browsable "sequence space": every sequence is embedded, reduced to 2D
coordinates, and plotted so you can filter, color, and inspect it live in a
local web app. Structure prediction (Boltz-2) and structural alignment (RMSD)
are wired into the same viewer for on-demand structural context.

Everything revolves around one file format, the **`.sse.tsv` datafile** — a
plain TSV with a header row, a `Type` row (`id` / `label` / `coordinate` per
column), and one row per sequence. Every tool in this repo reads and/or
additively writes that one file; nothing is ever destroyed or reordered, so
you can always see exactly which tool contributed which column.

## Installation

SSE is a Python project. The instructions below assume no prior setup — if
you already have a Python/conda environment you like to manage yourself, you
just need the packages listed in [`env/requirement.yaml`](env/requirement.yaml).

### 1. Get a copy of this repository

If you have `git` installed:

```bash
git clone https://github.com/<your-org-or-user>/SequenceSpaceExplorer.git
cd SequenceSpaceExplorer
```

Otherwise, use your Git host's "Download ZIP" button and unzip it, then open a
terminal in the unzipped folder.

### 2. Install Conda

SSE depends on several scientific Python packages (PyTorch, scikit-learn,
Biopython, ...) that need to agree on compatible versions. **Conda** is a tool
that creates a self-contained, isolated Python environment with exactly the
right package versions, without touching any Python you may already have
installed. If you don't already have `conda` or `mamba` on your machine:

1. Install **Miniconda** (a small, no-frills installer for conda) — see
   [docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
   and pick the installer for your operating system (Windows/macOS/Linux).
2. Run the installer, accepting the defaults is fine for everyone except
   advanced users with a reason not to.
3. Open a **new** terminal window (the installer needs a fresh terminal to
   take effect) — on Windows, use the "Anaconda Prompt" that the installer
   added to your Start Menu; on macOS/Linux, use your regular terminal app.
4. Confirm it worked:

   ```bash
   conda --version
   ```

   This should print something like `conda 24.x.x`. If it prints
   "command not found" instead, the installer likely needs that fresh
   terminal window from step 3, or a computer restart.

### 3. Create the SSE environment

From the repository folder (the one containing this README), run:

```bash
conda env create -f env/requirement.yaml
```

This reads [`env/requirement.yaml`](env/requirement.yaml) and downloads and
installs everything SSE needs — Python itself, plus every package listed —
into a new, isolated environment named `sse`. This step downloads several
gigabytes (PyTorch and friends are large) and can take a while depending on
your connection; that's expected.

Whenever you want to run SSE, activate the environment first:

```bash
conda activate sse
```

You'll see your terminal prompt change to show `(sse)` at the start — that's
your confirmation the right environment is active. Do this once per terminal
session, before running any `python scripts/...` command from the Quickstart
below.

### 4. (Optional) GPU acceleration

`env/requirement.yaml` installs a CPU-only build of PyTorch, which works on
any machine but embeds sequences slowly (step 3 of the Quickstart below). If
your machine has an NVIDIA GPU (or you're on Apple Silicon), you can swap in
a much faster build after creating the environment:

```bash
conda activate sse
pip uninstall torch
```

Then follow the official installer picker at
[pytorch.org/get-started/locally](https://pytorch.org/get-started/locally/) —
choose your OS, "Pip", and the CUDA version matching your GPU driver (or
"Default" on Apple Silicon) — and run the command it gives you. The
coordinate-computation script auto-detects the best available device
(`--device auto`: CUDA > MPS > CPU), so no further configuration is needed
once the right PyTorch build is installed.

### 5. Verify it worked

```bash
python scripts/sse_visualizer.py --help
```

If this prints a usage message (rather than an error about a missing
package), the environment is set up correctly. See the Quickstart below for
what to run next.

### A note on structure prediction

The visualizer's Boltz-2 structure-prediction feature calls a hosted NVIDIA
API rather than running a local model, so it needs no extra install step —
just an API key, which you paste into the visualizer's UI when you use that
feature for the first time. Everything else in SSE runs fully locally.

### Keeping the environment up to date

If `env/requirement.yaml` changes later (e.g. a new dependency is added),
update your existing environment instead of recreating it:

```bash
conda env update -f env/requirement.yaml --prune
```

## Quickstart

These four steps take a raw source file to an interactive plot.

**1. Drop your source file into `initial_files/`** (or pass an absolute/relative
path directly — bare filenames are looked up there).

**2. Create the entry.** This is a one-time, bootstrap-only step per source file:

```bash
python scripts/sse_initialization.py my_hits.tsv --source em
```

`--source` selects the reader: `em` (EnzymeMiner-style tabular export), `fs`
(Foldseek JSON), or `fasta`. This creates `entries/my_hits/` containing
`my_hits.sse.tsv` plus a provenance manifest and empty subfolders
(`external_data/`, `structures/`, `figures/`, `msa_cache/`, `logs/`).

**3. Compute coordinates.** Embeds every sequence and reduces it to 2D:

```bash
python scripts/sse_coordinates.py my_hits --embedder esmc --reducer pca
```

This adds `coordinate`-typed columns (e.g. `esmc600m_mean_PC1/PC2`) to the
datafile. Re-run with a different `--embedder`/`--reducer` to add another,
independent coordinate system to the same entry — the viewer lets you switch
between them.

**4. Launch the visualizer:**

```bash
python scripts/sse_visualizer.py my_hits
```

Open the printed `http://127.0.0.1:8051` URL. From there you can filter by any
label column, build and save colored layers, click a point for its full row,
and (if configured) predict structures with Boltz-2 or run RMSD against a
query. For a full tour of every panel and option in the app, see the
[visualizer user guide](docs/visualizer_guide.md).

Optional enrichment steps — run any time after step 2, in any order, before or
after step 3:

```bash
# Pull NCBI taxonomy lineages onto every row
python scripts/fetch_taxonomy.py my_hits --email you@example.com

# Merge in your own external data (e.g. experimental measurements)
python scripts/merge_external.py my_hits my_measurements.csv
```

Every script accepts an entry **stem** (`my_hits`), an entry **directory**, or
a direct path to the `.sse.tsv` file. Run any script with `--help` for its
full flag reference.

## Repository layout

| Path | Contents |
|---|---|
| [`scripts/`](scripts/) | The CLI entry points you actually run: create an entry, compute coordinates, fetch taxonomy, merge external data, launch the visualizer. See [scripts/README.md](scripts/README.md). |
| [`sse_tools/`](sse_tools/) | The importable library backing the scripts: datafile I/O, readers, embedders, reducers, taxonomy strategies, and the visualizer's Dash backend (layers, jobs, Boltz, RMSD). See [sse_tools/README.md](sse_tools/README.md). |
| [`initial_files/`](initial_files/) | Raw source files you feed into `sse_initialization.py` (bare filenames given to scripts are looked up here). See [initial_files/README.md](initial_files/README.md). |
| [`entries/`](entries/) | Generated output: one subfolder per entry, each holding its `.sse.tsv` datafile, provenance manifest, and caches (structures, MSAs, figures, logs). See [entries/README.md](entries/README.md). |
| [`env/`](env/) | [`requirement.yaml`](env/requirement.yaml), the conda environment definition — see Installation above. |
| [`docs/`](docs/) | Longer-form documentation, currently the visualizer's user guide (Markdown and Word versions). See [docs/README.md](docs/README.md). |

## Design principles

- **One datafile per entry, additive only.** `sse_initialization.py` is the
  only tool that creates a `.sse.tsv`; every other tool merges new columns
  into it via `sse_tools.common.merge_columns`, which never reorders or drops
  rows and refuses to silently overwrite a column (`--force` required).
- **The `Type` row is the contract.** Every column is `id`, `label`, or
  `coordinate`. The visualizer derives its whole UI (which columns are
  filterable, which pairs are plottable axes) from that row — no separate
  config file to keep in sync.
- **Plug-in registries, not branching.** Source formats (readers), embedding
  models, dimensionality reducers, and taxonomy-resolution strategies are each
  a small registry of interchangeable modules. Adding support for a new one
  means adding a file and a registry entry, not editing a script.
- **Provenance is tracked.** Each entry's `logs/*.manifest.json` records which
  tool/version/params produced every column and its coverage, so a datafile
  documents its own history.

## License

MIT — see [LICENSE](LICENSE).
