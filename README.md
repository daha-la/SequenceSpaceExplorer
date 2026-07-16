# Sequence Space Explorer

Sequence Space Explorer builds an entry-centered, extensible protein sequence-space dataset. It can initialize data from EnzymeMiner/generic TSV, FASTA, or Foldseek JSON; generate sequence- or structure-aware embeddings; reduce and cluster the embedding space; merge taxonomy and external measurements; predict structures and binding with Boltz-2 and compare them by RMSD; and explore the result in an interactive Dash visualizer. Structure prediction and RMSD run as a pipeline module (`scripts/sse_boltz.py`) over sequences you select and export from the visualizer.

## Documentation

- [Practical user guide](docs/USER_GUIDE.md) — installation, recommended workflow, copyable commands, visualizer use, and safe reruns.
- [Detailed pipeline reference](docs/PIPELINE_REFERENCE.md) — data contracts, every script and major option, outputs, caches, decision guidance, extension points, and troubleshooting.

## Pipeline control center

[`pipeline-ui-poc/`](pipeline-ui-poc/) contains the browser-based pipeline control center. It represents all eight pipeline entry points, discovers real entries, validates and executes commands through a local runner, streams logs and status, supports cancellation and uploads, persists job history, and serializes data-writing commands for each entry.

On macOS, double-click [`Start Sequence Space Explorer.command`](pipeline-ui-poc/Start%20Sequence%20Space%20Explorer.command) to launch the complete local app. Use the **Shut down** button in the app to stop the interface, runner, and active pipeline process.

## Quick start

Run commands from the repository root:

```bash
# Create an entry from a TSV in initial_files/
python scripts/sse_initialization.py my_enzymes.tsv \
  --source em \
  --name my_enzymes

# Generate ESM-C embeddings and PCA coordinates
python scripts/sse_coordinates.py my_enzymes \
  --embedder esmc \
  --reducer pca \
  --n-components 10

# Add a two-dimensional UMAP from the same embedding cache
python scripts/sse_coordinates.py my_enzymes \
  --embedder esmc \
  --reducer umap \
  --n-components 2

# Open the interactive explorer
python scripts/sse_visualizer.py my_enzymes
```

See the [practical user guide](docs/USER_GUIDE.md) before the first run, particularly for environment setup, input-specific query handling, and the recommended order for annotations and clustering.
