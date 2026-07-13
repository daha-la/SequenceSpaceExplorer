# sse_tools/

The importable library behind every script in [`scripts/`](../scripts/). This
folder has no CLI of its own — nothing here is meant to be run directly. It
exists so the scripts stay thin (argument parsing + orchestration) while the
actual logic (datafile I/O, readers, embedders, reducers, taxonomy lookups,
the visualizer's backend) lives in one place, is unit-testable, and is
reused instead of duplicated.

`sse_tools/common.py` is the foundation: every other module, and every
script, imports from it rather than from each other. It owns the datafile
format itself — the `Type` row contract, reading/writing `.sse.tsv`,
resolving an `ENTRY` argument to a path, and `merge_columns`, the one
function every tool uses to additively write into a datafile (left-join by
id, never reorders/drops rows, refuses silent column collisions).

Comments throughout this package cite `spec §X` — that's
[`docs/SSE_datafile_spec.md`](../docs/SSE_datafile_spec.md), the
formal contract for the datafile format and the guarantees `merge_columns`
and the readers/embedders/reducers/taxonomy plug-ins rely on.

## The plug-in pattern

Four areas of the codebase are deliberately built as small plug-in
registries rather than branching logic, so adding support for something new
never means editing an existing script:

```
readers/     source format  -> a flat table           (em, fasta, fs)
embedders/   sequence/struct -> an embedding vector   (esmc, prostt5, saprot)
reducers/    embedding matrix -> 2D/nD coordinates    (pca, umap, tsne)
taxonomy/    datafile row -> an NCBI taxId            (em, foldseek)
```

Each follows the same shape: a `base.py` defining the contract (a dataclass
or a class to subclass), one file per implementation, and an explicit
`REGISTRY` dict in `__init__.py` mapping a short CLI-facing name to it.
**To add a new one, write one file implementing the contract and add one
registry line — nothing else in the codebase changes.** `readers/foldseek.py`
and `embedders/esmc.py` are the best worked examples to copy from for a new
reader/embedder respectively.

## Layout

| Path | Role |
|---|---|
| [`common.py`](common.py) | Foundation: errors, reserved column/Type-row names, datafile read/write, path resolution, `merge_columns`. Everything else depends on this; it depends on nothing else in the package. |
| [`compute_seq_features.py`](compute_seq_features.py) | Sequence-derived feature columns (length, charge, pI, GRAVY, ...) computed once at entry creation. |
| [`readers/`](readers/) | One module per source file format, turning it into a flat table. `em` (tabular/EnzymeMiner), `fasta`, `fs` (Foldseek webserver JSON). Used only by `sse_initialization.py`. |
| [`embedders/`](embedders/) | One module per embedding model, plus the shared streaming/resume/device-selection machinery (`base.py`). `esmc` (sequence-based); `prostt5`/`saprot` (structure-based, via `ca_reconstruct.py` + `structure_input.py`, which rebuild a backbone from Foldseek Cα traces and convert it to a 3Di sequence). Used only by `sse_coordinates.py`. |
| [`reducers/`](reducers/) | One module per dimensionality-reduction method: `pca`, `umap`, `tsne`. Used only by `sse_coordinates.py`. |
| [`taxonomy/`](taxonomy/) | One module per taxId-resolution strategy (`em`: NCBI efetch by accession; `foldseek`: reuses source labels already on the datafile), plus shared NCBI E-utilities/lineage-expansion code (`base.py`). Used only by `fetch_taxonomy.py`. |
| [`visualizer_state.py`](visualizer_state.py) | Turns one loaded datafile into the visualizer's state: column classification (continuous/boolean/categorical/tag-split) from the `Type` row and value shape, coordinate-system grouping, default axes. Pure logic, no Dash. |
| [`layers.py`](layers.py) | Persistence for saved filter "layers" (`logs/layers.json`) and revalidating them against the current datafile. |
| [`jobs.py`](jobs.py) | Persistence for Boltz/RMSD job records (`logs/jobs.json`) launched from the visualizer. |
| [`boltz.py`](boltz.py) | Boltz-2 structure-prediction backend: MSA generation, submitting/polling the hosted NVIDIA API, on-disk structure caching, writing prediction scalars (pTM, pLDDT) back into the datafile. |
| [`rmsd.py`](rmsd.py) | Structural-alignment backend: parses predicted `.cif` structures, Kabsch and CE-based RMSD, results log, writing RMSD columns back into the datafile. |
