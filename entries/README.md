# entries/

Generated output. Every "entry" — one source file turned into a sequence
space — gets its own subfolder here, named after the entry's stem (e.g.
`entries/my_hits/` for an entry created as `my_hits`). Nothing under
`entries/` is meant to be hand-edited; every file in it is created and
updated by the scripts in [`scripts/`](../scripts/).

## Layout of one entry folder

```
entries/<stem>/
├── <stem>.sse.tsv          the datafile — the only file that matters long-term
├── external_data/          drop files here for merge_external.py to merge in
├── embeddings/             cached embedding matrices + 3Di sequence cache
├── structures/             Boltz-2 predicted structures (apo/ and holo/)
├── msa_cache/              cached MSAs used for Boltz-2 predictions
├── figures/                diagnostic plots + manual exports from the visualizer
└── logs/                   manifest, job state, saved layers, run logs
```

### `<stem>.sse.tsv`

The datafile: a header row, a `Type` row (`id`/`label`/`coordinate` per
column), then one row per sequence. This is the one file every tool reads
and additively writes to, and the only thing you'd actually need to hand off
or back up — everything else in the entry folder is a cache, a log, or an
input staging area that exists to make regenerating that file faster.

### `external_data/`

Where `merge_external.py` looks up a bare filename for the dataset (and
optional id-translator table) you're merging into this entry — the
per-entry equivalent of the top-level `initial_files/`. Files here are
never modified.

### `embeddings/`

Raw embedding matrices, one `<tag>.emb.tsv` per embedder/pooling
combination `sse_coordinates.py` has run (e.g. `esmc600m_mean.emb.tsv`).
Kept so re-running with a different `--reducer` reuses the existing
embeddings instead of re-embedding from scratch, and so an interrupted run
can resume. `embeddings/3di_cache/` additionally caches the 3Di sequences
computed from Foldseek structures for the `prostt5`/`saprot` embedders.

### `structures/`

Boltz-2 predictions launched from the visualizer, as `.cif` files:
`apo/<sequence_id>/<sequence_id>_Rank_N.cif` for a plain structure
prediction, `holo/<sequence_id>__<ligand_label>__<hash>/..._Rank_N.cif` for
a prediction with a bound ligand (SMILES), where `<hash>` distinguishes
different ligands run against the same sequence.

### `msa_cache/`

Multiple sequence alignments fetched for Boltz-2 predictions, cached as
`<sequence_hash>.a3m` so the same sequence is never re-aligned twice.

### `figures/`

PNGs written to automatically: `<coordinate-tag>_<reducer>.png`, the
diagnostic plot each reducer produces for a coordinate system (e.g.
explained-variance for PCA). Also holds `sse_figure_<timestamp>.<ext>`
files exported manually from the visualizer's plot (PNG/SVG/PDF).

### `logs/`

Everything else:

| File | Contents |
|---|---|
| `<stem>.sse.manifest.json` | Column-level provenance: which tool/version/params produced each column and its coverage. Written at creation, appended to by every merge. |
| `jobs.json` | State of Boltz/RMSD jobs launched from the visualizer (used to resume/poll in-progress work). |
| `layers.json` | Saved filter "layers" from the visualizer's sidebar. |
| `boltz_log.csv` | One row per Boltz-2 prediction run (status, pTM, pLDDT, cif paths, ...). |
| `rmsd_log.csv` | One row per RMSD comparison run. |
| `<tag>_<reducer>.log` | Text log from a `sse_coordinates.py` run for that coordinate system. |
