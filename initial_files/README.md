# initial_files/

Staging area for the raw source files you turn into entries. This is where
`sse_initialization.py` looks for a bare filename (see the Quickstart in the
[top-level README](../README.md)) — passing a full/relative path elsewhere on
disk works too, but keeping sources here means every entry's provenance
(`source_file` in its manifest) stays resolvable later.

Nothing in this folder is ever modified by SSE. `sse_initialization.py` reads
a file here once to create an entry; it's never re-read after that — **except
for Foldseek JSON**, which `sse_coordinates.py` reads again later if you use a
structure-based embedder (`prostt5`/`saprot`), since the Cα coordinates it
needs are only in the original JSON, not in the datafile. Don't delete or
move a Foldseek source file out of this folder once you've created its entry
if you plan to use those embedders.

The file's stem becomes the entry's default name (`entries/<stem>/`, override
with `sse_initialization.py ... --name`), so a distinctive filename is worth
choosing up front.

## What to put here

One file per `--source` you intend to use:

| `--source` | Format | Requirements |
|---|---|---|
| `em` | Tabular TSV (EnzymeMiner selection-table export, or any similarly-shaped tabular search result) | A header row; an ID column (default `Accession`, override with `--id_col`) and a sequence column (default `Sequence`, override with `--seq_col`). |
| `fasta` | Standard FASTA (`.fasta`/`.faa`/`.fa`/`.txt`) | One `>header` + sequence per record. The ID is taken from the first whitespace-delimited token of the header (UniProt-style `db\|accession\|name` headers are unwrapped to the accession). |
| `fs` | Foldseek webserver JSON (the raw JSON downloadable from a Foldseek search result page) | Must contain the search's `queries` and per-database `results`/`alignments` blocks, as produced by the Foldseek webserver — not the plain-text `.m8` format. |

Each source also auto-detects which row(s) get flagged as the entry's
`query` (reference) sequences differently — see ["What 'query'
means"](../scripts/README.md#what-query-means) in `scripts/README.md`.

Both example files currently in this folder are real inputs for the two
`entries/` examples in this repo, and are a good reference for the expected
shape of each format:

- `EnzymeMiner_Selection_Table_ri4plk.tsv` — `em`
- `Foldseek_2026_05_21_13_42_55.json` — `fs`

To support a source format not listed here, see the plug-in pattern for
`readers/` in [`sse_tools/README.md`](../sse_tools/README.md).
