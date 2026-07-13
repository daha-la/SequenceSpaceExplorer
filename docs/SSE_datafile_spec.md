# SSE Datafile & Contract Specification

**Audience:** SSE tool authors.

---

## 0. Terminology

| Term | Meaning |
|------|---------|
| **SSE** | The whole project: repository, the citable thing. |
| **SSE visualizer** | The Dash app. A thin viewer over one SSE datafile. |
| **SSE tools / scripts** | Everything that is not the visualizer (creation, embedding, feature, RMSD, merge, etc.). |
| **SSE datafile** | The canonical, per-dataset, text-tabular file that is the single source of truth for visualization. One per source dataset. |
| **entry** | The on-disk folder `entries/<stem>/` holding one datafile and everything derived from it. |
| **source file** | A raw input: an EnzymeMiner TSV, a Foldseek JSON, or a FASTA file. |

---

## 1. Design principles

1. **The datafile is the single source of truth.** All information to be shown or used lives in one per-dataset file: original source data, post-processed features, coordinates, and any new labels (RMSD, third-party predictions).
2. **The visualizer renders from the datafile.** The Dash app reads the SSE datafile + Type row as its rendering source of truth. It does not infer plot labels from logs or the manifest. Some app-side workflows (Boltz and RMSD) do write results back, but only through shared backend modules and `common.merge_columns`; their logs/jobs remain operational history, not rendering inputs.
3. **Every new capability is a tool that adds columns**, not an app feature. New labels arrive by pointing a tool at the datafile (or by ID-matched merge from a third-party tool). This minimizes hardcoding and lets users extend SSE themselves.
4. **Text-tabular format** so advanced users can inspect/use the datafile outside the visualizer. (This choice is *forced* by the Type row living as a data row — a typed binary schema like Parquet cannot carry it. Accepted cost: float precision on coordinates, and in-app `to_numeric` coercion remain.)
5. **Creation is bootstrap-only.** An entry is created once and never re-created to update; all post-creation changes go through merge.

---

## 2. Repository layout

```
SSE/
  entries/                         # one folder per dataset
    <stem>/
      <stem>.sse.tsv               # THE datafile (entry root)
      external_data/               # user drop-zone for third-party outputs (SoluProt, etc.)
      structures/                  # Boltz-2 outputs
        apo/<sequence_id>/...
        holo/<sequence_id>__<label>__<hash>/...
      figures/                     # exported figures + reducer diagnostics
      msa_cache/                   # cached ColabFold MSAs
      embeddings/                  # embedding matrices (resume) + 3di_cache/ (struct)
      logs/                        # all operational history + manifest + app state
        boltz_log.csv
        rmsd_log.csv
        <tag>_<reducer>.log        # one per sse_coordinates.py run, e.g. esmc600m_mean_pca.log
        layers.json                # persisted saved layers
        jobs.json                  # persisted Boltz/RMSD job state
        <stem>.sse.manifest.json   # provenance manifest (see §11)
      tmp_taxonomy/                # fetch_taxonomy.py's resumable batch cache (§11B.3); deleted on success
  initial_files/                   # raw EM TSVs, Foldseek JSONs, FASTAs
  scripts/                         # USER-FACING tools (run by hand) + the visualizer
    sse_initialization.py          # the creation tool (§6)
    sse_coordinates.py             # coordinate-generation tool (§11A)
    fetch_taxonomy.py              # taxonomy enrichment tool (§11B)
    merge_external.py              # external dataset merge tool (§10.8)
    sse_visualizer.py              # Dash visualizer (§14)
  sse_tools/                       # importable library; users never run these directly
    common.py                      # SSEError/abort, reserved names + Type tokens,
                                   # read/write/merge datafile, path helpers (incl.
                                   # resolve_entry_path, shared by every ENTRY-taking
                                   # script), manifest I/O, process-local write lock
    compute_seq_features.py        # the feature set (§9), imported by creation
    visualizer_state.py            # entry/datafile load context for the visualizer (§14.1);
                                   # wraps common.resolve_entry_path with its own
                                   # working-directory setup
    layers.py                      # persistent saved-layer I/O (§14.4)
    jobs.py                        # persistent job-state I/O (§14.6)
    boltz.py                       # Boltz-2 backend (§14.7)
    rmsd.py                        # RMSD backend (§14.8)
    readers/
      base.py                      # the ReaderResult contract (only)
      __init__.py                  # explicit REGISTRY of sources
      em.py  fasta.py  foldseek.py # one reader per source
    embedders/
      base.py  esmc.py  prostt5.py  saprot.py  structure_input.py  ca_reconstruct.py
    reducers/
      base.py  pca.py  umap.py  tsne.py
    taxonomy/                      # taxId-resolution strategies (§11B), mirrors readers/
      base.py                      # shared eutils/lineage-expansion primitives (strategy-agnostic)
      __init__.py                  # explicit REGISTRY + detect_strategy() auto-detection
      em.py  foldseek.py           # one strategy per source
```

- `<stem>` = the filename stem of the **initial source file used to initialize the datafile**, taken **verbatim** (`my_enzymes.tsv` -> `my_enzymes`). (Every source is a file, so a stem always exists; `--name` overrides it.)
- `external_data/` is created empty at initialization. It is the **recommended drop-zone** for raw third-party outputs (e.g. `soluprot_out.tsv`) so they stay co-located with the entry (one archivable, reproducible folder). `merge_external.py`'s `external_file` and `--translator` arguments both resolve via `common.resolve_input_path(name, external_data_dir)` — a full/relative path is used as-is, a bare filename is looked up here first. Still not a hard constraint: a full path elsewhere on disk always works too (§10.8).
- **Tool/library split:** anything a user runs by hand lives in `scripts/`; everything those scripts import lives in `sse_tools/`. Import direction is one-way (`scripts/` -> `sse_tools/`, never the reverse). This is what makes adding a source, embedder, reducer, or backend a library change rather than a user-facing script rewrite.
- **Shared foundation `common.py`:** `common.py` is the bottom of the dependency graph (imports nothing else from `sse_tools/`). It holds tool-agnostic primitives: the error type + `abort`, reserved column names + Type-token constants, strict `read_datafile`/`write_datafile`, atomic `merge_columns`, path resolution, and manifest I/O (§11). `readers/base.py` holds **only** the `ReaderResult` contract; readers import shared helpers from `common`.
- The old standalone `output/` tree is removed. Runtime artifacts are co-located under `entries/<stem>/`.

**Stem collision (handled):** two different sources with the same stem (`results.tsv`, `results.json`) both want `entries/results/`. This is a collision -> **abort** (the existing-entry check, §6.3); use `--name` to disambiguate. Verbatim names; no source-type qualifier.

---

## 3. The SSE datafile format

### 3.1 Structure

```
<header row>           # column names, tab-separated
<Type row>             # one type token per column (this is data row 0)
<data row 1>
<data row 2>
...
```

- **Type row** is the first data row. Its cell under the ID column reads literally `id`. Every column has a non-blank token.
- **Type tokens:** `id` | `label` | `coordinate`.
  - `id` — the ID column (exactly one).
  - `label` — anything used for filtering/colouring (or hidden via classification). Includes Sequence, features, source annotations, `query`, RMSD, third-party labels.
  - `coordinate` — plot axes (e.g. `PC1`…`PC10`, `UMAP1`, `UMAP2`).
- **Detection:** row 0's first cell, stripped + lowercased, `== "id"` (i.e. the ID column's Type cell). A datafile **missing the Type row → hard abort** in the visualizer and every tool. (Justified: SSE's own scripts always write it.)

### 3.2 Read/write mechanics

- **Write order:** header, then Type row, then data.
- **Read:** parse normally, then peel `df.iloc[0]` as the type map (a Series keyed by column name), drop it, `reset_index(drop=True)`. Dtypes are established *after* the peel; coordinates/continuous labels are re-coerced with `pd.to_numeric(errors="coerce")` in the visualizer.
- **Strict validation:** duplicate physical header names abort before pandas can mangle them. Every column must have a non-blank Type token. Type tokens must be exactly `id`, `label`, or `coordinate`. Exactly one column must have Type `id`, and the **first column's Type token must be `id`**. A missing or malformed Type row aborts in every tool.
- **Unique ID rows required for merges.** `merge_columns` aborts if the existing datafile ID column or the incoming merge table has duplicate IDs. This prevents left-join fan-out and preserves the never-add-rows invariant.
- **TSV-safe scalar values:** first-party writebacks that are numeric in memory (Boltz scores, RMSD values) stringify before writing. Numeric-looking strings are later detected as continuous by the visualizer.

### 3.3 Format details

| Aspect | Rule |
|--------|------|
| Delimiter | **Tab** (TSV). Forced by tag-split values containing commas; EnzymeMiner is already TSV. |
| Null | **Empty field only.** Read with `keep_default_na=False`; do **not** treat `NA`/`NaN`/`null`/`N/A` strings as null (biological tokens legitimately look like that). Numeric columns still become NaN via `to_numeric(errors="coerce")`. |
| Encoding | Write **UTF-8, no BOM**. Read with **`utf-8-sig`** (tolerates a BOM added by Excel round-trips). |
| Quoting | **Lossless.** `QUOTE_MINIMAL`, `"` as quote char, doubled `""` for an embedded quote, **no** escapechar. This is the pandas + Excel default, so `to_csv(sep="\t")` / `read_csv(sep="\t")` round-trip with no extra options. A third-party writer using a different dialect (escapechar, single quotes) is a **contract violation**. |

**Consequence of lossless quoting (note in user docs):** a field with an embedded newline becomes a *physical multi-line record*. pandas handles it; line-oriented tools (`wc -l`, `cut`, `awk`, naive `head`) and Excel's wrapped cells do not. "Inspectable outside the visualizer" therefore means *with a CSV-aware reader* for the rare rows with embedded newlines.

**Writer:** the write side is `common.write_datafile(df, types, path)`, shared by every tool that emits a datafile. It takes an explicit `types` mapping (column → token) rather than inferring, so creation passes `id`/`label` and the embedding tool passes `coordinate` for its PC columns through the same writer. It aborts if any column lacks a token.

**Reader:** `common.read_datafile(path) -> (df, types)` is the pair — it does the Type-row peel (read, split off row 0 as the type map, reset index) and enforces the strict validation above. `common.merge_columns(...)` builds on `read_datafile`/`write_datafile` to additively left-join new columns into a datafile (see §10 and §12).

### 3.4 Naming

- Datafile: `<stem>.sse.tsv`.
- Manifest: `<stem>.sse.manifest.json` (in `logs/`, §11).

---

## 4. Reserved / structural columns

Guaranteed present in every datafile, produced by creation, depended on downstream. They are **overwritable only with `--force`** — the same rule as every other column (§10.6). There is **no special "never" tier**: overwriting *anything* requires `--force`; reserved columns are not exempt from that, they simply share it. A bare (non-`--force`) write that collides with a reserved column aborts like any collision.

> **Implementation bug, found in code review, not yet fixed:** `common.merge_columns` currently aborts on a reserved-column (`id`/`Sequence`/`query`) name collision **unconditionally**, before even checking `force` — so the `--force` override this section requires is currently unreachable for reserved columns specifically (non-reserved collisions correctly honor `force=True`). This contradicts the "no special never tier" rule above. Needs a one-line fix in `merge_columns` (gate the reserved-name abort on `not force`) before any tool relies on force-replacing a reserved column.

| Column | Type | Guarantee |
|--------|------|-----------|
| `id` (the ID column) | `id` | Exactly one. Unique. The join key for all merges. |
| `Sequence` | `label` | Present **and populated on every row** (see §7). Enables embedding, Boltz, features. |
| `query` | `label` (boolean) | Always emitted, uniform across all sources, default `False`. Drives star rendering and is filterable (§8). |

`Closest query` is **not** reserved — it is an EnzymeMiner-only optional source label (§8).

> **First-party vs external authority.** "External merge" (§10) and "first-party tool" (creation, feature/query §9, embedding, RMSD) are the same write path and obey the same `--force` rule. The only difference is that first-party tools *own* the columns they emit, so re-running them to update those columns is expected (with `--force`); an external merge overwriting a first-party column is possible but unusual (also `--force`). This removes the apparent §4/§9 contradiction: the feature/query tool rewriting `query` on a re-run is allowed because it owns `query` and the user passed `--force`.

---

## 5. ID resolution per source (creation)

**Universal rule:** duplicate IDs **abort** creation with a report. Uniqueness is checked *after* any source-specific dedup/resolution below.

### 5.0 Reader architecture
Source handling is split into **readers** (one per source, in `sse_tools/readers/`) and a source-agnostic **pipeline** (in `scripts/sse_initialization.py`, §6). A reader knows everything source-specific; the pipeline knows nothing about sources. A reader is `reader(path, args) -> ReaderResult`, returning a table plus:

| field | meaning |
|-------|---------|
| `table` | resolved DataFrame (all columns read as `str`), containing the ID and Sequence columns plus any extra → `label` columns. Columns starting with `_` are internal, dropped on write. |
| `id_col`, `seq_col` | which columns hold the resolved ID and the sequence (pipeline renames them to `id` / `Sequence`). |
| `source` | manifest tag (`em`/`fs`/`fasta`). |
| `match_col`, `match_label` | the column `--query` matches against (raw pre-resolution names) and its human name for error messages. Defaults to `id_col`. |
| `auto_query` | boolean mask of source-derived queries, or `None` if the source has no native query concept (§8). |
| `notes` | source facts for the manifest header. |

**`--source` is mandatory** (no auto-detection): it selects the reader. An unknown value tells the user to fix `--source` or add a reader.

**Output ID column is named `id` uniformly** regardless of source — the raw `Accession`/`target`/FASTA-token value goes into a column literally named `id` (Type `id`). The visualizer keys off the Type token, not the name, so this is cosmetic but uniform across sources.

### 5.1 EnzymeMiner (reader `em`, also the generic tabular reader)
- ID = `--id_col` (default `Accession`); sequence = `--seq_col` (default `Sequence`). Defaulting these makes EnzymeMiner the zero-flag case while the same reader handles any other table.
- `--query` matches the raw `--id_col` value.

### 5.2 Foldseek (reader `fs`)
- ID = `target` field of each hit.
- **Dedup before the uniqueness check:** the JSON has multiple database blocks (the reference file: 9). The same `target` can be a hit in several. Reduce to **one row per unique target**, aggregating the databases it appeared in into a **`Databases`** tag-split label column. (Without dedup-then-check, the abort rule misfires on normal files.) Source labels carried through: `taxName`, `taxId`, `Description`, `Databases`.
- **Keep the `_N` suffix verbatim.** `A0A978RW52`, `A0A978RW52_2`, `A0A978RW52_3` are distinct targets. Never suffix-normalize (would collapse them and abort).
- **Query row:** see §8.2. `--query` matches the raw `target`.

### 5.3 FASTA (reader `fasta`)
Per header, in order:
1. Drop leading `>`, cut at first whitespace → take the **first token**. The post-whitespace remainder of the header is **kept as a `Description` label**. It is heterogeneous free text, so it will typically classify as high-cardinality → `skip` (visible in the details panel, not as a filter), which is the right place for it.
2. Split the token on `|`.
3. If `field[0]` is a **known db tag** → candidate ID = `field[1]`. Else → candidate ID = `field[0]`.
4. Whitespace → `_` (guard).

`--query` matches the **full header** (the header line with the leading `>` stripped) — so the user supplies what they see in the file, not the resolved token.

**Collision handling (two-pass, set-wide):**
- Pass 1 resolves every header to its candidate ID and records the *next unused pipe-field* as the fallback.
- Pass 2 finds candidate IDs appearing more than once; for every row in a colliding group, **append `_` + the next unused pipe-field** (the one not already used as the candidate: for the db-tag branch that's `field[2]`, for the plain branch that's `field[1]`).
- If no unused pipe-field exists (e.g. bare `>P12345` duplicated), the collision is **unresolvable**.
- After appends, **re-check uniqueness**; any surviving duplicates **abort** with a report.

**Known db tags:** `{sp, tr, up, ur, ur100, ur90, ur50}`. Extend in `sse_tools/readers/fasta.py` if other tags are needed. Confirm this set covers your FASTA sources.

### 5.4 Adding a source
Write a reader module in `sse_tools/readers/` (copy `foldseek.py` for a non-tabular source), import it in `readers/__init__.py`, add one line to `REGISTRY`. No user-facing script changes; the pipeline is untouched.

---

## 6. Creation contract

Creation turns one source file into a fresh entry. **Bootstrap-only.**

### 6.1 Pipeline order
1. Resolve IDs per §5.
2. **Determine query membership** (§8) — *before* any row drop, because a query row gets a stricter drop rule (step 4). Query membership is computable pre-drop for every source: Foldseek from `queries[0]` / the `seqId==100` resolution; EM from the set of IDs appearing as values in `Closest query`.
3. **Structural sequence check:** if the source has *no sequence information at all* (no EM sequence column / no Foldseek `tSeq` / FASTA headers without sequence lines) → **abort**.
4. **Per-row sequence resolution, with a split drop rule:**
   - **Non-query** row with empty/unusable sequence → **drop-with-report** (sequenceless Foldseek hits are normal sparsity).
   - **Query** row with empty/unusable sequence → **abort**, reporting the query ID. Losing the reference protein is a large error, not normal sparsity.
   - Report the dropped (non-query) count and IDs.
5. Compute features via `sse_tools/compute_seq_features.py`; finalize the `query` column (§9, §8).
6. Build header + Type row + data; write the datafile; write the initial manifest (§11).

**Dangling `Closest query` note (EM):** because non-query rows can be dropped at step 4, a surviving row's `Closest query` value may point at a dropped ID. Per pressure-test item 5 this is treated as a **report-and-abort** condition: if dropping a row removes an ID that another surviving row references as its `Closest query`, that is a missing-query-sequence error and creation aborts with a report. (A query that any row depends on must have a usable sequence.)

**Result guarantee:** `Sequence` present as a column and populated on every surviving row. Coordinates are **absent** (added later by the embedding tool, §13).

The pipeline (steps 1–6) lives in `scripts/sse_initialization.py`: `build_entry` (steps 1–4, 6-assembly) and `build_manifest` (header assembly) are creation-specific; the readers (step 1), feature computation (step 5), and the shared `write_datafile`/manifest I/O are imported from `sse_tools/`.

### 6.0 CLI surface
```
sse_initialization.py INPUT_FILE
    --source {em, fasta, fs}         # MANDATORY; selects the reader
    [--id_col NAME]                  # em/tabular: ID column (default "Accession")
    [--seq_col NAME]                 # em/tabular: sequence column (default "Sequence")
    [--query VALUE [VALUE ...]]      # mark queries; matches per-source field (§8); OVERRIDES auto-detection
    [--name STEM]                    # override the entry stem
    [--entries-dir DIR]              # override the entries/ root
    [--initial-files-dir DIR]        # override the initial_files/ lookup dir
    [--force]                        # delete + rebuild an existing entry
    [-h/--help]
```
- `INPUT_FILE` resolution: if the string is an existing path (absolute or relative), use it; else look for `initial_files/<INPUT_FILE>`; else abort naming both locations. The stem derives from the resolved file's basename, decoupled from the path typed.

### 6.2 "Usable" sequence policy
After strip + uppercase, a sequence is **usable** iff: non-empty **and** contains only the 20 standard residues plus `U`, `B`, `Z`, `X` (the codes ESM-C tokenizes). Anything with gaps `-`, stops `*`, or other characters is **unusable** → dropped.

### 6.3 Collision / atomicity
- If `entries/<stem>/` already exists → **abort**, naming the path (protects expensive derived data).
- `--force` flag (default off): delete the existing entry and rebuild, with a warning.
- Build in a temp location; move into `entries/<stem>/` **only on success** (a crash leaves no half-entry to block a retry).

### 6.4 Validation that aborts (creation)
All aborts raise `common.SSEError` (one project-wide error type) via `common.abort(msg)`; the CLI catches it, prints `ERROR: <msg>` to stderr, and exits 1. Nothing is written on abort. Cases:
- No Type row writable (internal error).
- No sequence information in source (§6.1.3).
- Duplicate IDs after resolution (§5).
- (FASTA) unresolvable ID collisions (§5.3).
- A **query** row has no usable sequence (§6.1.4).
- A surviving row's `Closest query` references a dropped ID (§6.1, dangling-query note).
- (`--query`) a value matches no row (§8.5).

---

## 7. Sequence column
- Required as a column **and** populated on every row (can't embed without it).
- Rows with no usable sequence are **dropped at creation** (§6.1.4, non-query) or **abort creation** (§6.1.4, query), never carried as blanks. (So the §6.2 policy is enforced at creation; downstream tools may assume every row has a usable sequence.)
- Visualizer consequence: drop the load-time "find a Sequence column" hunt and the "Boltz disabled if no Sequence column" branch — the column is guaranteed (§14).

---

## 8. Query identification

A boolean **`query`** column, Type `label` (so it is filterable: All/True/False). Star rendering in the visualizer is keyed off **`query == True`** (no longer off `Closest query`).

### 8.1 EnzymeMiner
- `query = True` for rows whose **ID appears as a value** in the `Closest query` column (those rows *are* designated queries).
- `Closest query` **survives** as an EM-only categorical label column (it is consumed to derive `query`, but retained because it is meaningful for EM results: which query each row is nearest to).
  - Note: with a single query it is constant → classifier → `skip` → usually invisible in the filter panel. Becomes a useful categorical only with multiple queries. "Survives" ≠ "always shown."

### 8.2 Foldseek
- The query AA sequence is available directly from `queries[0].sequence` (header in `queries[0].header`, e.g. `job_A …`). The `query` block holds only `pdb`/`qCa`.
- **Query ID resolution:** if a hit has `seqId == 100` **and** its `tSeq` equals the query sequence, the query **takes that target's ID** and the duplicate hit row is **dropped** (prevents two rows — `job_A` and the self-hit — for one physical protein). Otherwise the query keeps its header token (e.g. `job_A`).
- `query = True` on that single resolved row (single-query case).

### 8.3 FASTA
- No query concept. `query` column is **present but all-`False`** (uniform schema; avoids "does this column exist" checks downstream).

### 8.4 Uniformity
`query` is **always emitted**, every source, default `False`. (Contrast `Closest query`: present only when EM provided it.)

### 8.5 The `--query` override
The auto-detections above (8.1 EM, 8.2 Foldseek) are **retained**; FASTA has no auto-detection. `--query` is a manual override that, when supplied, **replaces** the auto-detected set entirely (a user who passes `--query` knows their data) — it is not a union. It is the only way to flag queries for FASTA.
- **Match field per source (the raw, pre-resolution name):** EM → raw `--id_col` value; FASTA → full header (leading `>` stripped); Foldseek → raw `target`. Exposed per reader as `match_col`/`match_label` (§5.0).
- **Unmatched `--query` value → abort** (loudly, naming the value), before anything is written. A value matching multiple rows is allowed and reported.

---

## 9. Sequence features (`sse_tools/compute_seq_features.py`)

Sequence-derived features computed at creation (step 5 of §6.1). The module is importable and the feature set is kept separate from creation so it can grow over time, and so a future standalone "recompute features on an existing datafile" tool can reuse it. No GPU coupling (pure CPU). New user-invented features arrive via the merge contract.

**Query identification is NOT here** — it moved into the readers (`auto_query`, §5.0/§8). This module is features only. (Earlier drafts bundled the two; the reader split separated them.)

> **Re-running / standalone use.** A future standalone feature tool that writes onto an existing datafile is a first-party write: it overwrites the feature columns it owns, but only with `--force`, like every write (§4). A `--force` rewrite is also refused if a saved layer references a column being replaced (§10.6). Not built yet — for now features run only at creation.

### 9.1 Feature set (column name in the datafile → definition)

**Always computed (count-based, robust to ambiguous residues `U`/`B`/`Z`/`X`):**

| Column | Definition |
|--------|------------|
| `length` | residue count (rendered as integer) |
| `acidic_count` | count(D) + count(E) (integer) |
| `basic_count` | count(K) + count(R) (integer) |
| `acidic_ratio` | `acidic_count` / `length` |
| `basic_ratio` | `basic_count` / `length` |
| `ED_RK_ratio` | (D+E) / (K+R); **null** when K+R = 0. (His excluded.) |
| `ED_IK_ratio` | (D+E) / (I+K); its **own** formula — I is not basic; **null** when I+K = 0. |

**Conditional (pKa/ProtParam; `null` when the sequence contains `X`/`B`/`Z`/`U`):**

| Column | Definition |
|--------|------------|
| `net_charge_pH7` | **pKa-weighted net charge at pH 7.0** (Biopython `ProteinAnalysis.charge_at_pH(7.0)`, His included). Nullable: the pKa model is undefined for ambiguous residues. NB `charge_at_pH` does *not* raise on ambiguous residues (it silently ignores them), so the conditional group is gated explicitly on `set(seq) ⊆ standard 20`, not on a ProtParam exception. |
| `MW` | molecular weight |
| `pI` | isoelectric point |
| `aromaticity` | aromaticity |
| `instability_index` | instability index |
| `GRAVY` | grand average of hydropathy |

- The count-based group is populated for every row; the conditional group is `null` on ambiguous-residue rows (a continuous column with gaps — visualizer greys valueless points). Integer counts (`length`, `acidic_count`, `basic_count`) render without a trailing `.0` (stored as nullable `Int64`).
- biopython is already a dependency.
- **Redundancy note:** EnzymeMiner ships its own `Sequence length` column (a source label); the computed `length` feature coexists with it. Kept deliberately — `length` is the uniform, guaranteed-present name across all sources; downstream tools should key on `length`.

---

## 10. Merge contract

Merge adds columns to an **existing** entry. The matched pair to creation. The full external merge CLI is still a future wrapper, but the shared core, `common.merge_columns(...)`, is implemented and used by coordinates, Boltz, and RMSD.

### 10.1 Invariants (implemented in `common.merge_columns`)
- **Additive left-join on ID.** Adds columns; **never reorders or removes rows** except for explicit column drops requested by the caller. Row count must be unchanged after merge.
- **One-to-one merge validation.** Existing datafile IDs must be unique; incoming IDs must be unique. Duplicate incoming IDs abort before pandas can fan out rows.
- **Atomic write:** `merge_columns` performs read -> validate -> optional column drop -> merge -> write temp -> `os.replace` as one critical section.
- **Process-local write serialization:** the full critical section is wrapped by a process-wide `threading.Lock` in `common.py`. This prevents lost updates between Dash callbacks / Boltz worker threads / RMSD callbacks inside one running visualizer process. Cross-process writes still require operational discipline or a future OS-level file lock.
- Entry must already exist (merge != create).
- Merge is the chokepoint that updates the **manifest** (§11).
- **Merge writes the Type token only.** Label sub-classification (continuous/boolean/categorical/tag_split/skip) stays in the visualizer at load.

### 10.2 ID matching
- **Exact string match only.** No fuzzy matching, no built-in normalization, ever.
- Practical because the datafile owns the ID space and carries ID+Sequence: the intended loop is *generate external labels from datafile-exported IDs* via the **`export-ids` helper** (not yet built), so returned IDs match by construction.

### 10.3 Optional translator table (`merge_external.py --translator`)
**Supersedes the earlier draft of this section**, which specified fixed header names `foreign_id`/`canonical_id`. What's actually built is simpler and positional: **column 0 = datafile id, column 1 = external id**, header text ignored entirely (the two sources rarely agree on what to call their id columns, so requiring specific header text would just be friction).
- **Pre-pass:** every external id is looked up in the translator; a hit is remapped onto the datafile id before matching, a miss is **warned about and dropped** (`--translator` id(s) not found), not aborted — an incomplete translator is normal, not an error.
- **Validation before any write:** a translator row with a blank id in either column aborts. A **duplicate external id** in the translator (two rows disagreeing about where one external id maps) aborts, naming the id(s).
- **Many-to-one is not specially validated here.** If the translator legitimately maps two different external ids onto the same datafile id, that surfaces downstream as an ordinary duplicate-incoming-id abort in `merge_columns` (§10.1) — there is no value-agreement leniency built (the original draft's "legitimate only if the collapsing rows agree" rule was not implemented; a collision aborts unconditionally, regardless of whether the two rows' values would have agreed).

### 10.4 Sequence verification rail (not built in `merge_external.py` v1.0)
- **Design intent, unchanged from the original draft:** on by default whenever the incoming file carries sequences; compare each matched id's incoming sequence to the datafile sequence; mismatch aborts by default, `--allow-sequence-mismatch` to override.
- **As shipped, this rail does not exist.** `merge_external.py` has no sequence-column detection or comparison at all — an external file with a sequence column that happens to match a datafile id by coincidence (not by referring to the same protein) would merge silently. Real risk when translator-mapped ids or short/generic external ids are involved. Left as follow-up work; do not rely on this protection until it is built.

### 10.5 Unmatched handling (`merge_external.py`)
- **Datafile rows with no incoming value -> empty cell** (greyed/clickable in the visualizer, like any missing continuous value). The manifest records the column's coverage (matched N of total).
- **Incoming IDs matching nothing -> warn-and-drop.** `merge_external.py` prints a warning naming (up to 20) unmatched external ids and drops those rows before merging; there is no separate structured match-report file, the warning is the report.
- **0%-match -> abort** (a column matching nothing means "wrong file/entry"; never write an all-empty column). Also enforced after translator drops: if nothing survives id-matching, `merge_external.py` aborts before calling `merge_columns`.

### 10.6 Column-name collision and replacement
- Exact-name collision -> **flag and abort**, unless `force=True` at the core or `--force` at the CLI layer.
- `force=True` = **full-column replace**: drop the existing column, write the incoming one (Type token + manifest entry rewritten; unmatched-rows-empty applies fresh). No half-and-half hybrid.
- `drop_columns=` is available in `merge_columns` for atomic full-system replacement. It is used by `sse_coordinates.py` to replace a coordinate system while removing old orphan columns (e.g. replacing 10 PCs with 2 PCs drops old PC3..PC10 in the same atomic write as the new PC1..PC2).
- **"`--force` at the CLI layer" now has two meanings in practice, worth distinguishing.** `merge_external.py`'s `--force` is exactly this section's flag: permission to overwrite. `fetch_taxonomy.py`'s `--force` means something broader ("wipe cache, ignore the datafile's existing taxonomy columns, refetch everything from scratch") — whether the *merge* is permitted to overwrite (`merge_force` internally) is computed separately and is `True` automatically whenever the datafile already carries taxonomy columns, since a plain rerun always reconstructs those columns' full prior contents before merging (§11B). Don't assume every tool's `--force` maps 1:1 onto `merge_columns(force=...)`.
- **Refuse force-replace of a column referenced by any saved layer is not built.** A saved layer can reference a column in its colour state (`cont_col`) or its filter state (`filter_state`); replacing such a column would invalidate the layer's stored colour/filter meaning. Still true for `merge_external.py` v1.0 as well as the original external-merge design: it does not read `logs/layers.json` before allowing a force-replace.
- **Reserved columns (`id`, `Sequence`, `query`) follow the same `--force` rule as every column** (§4) — overwritable only with `--force`, never silently. They are not a special "never" tier. A force-replace of `query` additionally hits the saved-layer-reference refusal whenever a layer filters on `query`. **See the §4 note: this is currently broken in `merge_columns` for reserved columns specifically.**

### 10.7 Type assignment for merged columns
- The user supplies the Type token for external labels, validated to **`label` | `coordinate`** (`id` is never valid for a merged column). Invalid token -> abort.
- SSE's own tools hardcode the correct token: embedding -> `coordinate`; Boltz/RMSD -> `label`; `fetch_taxonomy.py` -> `label` (§11B). `merge_external.py` is the one tool where the user chooses, via `--type {label, coordinate}` (default `label`), since it doesn't know in advance what an arbitrary external file contains.

### 10.8 CLI surface (`scripts/merge_external.py`)
```
merge_external.py ENTRY EXTERNAL_FILE
    [--id-col NAME]        # column in EXTERNAL_FILE to use as its id (default: first column)
    [--columns C1,C2,...]  # subset of external columns to merge (default: all except id)
    [--translator PATH]    # id-translation table, §10.3
    [--type {label, coordinate}]   # Type token for every merged column (default: label)
    [--delimiter CHAR]     # override delimiter auto-detection (both EXTERNAL_FILE and --translator)
    [--force]              # permit overwriting colliding column(s), §10.6
    [--entries-dir DIR]
```
- `ENTRY` resolves via `common.resolve_entry_path` (§2).
- `EXTERNAL_FILE`/`--translator` resolve via `common.resolve_input_path`: a full/relative path is used as-is; a bare filename is looked up in `entries/<stem>/external_data/` (§2).
- Delimiter is inferred from the file extension (`.tsv` -> tab, `.csv` -> comma); `--delimiter` overrides for files that don't follow that convention.
- Provenance is written as `provenance_source="external"`, `tool="merge_external"`, with `params` recording the source filename, resolved id column, whether a translator was used, and the Type token — so the manifest (§11.2) fully answers "where did this column come from" without needing the original command line.

---

## 11. Provenance: manifest + logs

Two distinct artifacts, both in `entries/<stem>/logs/`.

### 11.1 Logs (per-run operational history)
`boltz_log.csv`, `rmsd_log.csv`, embedding logs, and app state files (`layers.json`, `jobs.json`). Tool-specific operational records answer "what did this tool do, and when." The visualizer may use `layers.json`/`jobs.json` for app state, but it does not use logs or the manifest as a substitute for datafile columns.

### 11.2 Manifest (per-column current-state provenance)
`logs/<stem>.sse.manifest.json`. (In `logs/` so non-advanced users aren't confused by it at the entry root; its audience is tools + advanced users.)

- **Written at the creation chokepoint, extended at the merge chokepoint.** Provenance is the merge contract's responsibility, not each tool's — so a user-written or third-party tool doesn't have to maintain provenance by hand.
- **The visualizer does not depend on the manifest.** A missing/stale manifest = degraded provenance, never a load failure. (Datafile + Type row = source of truth for rendering; manifest = source of truth for provenance; logs = operational trail.)

**Per-column entry** (uniform schema; every key always present, blanks = unknown). Built by `common.make_column_entry(...)`:
- `name`, `type` (`id`/`label`/`coordinate`)
- `provenance_source`: `"sse"` | `"external"` (so auto-filled vs user-supplied is queryable, not inferred)
- `tool`, `version`, `params`, `notes`
- `coverage`: `N/total` populated-cell count, from `common.column_coverage(df)` (esp. for partially-populated merged columns). A reference to the relevant log (e.g. `logs/rmsd_log.csv`) can go in `notes` — referenced, not duplicated.

**Provenance tiers:**
- SSE's own tools (creation, feature, embedding, RMSD, taxonomy §11B) **auto-fill** provenance reliably (`provenance_source="sse"`).
- **External merges (`merge_external.py`, built differently than originally drafted):** provenance is **auto-filled, not user-supplied**. There are no `--source-tool`/`--source-version`/`--notes` flags. `tool="merge_external"`, `version` is the script's own version string, and `notes` is auto-generated as `"merged from <source filename>"`; `params` records the source filename, resolved id column, whether a translator was used, and the Type token. `provenance_source="external"` is still set correctly. If per-merge user-supplied notes/versioning turn out to matter in practice, that's an open flag addition, not something already wired up.

**Dataset header (manifest):**
- `schema_version` — integer label of which datafile-rules version this file was built under. **Written for the record, not checked on load** (item 4): this redesign assumes every datafile is created from scratch under the current ruleset, so the visualizer does **not** gate loading on it and there is no migration path. It exists purely as a provenance stamp on a thesis-era file (which ruleset produced it). One line, bumped only if the format deliberately changes.
- `source_file`, `source_type`, `creation_tool` (+ version), `created_utc`, `id_resolution` (reader notes), `row_counts` (kept / dropped / dropped_reason / queries), `dropped_ids`.
- `last_modified_utc`, `last_tool` — **auto-stamped on every manifest write.** At creation these equal the creation timestamp/tool; after a merge they show the merge tool, so the header answers "when was this datafile last touched, by what" at a glance. The rest of the header is **creation-only and immutable**; merge only appends to `columns` and re-stamps these two fields.

### 11.3 Manifest API (in `common.py`)
Shared so creation and merge use one implementation (neither reaches into the other):
- `make_column_entry(name, type, *, provenance_source="sse", tool, version, params, notes, coverage)` → one uniform column dict.
- `column_coverage(df)` → `{col: "N/total"}`.
- `read_manifest(path)` → load an existing manifest (aborts if missing; used by merge to extend).
- `write_manifest(path, manifest, tool)` → stamp `last_modified_utc`/`last_tool`, then write JSON.

Creation assembles the immutable header itself (its creation-specific facts) and emits one `make_column_entry(..., provenance_source="sse")` per column. Merge will `read_manifest` → append external-provenance entries → `write_manifest(..., tool="sse_merge")`.

---

## 11A. Coordinate generation

Coordinates are added post-creation by a separate tool, the first consumer of `TYPE_COORDINATE`, the typed `write_datafile`, and `merge_columns`. Two independent plug-in axes, each an explicit registry mirroring the readers:

- **Embedders** (`sse_tools/embedders/`, `--embedder`, default `esmc`): turn the entry's rows into an embedding matrix. `esmc` is universal (uses the `Sequence` column); `prostt5` and `saprot` are **Foldseek-only** (need a 3Di sequence) and abort unless the entry's `source_type == "fs"`. The streaming-to-disk loop with ID-based resume, pooling (mean/max/min), and device selection are shared in `embedders/base.py`; each embedder supplies only `tag`, `prepare`, `load_model`, `encode_batch`. Heavy imports (`torch`, model libraries) are lazy and are not required when an existing completed embedding cache fully covers the requested IDs.
- **Reducers** (`sse_tools/reducers/`, `--reducer`, default `pca`): matrix -> coordinates + metadata + an optional diagnostic figure. Each owns its component **label** (`pca` -> `PC`, scree figure; `umap` -> `UMAP`; `tsne` -> `TSNE`, both 2D-projection figures). UMAP-specific knobs (`--umap-neighbors`, `--umap-min-dist`, `--umap-metric`) live in the shared CLI.

### 11A.1 Coordinate naming
Columns are named **`<tag>_<LABEL><n>`**, where:
- `tag` encodes everything that defines the embedding identity: model variant **and** pooling — e.g. `esmc600m_mean`, `esmc300m_max`, `saprot_mean`. `--label` overrides it. `n_components` is **not** in the tag (it changes reduction depth, not embedding identity).
- `LABEL` is the reducer's component label, `n` the 1-based index: `esmc600m_mean_PC1`, `saprot_mean_UMAP1`.

This lets multiple coordinate systems coexist in one datafile without collision. The visualizer groups them by "prefix up to the trailing number" (`esmc600m_mean_PC`, `esmc600m_mean_UMAP`, `saprot_mean_PC` are three systems). Default mode pairs axes only inside one coordinate system; advanced free-axis mode allows cross-system comparisons with a warning (§14.2).

### 11A.2 Embedding cache and rerun semantics
The embedding matrix persists in `embeddings/<tag>.emb.tsv`; interrupted runs use `<tag>.emb.tsv.part`.

- If a completed cache exists and covers all requested IDs, it is reused directly and the model is not loaded.
- If a completed cache exists but the entry has added rows, the completed cache is copied to `.part` and only missing IDs are embedded.
- If only `.part` exists, the interrupted run resumes from it.
- `--rereduce`: reuse cached embeddings, rerun the reducer, and replace the existing coordinate system.
- `--force`: same as `--rereduce`; it is the cheap overwrite path and does **not** recompute embeddings.
- `--reembed`: delete cached embeddings for this tag, recompute embeddings, rerun the reducer, and replace the coordinate system.

Rationale: users changing `--n-components`, UMAP settings, or PCA depth should not accidentally recompute tens of thousands of expensive embeddings. Only `--reembed` means "redo upstream embedding".

### 11A.3 Coordinate-system replacement
When replacing an existing coordinate system, the tool detects all existing columns matching the same `<tag>_<LABEL><n>` prefix and passes them to `merge_columns(..., drop_columns=...)`. This removes stale orphan columns and writes the replacement in one atomic `merge_columns` operation.

Example: replacing `esmc600m_mean_PC1..PC10` with `--n-components 2 --force` leaves exactly `PC1..PC2`; old `PC3..PC10` do not survive.

### 11A.4 Unembeddable rows: abort by default
A coordinate is a *position*; a row without one is **absent from that scatter**, not greyed. If an embedder cannot process a row, the tool **aborts loudly** by default (pre-flighting structure skips before any GPU time, naming counts and reasons). `--include_empty` proceeds and leaves those rows with **empty coordinate cells** in that system (the manifest `coverage` records the partial fill). Consequence: ESM-C covers every row (creation guaranteed usable sequences), while a structure system legitimately covers only the hits with usable tCa.

### 11A.5 Persistence & I/O
- Inputs come from the datafile via `read_datafile` (`id` + `Sequence`); tCa for structure embedders comes from the **original Foldseek JSON** in `initial_files/` (located via the manifest `source_file`, `--foldseek-json` override), keyed on `id` — coordinates are never stored in the datafile.
- Outputs merge into the datafile as `coordinate` columns via `merge_columns` (provenance `sse`/`sse_coordinates`/params=embedder·model·pooling·reducer·n). The reducer figure is written to `figures/`, and a run log is written to `logs/`.

### 11A.6 CLI
```
sse_coordinates.py ENTRY            # entry stem (in entries/) or path to .sse.tsv
    [--embedder {esmc, prostt5, saprot}]   # default esmc
    [--reducer {pca, umap, tsne}]          # default pca
    [--pooling {mean, max, min}]           # default mean
    [--device {auto, cuda, mps, cpu}]      # default auto (cuda > mps > cpu)
    [--esmc-model {esmc_300m, esmc_600m}]  # esmc variant; default esmc_600m
    [--n-components N]                     # default 10 (use 2 for a UMAP landscape)
    [--umap-neighbors N] [--umap-min-dist F] [--umap-metric M]
    [--tsne-perplexity F]
    [--batch-size N] [--write-every N] [--max-residues N]
    [--foldseek-json PATH]
    [--label TAG]
    [--include_empty]
    [--rereduce]                           # reuse embeddings, rerun reducer, replace system
    [--force]                              # alias for --rereduce
    [--reembed]                            # recompute embeddings too
    [--entries-dir DIR] [--initial-files-dir DIR]
```

---

## 11B. Taxonomy enrichment

A second post-creation add-on tool, structurally parallel to coordinate generation (§11A): resolves each row's NCBI taxonomy and merges it in as `label` columns via `merge_columns`. taxId resolution is source-dependent, so it is split into an explicit **strategy registry** (`sse_tools/taxonomy/`), mirroring `readers/` (§5.0) rather than branching inside one script.

### 11B.1 Strategies
- **`em`** (universal fallback): `efetch protein` on the datafile's `id` column, indexed with and without a version suffix. `detect()` always returns `False` — it has no positive signal of its own and is only chosen when nothing else matches.
- **`foldseek`**: uses the `Databases`/`taxId` source labels already carried into the datafile by the Foldseek reader (§5.2), rather than looking anything up from scratch for most rows. Three resolution paths, dispatched per row by `Databases` tag membership:
  - **A.** an existing `taxId` value (`afdb*`, `pdb100`, `BFVD`, `bfmd`, `cath50`) is used directly.
  - **B.** `gmgcl_id` -> the row's id (stripped of Foldseek's `_trun_<n>` suffix) is resolved to a taxId via the GMGC unigene API.
  - **C.** `mgnify_esm30` -> no per-protein taxId exists; flagged via a `NO_TAXONOMY` sentinel, which becomes `tax_status = "no_taxonomy"` (a definitive fact, not a failure — never retried by `--retry-failed`, §11B.3).
  - `detect()` fires when the datafile carries both `Databases` and `taxId` columns.
  - **Judgment call, not dictated by any prior script:** because §5.2's Foldseek dedup aggregates every database a target appeared in into one comma-joined `Databases` tag_split value, a single row can carry more than one of the tags above. Resolution priority is **A > B > C** when tags overlap on one row. This was never validated against real multi-tag data; confirm it matches intent before relying on it at scale.
- **`auto`** (default, `--strategy auto`): `taxonomy.detect_strategy(df, types)` tries `foldseek` first, falls back to `em`. `--strategy {auto, em, foldseek}` overrides.
- Lineage expansion (`efetch taxonomy`, `superkingdom` .. `species`, with the NCBI `domain`->`superkingdom` rank-rename tolerated) is shared across strategies in `taxonomy/base.py`, since it no longer depends on how the taxId was obtained.

### 11B.2 Merged columns
`tax_status`, `taxid`, `tax_organism`, `superkingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species` — all Type `label` (§10.7). `tax_status` values: `ok`, `no_taxonomy` (terminal, definitive), `taxid_unresolved`, `lineage_unresolved` (both retryable, §11B.3).

### 11B.3 Cache, resume, and rerun contract
- Progress is cached to `entries/<stem>/tmp_taxonomy/taxonomy_cache.tsv`, saved after every batch, so an interrupted run resumes from where it left off. The cache directory is deleted only after a **successful** merge into the datafile.
- **A plain rerun after a successful prior run does not refetch already-resolved rows**, even though the tmp cache is gone by then: already-resolved rows are re-seeded directly from the datafile's own existing taxonomy columns (if present), not only from the tmp cache. A rerun with nothing left to do is a no-op merge.
- **`--force`** wipes the tmp cache, skips seeding from the datafile entirely, and refetches every row from scratch.
- **`merge_force` (whether the final merge may overwrite existing taxonomy columns) is decoupled from `--force`'s "wipe everything" meaning** (§10.6): it is automatically `True` whenever the datafile already carries taxonomy columns, because the reconstructed column always carries every previously-resolved value forward untouched plus whatever's newly resolved — never a destructive overwrite, so it doesn't need the user to separately opt in with `--force`.
- **`--retry-failed`** re-attempts only rows with a retryable status (`taxid_unresolved`/`lineage_unresolved`), leaving `ok`/`no_taxonomy` rows untouched and without wiping the cache. Rows that already have a taxId but failed only at lineage expansion skip the taxId-resolution step entirely on retry (no redundant `efetch protein`/GMGC calls).
- Every NCBI/GMGC batch call prints progress (`  em batch N/M (...)`, etc.) to stderr and respects NCBI rate limits (0.4s / 0.12s-with-API-key delay between batches) — a long-running fetch with many batches is expected to look slow, not to look silent.

### 11B.4 CLI
```
fetch_taxonomy.py ENTRY
    --email EMAIL                          # required by NCBI
    [--api-key KEY]                        # optional, raises NCBI rate limit to 10/s
    [--strategy {auto, em, foldseek}]       # default auto
    [--batch N]                            # ids/taxIds per NCBI request batch, default 100
    [--gmgc-batch N]                       # unigene ids per GMGC request, default 50 (foldseek strategy)
    [--force]                              # wipe cache, ignore existing columns, refetch everything
    [--retry-failed]                       # re-attempt only taxid_unresolved/lineage_unresolved rows
    [--entries-dir DIR]
```
`ENTRY` resolves via `common.resolve_entry_path` (§2).

---

## 12. Writer model

The original design called for a single designated writer thread. The implemented visualizer uses the same safety principle at the shared write chokepoint: **every datafile write goes through `common.merge_columns`, and `merge_columns` serializes its full critical section with a process-local lock.**

- **Protected within one running visualizer process:** concurrent Dash callbacks, Boltz worker threads, and RMSD callbacks cannot interleave stale read -> stale write updates, because they all acquire the same `common._DATAFILE_WRITE_LOCK` through `merge_columns`.
- **Atomic flush:** every datafile update writes a temp file and uses `os.replace`; a crash cannot leave a half-written datafile.
- **No direct app-side TSV writes:** Boltz and RMSD backend modules write scalar columns through `merge_columns`. CIF structure files, MSA caches, logs, and job/layer JSON files are separate operational artifacts.
- **Cross-process caveat:** a CLI tool and the visualizer writing the same entry at exactly the same time are not protected by the process-local lock. Avoid external writes while the visualizer is running, or add an OS-level file lock in a future hardening pass.

---

## 13. Reload contract

"Reload" re-reads the datafile and refreshes datafile-derived state **without resetting the session**.

- **Refreshed (server-side, datafile-derived):** dataframe, Type map, column metadata, query IDs, name columns, coordinate-column list, coordinate systems, axis ranges, continuous/categorical/boolean filter metadata.
- **Preserved (browser stores):** saved layers, selection, active filters, colour mode, appearance settings, current axis choice, free-axis/system mode, and app workflow state unless the relevant columns/IDs vanished.
- **Invariant:** reload may **add columns** but should not reorder rows. `merge_columns` preserves row order and row count; saved layers are stored by ID and therefore survive reloads even if a future tool changes row order.
- **Layer validation:** reload validates saved layer IDs against the newly loaded datafile. Vanished IDs are removed from layers and `layers.json` is rewritten. **Implemented centrally inside `reload_state()` itself**, not duplicated per call site — this was previously implemented separately in `build_app` (startup) and the manual-reload callback, while the Boltz-poll and RMSD-completion reload paths called `reload_state()` bare and skipped it entirely. Since `merge_columns` never removes rows (§10.1), those two paths couldn't actually have produced a stale layer ID before this fix — but centralizing means any future reload trigger gets validation for free rather than needing to remember to wire it in.
- **Two-phase dynamic-panel rebuild:** Dash pattern-matched filter/column components are cleared first, then rebuilt in a second callback. This avoids the Dash frontend `Cannot read properties of undefined (reading 'props')` error that occurs when replacing pattern-matched component trees in one update.
- **State locking:** server-side visualizer globals are swapped under a lock after a full new state is built, so render callbacks do not see half-updated combinations of dataframe/axes/metadata.
- **Bug fixed:** `_name_cols` is part of the reload/global state; name-column search works after load/reload.

---

## 14. Visualizer setup

The visualizer rewrite is implemented as a Dash app that opens one entry and treats the `.sse.tsv` as its rendering source of truth. It no longer performs a two-file annotation/coordinate upload or a load-time annotation-coordinate merge. It also no longer merges Boltz/RMSD logs at load time; Boltz/RMSD write scalar results into the datafile through `merge_columns`.

### 14.1 Module split

| File | Role |
|------|------|
| `scripts/sse_visualizer.py` | Dash layout + callbacks. Thin UI layer over the datafile and backend modules. |
| `sse_tools/visualizer_state.py` | Resolves entry paths and loads `VisualizerState`: dataframe, Type map, ID column, label columns, coordinate columns, coordinate systems, query IDs, name columns, column classifications, ranges. |
| `sse_tools/layers.py` | Reads/writes/validates `logs/layers.json`. Saved layers are ID-based and durable across app restarts. |
| `sse_tools/jobs.py` | Reads/writes `logs/jobs.json`; stale non-terminal jobs become `interrupted` at app startup. |
| `sse_tools/boltz.py` | Boltz-2 backend: API validation, apo/holo prediction, MSA cache, CIF writing, log/job updates, datafile score writeback. |
| `sse_tools/rmsd.py` | RMSD backend: structure discovery, sequence-guided and CE RMSD, rank handling, cache/log use, datafile writeback. |

Import direction remains `scripts/ -> sse_tools/`. Backends do not reach into the Dash callback layer.

### 14.2 Coordinate UI and plotting

- **Two axis modes:**
  - Coordinate-system mode (default): choose one coordinate system and then choose X/Y axes within it.
  - Advanced free-axis mode: choose any coordinate column for X and any coordinate column for Y. Cross-system axis pairings are allowed for expert comparisons and should be visually/warning-labeled as such.
- **NaN-awareness is mandatory and implemented:** every trace is gated by coordinate coverage. Background, query markers, working layer, saved layers, and selection overlays only render rows covered by the active axes. Axis ranges ignore NaNs. The point-count text reports how many filtered/layer points have coordinates in the active view.
- **No-coordinate empty state:** a valid datafile with no coordinate columns renders an informative message rather than crashing.
- **Details panel:** hides **all** `coordinate`-typed columns, not only the active X/Y columns.
- **Query markers:** rendered from `query == True`; query marker appearance is controlled independently from filtered points.
- **Point selection (`selection-store`):** clicking a plotted point toggles its id in/out of a client-side selection list (a separate mechanism from filters/layers — a click doesn't touch the working filter). Selected points render with a distinct highlight ring in a user-chosen colour. This selection list is the mechanism two other features build on: it's the query scope for RMSD's "selected sequences" (§14.8), and the most-recently-clicked id is separately tracked (`boltz-clicked-id-store`) as the single target for a Boltz-2 submission (§14.7) — clicking never adds to a Boltz-2 "batch," there is no batch mode for Boltz-2. A **→ Working layer** action copies the current selection's ids into the ID-search filter, folding a selection into the working filter by explicit action rather than automatically.

### 14.3 Filters, colour, and working layer

- Column classification happens at load/reload from Type `label` columns: continuous, boolean, categorical, tag-split, or skip/high-cardinality. Type `coordinate` columns are not shown as metadata/filter labels.
- **High-cardinality categorical narrowing (`visualizer_state.py` + `sse_visualizer.py`).** Two mechanisms, previously undocumented here, work together on any label column — not taxonomy-specific, despite being what makes taxonomy ranks usable:
  - **Rescue pass at classification (`_rescue_high_card_categoricals`, in `build_col_meta`).** A label column above `MAX_CAT_UNIQUE` (200) is normally classified `skip` — a flat dropdown with hundreds of entries is unusable. It is promoted back to `categorical` if it **nests** with an already-categorical column: another categorical where >=98% of the child's values map to exactly one parent value (order-agnostic; a column nests either as parent or child). The pass iterates, so a chain like `family -> genus -> species` fully recovers even though `genus`/`species` individually have thousands of unique values. Column names are never inspected — a column is rescued purely because it happens to nest cleanly, whatever it's called.
  - **Live option narrowing (`narrow_cat_options` in `sse_visualizer.py`).** Every categorical dropdown is rebuilt on every filter change to show only values that still yield rows given all *other* active filters (leave-one-out — a column's own selection doesn't narrow itself), with inline counts (`"Bacilli (412)"`). A value the user already selected that goes dead is kept, labeled `— 0`, so it's visible/removable rather than silently disappearing.
  - Together: a rescued column like `species` is unusable as a flat list of thousands of entries, but becomes a short, live-narrowed list the moment a coarser related column (`family`, `genus`) is filtered first.
- Continuous filters have two-way-linked sliders and min/max number inputs. Editing number boxes moves the slider and affects the working layer immediately. Decimal commas are tolerated (`0,059` -> `0.059`).
- Missing continuous values are rendered as grey/no-value traces when relevant; they remain clickable/selectable.
- Search by ID/name is preserved and can append selection IDs into the working layer.
- The working layer is independent from saved layers. Reordering, toggling, saving, deleting, or restarting must **not** load a saved layer's filters into the working layer. Loading a saved layer into the working layer only happens after an actual click on the load button, guarded against pattern-matched component creation events.

### 14.4 Saved layers

- Saved layers persist in `entries/<stem>/logs/layers.json` and are loaded at startup.
- Layers are ID-based, not positional. This is more robust than the old base-index-only approach while still preserving display state across reloads.
- Each layer stores name, visibility, colour mode, fixed/continuous colour settings, alpha, point size, selected IDs, and `filter_state` for restoring the working layer by explicit user action.
- Save-layer write path treats `layers.json` as canonical and merges with browser store as a fallback. This prevents stale browser state from overwriting previous layers.
- Reorder/toggle/delete/load callbacks check for a real button click (`n_clicks > 0` on the triggered button) before mutating state. This prevents dynamic sidebar re-renders from accidentally deleting layers or loading filters.
- Clear-all writes an empty layer list.

### 14.5 Extraction and figure export

- Visible saved layers can be extracted to CSV: unique ids across all visible layers, coordinate columns stripped. Offered as a browser download (`dcc.send_data_frame`) **and** simultaneously saved to `entries/<stem>/logs/extraction_<run_id>.csv`, so every extraction has a durable record even if the download itself isn't kept.
- Figure export is integrated into the right panel and uses Kaleido (`pio.to_image`, run off the main thread with a timeout) so a slow/failed render can't hang the app.
- Export controls actually present: format (PNG/SVG/PDF), resolution (150/300/600 dpi — ignored for SVG), legend on/off, transparent-background toggle, axis line colour, axis-label/tick colour, marker-edge colour, and background colour (used unless transparent is on) — all as hex inputs. There is **no** separate gridline control; gridlines are hardcoded off (`showgrid=False`) both on the live plot and in every export, not user-configurable. Width × height in pixels, with a one-click reset to the 1200×800 default. Destination: browser download, or saved directly to `entries/<stem>/figures/`.
- Export operates on a copy of the current Plotly figure (legend/colour/axis overrides applied to that copy); it does not mutate the live app figure or its stored state.

### 14.6 Job state

- Job state persists in `logs/jobs.json`.
- On app startup, non-terminal jobs from a previous interrupted session are marked `interrupted` rather than left as `queued`/`running` forever.
- Polling callbacks update job tables and trigger datafile reloads after completed writebacks.

### 14.7 Boltz-2 integration

- UI supports API-key validation, MSA on/off, apo/holo predictions, SMILES input, ligand labels, force rerun/cache behavior, and advanced prediction parameters (`recycling_steps`, `sampling_steps`, `diffusion_samples`, `step_scale`).
- Prediction parameter layout is a compact two-column label/input layout. Numeric parsing tolerates comma decimals for `step_scale`.
- Outputs:
  - CIF structures are written under `structures/apo/<sequence_id>/...` or `structures/holo/<sequence_id>__<label>__<hash>/...`.
  - Operational details are written to `logs/boltz_log.csv` and `logs/jobs.json`.
  - The SSE datafile receives **only score columns**:
    - `boltz_apo_ptm`
    - `boltz_apo_plddt`
    - `boltz_holo_<holotag>_ptm`
    - `boltz_holo_<holotag>_plddt`
  - The SSE datafile does **not** keep status/MSA columns (`boltz_apo_status`, `boltz_apo_msa_used`, `boltz_holo_<holotag>_status`, `boltz_holo_<holotag>_msa_used`). Status and MSA usage belong in logs/jobs, not in the main datafile. If any of those obsolete columns exist on the datafile (written by an older app version), `write_boltz_scalars` detects and drops them in the same atomic `merge_columns` write that adds the current score columns — self-healing, not a one-time migration a user has to run.
- Cached jobs write the same score-only columns as fresh jobs. Numeric values are written as TSV-safe numeric strings so they classify as continuous after reload.
- When score columns are written, the app reloads and rebuilds filter/colour options so pTM/pLDDT become available for filtering and colouring.
- **Third-party network calls, for the record:** prediction goes to NVIDIA's hosted Boltz-2 endpoint (`health.api.nvidia.com`); MSA generation (when **Use MSA** is on) goes to the public ColabFold server (`api.colabfold.com`). Both are outside SSE's control and outside this repo — relevant if an entry's sequences are not meant to leave the local machine.

### 14.8 RMSD integration

- RMSD analysis discovers completed apo structures and supports sequence-guided (`seq`), structure-based CE (`ce`), or both.
- Reference structure and reference rank are user-selectable.
- Per-sequence query-rank overrides are stored as real Dash state with pattern IDs and preserved across structure refreshes. RMSD calculation reads the live input values at click time, so typing a rank and immediately calculating works.
- Results are cached in `logs/rmsd_log.csv` and written to the datafile as label columns:
  - `RMSD_vs_<reference_id>_r<reference_rank>_seq`
  - `RMSD_vs_<reference_id>_r<reference_rank>_ce`
- RMSD values are stringified before writeback and classify as continuous after reload.
- Completed RMSD writeback triggers reload/rebuild so new RMSD columns are available for filtering/colouring.

### 14.9 Theming

Superseded: the theme/CSS overhaul this section used to describe as deferred has since been built. Four built-in themes ship in `scripts/sse_visualizer.py` — Clean Lab and Dark Lab (light/dark, low-decoration), Rose Quartz and Deep Canopy (light/dark, "Happy Hues"-derived palettes) — selectable from a dropdown in the app header.

- **Mechanism:** every theme is a block of CSS custom properties (`--page`, `--surface`, `--surface-2`, `--border`, `--border-soft`, `--text`, `--text-muted`, `--text-faint`, `--accent`, `--accent-weak`, `--on-accent`, `--success`, `--warning`, `--danger`, `--boltz`, `--shadow-card`, plus `color-scheme`), scoped under `html[data-theme="<name>"]`. All of it lives inline in one `<style>` block inside `SSE_INDEX_STRING` (the app's Dash `index_string`) — **not** a Dash `assets/` folder; no `assets/` directory exists in this repo. Component-level rules on top of the variables re-theme Dash's own bundled CSS for controls that don't otherwise pick up custom properties (`dcc.Input`/`Dropdown` containers, `rc-slider`/newer Radix slider internals, checkbox/radio accent colour), each pinned with `!important` to beat the bundled styles.
- **Switching:** a `dcc.Dropdown` (`theme-select`) drives a **clientside** callback (`document.documentElement.setAttribute('data-theme', v)`) — the swap is instant, no server round-trip and no page reload.
- **Persistence:** the dropdown uses Dash's built-in `persistence=True` (browser-local, keyed by component id), so the last-picked theme is remembered across restarts of the app in the same browser. This is separate from, and unrelated to, any datafile or `jobs.json`/`layers.json` state — a theme choice is purely a client-side UI preference.
- Purely cosmetic: no theme affects plot data, filter behaviour, or anything written to the datafile.

---

## 15. Deferred / open items

| Item | Status | Notes |
|------|--------|-------|
| **Coordinate grouping** | BUILT (§11A, §14.2) | Default system mode groups by prefix-up-to-trailing-number; advanced free-axis mode allows cross-system pairings. |
| **Embedding tool** | BUILT (§11A) | `sse_coordinates.py` + embedders/reducers. Cache reuse, `--rereduce`, cheap `--force`, and `--reembed` semantics are implemented. |
| **Visualizer rewrite** | BUILT (§14) | Single-datafile Dash app, persisted layers/jobs, Boltz/RMSD backends, NaN-aware plotting, manual reload, advanced export, point selection (§14.2), theming (§14.9). |
| **Taxonomy enrichment tool** | BUILT (§11B) | `fetch_taxonomy.py` + `sse_tools/taxonomy/` (em/foldseek strategies, auto-detect). Resumable cache, `--retry-failed`. Foldseek's A>B>C multi-tag priority (§11B.1) is an unvalidated judgment call. |
| **External merge CLI** | BUILT v1 (§10.3-§10.8) | `merge_external.py`: translator table (positional, not the original fixed-header design), `--columns` subset, `--type`, warn-and-drop, `external_data/` auto-resolution. **NOT built:** sequence-verification rail (§10.4), saved-layer force-replace refusal (§10.6). |
| **Shared entry-path resolution** | BUILT (§2) | `common.resolve_entry_path` unifies three previously near-duplicate implementations (`fetch_taxonomy.py`, `merge_external.py`, and the visualizer's richer `EntryContext` wrapper); also fixed a bug where the CLI tools couldn't accept an entry *directory* path, only a bare stem or a literal `.sse.tsv` path. |
| **Reserved-column force-overwrite bug** | OPEN | `merge_columns` aborts on an `id`/`Sequence`/`query` name collision unconditionally, before checking `force` — so `--force` is currently unreachable for reserved columns specifically, contradicting §4's "no special never tier" rule. One-line fix identified, not yet applied. |
| **`export-ids` helper** | DECIDED, NOT BUILT | Pull ID+Sequence from the datafile to feed third-party tools so returned IDs match exactly. |
| **OS-level cross-process file lock** | DEFERRED | Current write lock is process-local. Needed only if users commonly run CLI writes while the visualizer is active. |
| **Theme system** | BUILT (§14.9) | Four built-in themes, header dropdown, instant clientside switch, browser-persisted choice. Superseded the earlier "deferred, rolled back" entry this row used to carry. |
| **FASTA known-db-tag set** | IMPLEMENTED (confirm) | `{sp, tr, up, ur, ur100, ur90, ur50}` in `fasta.py`; confirm it covers your sources. |
| **Old-datafile compatibility** | DECIDED (none) | Everything is created from scratch under the current ruleset. `schema_version` is written but not checked on load (§11). |
