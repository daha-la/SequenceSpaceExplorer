# SSE Visualizer — User Guide

This guide covers every panel, tool, and option in the Sequence Space
Explorer visualizer — the interactive plot you get from running
`sse_visualizer.py`. For how to install SSE and produce an entry to open,
see the [top-level README](../README.md); for the launch command itself,
see [`scripts/README.md`](../scripts/README.md#sse_visualizerpy).

## Contents

1. [Launching](#1-launching)
2. [Layout at a glance](#2-layout-at-a-glance)
3. [Left panel: plot controls](#3-left-panel-plot-controls)
4. [Center: the plot](#4-center-the-plot)
5. [Right panel: layers and structural analysis](#5-right-panel-layers-and-structural-analysis)
6. [Header controls](#6-header-controls)
7. [Concepts worth understanding](#7-concepts-worth-understanding)

## 1. Launching

```bash
python scripts/sse_visualizer.py <entry> [--port 8051]
```

`<entry>` is an entry stem, an entry directory, or a direct path to a
`.sse.tsv` file. Once it starts, open the printed address
(`http://127.0.0.1:8051` by default) in a browser. The app runs entirely on
your machine — nothing about the plot itself is sent anywhere — the only
features that reach the network are NCBI taxonomy fetching (a separate
script) and, if you use them, Boltz-2 structure prediction and MSA lookup
(§5).

If the entry has no coordinate columns yet (you haven't run
`sse_coordinates.py`), the plot area shows a message telling you so instead
of a blank plot — that's expected, not an error.

## 2. Layout at a glance

The app is three columns:

- **Left** — controls for what's plotted and how it looks: coordinate
  system, appearance, colour, filters, ID search, saving layers, and column
  type overrides.
- **Center** — the plot itself, a selection toolbar above it, and a details
  panel below it that fills in when you click a point.
- **Right** — saved layers, data extraction, figure export, and the two
  structural-analysis tools (Boltz-2 prediction and RMSD).

Every panel on the left and right is a collapsible section (click its
heading to expand/collapse) — collapse the ones you're not using to cut down
on scrolling.

## 3. Left panel: plot controls

### Coordinates

Chooses what's on the X and Y axes.

- **Coordinate system mode** (default) — pick one coordinate system from the
  dropdown (e.g. `esmc600m_mean_PC`), then pick which two of its columns go
  on X and Y (e.g. `PC1`/`PC2`, or `PC1`/`PC3` to look at a different plane
  of the same PCA). This is what you want almost all the time.
- **Advanced free-axis mode** — pick literally any coordinate column for X
  and any for Y, including from two *different* coordinate systems. A
  warning appears when you do this, since plotting, say, a PCA axis against
  a UMAP axis is a valid thing to do but not something with an obvious
  geometric meaning — the warning is a reminder, not an error.

### Appearance

Cosmetic controls, grouped by what they affect:

- **Filtered points** — opacity, size, and marker shape of the points
  currently passing your working filter (see [§7](#7-concepts-worth-understanding)
  for what "working filter" means).
- **Background points** — size of the light-grey points behind everything
  else, representing every sequence with coordinates regardless of filter
  state, so you always have visual context for where the filtered subset
  sits in the whole space.
- **Query markers** — opacity, size, and stacking position (on top of
  everything / below saved layers / hidden entirely) of the star/diamond/etc.
  markers used for rows flagged as `query` (reference sequences).
- **Working filter position** — whether the working filter's points draw on
  top of or underneath your saved layers.

### Colour

How the working filter's points are coloured:

- **Fixed colour** — pick one colour from the swatches for every point.
- **Continuous** — colour by a numeric column, with a colormap, an optional
  reverse toggle, and a choice of colour range: **Subset** scales the
  colormap to the min/max of only the currently-filtered points (more
  visual contrast within your subset); **Global** scales it to the min/max
  across the whole entry (lets you compare colour intensity meaningfully
  between different filters/layers).

### Filters

One control per filterable label column, grouped into four kinds — SSE
infers each column's kind automatically from its values (see
**Column settings** below for exactly how), and you can override a wrong
guess there too.

- **Continuous** — a checkbox to enable it, a range slider, and matching
  min/max number boxes (typing a number and the slider stay in sync either
  direction).
- **Boolean** — All / True / False.
- **Categorical** — a multi-select dropdown. Its option list live-updates to
  show only values that actually occur among rows passing your *other*
  enabled filters, with a live count per value — so as you narrow down with
  one filter, the next one's dropdown automatically reflects what's still
  reachable.
- **Tag-split** — for columns whose cells hold comma-separated tags (e.g. a
  `Databases` column listing several source databases per row); the
  dropdown lists individual tags rather than whole cell values, and a row
  matches if it contains any of the tags you pick.

Enabled filters combine with AND logic. A filter you haven't checked the
enable-box for has no effect, even if you've touched its slider/dropdown —
this lets you pre-stage a filter's value without it affecting the plot yet.

### Search by ID

Enable it, then paste a comma-separated list of IDs or names into the box.
It matches against the datafile's id column and any column that looks like
a name/label field (column names containing "name", "label", "alias",
"accession", or "gene"), and folds into the working filter alongside
whatever else is enabled. A status line reports how many of your pasted
values actually matched something.

### Save layer

Snapshots the *current* working filter — its filter conditions and its
colour/style settings — as a permanent, independently-toggleable **saved
layer**, listed in the right panel. Give it a name, or leave it blank for an
automatically generated one summarizing the active filters. See
[§7](#7-concepts-worth-understanding) for how layers relate to the working
filter.

### Column settings

Lets you override how a column was auto-classified (Continuous / Boolean /
Categorical / Tag split / Skip) — useful when the automatic guess gets it
wrong, e.g. a numeric-looking ID column you don't want treated as
continuous. Change a dropdown, then click **Rebuild filter panel** to apply
your overrides (they don't take effect live as you change them, so you can
adjust several before rebuilding).

**How the automatic guess works.** Every `label` column is sorted into one
of the five kinds above the first time it's loaded, by this rule:

1. **Empty or constant columns are skipped.** No values, or only one
   distinct value, isn't filterable on anything.
2. **Numeric columns are Continuous** if at least 80% of their non-empty
   values parse as a number.
3. **Boolean-looking columns are Boolean** if every non-empty value is one
   of `true`/`false`/`yes`/`no`/`1`/`0` (case-insensitive).
4. Otherwise, for text columns:
   - Date-shaped values (`YYYY-MM-DD...` in most rows) are **skipped** —
     dates aren't currently a supported filter kind.
   - A column that's mostly comma-separated *numbers* (not tags) is
     **skipped** rather than treated as Tag-split.
   - A column is **skipped** as free text if it's both high-cardinality
     (≥80% of values are unique) and long (≥30 characters average) —
     e.g. a description or notes field.
   - A column with **more than 200 distinct values** is skipped — a flat
     dropdown with hundreds of entries isn't usable — *unless* it's rescued
     (see below).
   - Otherwise: if ≥30% of values contain a comma, or ≥5% do *and* the
     individual comma-separated tokens repeat across rows more than the
     whole cell values do, it's **Tag-split** (e.g. a `Databases` column
     listing several source databases per row). Everything else is plain
     **Categorical**.

**The 200-unique-value cutoff can be rescued.** A column skipped only for
being high-cardinality (not for being numeric/boolean/date-like/free-text)
is promoted back to Categorical if it "nests" under — or has nesting under
it — another column that's already Categorical or Tag-split: at least 40%
row coverage in both columns, and at least 98% of the time each value in
the high-cardinality column maps to exactly one value in the partner
column (a functional dependency, checked in both directions, so it doesn't
matter which of the two columns is "wider"). This is how a `species` column
with hundreds of distinct values becomes filterable once it's paired with a
coarser `genus` or `family` column — narrowing the coarser filter first
shrinks the rescued column's dropdown down to a handful of live values.
It's a general nesting check, not special-cased to taxonomy columns by
name, so it applies to any two columns with that kind of relationship. A
high-cardinality column with no such partner stays skipped, and you can
still promote it manually here if you want it anyway.

## 4. Center: the plot

Above the plot:

- A point-count summary (`Showing X / Y filtered sequences · Z have
  coordinates here · N saved layer(s) visible`) and an eye icon to
  show/hide the working filter's points without disabling its filters.
- A selection bar: how many points are currently selected, a colour picker
  for the selection highlight ring, **Clear**, and **→ Working layer**,
  which copies your current point selection into the ID Search box.

The plot itself is a standard pannable/zoomable Plotly scatter plot.
**Click a point to select or deselect it** (clicking again toggles it back
off) — this is how you build up a multi-point selection for RMSD's "selected
sequences" scope, or to move into the ID search box. Clicking also fills in
the **details panel** below the plot with every column's value for that
row, and marks that sequence as the current Boltz-2 target (§5).

## 5. Right panel: layers and structural analysis

### Saved layers

Every layer you've saved, in stacking order (top of the list draws last,
i.e. on top). Each has: reorder arrows, a **⤴** button to load its filter
state back into the working filter (so you can pick up editing it), a
visibility toggle, and a delete button. **Clear all** removes every saved
layer. Layers persist to `entries/<stem>/logs/layers.json` and survive
restarting the app.

### Extract

Downloads a CSV of every unique sequence across your currently *visible*
saved layers (coordinate columns stripped, since those aren't meaningful
outside the plot), and keeps a copy in `entries/<stem>/logs/` for your
records. This is the main way to get a filtered subset of an entry out of
the visualizer and into something else (a new FASTA, a spreadsheet, ...).

### Export figure

Exports the plot exactly as currently displayed (same filters, colours,
layers) as an image:

- **Format**: PNG, SVG, or PDF.
- **Resolution**: 150/300/600 dpi (rasterized formats only; ignored for SVG,
  which is already resolution-independent).
- Toggle the legend and a transparent background.
- Override axis line, axis-label, and marker-edge colours (hex), and a
  background colour (used unless transparent is on).
- Custom pixel width/height, with a one-click reset to the default
  1200×800.
- **Save to**: a browser download, or directly into the entry's
  `figures/` folder.

Exporting requires the `kaleido` package (already in
[`env/requirement.yaml`](../env/requirement.yaml)); if it's missing or the
render times out, the status line explains what went wrong.

### Boltz-2 structure prediction

Predicts a 3D structure for **one sequence at a time** via the hosted
Boltz-2 API, and optionally with a bound small-molecule ligand. Unlike
RMSD's scope option below, Boltz-2 has no bulk mode — if you've built a
multi-point selection (e.g. for RMSD), submitting to Boltz-2 ignores it and
predicts only the single sequence you most recently clicked.

**What leaves your machine.** This feature calls two external, third-party
services over the network, so don't use it on sequences you need to keep
private: your sequence (and, if you add one, the ligand SMILES) is sent to
NVIDIA's hosted Boltz-2 API (`health.api.nvidia.com`) to run the prediction,
and — unless you turn off **Use MSA** — is separately sent to the public
ColabFold MSA server (`api.colabfold.com`) to build an alignment first.
Nothing else in SSE sends sequence data anywhere.

**Getting an API key.** You need an NVIDIA NGC personal API key (it starts
with `nvapi-`):

1. Sign in to (or create) an NVIDIA account at
   [org.ngc.nvidia.com/account/api-keys](https://org.ngc.nvidia.com/account/api-keys)
   and click **Generate Personal Key**.
2. In the **Services Included** dropdown, make sure **NGC Catalog** is
   selected (you can add other services too if you plan to reuse the same
   key elsewhere).
3. Generate the key and copy it immediately — NGC shows it to you only
   once. Personal keys can be given an expiration date and revoked/rotated
   later from the same page if you lose it or want to retire it.

See NVIDIA's own [Boltz-2 getting-started
guide](https://docs.nvidia.com/nim/bionemo/boltz2/latest/getting-started.html)
for more background on the underlying service.

**The key is not saved anywhere by SSE.** It lives only in the browser page
for the current session — it isn't written to the datafile, to any file on
disk, or persisted across a restart of the visualizer, so you'll need to
paste it in again each time you relaunch the app.

1. Paste your API key and click **Check API key** to validate it — the
   submit button stays disabled until validation succeeds (changing the
   key text invalidates it again, so you always re-check after editing it).
2. **Click a sequence point in the plot** to select the prediction target —
   the submit button also stays disabled until you've done this.
3. Leave **Use MSA** on (recommended: an alignment-informed prediction is
   more accurate) unless you have a specific reason to skip it — e.g. to
   avoid sending the sequence to the ColabFold server, or for speed.
4. Optionally paste one or more ligand SMILES strings (one per line) and a
   short label, to predict a **holo** (ligand-bound) structure instead of
   **apo**. Different ligands against the same sequence are cached and kept
   separately.
5. Optionally open **Prediction parameters** to adjust recycling steps,
   sampling steps, diffusion samples, or step scale — the Boltz-2 defaults
   are reasonable starting points for most proteins.
6. **Force re-run** ignores an existing cached prediction for this exact
   sequence/ligand combination and predicts again from scratch.
7. Click **Send to Boltz-2**. A job table below tracks status (`queued` →
   `msa` → `predicting` → `done`/`error`) and refreshes automatically every
   few seconds while anything is in progress — you don't need to keep the
   panel open or reload manually.

When a prediction finishes, its pTM/pLDDT scores are written back into the
datafile as new columns automatically, and the app reloads its filter/colour
controls so those columns are immediately usable — no manual "Reload
datafile" needed. Predicted structures land in
`entries/<stem>/structures/apo/` or `.../holo/` as `.cif` files (see
[`entries/README.md`](../entries/README.md)).

### RMSD analysis

Structurally aligns predicted structures against a reference and reports
RMSD (root-mean-square deviation), once you have at least one completed
Boltz-2 apo prediction. **RMSD only ever considers `apo` predictions** —
a `holo` (ligand-bound) structure never appears as a reference or query
option here, even if you've predicted one for the same sequence. If you
need to compare a sequence structurally, predict (or also predict) it as
apo.

- **Reference structure**: which predicted sequence to align everything
  else against, and which ranked model of it (`Reference rank` — Boltz-2
  produces multiple ranked models per prediction; rank 0 is the top-scoring
  one).
- **Scope**: compare against **all** completed apo structures, or only your
  currently **selected** sequences (built by clicking points in the plot).
- **Advanced options — per-sequence rank**: override which ranked model to
  use for specific query sequences, instead of rank 0 for all of them.
- **Alignment method**: **sequence-guided** (Kabsch superposition using a
  pairwise sequence alignment to match up residues — the right default when
  structures are homologous but not identical), **structure-based** (CE —
  purely geometric alignment, useful when sequence identity is too low for a
  reliable sequence alignment), or **both**.

Click **Calculate RMSDs**. Results (per query: rank, aligned-residue count,
RMSD in Å, method, and whether it was served from cache) appear in a table,
and — like Boltz-2 — are merged back into the datafile as new columns
automatically.

## 6. Header controls

- **Theme**: four colour themes (Clean Lab, Dark Lab, Rose Quartz, Deep
  Canopy) — purely cosmetic, your choice persists across sessions.
- **Reload datafile**: re-reads the `.sse.tsv` from disk. Use this after
  running another script (`sse_coordinates.py`, `fetch_taxonomy.py`,
  `merge_external.py`) against the same entry in a separate terminal while
  the visualizer is open, so the new columns show up without restarting the
  app. (Boltz-2 and RMSD reload automatically on completion, as noted
  above — this button is for changes made *outside* the visualizer.)

## 7. Concepts worth understanding

**Working filter vs. saved layers.** The "working filter" is whatever the
Filters/Colour/Search panels currently say — a live, editable, single
selection. A "saved layer" is a frozen snapshot of a working filter you
liked, kept around so you can compare several filtered subsets on the plot
at once (each with its own colour/style), independently of whatever you're
currently editing in the working filter. You can have any number of saved
layers visible simultaneously, but only one working filter at a time.

**Everything the visualizer writes lands in the datafile.** Boltz-2 scores
and RMSD results aren't kept separately from the rest of your data — they
become ordinary columns in the `.sse.tsv`, filterable and colourable like
any other column, and visible to every other SSE tool from then on.

**The plot only shows what has coordinates.** A row can exist in the
datafile, pass every filter, and still not appear on the plot if it has no
value in the current coordinate system's columns (e.g. it was skipped
during embedding — see [`sse_coordinates.py`](../scripts/README.md#sse_coordinatespy)'s
`--include_empty`). The point-count line above the plot reports both
numbers so this is never silently invisible.
