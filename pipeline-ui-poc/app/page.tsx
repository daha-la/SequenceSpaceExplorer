"use client";

import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import {
  cancelPipelineJob,
  getEntryFiles,
  getJob,
  getJobs,
  getRunnerState,
  requestPipelineShutdown,
  revealEntryPath,
  submitJob,
  uploadRunnerFile,
  type EntryFile,
  type EntrySummary,
  type PipelineJob,
} from "./runner-client";

type ToolId =
  | "initialize"
  | "taxonomy"
  | "external"
  | "coordinates"
  | "distance"
  | "cluster"
  | "boltz"
  | "visualizer";

type FieldState = Record<string, string | boolean | number>;

// `guide` is the anchor id of this step's section in the visual pipeline guide
// (served at /pipeline-guide.html; source: docs/PIPELINE_WALKTHROUGH.html).
const GUIDE_URL = "/pipeline-guide.html";

const tools: Array<{
  id: ToolId;
  label: string;
  script: string;
  stage: string;
  description: string;
  mark: string;
  guide: string;
}> = [
  {
    id: "initialize",
    label: "Create entry",
    script: "sse_initialization.py",
    stage: "1 · Entry",
    description: "Start from TSV, FASTA, or Foldseek JSON",
    mark: "01",
    guide: "stage-initialization",
  },
  {
    id: "taxonomy",
    label: "Fetch taxonomy",
    script: "fetch_taxonomy.py",
    stage: "2 · Enrich",
    description: "Resolve NCBI taxIds and lineages",
    mark: "02",
    guide: "stage-taxonomy",
  },
  {
    id: "external",
    label: "Merge external data",
    script: "merge_external.py",
    stage: "2 · Enrich",
    description: "Join assay, annotation, or coordinate columns",
    mark: "03",
    guide: "stage-external",
  },
  {
    id: "coordinates",
    label: "Build coordinates",
    script: "sse_coordinates.py",
    stage: "3 · Embed",
    description: "Embed sequences and reduce their dimensions",
    mark: "04",
    guide: "stage-embed",
  },
  {
    id: "distance",
    label: "Query distances",
    script: "sse_esmc_distance.py",
    stage: "4 · Analyze",
    description: "Measure embedding distance to references",
    mark: "05",
    guide: "stage-distance",
  },
  {
    id: "cluster",
    label: "Cluster space",
    script: "sse_cluster.py",
    stage: "4 · Analyze",
    description: "Run k-means or HDBSCAN with diagnostics",
    mark: "06",
    guide: "stage-cluster",
  },
  {
    id: "boltz",
    label: "Structure & binding",
    script: "sse_boltz.py",
    stage: "4 · Analyze",
    description: "Boltz-2 structure prediction + RMSD for a selection",
    mark: "07",
    guide: "stage-structure",
  },
  {
    id: "visualizer",
    label: "Open explorer",
    script: "sse_visualizer.py",
    stage: "5 · Explore",
    description: "Launch the interactive Dash visualization",
    mark: "08",
    guide: "stage-explore",
  },
];

const initialValues: Record<ToolId, FieldState> = {
  initialize: {
    input: "my_enzymes.tsv",
    source: "em",
    name: "my_enzymes",
    idCol: "Accession",
    seqCol: "Sequence",
    query: "",
    force: false,
  },
  taxonomy: {
    email: "researcher@example.org",
    strategy: "auto",
    apiKey: "",
    batch: 100,
    gmgcBatch: 50,
    rerun: "resume",
  },
  external: {
    file: "measurements.csv",
    idCol: "",
    columns: "pI,Melting_temperature",
    translator: "",
    type: "label",
    delimiter: "auto",
    force: false,
  },
  coordinates: {
    embedder: "esmc",
    reducer: "umap",
    pooling: "mean",
    normalize: true,
    device: "auto",
    esmcModel: "esmc_600m",
    prostCheckpoint: "Rostlab/ProstT5",
    saprotCheckpoint: "westlake-repl/SaProt_650M_AF2",
    components: 2,
    neighbors: 15,
    minDist: 0.1,
    metric: "euclidean",
    perplexity: 30,
    tsnePca: 50,
    batchSize: 32,
    writeEvery: 1000,
    maxResidues: 1500,
    foldseekJson: "",
    label: "",
    includeEmpty: false,
    rerun: "new",
  },
  distance: {
    embedding: "esmc600m_mean",
    queryMode: "marked",
    queryIds: "OleD_S1, AgepGT_S3",
    raw: false,
    force: false,
  },
  cluster: {
    embedding: "esmc600m_mean",
    raw: false,
    space: "pca",
    pcaMode: "dims",
    pcaDims: 50,
    pcaVariance: 0.95,
    clusterer: "kmeans",
    kMode: "auto",
    k: 12,
    kMin: 2,
    kMax: 20,
    minClusterSize: 50,
    minSamples: "",
    label: "",
    force: false,
    analysis: true,
    topN: 5,
    fdr: 0.05,
  },
  boltz: {
    selection: "",
    apiKey: "",
    smiles: "",
    smilesLabel: "",
    useMsa: true,
    recyclingSteps: 3,
    samplingSteps: 200,
    diffusionSamples: 5,
    stepScale: 1.638,
    force: false,
    rmsd: false,
    rmsdReference: "",
    rmsdRefRank: 0,
    rmsdMethod: "seq",
    rmsdScope: "all",
  },
  visualizer: {
    port: 8051,
  },
};

const landingClusterShapes = [
  { color: "#c7d2d0", lobes: [[49, 25, 8, 5, 14], [43, 39, 5, 8, 13], [56, 42, 7, 4, 12], [48, 53, 4, 8, 10], [35, 45, 7, 4, 10], [60, 26, 5, 3, 8], [39, 17, 3, 3, 5], [34, 61, 3, 4, 6], [62, 53, 3, 6, 6]] },
  { color: "#96b5e7", lobes: [[58, 71, 14, 12, 38], [70, 67, 12, 9, 28], [79, 73, 11, 7, 23], [63, 84, 8, 5, 15], [52, 64, 6, 7, 13], [84, 57, 8, 5, 15]] },
  { color: "#ff169b", lobes: [[76, 25, 10, 6, 24], [84, 32, 8, 6, 20], [88, 44, 4, 8, 15], [73, 40, 5, 4, 12], [91, 23, 3, 4, 8], [69, 18, 5, 4, 10]] },
  { color: "#20dfca", lobes: [[22, 49, 6, 7, 20], [17, 42, 3, 3, 8], [26, 57, 4, 3, 10]] },
  { color: "#ae36c2", lobes: [[35, 76, 5, 4, 17], [42, 85, 4, 4, 11], [38, 91, 3, 3, 7]] },
  { color: "#7560ff", lobes: [[31, 40, 4, 5, 18], [28, 34, 3, 3, 9]] },
  { color: "#a9bd4a", lobes: [[67, 47, 5, 8, 20], [72, 53, 4, 5, 12]] },
  { color: "#3abfff", lobes: [[52, 56, 4, 3, 14]] },
  { color: "#ff4c44", lobes: [[39, 59, 3, 3, 10]] },
  { color: "#ff875f", lobes: [[25, 68, 3, 2, 8]] },
  { color: "#e75bff", lobes: [[35, 28, 3, 3, 8]] },
  { color: "#f3a6b7", lobes: [[51, 20, 4, 3, 10]] },
] as const;

const landingPoints = landingClusterShapes.flatMap((shape, cluster) =>
  shape.lobes.flatMap(([cx, cy, rx, ry, count], lobe) => {
    const renderedCount = count * 7;
    return Array.from({ length: renderedCount }, (_, index) => {
      const angle = index * 2.399963 + lobe * .81 + cluster * .37;
      const radius = Math.sqrt((index + .45) / renderedCount);
      const wobble = Math.sin((index + 1) * 4.73 + cluster) * .55;
      return {
        cluster,
        color: shape.color,
        x: cx + Math.cos(angle) * rx * radius + wobble,
        y: cy + Math.sin(angle) * ry * radius + Math.cos(index * 3.17) * .4,
        size: index % 43 === 0 ? 4.5 : index % 10 === 0 ? 3.2 : 2.15,
        opacity: cluster === 0 ? .88 + (index % 4) * .035 : .94 + (index % 3) * .025,
      };
    });
  })
);

function flag(value: unknown, name: string) {
  return value ? ` ${name}` : "";
}

function quoted(value: unknown) {
  return String(value).includes(" ") ? `"${value}"` : String(value);
}

// Mirrors Embedder.tag() in sse_tools/embedders/, and --label overriding it in
// sse_coordinates.py, so the panel can name the cache files a run will touch.
function embeddingTag(v: FieldState) {
  if (v.label) return String(v.label);
  if (v.embedder === "esmc") return `esmc${String(v.esmcModel).replace("esmc_", "")}_${v.pooling}`;
  return `${v.embedder}_${v.pooling}`;
}

function buildCommand(tool: ToolId, v: FieldState, entry: string) {
  const root = "python scripts/";
  if (tool === "initialize") {
    let cmd = `${root}sse_initialization.py ${quoted(v.input)} --source ${v.source}`;
    if (v.name) cmd += ` --name ${quoted(v.name)}`;
    if (v.source === "em") {
      cmd += ` --id_col ${quoted(v.idCol)} --seq_col ${quoted(v.seqCol)}`;
    }
    if (v.query) cmd += ` --query ${v.query}`;
    return cmd + flag(v.force, "--force");
  }
  if (tool === "taxonomy") {
    let cmd = `${root}fetch_taxonomy.py ${entry} --email ${quoted(v.email)} --strategy ${v.strategy}`;
    if (v.apiKey) cmd += " --api-key ••••••••";
    cmd += ` --batch ${v.batch} --gmgc-batch ${v.gmgcBatch}`;
    if (v.rerun === "retry") cmd += " --retry-failed";
    if (v.rerun === "force") cmd += " --force";
    return cmd;
  }
  if (tool === "external") {
    let cmd = `${root}merge_external.py ${entry} ${quoted(v.file)}`;
    if (v.idCol) cmd += ` --id-col ${quoted(v.idCol)}`;
    if (v.columns) cmd += ` --columns ${quoted(v.columns)}`;
    if (v.translator) cmd += ` --translator ${quoted(v.translator)}`;
    cmd += ` --type ${v.type}`;
    if (v.delimiter !== "auto") cmd += ` --delimiter ${quoted(v.delimiter)}`;
    return cmd + flag(v.force, "--force");
  }
  if (tool === "coordinates") {
    let cmd = `${root}sse_coordinates.py ${entry} --embedder ${v.embedder} --reducer ${v.reducer}`;
    cmd += ` --pooling ${v.pooling} --device ${v.device} --n-components ${v.components}`;
    cmd += v.normalize ? " --normalize" : " --no-normalize";
    if (v.embedder === "esmc") cmd += ` --esmc-model ${v.esmcModel}`;
    if (v.embedder === "prostt5") cmd += ` --prostt5-checkpoint ${quoted(v.prostCheckpoint)}`;
    if (v.embedder === "saprot") cmd += ` --saprot-checkpoint ${quoted(v.saprotCheckpoint)}`;
    if (v.reducer === "umap") {
      cmd += ` --umap-neighbors ${v.neighbors} --umap-min-dist ${v.minDist} --umap-metric ${v.metric}`;
    }
    if (v.reducer === "tsne") {
      cmd += ` --tsne-perplexity ${v.perplexity} --tsne-pca ${v.tsnePca}`;
    }
    cmd += ` --batch-size ${v.batchSize} --write-every ${v.writeEvery}`;
    if (v.embedder !== "esmc") {
      cmd += ` --max-residues ${v.maxResidues}`;
      if (v.foldseekJson) cmd += ` --foldseek-json ${quoted(v.foldseekJson)}`;
      cmd += flag(v.includeEmpty, "--include_empty");
    }
    if (v.label) cmd += ` --label ${quoted(v.label)}`;
    if (v.rerun === "rereduce") cmd += " --rereduce";
    if (v.rerun === "reembed") cmd += " --reembed";
    return cmd;
  }
  if (tool === "distance") {
    let cmd = `${root}sse_esmc_distance.py ${entry} --embedding ${quoted(v.embedding)}`;
    if (v.queryMode === "explicit" && v.queryIds) {
      cmd += ` --query-id ${String(v.queryIds).replaceAll(",", " ")}`;
    }
    return cmd + flag(v.raw, "--raw") + flag(v.force, "--force");
  }
  if (tool === "cluster") {
    let cmd = `${root}sse_cluster.py ${entry} --embedding ${quoted(v.embedding)} --clusterer ${v.clusterer}`;
    cmd += ` --space ${v.space}`;
    if (v.space === "pca") {
      cmd += v.pcaMode === "variance" ? ` --pca-variance ${v.pcaVariance}` : ` --pca-dims ${v.pcaDims}`;
    }
    if (v.clusterer === "kmeans") {
      cmd += v.kMode === "fixed" ? ` --k ${v.k}` : ` --k-min ${v.kMin} --k-max ${v.kMax}`;
    } else {
      cmd += ` --min-cluster-size ${v.minClusterSize}`;
      if (v.minSamples) cmd += ` --min-samples ${v.minSamples}`;
    }
    if (v.label) cmd += ` --label ${quoted(v.label)}`;
    cmd += flag(v.raw, "--raw") + flag(v.force, "--force");
    if (!v.analysis) cmd += " --no-analysis";
    else cmd += ` --analysis-top-n ${v.topN} --fdr ${v.fdr}`;
    return cmd;
  }
  if (tool === "boltz") {
    let cmd = `${root}sse_boltz.py ${entry}`;
    if (v.selection) cmd += ` --selection ${quoted(v.selection)}`;
    if (v.apiKey) cmd += " --api-key ••••••••";
    if (v.smiles) cmd += ` --smiles ${quoted(v.smiles)}`;
    if (v.smilesLabel) cmd += ` --smiles-label ${quoted(v.smilesLabel)}`;
    if (!v.useMsa) cmd += " --no-msa";
    cmd += ` --recycling-steps ${v.recyclingSteps} --sampling-steps ${v.samplingSteps}`;
    cmd += ` --diffusion-samples ${v.diffusionSamples} --step-scale ${v.stepScale}`;
    cmd += flag(v.force, "--force");
    if (v.rmsd) {
      cmd += ` --rmsd --rmsd-reference ${quoted(v.rmsdReference)} --rmsd-ref-rank ${v.rmsdRefRank}`;
      cmd += ` --rmsd-method ${v.rmsdMethod} --rmsd-scope ${v.rmsdScope}`;
    }
    return cmd;
  }
  return `${root}sse_visualizer.py ${entry} --port ${v.port}`;
}

function Field({
  label,
  hint,
  children,
}: {
  label: string;
  hint?: string;
  children: React.ReactNode;
}) {
  return (
    <label className="field">
      <span className="field-label">{label}</span>
      {children}
      {hint ? <span className="field-hint">{hint}</span> : null}
    </label>
  );
}

function FilePathInput({
  value,
  onChange,
  accept,
  runnerOnline,
}: {
  value: string;
  onChange: (value: string) => void;
  accept?: string;
  runnerOnline: boolean;
}) {
  const [uploading, setUploading] = useState<number | null>(null);
  const [error, setError] = useState("");
  const picker = useRef<HTMLInputElement>(null);

  const chooseFile = async (file?: File) => {
    if (!file) return;
    setError("");
    setUploading(0);
    try {
      const uploaded = await uploadRunnerFile(file, setUploading);
      onChange(uploaded.path);
    } catch (reason) {
      setError(reason instanceof Error ? reason.message : "Upload failed");
    } finally {
      setUploading(null);
    }
  };

  return (
    <span className="file-path-input">
      <input value={value} onChange={(event) => onChange(event.target.value)} />
      <span className="file-upload-row">
        <span>{uploading === null ? "Use a path available to the runner, or upload a file." : `Uploading… ${uploading}%`}</span>
        <button type="button" disabled={!runnerOnline || uploading !== null} onClick={() => picker.current?.click()}>Browse</button>
        <input ref={picker} className="visually-hidden" type="file" accept={accept} disabled={!runnerOnline || uploading !== null} onChange={(event) => void chooseFile(event.target.files?.[0])} />
      </span>
      {error ? <span className="field-error">{error}</span> : null}
    </span>
  );
}

function Toggle({
  checked,
  onChange,
  label,
  description,
}: {
  checked: boolean;
  onChange: (value: boolean) => void;
  label: string;
  description?: string;
}) {
  return (
    <label className="toggle-row">
      <input type="checkbox" checked={checked} onChange={(e) => onChange(e.target.checked)} />
      <span className="toggle-control" aria-hidden="true"><span /></span>
      <span>
        <strong>{label}</strong>
        {description ? <small>{description}</small> : null}
      </span>
    </label>
  );
}

// The pipeline picks a vector geometry exactly once, in sse_coordinates.py: it
// L2-normalizes before reducing and writes embeddings/normalized/<tag>.emb.tsv.
// No later tool normalizes anything - they just resolve that cache. This panel
// is where that contract has to be legible, because nowhere else can show it.
function GeometrySection({
  values,
  setValue,
  entry,
}: {
  values: FieldState;
  setValue: (key: string, value: string | number | boolean) => void;
  entry?: EntrySummary;
}) {
  const normalize = Boolean(values.normalize);
  const tag = embeddingTag(values);
  const caches = entry?.embeddings ?? [];
  const stranded = !normalize && caches.some((item) => item.tag === tag && item.normalized);

  return (
    <div className="subsection geometry">
      <div className="subsection-title">Vector geometry</div>
      <p className="geometry-lede">
        The pipeline chooses a geometry <b>once, here</b>. This run writes the vectors it reduced,
        and <b>Query distances</b> and <b>Cluster space</b> read those same vectors back from disk.
        Neither of them normalizes anything on its own, so whatever you pick below is the geometry
        every later analysis inherits.
      </p>

      <Toggle
        checked={normalize}
        onChange={(value) => setValue("normalize", value)}
        label="L2-normalize embeddings"
        description={normalize
          ? "Vectors are scaled to unit length, so distance and clustering are explicitly cosine geometry."
          : "Coordinates, distances, and clustering all fall back to raw, unscaled model output."}
      />

      <dl className="geometry-contract">
        <div>
          <dt>Reduces</dt>
          <dd>{normalize ? "Unit-length vectors" : "Raw model output, unscaled"}</dd>
        </div>
        <div>
          <dt>Writes</dt>
          <dd>
            <code>embeddings/{tag}.emb.tsv</code>
            {normalize
              ? <> plus its normalized sibling <code>embeddings/normalized/{tag}.emb.tsv</code></>
              : <> only. No normalized sibling is written.</>}
          </dd>
        </div>
        <div>
          <dt>Downstream</dt>
          <dd>{normalize
            ? "Query distances and Cluster space prefer the normalized sibling automatically."
            : "Query distances and Cluster space read the raw cache."}</dd>
        </div>
        <div>
          <dt>Recorded as</dt>
          <dd><code>normalize={normalize ? "l2" : "none"}</code> in the provenance of every column this run writes</dd>
        </div>
      </dl>

      <div className="geometry-state">
        <div className="geometry-state-title">Embedding caches in this entry</div>
        {caches.length ? (
          <ul>
            {caches.map((item) => (
              <li key={item.tag} className={item.normalized ? "l2" : "raw"}>
                <code>{item.tag}</code>
                <span>{item.normalized ? "L2-normalized sibling on disk" : "Raw only"}</span>
                {item.tag === tag ? <b>this run</b> : null}
              </li>
            ))}
          </ul>
        ) : (
          <p>None yet. This run creates the first one and sets the geometry for everything after it.</p>
        )}
      </div>

      {stranded ? (
        <div className="notice warning geometry-alert">
          <span>Mixed geometry</span>
          <span className="notice-body">
            A normalized cache already exists at <code>embeddings/normalized/{tag}.emb.tsv</code>. These
            coordinates would be reduced from raw vectors while Query distances and Cluster space keep
            preferring that normalized file, leaving one entry with two geometries. Keep normalization
            on, or delete that file first.
          </span>
        </div>
      ) : null}
    </div>
  );
}

// Distance and clustering never re-normalize; they resolve whichever cache
// Build coordinates left behind. Say so, rather than implying a free choice.
function EmbeddingCacheField({
  values,
  setValue,
  entry,
}: {
  values: FieldState;
  setValue: (key: string, value: string | number | boolean) => void;
  entry?: EntrySummary;
}) {
  const caches = entry?.embeddings ?? [];
  const selected = String(values.embedding ?? "");
  const cache = caches.find((item) => item.tag === selected);
  const usingNormalized = Boolean(cache?.normalized) && !values.raw;

  return (
    <>
      <Field label="Embedding cache" hint="Built by Build coordinates. Its geometry is fixed on disk; this tool only reads it.">
        <select value={selected} onChange={(event) => setValue("embedding", event.target.value)}>
          {caches.length
            ? caches.map((item) => (
              <option key={item.tag} value={item.tag}>{item.tag} · {item.normalized ? "normalized available" : "raw only"}</option>
            ))
            : <option value={selected}>{selected || "No embedding cache found"}</option>}
        </select>
      </Field>
      <div className={`notice ${usingNormalized ? "success" : "warning"}`}>
        <span>{usingNormalized ? "L2 geometry" : "Raw geometry"}</span>
        <span className="notice-body">{usingNormalized
          ? `Reads embeddings/normalized/${selected}.emb.tsv and records normalize=l2.`
          : cache && !cache.normalized
            ? `${selected} has no normalized sibling, so this reads raw vectors and records normalize=none. Re-run Build coordinates with normalization on to change that.`
            : `Raw vectors requested, so this reads embeddings/${selected}.emb.tsv and records normalize=none.`}</span>
      </div>
    </>
  );
}

function ConfigurationForm({
  tool,
  values,
  setValue,
  runnerOnline,
  entry,
}: {
  tool: ToolId;
  values: FieldState;
  setValue: (key: string, value: string | number | boolean) => void;
  runnerOnline: boolean;
  entry?: EntrySummary;
}) {
  const input = (key: string, type = "text") => ({
    type,
    value: String(values[key] ?? ""),
    onChange: (e: React.ChangeEvent<HTMLInputElement>) =>
      setValue(key, type === "number" ? Number(e.target.value) : e.target.value),
  });
  const select = (key: string) => ({
    value: String(values[key] ?? ""),
    onChange: (e: React.ChangeEvent<HTMLSelectElement>) => setValue(key, e.target.value),
  });

  if (tool === "initialize") {
    return (
      <>
        <div className="field-grid two">
          <Field label="Source file" hint="A path or a filename in initial_files/">
            <FilePathInput value={String(values.input || "")} onChange={(value) => setValue("input", value)} accept=".tsv,.csv,.fasta,.fa,.faa,.json" runnerOnline={runnerOnline} />
          </Field>
          <Field label="Source format">
            <select {...select("source")}>
              <option value="em">EnzymeMiner / generic TSV</option>
              <option value="fasta">FASTA</option>
              <option value="fs">Foldseek JSON</option>
            </select>
          </Field>
        </div>
        <Field label="Entry name" hint="Defaults to the source filename stem">
          <input {...input("name")} />
        </Field>
        {values.source === "em" ? (
          <div className="field-grid two">
            <Field label="ID column"><input {...input("idCol")} /></Field>
            <Field label="Sequence column"><input {...input("seqCol")} /></Field>
          </div>
        ) : null}
        <Field label="Query values" hint={values.source === "fasta" ? "Use complete FASTA headers; leave blank for no marked query" : "Space-separated values; overrides automatic query detection"}>
          <input {...input("query")} placeholder="Optional" />
        </Field>
        <Toggle checked={Boolean(values.force)} onChange={(v) => setValue("force", v)} label="Delete and rebuild an existing entry" description="Destructive: this removes the complete entry directory." />
      </>
    );
  }

  if (tool === "taxonomy") {
    return (
      <>
        <div className="field-grid two">
          <Field label="NCBI contact email"><input {...input("email", "email")} /></Field>
          <Field label="Resolution strategy">
            <select {...select("strategy")}>
              <option value="auto">Auto-detect</option>
              <option value="em">NCBI protein accession</option>
              <option value="foldseek">Foldseek metadata</option>
            </select>
          </Field>
        </div>
        <Field label="NCBI API key" hint="Optional; passed to this job without being written to job history">
          <input {...input("apiKey", "password")} placeholder="Optional" />
        </Field>
        <div className="field-grid two">
          <Field label="NCBI batch size"><input {...input("batch", "number")} min="1" /></Field>
          <Field label="GMGC batch size"><input {...input("gmgcBatch", "number")} min="1" /></Field>
        </div>
        <Field label="Run mode">
          <select {...select("rerun")}>
            <option value="resume">Resume / preserve completed rows</option>
            <option value="retry">Retry unresolved rows</option>
            <option value="force">Refetch every row</option>
          </select>
        </Field>
      </>
    );
  }

  if (tool === "external") {
    return (
      <>
        <Field label="External CSV or TSV"><FilePathInput value={String(values.file || "")} onChange={(value) => setValue("file", value)} accept=".csv,.tsv,.txt" runnerOnline={runnerOnline} /></Field>
        <div className="field-grid two">
          <Field label="External ID column" hint="Blank uses the first column"><input {...input("idCol")} placeholder="First column" /></Field>
          <Field label="Column type">
            <select {...select("type")}><option value="label">Label / metadata</option><option value="coordinate">Coordinate / plot axis</option></select>
          </Field>
        </div>
        <Field label="Columns to merge" hint="Comma-separated; blank merges every non-ID column"><input {...input("columns")} /></Field>
        <Field label="ID translator" hint="Optional two-column table: SSE ID, then external ID"><FilePathInput value={String(values.translator || "")} onChange={(value) => setValue("translator", value)} accept=".csv,.tsv,.txt" runnerOnline={runnerOnline} /></Field>
        <Field label="Delimiter">
          <select {...select("delimiter")}><option value="auto">Infer from extension</option><option value=",">Comma</option><option value="\t">Tab</option><option value=";">Semicolon</option></select>
        </Field>
        <Toggle checked={Boolean(values.force)} onChange={(v) => setValue("force", v)} label="Replace colliding columns" />
      </>
    );
  }

  if (tool === "coordinates") {
    const structure = values.embedder !== "esmc";
    return (
      <>
        <div className="choice-grid three" role="radiogroup" aria-label="Embedding model">
          {[{id:"esmc",name:"ESM-C",sub:"Sequence"},{id:"prostt5",name:"ProstT5",sub:"3Di structure"},{id:"saprot",name:"SaProt",sub:"Sequence + 3Di"}].map((item) => (
            <button key={item.id} className={`choice-card ${values.embedder === item.id ? "selected" : ""}`} onClick={() => setValue("embedder", item.id)} type="button" aria-pressed={values.embedder === item.id}>
              <span className="choice-dot" /><strong>{item.name}</strong><small>{item.sub}</small>
            </button>
          ))}
        </div>
        {structure ? <div className="notice warning"><span>Foldseek required</span> This model needs a Foldseek entry and source JSON with C-alpha coordinates.</div> : null}
        <div className="field-grid three">
          <Field label="Reduction">
            <select {...select("reducer")}><option value="pca">PCA</option><option value="umap">UMAP</option><option value="tsne">t-SNE</option></select>
          </Field>
          <Field label="Components"><input {...input("components", "number")} min="1" /></Field>
          <Field label="Pooling"><select {...select("pooling")}><option>mean</option><option>max</option><option>min</option></select></Field>
        </div>
        <div className="field-grid two">
          <Field label="Compute device"><select {...select("device")}><option value="auto">Auto: CUDA → MPS → CPU</option><option value="cuda">CUDA</option><option value="mps">Apple MPS</option><option value="cpu">CPU</option></select></Field>
          {values.embedder === "esmc" ? (
            <Field label="ESM-C model"><select {...select("esmcModel")}><option value="esmc_600m">ESM-C 600M</option><option value="esmc_300m">ESM-C 300M</option></select></Field>
          ) : (
            <Field label="Model checkpoint"><input {...input(values.embedder === "prostt5" ? "prostCheckpoint" : "saprotCheckpoint")} /></Field>
          )}
        </div>
        <GeometrySection values={values} setValue={setValue} entry={entry} />

        {values.reducer === "umap" ? (
          <div className="subsection">
            <div className="subsection-title">UMAP neighborhood</div>
            <div className="field-grid three">
              <Field label="Neighbors"><input {...input("neighbors", "number")} min="2" /></Field>
              <Field label="Minimum distance"><input {...input("minDist", "number")} min="0" step="0.05" /></Field>
              <Field label="Metric"><input {...input("metric")} /></Field>
            </div>
          </div>
        ) : null}
        {values.reducer === "tsne" ? (
          <div className="subsection">
            <div className="subsection-title">t-SNE neighborhood</div>
            <div className="field-grid two">
              <Field label="Perplexity"><input {...input("perplexity", "number")} min="1" /></Field>
              <Field label="PCA pre-reduction"><input {...input("tsnePca", "number")} min="0" /></Field>
            </div>
          </div>
        ) : null}
        {structure ? (
          <div className="subsection">
            <div className="subsection-title">Structure input</div>
            <div className="field-grid two">
              <Field label="Foldseek JSON" hint="Blank resolves it from the manifest"><FilePathInput value={String(values.foldseekJson || "")} onChange={(value) => setValue("foldseekJson", value)} accept=".json" runnerOnline={runnerOnline} /></Field>
              <Field label="Maximum residues"><input {...input("maxResidues", "number")} min="3" /></Field>
            </div>
            <Toggle checked={Boolean(values.includeEmpty)} onChange={(v) => setValue("includeEmpty", v)} label="Keep rows without usable structures" description="Their coordinate cells remain empty." />
          </div>
        ) : null}
        <details className="advanced">
          <summary>Advanced run and cache options</summary>
          <div className="field-grid three advanced-body">
            <Field label="Batch size"><input {...input("batchSize", "number")} min="1" /></Field>
            <Field label="Write every"><input {...input("writeEvery", "number")} min="1" /></Field>
            <Field label="Custom tag"><input {...input("label")} placeholder="Automatic" /></Field>
          </div>
          <Field label="If this coordinate system already exists">
            <select {...select("rerun")}><option value="new">Stop and ask</option><option value="rereduce">Reuse embeddings and re-reduce</option><option value="reembed">Recompute embeddings and reduce</option></select>
          </Field>
        </details>
      </>
    );
  }

  if (tool === "distance") {
    return (
      <>
        <EmbeddingCacheField values={values} setValue={setValue} entry={entry} />
        <Field label="Reference sequences">
          <select {...select("queryMode")}><option value="marked">Use rows marked query=True</option><option value="explicit">Choose explicit IDs</option></select>
        </Field>
        {values.queryMode === "explicit" ? <Field label="Query IDs" hint="Comma-separated SSE IDs"><input {...input("queryIds")} /></Field> : <div className="notice success"><span>{formatCount(entry?.queries)} queries detected</span> Rows marked query=True in this entry.</div>}
        <Toggle checked={Boolean(values.raw)} onChange={(v) => setValue("raw", v)} label="Use raw embeddings" description="Ignore the normalized sibling and measure distances on raw vectors instead." />
        <Toggle checked={Boolean(values.force)} onChange={(v) => setValue("force", v)} label="Replace existing distance columns" />
      </>
    );
  }

  if (tool === "cluster") {
    return (
      <>
        <EmbeddingCacheField values={values} setValue={setValue} entry={entry} />
        <div className="choice-grid two" role="radiogroup" aria-label="Clustering technique">
          {[{id:"kmeans",name:"K-means",sub:"Assign every sequence"},{id:"hdbscan",name:"HDBSCAN",sub:"Discover groups + noise"}].map((item) => (
            <button key={item.id} className={`choice-card ${values.clusterer === item.id ? "selected" : ""}`} onClick={() => setValue("clusterer", item.id)} type="button" aria-pressed={values.clusterer === item.id}>
              <span className="choice-dot" /><strong>{item.name}</strong><small>{item.sub}</small>
            </button>
          ))}
        </div>
        <div className="field-grid two">
          <Field label="Clustering space"><select {...select("space")}><option value="pca">PCA-reduced embeddings</option><option value="full">Full embedding dimensions</option></select></Field>
          <Field label="Output tag"><input {...input("label")} placeholder="Use embedding tag" /></Field>
        </div>
        {values.space === "pca" ? (
          <div className="subsection">
            <div className="segmented">
              <button type="button" className={values.pcaMode === "dims" ? "active" : ""} onClick={() => setValue("pcaMode", "dims")}>Fixed dimensions</button>
              <button type="button" className={values.pcaMode === "variance" ? "active" : ""} onClick={() => setValue("pcaMode", "variance")}>Variance target</button>
            </div>
            <Field label={values.pcaMode === "dims" ? "PCA dimensions" : "Variance fraction"}>
              <input {...input(values.pcaMode === "dims" ? "pcaDims" : "pcaVariance", "number")} step={values.pcaMode === "dims" ? 1 : 0.01} />
            </Field>
          </div>
        ) : null}
        {values.clusterer === "kmeans" ? (
          <div className="subsection">
            <div className="segmented">
              <button type="button" className={values.kMode === "auto" ? "active" : ""} onClick={() => setValue("kMode", "auto")}>Auto-select k</button>
              <button type="button" className={values.kMode === "fixed" ? "active" : ""} onClick={() => setValue("kMode", "fixed")}>Fixed k</button>
            </div>
            {values.kMode === "auto" ? <div className="field-grid two"><Field label="Minimum k"><input {...input("kMin", "number")} /></Field><Field label="Maximum k"><input {...input("kMax", "number")} /></Field></div> : <Field label="Number of clusters"><input {...input("k", "number")} /></Field>}
          </div>
        ) : (
          <div className="field-grid two subsection">
            <Field label="Minimum cluster size"><input {...input("minClusterSize", "number")} /></Field>
            <Field label="Minimum samples" hint="Blank follows cluster size"><input {...input("minSamples", "number")} placeholder="Automatic" /></Field>
          </div>
        )}
        <Toggle checked={Boolean(values.analysis)} onChange={(v) => setValue("analysis", v)} label="Generate Tier-2 analysis" description="Profiles, enrichment, and representative sequences." />
        {values.analysis ? <div className="field-grid two"><Field label="Representatives per cluster"><input {...input("topN", "number")} /></Field><Field label="FDR threshold"><input {...input("fdr", "number")} step="0.01" /></Field></div> : null}
        <Toggle checked={Boolean(values.raw)} onChange={(v) => setValue("raw", v)} label="Use raw embedding geometry" description="Ignore the normalized sibling and cluster on raw vectors instead." />
        <Toggle checked={Boolean(values.force)} onChange={(v) => setValue("force", v)} label="Replace an existing matching clustering" />
      </>
    );
  }

  if (tool === "boltz") {
    const selections = entry?.selections ?? [];
    const activeSelection = (values.selection ? selections.find((sel) => sel.name === values.selection) : selections[0]) ?? selections[0];
    const refOptions = activeSelection?.ids ?? [];
    const datafileName = `${entry?.name ?? "<entry>"}.sse.tsv`;
    return (
      <>
        <Field label="Selection to analyze" hint="Exported from the visualizer's “Export selection for Boltz” button">
          {selections.length ? (
            <select {...select("selection")}>
              <option value="">Most recent selection</option>
              {selections.map((sel) => (
                <option key={sel.name} value={sel.name}>{sel.name}{sel.count != null ? ` (${sel.count} seq)` : ""}</option>
              ))}
            </select>
          ) : (
            <div className="notice"><span>No selections found</span> Select points in the explorer and click “Export selection for Boltz”, then reload this page.</div>
          )}
        </Field>
        <Field label="NVIDIA API key" hint="Required; passed to this job without being written to job history">
          <input {...input("apiKey", "password")} placeholder="nvapi-…" />
        </Field>
        <Field label="Substrate SMILES" hint="Optional; one per line. Present = holo prediction, empty = apo">
          <textarea value={String(values.smiles ?? "")} onChange={(e) => setValue("smiles", e.target.value)} placeholder="One SMILES per line" rows={2} />
        </Field>
        <Field label="Ligand label" hint="Optional; names the holo output columns/folders">
          <input {...input("smilesLabel")} placeholder="e.g. UDP-Glc" />
        </Field>
        <Toggle checked={Boolean(values.useMsa)} onChange={(v) => setValue("useMsa", v)} label="Generate MSA (recommended)" description="Uses ColabFold; more accurate but slower. Off runs single-sequence." />
        <details className="advanced">
          <summary>Prediction parameters</summary>
          <div className="field-grid two advanced-body">
            <Field label="Recycling steps"><input {...input("recyclingSteps", "number")} min="1" max="10" /></Field>
            <Field label="Sampling steps"><input {...input("samplingSteps", "number")} min="10" max="500" /></Field>
            <Field label="Diffusion samples"><input {...input("diffusionSamples", "number")} min="1" max="10" /></Field>
            <Field label="Step scale"><input {...input("stepScale", "number")} step="0.001" min="0.1" max="5" /></Field>
          </div>
        </details>
        <Toggle checked={Boolean(values.rmsd)} onChange={(v) => setValue("rmsd", v)} label="Measure structural RMSD" description="After prediction, run the RMSD comparison (Kabsch superposition) over the predicted apo structures." />
        {values.rmsd ? (
          <div className="subsection">
            <div className="subsection-title">RMSD structural comparison</div>
            <div className="notice"><div className="notice-body">Every predicted apo structure is aligned to one <strong>reference structure</strong> and its RMSD is measured. Results are appended to <code>{datafileName}</code> as <code>RMSD_vs_&lt;reference&gt;_r&lt;rank&gt;_&lt;method&gt;</code> columns — reload the explorer to color by them.</div></div>
            <Field label="Reference structure" hint="The structure every other one is compared against — one of the analyzed sequences">
              {refOptions.length ? (
                <select {...select("rmsdReference")}>
                  <option value="">Select a reference…</option>
                  {refOptions.map((id) => <option key={id} value={id}>{id}</option>)}
                </select>
              ) : (
                <input {...input("rmsdReference")} placeholder="e.g. OleD_S1" />
              )}
            </Field>
            <div className="field-grid two">
              <Field label="Reference rank" hint="Which predicted rank of the reference"><input {...input("rmsdRefRank", "number")} min="0" /></Field>
              <Field label="Alignment method"><select {...select("rmsdMethod")}><option value="seq">Sequence-guided</option><option value="ce">Structure-based (CE)</option><option value="both">Both</option></select></Field>
            </div>
            <Field label="Compare against" hint="Which query structures to measure against the reference">
              <select {...select("rmsdScope")}><option value="all">All predicted apo structures</option><option value="selected">Only this selection’s sequences</option></select>
            </Field>
          </div>
        ) : null}
        <Toggle checked={Boolean(values.force)} onChange={(v) => setValue("force", v)} label="Force re-run (ignore cached predictions)" />
        <div className="notice"><div className="notice-body">Outputs: ranked <code>.cif</code> structures under <code>structures/</code>, binding scores in <code>logs/boltz_log.csv</code>, and pTM/pLDDT (plus any RMSD) columns appended to <code>{datafileName}</code>. Use “Open structures folder” on the finished job to browse them.</div></div>
      </>
    );
  }

  return (
    <>
      <div className="notice success"><span>Entry ready</span> 22 coordinate columns and 2 clustering systems detected.</div>
      <Field label="Dash server port"><input {...input("port", "number")} min="1024" max="65535" /></Field>
      <div className="explorer-preview">
        <div className="mini-plot" aria-label="Illustrative sequence-space scatter plot">
          {Array.from({ length: 34 }).map((_, index) => <i key={index} style={{ "--x": `${8 + ((index * 29) % 82)}%`, "--y": `${10 + ((index * 47) % 75)}%`, "--d": `${(index % 5) * 0.04}s` } as React.CSSProperties} />)}
        </div>
        <div><strong>Sequence Space Explorer</strong><span>Coordinates, filters, layers, Boltz-2, RMSD, and exports</span></div>
      </div>
    </>
  );
}

function LandingPage({ onEnter }: { onEnter: () => void }) {
  return (
    <section className="landing-page">
      <div className="landing-copy">
        <div className="landing-kicker"><span /> Protein sequence intelligence</div>
        <h1>Sequence Space<br /><em>Explorer</em></h1>
        <p>Build, enrich, map, and analyze protein sequence spaces through one guided visual workflow.</p>
        <div className="landing-actions">
          <button type="button" className="landing-primary" onClick={onEnter}>Enter pipeline <span>→</span></button>
          <span className="landing-note"><b>{tools.length}</b> integrated tools · terminal optional</span>
        </div>
        <div className="landing-capabilities" aria-label="Pipeline capabilities">
          <span>Embed</span><i />
          <span>Reduce</span><i />
          <span>Cluster</span><i />
          <span>Explore</span>
        </div>
      </div>

      <div className="landing-visual" aria-label="Stylized sequence-space visualization example">
        <div className="hero-plot">
          <span className="density-field density-slate" />
          <span className="density-field density-magenta" />
          <span className="density-field density-gray" />
          <span className="density-field density-teal" />
          <span className="density-field density-purple" />
          {landingPoints.map((point, index) => (
            <i
              key={index}
              className={`data-point population-${point.cluster}`}
              style={{
                "--x": `${point.x}%`,
                "--y": `${point.y}%`,
                "--s": `${point.size}px`,
                "--c": point.color,
                "--o": point.opacity,
              } as React.CSSProperties}
            />
          ))}
        </div>
      </div>
    </section>
  );
}

const terminalStatuses = new Set(["succeeded", "failed", "cancelled", "interrupted"]);

function jobStatusLabel(job: PipelineJob) {
  if (job.ready) return "Explorer ready";
  return ({ queued: "Queued", running: "Running", cancelling: "Stopping", succeeded: "Complete", failed: "Failed", cancelled: "Cancelled", interrupted: "Interrupted" } as const)[job.status];
}

function formatCount(value: number | null | undefined) {
  return typeof value === "number" ? value.toLocaleString() : "—";
}

export default function Home() {
  const [tool, setTool] = useState<ToolId>("coordinates");
  const [entry, setEntry] = useState("EnzymeMiner_Selection_Table_ri4plk");
  const [entries, setEntries] = useState<EntrySummary[]>([]);
  const [values, setValues] = useState<Record<ToolId, FieldState>>(initialValues);
  const [job, setJob] = useState<PipelineJob | null>(null);
  const [runnerState, setRunnerState] = useState<"connecting" | "online" | "offline">("connecting");
  const [environmentReady, setEnvironmentReady] = useState(false);
  const [environmentMissing, setEnvironmentMissing] = useState<string[]>([]);
  const [runnerError, setRunnerError] = useState("");
  const [submitting, setSubmitting] = useState(false);
  const [copied, setCopied] = useState(false);
  const [showTools, setShowTools] = useState(true);
  const [theme, setTheme] = useState<"modern" | "original">(() => {
    if (typeof window === "undefined") return "modern";
    const saved = window.localStorage.getItem("sse-theme");
    return saved === "original" ? "original" : "modern";
  });
  const [showLanding, setShowLanding] = useState(true);
  const [entryFiles, setEntryFiles] = useState<EntryFile[] | null>(null);
  const [showFullLog, setShowFullLog] = useState(false);
  const [shuttingDown, setShuttingDown] = useState(false);
  const explorerWindow = useRef<Window | null>(null);
  const active = tools.find((item) => item.id === tool)!;
  const activeEntry = entries.find((item) => item.name === entry);
  const previewCommand = useMemo(() => buildCommand(tool, values[tool], entry), [tool, values, entry]);
  const command = job?.tool === tool && job.entry === (tool === "initialize" ? String(values.initialize.name || entry) : entry) ? job.command : previewCommand;
  const running = job && ["queued", "running", "cancelling"].includes(job.status);
  // A launched visualizer is a long-lived Dash service that never exits on its
  // own, so its job stays "running" and would otherwise block every run button
  // until a hard restart. Surfacing an explicit shut-down turns the same launch
  // control into the way back out.
  const visualizerRunning = Boolean(job && job.tool === "visualizer" && ["queued", "running", "cancelling"].includes(job.status));
  const jobId = job?.id;
  const jobStatus = job?.status;

  const refreshEntries = useCallback(async () => {
    try {
      const state = await getRunnerState();
      setEntries(state.entries);
      setRunnerState("online");
      setEnvironmentReady(state.health.environment.ready);
      setEnvironmentMissing(state.health.environment.missing);
      setRunnerError(state.health.environment.ready ? "" : `The selected Python environment is missing: ${state.health.environment.missing.join(", ")}. Install requirements.txt in that environment, then restart the runner.`);
      if (state.entries.length && !state.entries.some((item) => item.name === entry)) setEntry(state.entries[0].name);
    } catch (reason) {
      setRunnerState("offline");
      setEnvironmentReady(false);
      setRunnerError(reason instanceof Error ? reason.message : "The local runner is unavailable");
    }
  }, [entry]);

  useEffect(() => {
    const timer = window.setTimeout(() => void refreshEntries(), 0);
    return () => window.clearTimeout(timer);
  }, [refreshEntries]);

  useEffect(() => {
    if (runnerState !== "online" || !entry || job) return;
    void getJobs(entry).then((history) => { if (history[0]) setJob(history[0]); }).catch(() => undefined);
  }, [entry, job, runnerState]);

  useEffect(() => {
    if (!jobId || !jobStatus || terminalStatuses.has(jobStatus)) return;
    const interval = window.setInterval(async () => {
      try {
        const updated = await getJob(jobId);
        setJob(updated);
        if (terminalStatuses.has(updated.status)) void refreshEntries();
      } catch (reason) {
        setRunnerError(reason instanceof Error ? reason.message : "Could not refresh the job");
      }
    }, 1000);
    return () => window.clearInterval(interval);
  }, [jobId, jobStatus, refreshEntries]);

  useEffect(() => {
    if (job?.tool !== "visualizer" || !explorerWindow.current) return;
    const pendingWindow = explorerWindow.current;
    if (pendingWindow.closed) {
      explorerWindow.current = null;
      return;
    }
    if (job.ready && job.url) {
      const destination = `${job.url.replace(/\/$/, "")}/`;
      pendingWindow.location.replace(destination);
      explorerWindow.current = null;
      return;
    }
    if (terminalStatuses.has(job.status)) {
      pendingWindow.document.title = "Explorer did not start";
      pendingWindow.document.body.innerHTML = "<main style='max-width:620px;margin:15vh auto;padding:32px;font:16px/1.55 system-ui;color:#dcebe8;background:#002a32;border-radius:16px'><h1 style='margin-top:0'>Explorer did not start</h1><p>Return to Sequence Space Explorer and review the job log for details.</p></main>";
      explorerWindow.current = null;
    }
  }, [job]);

  const setValue = (key: string, value: string | number | boolean) => {
    setValues((current) => ({ ...current, [tool]: { ...current[tool], [key]: value } }));
  };

  const chooseTheme = (nextTheme: "modern" | "original") => {
    setTheme(nextTheme);
    window.localStorage.setItem("sse-theme", nextTheme);
  };

  const startJob = async () => {
    setSubmitting(true);
    setRunnerError("");
    if (tool === "visualizer") {
      const pendingWindow = window.open("", "_blank");
      if (!pendingWindow) {
        setRunnerError("The browser blocked the explorer window. Allow pop-ups for this app and try again.");
        setSubmitting(false);
        return;
      }
      pendingWindow.document.title = "Starting Sequence Space Explorer";
      pendingWindow.document.body.innerHTML = "<main style='max-width:620px;margin:15vh auto;padding:32px;font:16px/1.55 system-ui;color:#dcebe8;background:#002a32;border-radius:16px'><h1 style='margin-top:0'>Starting explorer…</h1><p>The visualizer is loading its data. This window will open it automatically when it is ready.</p></main>";
      pendingWindow.document.body.style.background = "#001f28";
      explorerWindow.current = pendingWindow;
    }
    try {
      const submitted = await submitJob({ tool, entry, values: values[tool] });
      setJob(submitted);
      if (tool === "initialize") setEntry(submitted.entry);
    } catch (reason) {
      if (tool === "visualizer" && explorerWindow.current) {
        explorerWindow.current.close();
        explorerWindow.current = null;
      }
      setRunnerError(reason instanceof Error ? reason.message : "Could not submit the job");
    } finally {
      setSubmitting(false);
    }
  };

  const cancelJob = async () => {
    if (!job) return;
    try { setJob(await cancelPipelineJob(job.id)); }
    catch (reason) { setRunnerError(reason instanceof Error ? reason.message : "Could not stop the job"); }
  };

  const shutDownApp = async () => {
    const confirmed = window.confirm("Shut down Sequence Space Explorer? Any active pipeline job will be stopped.");
    if (!confirmed) return;
    setRunnerError("");
    try {
      await requestPipelineShutdown();
      setShuttingDown(true);
    } catch (reason) {
      setRunnerError(reason instanceof Error ? reason.message : "Could not shut down the app");
    }
  };

  const copyCommand = async () => {
    try { await navigator.clipboard.writeText(command); } catch { /* clipboard can be unavailable in embedded previews */ }
    setCopied(true);
    window.setTimeout(() => setCopied(false), 1500);
  };

  const viewFiles = async () => {
    setRunnerError("");
    try { setEntryFiles(await getEntryFiles(entry)); }
    catch (reason) { setRunnerError(reason instanceof Error ? reason.message : "Could not list entry files"); }
  };

  const revealFolder = async (subpath: string) => {
    setRunnerError("");
    try { await revealEntryPath(entry, subpath); }
    catch (reason) { setRunnerError(reason instanceof Error ? reason.message : "Could not open the folder"); }
  };

  const completedTool = (id: ToolId) => Boolean(activeEntry && (
    id === "initialize" ||
    (id === "taxonomy" && activeEntry.taxonomyFields > 0) ||
    (id === "coordinates" && activeEntry.coordinates > 0) ||
    (id === "distance" && activeEntry.annotations > 0) ||
    (id === "cluster" && activeEntry.analyses > 0)
  ));

  return (
    <main className={`app-shell ${showLanding ? "landing-open" : ""}`} data-theme={theme}>
      <div className="ambient-space" aria-hidden="true">
        <span className="ambient-orbit orbit-one" /><span className="ambient-orbit orbit-two" /><span className="ambient-glow glow-one" /><span className="ambient-glow glow-two" />
        {Array.from({ length: 76 }).map((_, index) => <i key={index} style={{ "--x": `${2 + ((index * 37) % 92)}%`, "--y": `${8 + ((index * 53) % 84)}%`, "--s": `${2 + (index % 4)}px`, "--o": `${0.16 + (index % 6) * 0.09}`, "--delay": `${(index % 11) * -0.43}s` } as React.CSSProperties} />)}
      </div>
      <header className="topbar">
        <button className="brand" type="button" onClick={() => setShowLanding(true)} aria-label="Return to the Sequence Space Explorer landing page">
          <div className="brand-mark"><span /><span /><span /></div><div><strong>Sequence Space Explorer</strong><span>Pipeline control center</span></div>
        </button>
        <div className="top-actions">
          <div className={`runner-badge ${runnerState === "online" && !environmentReady ? "offline" : runnerState}`}><span /> {runnerState === "online" ? environmentReady ? "Runner connected" : "Setup needed" : runnerState === "offline" ? "Runner offline" : "Connecting"}</div>
          <div className="theme-switch" role="group" aria-label="Interface theme">
            <button type="button" className={theme === "modern" ? "active" : ""} aria-pressed={theme === "modern"} onClick={() => chooseTheme("modern")}><span>◐</span> Modern</button>
            <button type="button" className={theme === "original" ? "active" : ""} aria-pressed={theme === "original"} onClick={() => chooseTheme("original")}><span>○</span> Original</button>
          </div>
          <button className="icon-button" type="button" aria-label="Reconnect to runner" onClick={() => void refreshEntries()}>↻</button>
          <button className="shutdown-button" type="button" onClick={() => void shutDownApp()} disabled={runnerState !== "online"}><span aria-hidden="true">⏻</span><b>Shut down</b></button>
          <div className="avatar">WG</div>
        </div>
      </header>

      {showLanding ? <LandingPage onEnter={() => setShowLanding(false)} /> : <div className="workspace">
        <aside className={`sidebar ${showTools ? "" : "collapsed"}`}>
          <div className="side-heading"><span>Workspace</span><button type="button" aria-label="Collapse tools" onClick={() => setShowTools(!showTools)}>{showTools ? "‹" : "›"}</button></div>
          <div className="entry-card">
            <span className="entry-icon">E</span><span><select aria-label="Selected entry" value={entry} onChange={(event) => { const nextEntry = entries.find((item) => item.name === event.target.value); const embedding = nextEntry?.embeddings[0]?.tag; setEntry(event.target.value); setJob(null); if (embedding) setValues((current) => ({ ...current, distance: { ...current.distance, embedding }, cluster: { ...current.cluster, embedding } })); }} disabled={!entries.length}>{entries.length ? entries.map((item) => <option key={item.name}>{item.name}</option>) : <option>No entries found</option>}</select><small>{formatCount(activeEntry?.sequences)} sequences</small></span><b>⌄</b>
          </div>
          <nav aria-label="Pipeline tools">
            <div className="nav-label">Pipeline</div>
            {tools.map((item) => <button key={item.id} type="button" onClick={() => setTool(item.id)} className={tool === item.id ? "active" : ""}><span className="nav-mark">{item.mark}</span><span><strong>{item.label}</strong><small>{item.stage}</small></span>{completedTool(item.id) ? <i className="complete">✓</i> : null}</button>)}
          </nav>
          <div className="side-footer">
            <div className={`health ${runnerState === "online" && !environmentReady ? "offline" : runnerState}`}><span /><div><strong>{runnerState === "online" ? environmentReady ? "Environment ready" : "Python setup needed" : "Runner unavailable"}</strong><small>{runnerState === "online" ? environmentReady ? "Local pipeline access" : `${environmentMissing.length} packages missing` : "Start the local runner"}</small></div></div>
            <button type="button" onClick={() => void refreshEntries()}><span>↻</span> Refresh environment</button>
          </div>
        </aside>

        <section className="content">
          {runnerError ? <div className="runner-alert" role="alert"><strong>{runnerState === "offline" ? "Local runner not connected." : "Action needed."}</strong> {runnerError}<button type="button" onClick={() => void refreshEntries()}>Retry</button></div> : null}
          <div className="entry-header">
            <div><div className="eyebrow">Selected entry</div><div className="entry-title-row"><h1>{entry || "No entry selected"}</h1><span className="ready-pill">{activeEntry ? "Ready" : "Unavailable"}</span></div><p>{activeEntry ? `${activeEntry.source.toUpperCase()} source · Updated ${new Date(activeEntry.updatedAt).toLocaleString()}` : "Connect the local runner to discover entries."}</p></div>
            <div className="entry-actions"><button className="secondary-button" type="button" onClick={() => void viewFiles()} disabled={!activeEntry}>View files</button><button className="secondary-button" type="button" onClick={() => setTool("visualizer")} disabled={!activeEntry}>Open explorer ↗</button></div>
          </div>

          <div className="stats-grid">
            <div><span>Sequences</span><strong>{formatCount(activeEntry?.sequences)}</strong><small>{formatCount(activeEntry?.queries)} marked queries</small></div>
            <div><span>Annotations</span><strong>{formatCount(activeEntry?.annotations)}</strong><small>{formatCount(activeEntry?.taxonomyFields)} taxonomy fields</small></div>
            <div><span>Coordinates</span><strong>{formatCount(activeEntry?.coordinates)}</strong><small>{activeEntry?.embeddings.map((item) => item.tag).join(" · ") || "No embedding cache"}</small></div>
            <div><span>Analyses</span><strong>{formatCount(activeEntry?.analyses)}</strong><small>{activeEntry?.analyses ? "Cluster reports available" : "No cluster reports"}</small></div>
          </div>

          <div className="stepper" aria-label="Pipeline progress">{["Entry", "Enrich", "Embed", "Analyze", "Explore"].map((name, index) => { const done = Boolean(activeEntry && (index === 0 || (index === 1 && activeEntry.annotations > 2) || (index === 2 && activeEntry.coordinates > 0) || (index === 3 && activeEntry.analyses > 0))); return <div key={name} className={done ? "done" : index === 4 ? "current" : ""}><span>{done ? "✓" : index + 1}</span><b>{name}</b></div>; })}</div>

          <div className="panel-grid">
            <section className="panel configure-panel">
              <div className="panel-header"><div className="tool-heading"><span>{active.mark}</span><div><div className="eyebrow">{active.stage}</div><h2>{active.label}</h2><p>{active.description}</p></div></div><div className="panel-header-meta"><a className="guide-link" href={`${GUIDE_URL}#${active.guide}`} target="_blank" rel="noopener noreferrer" title={`Read about “${active.label}” in the visual pipeline guide`}><span aria-hidden="true">?</span> Guide</a><div className="script-chip">{active.script}</div></div></div>
              <div className="panel-body"><ConfigurationForm tool={tool} values={values[tool]} setValue={setValue} runnerOnline={runnerState === "online"} entry={activeEntry} /></div>
            </section>

            <aside className="run-column">
              <section className="panel run-panel">
                <div className="panel-header compact"><div><div className="eyebrow">Review</div><h2>Run command</h2></div><span className={`simulation ${runnerState}`}>{runnerState === "online" ? "Validated locally" : "Runner required"}</span></div>
                <div className="command-box"><div><span>Generated command</span><button type="button" onClick={copyCommand}>{copied ? "Copied" : "Copy"}</button></div><code>{command}</code></div>
                <div className="run-summary"><div><span>Writes to</span><strong>{tool === "visualizer" ? "Local browser session" : tool === "initialize" ? `entries/${String(values.initialize.name || "<entry>")}/` : `${entry}.sse.tsv`}</strong></div><div><span>Execution</span><strong>{tool === "visualizer" ? "Long-lived local service" : tool === "coordinates" ? "Model and cache dependent" : "Serialized for this entry"}</strong></div></div>
                {tool === "visualizer" && visualizerRunning ? (
                  <button className="primary-button danger" type="button" onClick={() => void cancelJob()} disabled={job?.status === "cancelling"}>⏻ {job?.status === "cancelling" ? "Shutting down explorer…" : "Shut down explorer"}</button>
                ) : (
                  <button className="primary-button" type="button" onClick={() => void startJob()} disabled={runnerState !== "online" || !environmentReady || submitting || Boolean(running) || (tool !== "initialize" && !activeEntry)}>▶ {submitting ? "Submitting…" : tool === "visualizer" ? "Start explorer" : `Run ${active.label.toLowerCase()}`}</button>
                )}
                <p className="simulation-note">Commands run without a shell in the configured Python environment. Writers for the same entry are queued to prevent data corruption.</p>
              </section>

              <section className={`panel job-panel ${job ? "visible" : ""}`} aria-live="polite">
                <div className="panel-header compact"><div><div className="eyebrow">Latest job</div><h2>{job ? tools.find((item) => item.id === job.tool)?.label || job.tool : "No active run"}</h2></div>{job ? <span className={`job-status ${job.status}`}>{jobStatusLabel(job)}</span> : null}</div>
                {job ? <>
                  <div className="progress-meta"><span>{job.ready ? "Explorer is accepting connections" : job.status === "succeeded" ? "Finished successfully" : job.status === "failed" ? `Exited with code ${job.exitCode ?? "?"}` : job.status === "queued" ? "Waiting for this entry" : job.status === "cancelled" ? "Stopped safely" : "Processing entry"}</span><strong>{job.progress === null ? "Live" : `${job.progress}%`}</strong></div>
                  <div className={`progress-track ${job.progress === null ? "indeterminate" : ""}`}><span style={job.progress === null ? undefined : { width: `${job.progress}%` }} /></div>
                  <div className="log-window">{job.logs.slice(-8).map((line, index) => <div className={line.stream} key={`${line.at}-${index}`}><span>{index === job.logs.slice(-8).length - 1 && ["running", "queued"].includes(job.status) ? "›" : line.stream === "stderr" ? "!" : "✓"}</span>{line.text}</div>)}</div>
                  <div className="job-actions">{["running", "queued", "cancelling"].includes(job.status) ? <button type="button" onClick={() => void cancelJob()} disabled={job.status === "cancelling"}>Cancel job</button> : <button type="button" onClick={() => setJob(null)}>Dismiss</button>}{job.ready && job.url ? <button type="button" onClick={() => window.open(job.url, "_blank", "noopener,noreferrer")}>Open explorer ↗</button> : <button type="button" onClick={() => setShowFullLog(true)} disabled={!job.logs.length}>View full log</button>}{job.tool === "boltz" && job.status === "succeeded" ? <button type="button" onClick={() => void revealFolder("structures")}>Open structures folder ↗</button> : null}</div>
                </> : <div className="empty-job"><span>○</span><p>Configure a tool and run it to see its real progress, logs, and exit status here.</p></div>}
              </section>

              <section className="panel guidance-panel"><div className="guidance-icon">i</div><div><strong>Workflow guidance</strong><p>{tool === "cluster" ? "Add taxonomy and external annotations before clustering if you want them included in enrichment tests. Clustering reads the embedding cache as it was written; it does not normalize." : tool === "coordinates" ? "This is the only tool that normalizes. Distances and clustering re-read the cache it writes, so the geometry chosen here is the geometry the whole entry gets." : tool === "distance" ? "Distances are Euclidean over the cache as it was written. On normalized vectors that is cosine geometry; on raw vectors it is not." : tool === "initialize" ? "Creation is bootstrap-only. Later changes should use the additive enrichment tools." : tool === "boltz" ? "Select points in the explorer and export them first. Prediction saves .cif structures under structures/ and appends pTM/pLDDT (and any RMSD) columns to the datafile; reload the explorer to see them." : "Every completed tool records its parameters in the entry manifest or run logs."}</p></div></section>
            </aside>
          </div>
        </section>
      </div>}

      {entryFiles ? <div className="modal-backdrop" role="presentation" onMouseDown={() => setEntryFiles(null)}><section className="modal-panel" role="dialog" aria-modal="true" aria-label={`Files for ${entry}`} onMouseDown={(event) => event.stopPropagation()}><div className="modal-header"><div><div className="eyebrow">Entry contents</div><h2>{entry}</h2></div><button type="button" onClick={() => setEntryFiles(null)} aria-label="Close">×</button></div><div className="file-list">{entryFiles.map((file) => <div key={file.path}><code>{file.path}</code><span>{file.size < 1024 ? `${file.size} B` : file.size < 1024 * 1024 ? `${(file.size / 1024).toFixed(1)} KB` : `${(file.size / 1024 / 1024).toFixed(1)} MB`}</span></div>)}</div></section></div> : null}
      {showFullLog && job ? <div className="modal-backdrop" role="presentation" onMouseDown={() => setShowFullLog(false)}><section className="modal-panel log-modal" role="dialog" aria-modal="true" aria-label="Full job log" onMouseDown={(event) => event.stopPropagation()}><div className="modal-header"><div><div className="eyebrow">{jobStatusLabel(job)}</div><h2>{job.command}</h2></div><button type="button" onClick={() => setShowFullLog(false)} aria-label="Close">×</button></div><div className="full-log">{job.logs.map((line, index) => <div className={line.stream} key={`${line.at}-${index}`}><time>{new Date(line.at).toLocaleTimeString()}</time><span>{line.text}</span></div>)}</div></section></div> : null}
      {shuttingDown ? <div className="shutdown-screen" role="status"><div className="brand-mark" aria-hidden="true"><span /><span /><span /></div><h1>Shutting down…</h1><p>Pipeline jobs and local services are being stopped safely. You can close this browser tab.</p></div> : null}
    </main>
  );
}
