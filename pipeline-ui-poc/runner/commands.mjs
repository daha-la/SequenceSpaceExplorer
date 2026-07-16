import path from "node:path";

export const TOOL_IDS = new Set([
  "initialize",
  "taxonomy",
  "external",
  "coordinates",
  "distance",
  "cluster",
  "boltz",
  "visualizer",
]);

export const WRITER_TOOLS = new Set([
  "initialize",
  "taxonomy",
  "external",
  "coordinates",
  "distance",
  "cluster",
  "boltz",
]);

const scripts = {
  initialize: "sse_initialization.py",
  taxonomy: "fetch_taxonomy.py",
  external: "merge_external.py",
  coordinates: "sse_coordinates.py",
  distance: "sse_esmc_distance.py",
  cluster: "sse_cluster.py",
  boltz: "sse_boltz.py",
  visualizer: "sse_visualizer.py",
};

function text(value, fallback = "") {
  return value === undefined || value === null ? fallback : String(value).trim();
}

function choice(value, allowed, fallback) {
  const candidate = text(value, fallback);
  if (!allowed.includes(candidate)) throw new Error(`Unsupported value: ${candidate}`);
  return candidate;
}

function integer(value, fallback, minimum = Number.MIN_SAFE_INTEGER, maximum = Number.MAX_SAFE_INTEGER) {
  const candidate = value === "" || value === undefined || value === null ? fallback : Number(value);
  if (!Number.isInteger(candidate) || candidate < minimum || candidate > maximum) {
    throw new Error(`Expected a whole number between ${minimum} and ${maximum}`);
  }
  return candidate;
}

function number(value, fallback, minimum = -Infinity, maximum = Infinity) {
  const candidate = value === "" || value === undefined || value === null ? fallback : Number(value);
  if (!Number.isFinite(candidate) || candidate < minimum || candidate > maximum) {
    throw new Error(`Expected a number between ${minimum} and ${maximum}`);
  }
  return candidate;
}

function required(value, label) {
  const candidate = text(value);
  if (!candidate) throw new Error(`${label} is required`);
  return candidate;
}

function add(args, flag, value) {
  if (value !== "" && value !== undefined && value !== null) args.push(flag, String(value));
}

function addFlag(args, flag, enabled) {
  if (enabled) args.push(flag);
}

export function validateEntryName(value) {
  const entry = required(value, "Entry");
  if (!/^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$/.test(entry) || entry === "." || entry === "..") {
    throw new Error("Entry names may contain letters, numbers, dots, underscores, and hyphens");
  }
  return entry;
}

export function buildJobSpec(payload, rootDir) {
  const tool = text(payload?.tool);
  if (!TOOL_IDS.has(tool)) throw new Error("Unknown pipeline tool");
  const values = payload?.values && typeof payload.values === "object" ? payload.values : {};
  const args = [path.join(rootDir, "scripts", scripts[tool])];
  let entry = tool === "initialize" ? validateEntryName(values.name || path.parse(text(values.input)).name) : validateEntryName(payload.entry);
  const secrets = {};

  if (tool === "initialize") {
    args.push(required(values.input, "Source file"));
    add(args, "--source", choice(values.source, ["em", "fasta", "fs"], "em"));
    add(args, "--name", entry);
    if (values.source === "em") {
      add(args, "--id_col", required(values.idCol, "ID column"));
      add(args, "--seq_col", required(values.seqCol, "Sequence column"));
    }
    const queries = text(values.query).split(/\s+/).filter(Boolean);
    if (queries.length) args.push("--query", ...queries);
    addFlag(args, "--force", values.force === true);
  }

  if (tool === "taxonomy") {
    args.push(entry);
    add(args, "--email", required(values.email, "NCBI contact email"));
    add(args, "--strategy", choice(values.strategy, ["auto", "em", "foldseek"], "auto"));
    add(args, "--batch", integer(values.batch, 100, 1, 10_000));
    add(args, "--gmgc-batch", integer(values.gmgcBatch, 50, 1, 10_000));
    if (values.rerun === "retry") args.push("--retry-failed");
    if (values.rerun === "force") args.push("--force");
    if (text(values.apiKey)) secrets.NCBI_API_KEY = text(values.apiKey);
  }

  if (tool === "external") {
    args.push(entry, required(values.file, "External data file"));
    add(args, "--id-col", text(values.idCol));
    add(args, "--columns", text(values.columns));
    add(args, "--translator", text(values.translator));
    add(args, "--type", choice(values.type, ["label", "coordinate"], "label"));
    if (values.delimiter !== "auto") add(args, "--delimiter", values.delimiter);
    addFlag(args, "--force", values.force === true);
  }

  if (tool === "coordinates") {
    args.push(entry);
    const embedder = choice(values.embedder, ["esmc", "prostt5", "saprot"], "esmc");
    const reducer = choice(values.reducer, ["pca", "umap", "tsne"], "pca");
    add(args, "--embedder", embedder);
    add(args, "--reducer", reducer);
    add(args, "--pooling", choice(values.pooling, ["mean", "max", "min"], "mean"));
    add(args, "--device", choice(values.device, ["auto", "cuda", "mps", "cpu"], "auto"));
    add(args, "--n-components", integer(values.components, 2, 1, 4096));
    args.push(values.normalize === false ? "--no-normalize" : "--normalize");
    if (embedder === "esmc") add(args, "--esmc-model", choice(values.esmcModel, ["esmc_300m", "esmc_600m"], "esmc_600m"));
    if (embedder === "prostt5") add(args, "--prostt5-checkpoint", required(values.prostCheckpoint, "ProstT5 checkpoint"));
    if (embedder === "saprot") add(args, "--saprot-checkpoint", required(values.saprotCheckpoint, "SaProt checkpoint"));
    if (reducer === "umap") {
      add(args, "--umap-neighbors", integer(values.neighbors, 15, 2, 100_000));
      add(args, "--umap-min-dist", number(values.minDist, 0.1, 0, 1));
      add(args, "--umap-metric", required(values.metric, "UMAP metric"));
    }
    if (reducer === "tsne") {
      add(args, "--tsne-perplexity", number(values.perplexity, 30, 1));
      add(args, "--tsne-pca", integer(values.tsnePca, 50, 0));
    }
    add(args, "--batch-size", integer(values.batchSize, 32, 1));
    add(args, "--write-every", integer(values.writeEvery, 1000, 1));
    if (embedder !== "esmc") {
      add(args, "--max-residues", integer(values.maxResidues, 1500, 3));
      add(args, "--foldseek-json", text(values.foldseekJson));
      addFlag(args, "--include_empty", values.includeEmpty === true);
    }
    add(args, "--label", text(values.label));
    if (values.rerun === "rereduce") args.push("--rereduce");
    if (values.rerun === "reembed") args.push("--reembed");
  }

  if (tool === "distance") {
    args.push(entry);
    add(args, "--embedding", required(values.embedding, "Embedding cache"));
    if (values.queryMode === "explicit") {
      const ids = text(values.queryIds).split(/[\s,]+/).filter(Boolean);
      if (!ids.length) throw new Error("At least one query ID is required");
      args.push("--query-id", ...ids);
    }
    addFlag(args, "--raw", values.raw === true);
    addFlag(args, "--force", values.force === true);
  }

  if (tool === "cluster") {
    args.push(entry);
    add(args, "--embedding", required(values.embedding, "Embedding cache"));
    const clusterer = choice(values.clusterer, ["kmeans", "hdbscan"], "kmeans");
    const space = choice(values.space, ["pca", "full"], "pca");
    add(args, "--clusterer", clusterer);
    add(args, "--space", space);
    if (space === "pca") {
      if (values.pcaMode === "variance") add(args, "--pca-variance", number(values.pcaVariance, 0.95, 0.01, 1));
      else add(args, "--pca-dims", integer(values.pcaDims, 50, 1));
    }
    if (clusterer === "kmeans") {
      if (values.kMode === "fixed") add(args, "--k", integer(values.k, 12, 2));
      else {
        const minimum = integer(values.kMin, 2, 2);
        const maximum = integer(values.kMax, 20, minimum);
        add(args, "--k-min", minimum);
        add(args, "--k-max", maximum);
      }
    } else {
      add(args, "--min-cluster-size", integer(values.minClusterSize, 50, 2));
      if (text(values.minSamples)) add(args, "--min-samples", integer(values.minSamples, 50, 1));
    }
    add(args, "--label", text(values.label));
    addFlag(args, "--raw", values.raw === true);
    addFlag(args, "--force", values.force === true);
    if (values.analysis === false) args.push("--no-analysis");
    else {
      add(args, "--analysis-top-n", integer(values.topN, 5, 1));
      add(args, "--fdr", number(values.fdr, 0.05, 0, 1));
    }
  }

  if (tool === "boltz") {
    args.push(entry);
    add(args, "--selection", text(values.selection));
    add(args, "--smiles", text(values.smiles));
    add(args, "--smiles-label", text(values.smilesLabel));
    if (values.useMsa === false) args.push("--no-msa");
    add(args, "--recycling-steps", integer(values.recyclingSteps, 3, 1, 10));
    add(args, "--sampling-steps", integer(values.samplingSteps, 200, 10, 500));
    add(args, "--diffusion-samples", integer(values.diffusionSamples, 5, 1, 10));
    add(args, "--step-scale", number(values.stepScale, 1.638, 0.1, 5));
    addFlag(args, "--force", values.force === true);
    if (values.rmsd === true) {
      args.push("--rmsd");
      add(args, "--rmsd-reference", required(values.rmsdReference, "RMSD reference id"));
      add(args, "--rmsd-ref-rank", integer(values.rmsdRefRank, 0, 0, 100));
      add(args, "--rmsd-method", choice(values.rmsdMethod, ["seq", "ce", "both"], "seq"));
      add(args, "--rmsd-scope", choice(values.rmsdScope, ["all", "selected"], "all"));
    }
    if (text(values.apiKey)) secrets.BOLTZ_API_KEY = text(values.apiKey);
  }

  if (tool === "visualizer") {
    args.push(entry);
    add(args, "--port", integer(values.port, 8051, 1024, 65535));
  }

  return {
    tool,
    entry,
    args,
    secrets,
    command: ["python", path.relative(rootDir, args[0]), ...args.slice(1)].map(shellQuote).join(" "),
    writesEntry: WRITER_TOOLS.has(tool),
  };
}

export function shellQuote(value) {
  const candidate = String(value);
  return /^[A-Za-z0-9_./:@%+=,-]+$/.test(candidate) ? candidate : `'${candidate.replaceAll("'", `'\\''`)}'`;
}
