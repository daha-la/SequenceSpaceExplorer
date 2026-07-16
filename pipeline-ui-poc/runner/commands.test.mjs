import assert from "node:assert/strict";
import test from "node:test";
import { buildJobSpec, validateEntryName } from "./commands.mjs";

const root = "/workspace/sse";

test("builds an argv array without invoking a shell", () => {
  const job = buildJobSpec({
    tool: "cluster",
    entry: "example_entry",
    values: { embedding: "esmc600m_mean", clusterer: "kmeans", space: "pca", pcaMode: "dims", pcaDims: 25, kMode: "fixed", k: 8, analysis: true, topN: 5, fdr: 0.05 },
  }, root);
  assert.equal(job.entry, "example_entry");
  assert.deepEqual(job.args.slice(0, 3), ["/workspace/sse/scripts/sse_cluster.py", "example_entry", "--embedding"]);
  assert.ok(job.args.includes("--k"));
  assert.equal(job.writesEntry, true);
});

test("keeps API keys out of argv and the displayed command", () => {
  const job = buildJobSpec({
    tool: "taxonomy",
    entry: "example_entry",
    values: { email: "researcher@example.org", apiKey: "secret-value", strategy: "auto", batch: 100, gmgcBatch: 50, rerun: "resume" },
  }, root);
  assert.equal(job.secrets.NCBI_API_KEY, "secret-value");
  assert.doesNotMatch(job.command, /secret-value/);
  assert.ok(!job.args.includes("secret-value"));
});

test("builds a boltz spec with selection, params, RMSD, and a secret key", () => {
  const job = buildJobSpec({
    tool: "boltz",
    entry: "example_entry",
    values: {
      selection: "selection_20260715_101500.json",
      apiKey: "nvapi-secret",
      smiles: "OC1CCCCC1",
      smilesLabel: "UDP-Glc",
      useMsa: true,
      recyclingSteps: 3,
      samplingSteps: 200,
      diffusionSamples: 5,
      stepScale: 1.638,
      force: false,
      rmsd: true,
      rmsdReference: "OleD_S1",
      rmsdRefRank: 0,
      rmsdMethod: "both",
      rmsdScope: "selected",
    },
  }, root);
  assert.equal(job.entry, "example_entry");
  assert.deepEqual(job.args.slice(0, 2), ["/workspace/sse/scripts/sse_boltz.py", "example_entry"]);
  assert.ok(job.args.includes("--selection"));
  assert.ok(job.args.includes("selection_20260715_101500.json"));
  assert.ok(job.args.includes("--rmsd"));
  assert.deepEqual(job.args.slice(job.args.indexOf("--rmsd-reference"), job.args.indexOf("--rmsd-reference") + 2), ["--rmsd-reference", "OleD_S1"]);
  assert.equal(job.writesEntry, true);
  assert.equal(job.secrets.BOLTZ_API_KEY, "nvapi-secret");
  assert.doesNotMatch(job.command, /nvapi-secret/);
  assert.ok(!job.args.includes("nvapi-secret"));
});

test("boltz omits --no-msa when MSA is on and drops RMSD flags when disabled", () => {
  const job = buildJobSpec({
    tool: "boltz",
    entry: "example_entry",
    values: { useMsa: true, rmsd: false },
  }, root);
  assert.ok(!job.args.includes("--no-msa"));
  assert.ok(!job.args.includes("--rmsd"));
});

test("rejects traversal and unknown tools", () => {
  assert.throws(() => validateEntryName("../outside"), /Entry names/);
  assert.throws(() => buildJobSpec({ tool: "shell", values: {} }, root), /Unknown pipeline tool/);
});
