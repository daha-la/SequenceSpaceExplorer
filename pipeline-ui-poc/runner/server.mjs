import { createServer } from "node:http";
import { createWriteStream, existsSync } from "node:fs";
import { mkdir, readFile, readdir, rename, stat, writeFile } from "node:fs/promises";
import { randomUUID } from "node:crypto";
import path from "node:path";
import { spawn, spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";
import { buildJobSpec, validateEntryName } from "./commands.mjs";

const runnerDir = path.dirname(fileURLToPath(import.meta.url));
const appDir = path.dirname(runnerDir);
const rootDir = path.dirname(appDir);
const entriesDir = path.join(rootDir, "entries");
const workDir = path.join(appDir, "work");
const uploadsDir = path.join(workDir, "uploads");
const jobsPath = path.join(workDir, "jobs.json");
const requestedShutdownFile = process.env.SSE_SHUTDOWN_FILE ? path.resolve(process.env.SSE_SHUTDOWN_FILE) : null;
const shutdownFile = requestedShutdownFile && path.dirname(requestedShutdownFile) === workDir ? requestedShutdownFile : null;
const port = Number(process.env.SSE_RUNNER_PORT || 8788);
const host = process.env.SSE_RUNNER_HOST || "127.0.0.1";
const allowedOrigins = new Set((process.env.SSE_UI_ORIGINS || "http://localhost:3000,http://127.0.0.1:3000,http://localhost:5173,http://127.0.0.1:5173").split(","));
const active = new Map();
const jobSecrets = new Map();
let jobs = [];
let saveChain = Promise.resolve();
let stopping = false;

await mkdir(uploadsDir, { recursive: true });
try {
  jobs = JSON.parse(await readFile(jobsPath, "utf8"));
  jobs = jobs.map((job) => ["queued", "running"].includes(job.status)
    ? { ...job, status: "interrupted", finishedAt: new Date().toISOString(), error: "The local runner restarted before this job finished." }
    : job);
} catch {
  jobs = [];
}
await persistJobs();

const environment = probeEnvironment();

function corsHeaders(req) {
  const origin = req.headers.origin;
  const headers = {
    "Access-Control-Allow-Headers": "Content-Type, X-File-Name",
    "Access-Control-Allow-Methods": "GET, POST, DELETE, OPTIONS",
    "Cache-Control": "no-store",
  };
  if (origin && allowedOrigins.has(origin)) headers["Access-Control-Allow-Origin"] = origin;
  return headers;
}

function originAllowed(req) {
  return !req.headers.origin || allowedOrigins.has(req.headers.origin);
}

function send(req, res, statusCode, body) {
  res.writeHead(statusCode, { "Content-Type": "application/json; charset=utf-8", ...corsHeaders(req) });
  res.end(JSON.stringify(body));
}

// The visualizer (Dash/Flask) and some tools spawn their own child processes, so
// a bare child.kill() can leave those grandchildren orphaned and holding a port -
// which is exactly what forced a hard restart between visualizer launches. Jobs
// are spawned in their own process group (detached) so we can signal the whole
// tree with a negative pid.
function terminateChild(child, signal) {
  if (!child) return;
  try {
    if (child.pid && process.platform !== "win32") process.kill(-child.pid, signal);
    else child.kill(signal);
  } catch {
    try { child.kill(signal); } catch { /* already exited */ }
  }
}

function stopRunner() {
  if (stopping) return;
  stopping = true;
  for (const child of active.values()) terminateChild(child, "SIGTERM");
  server.close(() => process.exit(0));
  setTimeout(() => process.exit(0), 5000).unref();
}

async function readJson(req, limit = 2_000_000) {
  const chunks = [];
  let size = 0;
  for await (const chunk of req) {
    size += chunk.length;
    if (size > limit) throw new Error("Request is too large");
    chunks.push(chunk);
  }
  return JSON.parse(Buffer.concat(chunks).toString("utf8") || "{}");
}

function persistJobs() {
  saveChain = saveChain.then(async () => {
    const temporary = `${jobsPath}.tmp`;
    await writeFile(temporary, JSON.stringify(jobs.slice(0, 250), null, 2));
    await rename(temporary, jobsPath);
  });
  return saveChain;
}

function publicJob(job) {
  const safe = { ...job };
  delete safe.secrets;
  delete safe.args;
  return safe;
}

async function inspectEntry(name) {
  const directory = path.join(entriesDir, name);
  const files = await readdir(directory, { withFileTypes: true });
  const datafile = files.find((item) => item.isFile() && item.name.endsWith(".sse.tsv"));
  if (!datafile) return null;
  const logsDir = path.join(directory, "logs");
  let manifest = {};
  try {
    const manifests = (await readdir(logsDir)).filter((file) => file.endsWith(".manifest.json"));
    if (manifests[0]) manifest = JSON.parse(await readFile(path.join(logsDir, manifests[0]), "utf8"));
  } catch { /* an entry can predate manifests */ }
  let columns = Array.isArray(manifest.columns) ? manifest.columns : [];
  if (!columns.length) {
    const header = (await readFile(path.join(directory, datafile.name), "utf8")).split(/\r?\n/, 2);
    const names = header[0]?.split("\t") || [];
    const types = header[1]?.split("\t") || [];
    columns = names.map((column, index) => ({ name: column, type: types[index] }));
  }
  const coordinateColumns = columns.filter((column) => column.type === "coordinate").length;
  const annotationColumns = columns.filter((column) => column.type === "label").length;
  const taxonomyColumns = columns.filter((column) => column.tool === "fetch_taxonomy" || String(column.name).startsWith("tax_")).length;
  let analyses = 0;
  try { analyses = (await readdir(path.join(directory, "cluster_analysis"), { withFileTypes: true })).filter((item) => item.isDirectory()).length; } catch { /* optional */ }
  const embeddings = [];
  try {
    const raw = (await readdir(path.join(directory, "embeddings"))).filter((file) => file.endsWith(".emb.tsv"));
    let normalized = [];
    try { normalized = await readdir(path.join(directory, "embeddings", "normalized")); } catch { /* optional */ }
    for (const file of raw) embeddings.push({ tag: file.replace(/\.emb\.tsv$/, ""), normalized: normalized.includes(file) });
  } catch { /* optional */ }
  const selections = [];
  try {
    const selectionsDir = path.join(directory, "selections");
    const selFiles = (await readdir(selectionsDir))
      .filter((file) => file.startsWith("selection_") && file.endsWith(".json"));
    for (const file of selFiles) {
      let count = null;
      let created = null;
      let ids = [];
      try {
        const payload = JSON.parse(await readFile(path.join(selectionsDir, file), "utf8"));
        count = payload.count ?? (Array.isArray(payload.sequences) ? payload.sequences.length : null);
        created = payload.created_utc ?? null;
        // Sequence ids let the UI offer an RMSD reference picker. Cap the list so
        // a very large selection does not bloat the entries payload.
        if (Array.isArray(payload.sequences)) ids = payload.sequences.slice(0, 2000).map((s) => String(s.id));
      } catch { /* a malformed cache still lists by name */ }
      selections.push({ name: file, count, created, ids });
    }
    // Newest first: created_utc when present, else filename (timestamped) order.
    selections.sort((a, b) => String(b.created || b.name).localeCompare(String(a.created || a.name)));
  } catch { /* optional */ }
  const fileStats = await stat(path.join(directory, datafile.name));
  return {
    name,
    source: manifest.source_type || "unknown",
    sequences: manifest.row_counts?.kept ?? null,
    queries: manifest.row_counts?.queries ?? null,
    annotations: annotationColumns,
    taxonomyFields: taxonomyColumns,
    coordinates: coordinateColumns,
    analyses,
    embeddings,
    selections,
    updatedAt: fileStats.mtime.toISOString(),
  };
}

async function listEntries() {
  let directories = [];
  try { directories = await readdir(entriesDir, { withFileTypes: true }); } catch { return []; }
  const inspected = await Promise.all(directories.filter((item) => item.isDirectory() && !item.name.startsWith(".")).map((item) => inspectEntry(item.name)));
  return inspected.filter(Boolean).sort((a, b) => b.updatedAt.localeCompare(a.updatedAt));
}

async function listFiles(entry) {
  const name = validateEntryName(entry);
  const base = path.join(entriesDir, name);
  const output = [];
  async function walk(directory) {
    for (const item of await readdir(directory, { withFileTypes: true })) {
      if (item.name.startsWith(".")) continue;
      const absolute = path.join(directory, item.name);
      if (item.isDirectory()) await walk(absolute);
      else {
        const details = await stat(absolute);
        output.push({ path: path.relative(base, absolute), size: details.size, modifiedAt: details.mtime.toISOString() });
      }
      if (output.length >= 500) return;
    }
  }
  await walk(base);
  return output.sort((a, b) => a.path.localeCompare(b.path));
}

// Open a folder inside an entry in the OS file browser. This only makes sense on
// the machine that owns the datasets (the local runner), which is exactly where
// this process runs. The subpath is confined to the entry directory.
async function revealPath(entry, subpath) {
  const name = validateEntryName(entry);
  const base = path.join(entriesDir, name);
  const requested = subpath ? String(subpath) : "";
  const target = path.resolve(base, requested);
  const within = target === base || target.startsWith(base + path.sep);
  if (!within) throw new Error("Path escapes the entry directory");
  await stat(target); // throws if the folder does not exist yet
  const opener = process.platform === "darwin" ? "open"
    : process.platform === "win32" ? "explorer"
    : "xdg-open";
  const child = spawn(opener, [target], { detached: true, stdio: "ignore" });
  child.on("error", () => { /* surfaced to the caller via the missing-path check only */ });
  child.unref();
  return path.relative(base, target) || ".";
}

function appendLog(job, stream, text) {
  const lines = text.toString().split(/\r?\n/).filter(Boolean);
  for (const line of lines) {
    job.logs.push({ at: new Date().toISOString(), stream, text: line });
    if (job.logs.length > 2000) job.logs.shift();
    const match = line.match(/https?:\/\/(?:127\.0\.0\.1|localhost):\d+/);
    if (job.tool === "visualizer" && match) {
      job.url = match[0];
      job.ready = true;
      job.progress = 100;
    }
  }
  void persistJobs();
}

function pythonExecutable() {
  if (process.env.SSE_PYTHON) return process.env.SSE_PYTHON;
  for (const candidate of [path.join(rootDir, ".venv", "bin", "python"), path.join(rootDir, "venv", "bin", "python")]) {
    if (existsSync(candidate)) return candidate;
  }
  return "python3";
}

function probeEnvironment() {
  const modules = ["numpy", "pandas", "Bio", "scipy", "sklearn", "matplotlib", "umap", "dash", "plotly", "requests", "tqdm"];
  const script = `import importlib.util,json\nmods=${JSON.stringify(modules)}\nprint(json.dumps({'missing':[m for m in mods if importlib.util.find_spec(m) is None]}))`;
  const checked = spawnSync(pythonExecutable(), ["-c", script], { cwd: rootDir, encoding: "utf8", timeout: 10_000 });
  if (checked.error || checked.status !== 0) return { ready: false, missing: modules, error: checked.error?.message || checked.stderr?.trim() || "Python probe failed" };
  try {
    const result = JSON.parse(checked.stdout.trim());
    return { ready: result.missing.length === 0, missing: result.missing };
  } catch {
    return { ready: false, missing: modules, error: "Python probe returned an invalid response" };
  }
}

function canStart(job) {
  if (!job.writesEntry) return true;
  return !jobs.some((candidate) => candidate.id !== job.id && candidate.entry === job.entry && candidate.writesEntry && candidate.status === "running");
}

function pumpQueue() {
  for (const job of jobs.slice().reverse()) {
    if (job.status === "queued" && canStart(job)) runJob(job);
  }
}

function runJob(job) {
  job.status = "running";
  job.startedAt = new Date().toISOString();
  job.progress = null;
  const child = spawn(pythonExecutable(), job.args, {
    cwd: rootDir,
    env: { ...process.env, ...(jobSecrets.get(job.id) || {}), PYTHONUNBUFFERED: "1" },
    stdio: ["ignore", "pipe", "pipe"],
    detached: process.platform !== "win32",
  });
  active.set(job.id, child);
  jobSecrets.delete(job.id);
  job.pid = child.pid;
  appendLog(job, "system", `Started ${job.command}`);
  child.stdout.on("data", (chunk) => appendLog(job, "stdout", chunk));
  child.stderr.on("data", (chunk) => appendLog(job, "stderr", chunk));
  child.on("error", (error) => appendLog(job, "system", error.message));
  child.on("close", (code, signal) => {
    active.delete(job.id);
    job.exitCode = code;
    job.finishedAt = new Date().toISOString();
    job.progress = code === 0 ? 100 : job.progress;
    if (job.status === "cancelling") job.status = "cancelled";
    else job.status = code === 0 ? "succeeded" : "failed";
    appendLog(job, "system", signal ? `Process stopped by ${signal}.` : `Process exited with code ${code}.`);
    void persistJobs();
    pumpQueue();
  });
  void persistJobs();
}

async function upload(req, res) {
  const original = path.basename(String(req.headers["x-file-name"] || "upload.dat"));
  const safe = original.replace(/[^A-Za-z0-9._-]+/g, "_").slice(-180) || "upload.dat";
  const relative = path.join("pipeline-ui-poc", "work", "uploads", `${Date.now()}-${randomUUID().slice(0, 8)}-${safe}`);
  const absolute = path.join(rootDir, relative);
  const output = createWriteStream(absolute, { flags: "wx" });
  let size = 0;
  for await (const chunk of req) {
    size += chunk.length;
    if (size > 1_000_000_000) {
      output.destroy(new Error("Upload exceeds 1 GB"));
      throw new Error("Upload exceeds 1 GB");
    }
    if (!output.write(chunk)) await new Promise((resolve) => output.once("drain", resolve));
  }
  await new Promise((resolve, reject) => output.end((error) => error ? reject(error) : resolve()));
  send(req, res, 201, { path: relative, name: original, size });
}

const server = createServer(async (req, res) => {
  try {
    if (req.method === "OPTIONS") {
      if (!originAllowed(req)) return send(req, res, 403, { error: "Origin not allowed" });
      res.writeHead(204, corsHeaders(req));
      return res.end();
    }
    if (!originAllowed(req)) return send(req, res, 403, { error: "Origin not allowed" });
    const url = new URL(req.url || "/", `http://${req.headers.host || `${host}:${port}`}`);

    if (req.method === "GET" && url.pathname === "/api/health") {
      return send(req, res, 200, { ok: true, root: rootDir, python: pythonExecutable(), activeJobs: active.size, environment, managedShutdown: Boolean(shutdownFile) });
    }
    if (req.method === "GET" && url.pathname === "/api/entries") return send(req, res, 200, { entries: await listEntries() });
    const filesMatch = url.pathname.match(/^\/api\/entries\/([^/]+)\/files$/);
    if (req.method === "GET" && filesMatch) return send(req, res, 200, { files: await listFiles(decodeURIComponent(filesMatch[1])) });
    const revealMatch = url.pathname.match(/^\/api\/entries\/([^/]+)\/reveal$/);
    if (req.method === "POST" && revealMatch) {
      const body = await readJson(req).catch(() => ({}));
      try {
        const revealed = await revealPath(decodeURIComponent(revealMatch[1]), body?.subpath);
        return send(req, res, 200, { ok: true, path: revealed });
      } catch (reason) {
        return send(req, res, 400, { error: reason instanceof Error ? reason.message : "Could not open folder" });
      }
    }
    if (req.method === "GET" && url.pathname === "/api/jobs") {
      const filtered = url.searchParams.get("entry") ? jobs.filter((job) => job.entry === url.searchParams.get("entry")) : jobs;
      return send(req, res, 200, { jobs: filtered.slice(0, 50).map(publicJob) });
    }
    const jobMatch = url.pathname.match(/^\/api\/jobs\/([A-Za-z0-9-]+)$/);
    if (req.method === "GET" && jobMatch) {
      const job = jobs.find((candidate) => candidate.id === jobMatch[1]);
      return job ? send(req, res, 200, { job: publicJob(job) }) : send(req, res, 404, { error: "Job not found" });
    }
    if (req.method === "POST" && url.pathname === "/api/jobs") {
      const payload = await readJson(req);
      const spec = buildJobSpec(payload, rootDir);
      const now = new Date().toISOString();
      const id = randomUUID();
      const { secrets, ...persistedSpec } = spec;
      if (Object.keys(secrets).length) jobSecrets.set(id, secrets);
      const job = { id, ...persistedSpec, values: { ...payload.values, apiKey: undefined }, status: "queued", progress: 0, logs: [], createdAt: now };
      jobs.unshift(job);
      await persistJobs();
      pumpQueue();
      return send(req, res, 202, { job: publicJob(job) });
    }
    if (req.method === "DELETE" && jobMatch) {
      const job = jobs.find((candidate) => candidate.id === jobMatch[1]);
      if (!job) return send(req, res, 404, { error: "Job not found" });
      if (job.status === "queued") {
        job.status = "cancelled";
        job.finishedAt = new Date().toISOString();
      } else if (job.status === "running") {
        job.status = "cancelling";
        terminateChild(active.get(job.id), "SIGTERM");
        setTimeout(() => terminateChild(active.get(job.id), "SIGKILL"), 5000).unref();
      }
      await persistJobs();
      return send(req, res, 200, { job: publicJob(job) });
    }
    if (req.method === "POST" && url.pathname === "/api/uploads") return await upload(req, res);
    if (req.method === "POST" && url.pathname === "/api/shutdown") {
      if (shutdownFile) await writeFile(shutdownFile, `${new Date().toISOString()}\n`, { flag: "w" });
      send(req, res, 202, { ok: true, managed: Boolean(shutdownFile) });
      if (!shutdownFile) setTimeout(stopRunner, 100).unref();
      return;
    }
    return send(req, res, 404, { error: "Not found" });
  } catch (error) {
    return send(req, res, 400, { error: error instanceof Error ? error.message : "Request failed" });
  }
});

server.listen(port, host, () => {
  console.log(`Sequence Space Explorer runner listening at http://${host}:${port}`);
  console.log(`Pipeline root: ${rootDir}`);
});

for (const signal of ["SIGINT", "SIGTERM"]) {
  process.on(signal, stopRunner);
}
