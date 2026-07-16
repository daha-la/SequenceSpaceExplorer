import { a as require_react, o as __toESM, t as require_jsx_runtime } from "../index.js";
//#region app/runner-client.ts
var import_react = /* @__PURE__ */ __toESM(require_react(), 1);
var runnerUrl = (process.env.NEXT_PUBLIC_SSE_RUNNER_URL || "http://127.0.0.1:8788").replace(/\/$/, "");
async function request(path, init) {
	const response = await fetch(`${runnerUrl}${path}`, init);
	const payload = await response.json().catch(() => ({}));
	if (!response.ok) throw new Error(payload.error || `Runner request failed (${response.status})`);
	return payload;
}
async function getRunnerState() {
	const [health, entries] = await Promise.all([request("/api/health"), request("/api/entries")]);
	return {
		health,
		entries: entries.entries
	};
}
async function getJobs(entry) {
	return (await request(`/api/jobs${entry ? `?entry=${encodeURIComponent(entry)}` : ""}`)).jobs;
}
async function getJob(id) {
	return (await request(`/api/jobs/${encodeURIComponent(id)}`)).job;
}
async function submitJob(payload) {
	return (await request("/api/jobs", {
		method: "POST",
		headers: { "Content-Type": "application/json" },
		body: JSON.stringify(payload)
	})).job;
}
async function cancelPipelineJob(id) {
	return (await request(`/api/jobs/${encodeURIComponent(id)}`, { method: "DELETE" })).job;
}
async function requestPipelineShutdown() {
	return await request("/api/shutdown", { method: "POST" });
}
async function getEntryFiles(entry) {
	return (await request(`/api/entries/${encodeURIComponent(entry)}/files`)).files;
}
async function revealEntryPath(entry, subpath) {
	return await request(`/api/entries/${encodeURIComponent(entry)}/reveal`, {
		method: "POST",
		headers: { "Content-Type": "application/json" },
		body: JSON.stringify({ subpath })
	});
}
async function uploadRunnerFile(file, onProgress) {
	return await new Promise((resolve, reject) => {
		const xhr = new XMLHttpRequest();
		xhr.open("POST", `${runnerUrl}/api/uploads`);
		xhr.setRequestHeader("Content-Type", "application/octet-stream");
		xhr.setRequestHeader("X-File-Name", file.name);
		xhr.upload.onprogress = (event) => {
			if (event.lengthComputable) onProgress?.(Math.round(event.loaded / event.total * 100));
		};
		xhr.onerror = () => reject(/* @__PURE__ */ new Error("Could not reach the local runner"));
		xhr.onload = () => {
			let payload = {};
			try {
				payload = JSON.parse(xhr.responseText);
			} catch {}
			if (xhr.status < 200 || xhr.status >= 300 || !payload.path) return reject(new Error(payload.error || "Upload failed"));
			resolve({
				path: payload.path,
				name: payload.name || file.name,
				size: payload.size || file.size
			});
		};
		xhr.send(file);
	});
}
//#endregion
//#region app/page.tsx
var import_jsx_runtime = require_jsx_runtime();
var tools = [
	{
		id: "initialize",
		label: "Create entry",
		script: "sse_initialization.py",
		stage: "1 · Entry",
		description: "Start from TSV, FASTA, or Foldseek JSON",
		mark: "01"
	},
	{
		id: "taxonomy",
		label: "Fetch taxonomy",
		script: "fetch_taxonomy.py",
		stage: "2 · Enrich",
		description: "Resolve NCBI taxIds and lineages",
		mark: "02"
	},
	{
		id: "external",
		label: "Merge external data",
		script: "merge_external.py",
		stage: "2 · Enrich",
		description: "Join assay, annotation, or coordinate columns",
		mark: "03"
	},
	{
		id: "coordinates",
		label: "Build coordinates",
		script: "sse_coordinates.py",
		stage: "3 · Embed",
		description: "Embed sequences and reduce their dimensions",
		mark: "04"
	},
	{
		id: "distance",
		label: "Query distances",
		script: "sse_esmc_distance.py",
		stage: "4 · Analyze",
		description: "Measure embedding distance to references",
		mark: "05"
	},
	{
		id: "cluster",
		label: "Cluster space",
		script: "sse_cluster.py",
		stage: "4 · Analyze",
		description: "Run k-means or HDBSCAN with diagnostics",
		mark: "06"
	},
	{
		id: "boltz",
		label: "Structure & binding",
		script: "sse_boltz.py",
		stage: "4 · Analyze",
		description: "Boltz-2 structure prediction + RMSD for a selection",
		mark: "07"
	},
	{
		id: "visualizer",
		label: "Open explorer",
		script: "sse_visualizer.py",
		stage: "5 · Explore",
		description: "Launch the interactive Dash visualization",
		mark: "08"
	}
];
var initialValues = {
	initialize: {
		input: "my_enzymes.tsv",
		source: "em",
		name: "my_enzymes",
		idCol: "Accession",
		seqCol: "Sequence",
		query: "",
		force: false
	},
	taxonomy: {
		email: "researcher@example.org",
		strategy: "auto",
		apiKey: "",
		batch: 100,
		gmgcBatch: 50,
		rerun: "resume"
	},
	external: {
		file: "measurements.csv",
		idCol: "",
		columns: "pI,Melting_temperature",
		translator: "",
		type: "label",
		delimiter: "auto",
		force: false
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
		minDist: .1,
		metric: "euclidean",
		perplexity: 30,
		tsnePca: 50,
		batchSize: 32,
		writeEvery: 1e3,
		maxResidues: 1500,
		foldseekJson: "",
		label: "",
		includeEmpty: false,
		rerun: "new"
	},
	distance: {
		embedding: "esmc600m_mean",
		queryMode: "marked",
		queryIds: "OleD_S1, AgepGT_S3",
		raw: false,
		force: false
	},
	cluster: {
		embedding: "esmc600m_mean",
		raw: false,
		space: "pca",
		pcaMode: "dims",
		pcaDims: 50,
		pcaVariance: .95,
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
		fdr: .05
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
		rmsdScope: "all"
	},
	visualizer: { port: 8051 }
};
var landingPoints = [
	{
		color: "#c7d2d0",
		lobes: [
			[
				49,
				25,
				8,
				5,
				14
			],
			[
				43,
				39,
				5,
				8,
				13
			],
			[
				56,
				42,
				7,
				4,
				12
			],
			[
				48,
				53,
				4,
				8,
				10
			],
			[
				35,
				45,
				7,
				4,
				10
			],
			[
				60,
				26,
				5,
				3,
				8
			],
			[
				39,
				17,
				3,
				3,
				5
			],
			[
				34,
				61,
				3,
				4,
				6
			],
			[
				62,
				53,
				3,
				6,
				6
			]
		]
	},
	{
		color: "#96b5e7",
		lobes: [
			[
				58,
				71,
				14,
				12,
				38
			],
			[
				70,
				67,
				12,
				9,
				28
			],
			[
				79,
				73,
				11,
				7,
				23
			],
			[
				63,
				84,
				8,
				5,
				15
			],
			[
				52,
				64,
				6,
				7,
				13
			],
			[
				84,
				57,
				8,
				5,
				15
			]
		]
	},
	{
		color: "#ff169b",
		lobes: [
			[
				76,
				25,
				10,
				6,
				24
			],
			[
				84,
				32,
				8,
				6,
				20
			],
			[
				88,
				44,
				4,
				8,
				15
			],
			[
				73,
				40,
				5,
				4,
				12
			],
			[
				91,
				23,
				3,
				4,
				8
			],
			[
				69,
				18,
				5,
				4,
				10
			]
		]
	},
	{
		color: "#20dfca",
		lobes: [
			[
				22,
				49,
				6,
				7,
				20
			],
			[
				17,
				42,
				3,
				3,
				8
			],
			[
				26,
				57,
				4,
				3,
				10
			]
		]
	},
	{
		color: "#ae36c2",
		lobes: [
			[
				35,
				76,
				5,
				4,
				17
			],
			[
				42,
				85,
				4,
				4,
				11
			],
			[
				38,
				91,
				3,
				3,
				7
			]
		]
	},
	{
		color: "#7560ff",
		lobes: [[
			31,
			40,
			4,
			5,
			18
		], [
			28,
			34,
			3,
			3,
			9
		]]
	},
	{
		color: "#a9bd4a",
		lobes: [[
			67,
			47,
			5,
			8,
			20
		], [
			72,
			53,
			4,
			5,
			12
		]]
	},
	{
		color: "#3abfff",
		lobes: [[
			52,
			56,
			4,
			3,
			14
		]]
	},
	{
		color: "#ff4c44",
		lobes: [[
			39,
			59,
			3,
			3,
			10
		]]
	},
	{
		color: "#ff875f",
		lobes: [[
			25,
			68,
			3,
			2,
			8
		]]
	},
	{
		color: "#e75bff",
		lobes: [[
			35,
			28,
			3,
			3,
			8
		]]
	},
	{
		color: "#f3a6b7",
		lobes: [[
			51,
			20,
			4,
			3,
			10
		]]
	}
].flatMap((shape, cluster) => shape.lobes.flatMap(([cx, cy, rx, ry, count], lobe) => {
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
			opacity: cluster === 0 ? .88 + index % 4 * .035 : .94 + index % 3 * .025
		};
	});
}));
function flag(value, name) {
	return value ? ` ${name}` : "";
}
function quoted(value) {
	return String(value).includes(" ") ? `"${value}"` : String(value);
}
function embeddingTag(v) {
	if (v.label) return String(v.label);
	if (v.embedder === "esmc") return `esmc${String(v.esmcModel).replace("esmc_", "")}_${v.pooling}`;
	return `${v.embedder}_${v.pooling}`;
}
function buildCommand(tool, v, entry) {
	const root = "python scripts/";
	if (tool === "initialize") {
		let cmd = `${root}sse_initialization.py ${quoted(v.input)} --source ${v.source}`;
		if (v.name) cmd += ` --name ${quoted(v.name)}`;
		if (v.source === "em") cmd += ` --id_col ${quoted(v.idCol)} --seq_col ${quoted(v.seqCol)}`;
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
		if (v.reducer === "umap") cmd += ` --umap-neighbors ${v.neighbors} --umap-min-dist ${v.minDist} --umap-metric ${v.metric}`;
		if (v.reducer === "tsne") cmd += ` --tsne-perplexity ${v.perplexity} --tsne-pca ${v.tsnePca}`;
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
		if (v.queryMode === "explicit" && v.queryIds) cmd += ` --query-id ${String(v.queryIds).replaceAll(",", " ")}`;
		return cmd + flag(v.raw, "--raw") + flag(v.force, "--force");
	}
	if (tool === "cluster") {
		let cmd = `${root}sse_cluster.py ${entry} --embedding ${quoted(v.embedding)} --clusterer ${v.clusterer}`;
		cmd += ` --space ${v.space}`;
		if (v.space === "pca") cmd += v.pcaMode === "variance" ? ` --pca-variance ${v.pcaVariance}` : ` --pca-dims ${v.pcaDims}`;
		if (v.clusterer === "kmeans") cmd += v.kMode === "fixed" ? ` --k ${v.k}` : ` --k-min ${v.kMin} --k-max ${v.kMax}`;
		else {
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
function Field({ label, hint, children }) {
	return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("label", {
		className: "field",
		children: [
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
				className: "field-label",
				children: label
			}),
			children,
			hint ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
				className: "field-hint",
				children: hint
			}) : null
		]
	});
}
function FilePathInput({ value, onChange, accept, runnerOnline }) {
	const [uploading, setUploading] = (0, import_react.useState)(null);
	const [error, setError] = (0, import_react.useState)("");
	const picker = (0, import_react.useRef)(null);
	const chooseFile = async (file) => {
		if (!file) return;
		setError("");
		setUploading(0);
		try {
			onChange((await uploadRunnerFile(file, setUploading)).path);
		} catch (reason) {
			setError(reason instanceof Error ? reason.message : "Upload failed");
		} finally {
			setUploading(null);
		}
	};
	return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("span", {
		className: "file-path-input",
		children: [
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
				value,
				onChange: (event) => onChange(event.target.value)
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("span", {
				className: "file-upload-row",
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: uploading === null ? "Use a path available to the runner, or upload a file." : `Uploading… ${uploading}%` }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
						type: "button",
						disabled: !runnerOnline || uploading !== null,
						onClick: () => picker.current?.click(),
						children: "Browse"
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
						ref: picker,
						className: "visually-hidden",
						type: "file",
						accept,
						disabled: !runnerOnline || uploading !== null,
						onChange: (event) => void chooseFile(event.target.files?.[0])
					})
				]
			}),
			error ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
				className: "field-error",
				children: error
			}) : null
		]
	});
}
function Toggle({ checked, onChange, label, description }) {
	return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("label", {
		className: "toggle-row",
		children: [
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
				type: "checkbox",
				checked,
				onChange: (e) => onChange(e.target.checked)
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
				className: "toggle-control",
				"aria-hidden": "true",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {})
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("span", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: label }), description ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("small", { children: description }) : null] })
		]
	});
}
function GeometrySection({ values, setValue, entry }) {
	const normalize = Boolean(values.normalize);
	const tag = embeddingTag(values);
	const caches = entry?.embeddings ?? [];
	const stranded = !normalize && caches.some((item) => item.tag === tag && item.normalized);
	return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
		className: "subsection geometry",
		children: [
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
				className: "subsection-title",
				children: "Vector geometry"
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("p", {
				className: "geometry-lede",
				children: [
					"The pipeline chooses a geometry ",
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("b", { children: "once, here" }),
					". This run writes the vectors it reduced, and ",
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("b", { children: "Query distances" }),
					" and ",
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("b", { children: "Cluster space" }),
					" read those same vectors back from disk. Neither of them normalizes anything on its own, so whatever you pick below is the geometry every later analysis inherits."
				]
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
				checked: normalize,
				onChange: (value) => setValue("normalize", value),
				label: "L2-normalize embeddings",
				description: normalize ? "Vectors are scaled to unit length, so distance and clustering are explicitly cosine geometry." : "Coordinates, distances, and clustering all fall back to raw, unscaled model output."
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("dl", {
				className: "geometry-contract",
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("dt", { children: "Reduces" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dd", { children: normalize ? "Unit-length vectors" : "Raw model output, unscaled" })] }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("dt", { children: "Writes" }), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("dd", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("code", { children: [
						"embeddings/",
						tag,
						".emb.tsv"
					] }), normalize ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(import_jsx_runtime.Fragment, { children: [" plus its normalized sibling ", /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("code", { children: [
						"embeddings/normalized/",
						tag,
						".emb.tsv"
					] })] }) : /* @__PURE__ */ (0, import_jsx_runtime.jsx)(import_jsx_runtime.Fragment, { children: " only. No normalized sibling is written." })] })] }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("dt", { children: "Downstream" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("dd", { children: normalize ? "Query distances and Cluster space prefer the normalized sibling automatically." : "Query distances and Cluster space read the raw cache." })] }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("dt", { children: "Recorded as" }), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("dd", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("code", { children: ["normalize=", normalize ? "l2" : "none"] }), " in the provenance of every column this run writes"] })] })
				]
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "geometry-state",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
					className: "geometry-state-title",
					children: "Embedding caches in this entry"
				}), caches.length ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("ul", { children: caches.map((item) => /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("li", {
					className: item.normalized ? "l2" : "raw",
					children: [
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("code", { children: item.tag }),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: item.normalized ? "L2-normalized sibling on disk" : "Raw only" }),
						item.tag === tag ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("b", { children: "this run" }) : null
					]
				}, item.tag)) }) : /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", { children: "None yet. This run creates the first one and sets the geometry for everything after it." })]
			}),
			stranded ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "notice warning geometry-alert",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Mixed geometry" }), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("span", {
					className: "notice-body",
					children: [
						"A normalized cache already exists at ",
						/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("code", { children: [
							"embeddings/normalized/",
							tag,
							".emb.tsv"
						] }),
						". These coordinates would be reduced from raw vectors while Query distances and Cluster space keep preferring that normalized file, leaving one entry with two geometries. Keep normalization on, or delete that file first."
					]
				})]
			}) : null
		]
	});
}
function EmbeddingCacheField({ values, setValue, entry }) {
	const caches = entry?.embeddings ?? [];
	const selected = String(values.embedding ?? "");
	const cache = caches.find((item) => item.tag === selected);
	const usingNormalized = Boolean(cache?.normalized) && !values.raw;
	return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(import_jsx_runtime.Fragment, { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
		label: "Embedding cache",
		hint: "Built by Build coordinates. Its geometry is fixed on disk; this tool only reads it.",
		children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("select", {
			value: selected,
			onChange: (event) => setValue("embedding", event.target.value),
			children: caches.length ? caches.map((item) => /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("option", {
				value: item.tag,
				children: [
					item.tag,
					" · ",
					item.normalized ? "normalized available" : "raw only"
				]
			}, item.tag)) : /* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
				value: selected,
				children: selected || "No embedding cache found"
			})
		})
	}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
		className: `notice ${usingNormalized ? "success" : "warning"}`,
		children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: usingNormalized ? "L2 geometry" : "Raw geometry" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
			className: "notice-body",
			children: usingNormalized ? `Reads embeddings/normalized/${selected}.emb.tsv and records normalize=l2.` : cache && !cache.normalized ? `${selected} has no normalized sibling, so this reads raw vectors and records normalize=none. Re-run Build coordinates with normalization on to change that.` : `Raw vectors requested, so this reads embeddings/${selected}.emb.tsv and records normalize=none.`
		})]
	})] });
}
function ConfigurationForm({ tool, values, setValue, runnerOnline, entry }) {
	const input = (key, type = "text") => ({
		type,
		value: String(values[key] ?? ""),
		onChange: (e) => setValue(key, type === "number" ? Number(e.target.value) : e.target.value)
	});
	const select = (key) => ({
		value: String(values[key] ?? ""),
		onChange: (e) => setValue(key, e.target.value)
	});
	if (tool === "initialize") return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(import_jsx_runtime.Fragment, { children: [
		/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "field-grid two",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Source file",
				hint: "A path or a filename in initial_files/",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(FilePathInput, {
					value: String(values.input || ""),
					onChange: (value) => setValue("input", value),
					accept: ".tsv,.csv,.fasta,.fa,.faa,.json",
					runnerOnline
				})
			}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Source format",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
					...select("source"),
					children: [
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
							value: "em",
							children: "EnzymeMiner / generic TSV"
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
							value: "fasta",
							children: "FASTA"
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
							value: "fs",
							children: "Foldseek JSON"
						})
					]
				})
			})]
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
			label: "Entry name",
			hint: "Defaults to the source filename stem",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("name") })
		}),
		values.source === "em" ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "field-grid two",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "ID column",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("idCol") })
			}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Sequence column",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("seqCol") })
			})]
		}) : null,
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
			label: "Query values",
			hint: values.source === "fasta" ? "Use complete FASTA headers; leave blank for no marked query" : "Space-separated values; overrides automatic query detection",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
				...input("query"),
				placeholder: "Optional"
			})
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
			checked: Boolean(values.force),
			onChange: (v) => setValue("force", v),
			label: "Delete and rebuild an existing entry",
			description: "Destructive: this removes the complete entry directory."
		})
	] });
	if (tool === "taxonomy") return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(import_jsx_runtime.Fragment, { children: [
		/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "field-grid two",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "NCBI contact email",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("email", "email") })
			}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Resolution strategy",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
					...select("strategy"),
					children: [
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
							value: "auto",
							children: "Auto-detect"
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
							value: "em",
							children: "NCBI protein accession"
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
							value: "foldseek",
							children: "Foldseek metadata"
						})
					]
				})
			})]
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
			label: "NCBI API key",
			hint: "Optional; passed to this job without being written to job history",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
				...input("apiKey", "password"),
				placeholder: "Optional"
			})
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "field-grid two",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "NCBI batch size",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
					...input("batch", "number"),
					min: "1"
				})
			}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "GMGC batch size",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
					...input("gmgcBatch", "number"),
					min: "1"
				})
			})]
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
			label: "Run mode",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
				...select("rerun"),
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: "resume",
						children: "Resume / preserve completed rows"
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: "retry",
						children: "Retry unresolved rows"
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: "force",
						children: "Refetch every row"
					})
				]
			})
		})
	] });
	if (tool === "external") return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(import_jsx_runtime.Fragment, { children: [
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
			label: "External CSV or TSV",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(FilePathInput, {
				value: String(values.file || ""),
				onChange: (value) => setValue("file", value),
				accept: ".csv,.tsv,.txt",
				runnerOnline
			})
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "field-grid two",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "External ID column",
				hint: "Blank uses the first column",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
					...input("idCol"),
					placeholder: "First column"
				})
			}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Column type",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
					...select("type"),
					children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: "label",
						children: "Label / metadata"
					}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: "coordinate",
						children: "Coordinate / plot axis"
					})]
				})
			})]
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
			label: "Columns to merge",
			hint: "Comma-separated; blank merges every non-ID column",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("columns") })
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
			label: "ID translator",
			hint: "Optional two-column table: SSE ID, then external ID",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(FilePathInput, {
				value: String(values.translator || ""),
				onChange: (value) => setValue("translator", value),
				accept: ".csv,.tsv,.txt",
				runnerOnline
			})
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
			label: "Delimiter",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
				...select("delimiter"),
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: "auto",
						children: "Infer from extension"
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: ",",
						children: "Comma"
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: "\\t",
						children: "Tab"
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: ";",
						children: "Semicolon"
					})
				]
			})
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
			checked: Boolean(values.force),
			onChange: (v) => setValue("force", v),
			label: "Replace colliding columns"
		})
	] });
	if (tool === "coordinates") {
		const structure = values.embedder !== "esmc";
		return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(import_jsx_runtime.Fragment, { children: [
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
				className: "choice-grid three",
				role: "radiogroup",
				"aria-label": "Embedding model",
				children: [
					{
						id: "esmc",
						name: "ESM-C",
						sub: "Sequence"
					},
					{
						id: "prostt5",
						name: "ProstT5",
						sub: "3Di structure"
					},
					{
						id: "saprot",
						name: "SaProt",
						sub: "Sequence + 3Di"
					}
				].map((item) => /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
					className: `choice-card ${values.embedder === item.id ? "selected" : ""}`,
					onClick: () => setValue("embedder", item.id),
					type: "button",
					"aria-pressed": values.embedder === item.id,
					children: [
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { className: "choice-dot" }),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: item.name }),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("small", { children: item.sub })
					]
				}, item.id))
			}),
			structure ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "notice warning",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Foldseek required" }), " This model needs a Foldseek entry and source JSON with C-alpha coordinates."]
			}) : null,
			/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "field-grid three",
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
						label: "Reduction",
						children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
							...select("reducer"),
							children: [
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
									value: "pca",
									children: "PCA"
								}),
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
									value: "umap",
									children: "UMAP"
								}),
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
									value: "tsne",
									children: "t-SNE"
								})
							]
						})
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
						label: "Components",
						children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
							...input("components", "number"),
							min: "1"
						})
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
						label: "Pooling",
						children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
							...select("pooling"),
							children: [
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", { children: "mean" }),
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", { children: "max" }),
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", { children: "min" })
							]
						})
					})
				]
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "field-grid two",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
					label: "Compute device",
					children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
						...select("device"),
						children: [
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
								value: "auto",
								children: "Auto: CUDA → MPS → CPU"
							}),
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
								value: "cuda",
								children: "CUDA"
							}),
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
								value: "mps",
								children: "Apple MPS"
							}),
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
								value: "cpu",
								children: "CPU"
							})
						]
					})
				}), values.embedder === "esmc" ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
					label: "ESM-C model",
					children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
						...select("esmcModel"),
						children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
							value: "esmc_600m",
							children: "ESM-C 600M"
						}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
							value: "esmc_300m",
							children: "ESM-C 300M"
						})]
					})
				}) : /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
					label: "Model checkpoint",
					children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input(values.embedder === "prostt5" ? "prostCheckpoint" : "saprotCheckpoint") })
				})]
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)(GeometrySection, {
				values,
				setValue,
				entry
			}),
			values.reducer === "umap" ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "subsection",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
					className: "subsection-title",
					children: "UMAP neighborhood"
				}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
					className: "field-grid three",
					children: [
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
							label: "Neighbors",
							children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
								...input("neighbors", "number"),
								min: "2"
							})
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
							label: "Minimum distance",
							children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
								...input("minDist", "number"),
								min: "0",
								step: "0.05"
							})
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
							label: "Metric",
							children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("metric") })
						})
					]
				})]
			}) : null,
			values.reducer === "tsne" ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "subsection",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
					className: "subsection-title",
					children: "t-SNE neighborhood"
				}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
					className: "field-grid two",
					children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
						label: "Perplexity",
						children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
							...input("perplexity", "number"),
							min: "1"
						})
					}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
						label: "PCA pre-reduction",
						children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
							...input("tsnePca", "number"),
							min: "0"
						})
					})]
				})]
			}) : null,
			structure ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "subsection",
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
						className: "subsection-title",
						children: "Structure input"
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
						className: "field-grid two",
						children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
							label: "Foldseek JSON",
							hint: "Blank resolves it from the manifest",
							children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(FilePathInput, {
								value: String(values.foldseekJson || ""),
								onChange: (value) => setValue("foldseekJson", value),
								accept: ".json",
								runnerOnline
							})
						}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
							label: "Maximum residues",
							children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
								...input("maxResidues", "number"),
								min: "3"
							})
						})]
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
						checked: Boolean(values.includeEmpty),
						onChange: (v) => setValue("includeEmpty", v),
						label: "Keep rows without usable structures",
						description: "Their coordinate cells remain empty."
					})
				]
			}) : null,
			/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("details", {
				className: "advanced",
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("summary", { children: "Advanced run and cache options" }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
						className: "field-grid three advanced-body",
						children: [
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
								label: "Batch size",
								children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
									...input("batchSize", "number"),
									min: "1"
								})
							}),
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
								label: "Write every",
								children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
									...input("writeEvery", "number"),
									min: "1"
								})
							}),
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
								label: "Custom tag",
								children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
									...input("label"),
									placeholder: "Automatic"
								})
							})
						]
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
						label: "If this coordinate system already exists",
						children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
							...select("rerun"),
							children: [
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
									value: "new",
									children: "Stop and ask"
								}),
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
									value: "rereduce",
									children: "Reuse embeddings and re-reduce"
								}),
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
									value: "reembed",
									children: "Recompute embeddings and reduce"
								})
							]
						})
					})
				]
			})
		] });
	}
	if (tool === "distance") return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(import_jsx_runtime.Fragment, { children: [
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(EmbeddingCacheField, {
			values,
			setValue,
			entry
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
			label: "Reference sequences",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
				...select("queryMode"),
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
					value: "marked",
					children: "Use rows marked query=True"
				}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
					value: "explicit",
					children: "Choose explicit IDs"
				})]
			})
		}),
		values.queryMode === "explicit" ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
			label: "Query IDs",
			hint: "Comma-separated SSE IDs",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("queryIds") })
		}) : /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "notice success",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("span", { children: [formatCount(entry?.queries), " queries detected"] }), " Rows marked query=True in this entry."]
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
			checked: Boolean(values.raw),
			onChange: (v) => setValue("raw", v),
			label: "Use raw embeddings",
			description: "Ignore the normalized sibling and measure distances on raw vectors instead."
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
			checked: Boolean(values.force),
			onChange: (v) => setValue("force", v),
			label: "Replace existing distance columns"
		})
	] });
	if (tool === "cluster") return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(import_jsx_runtime.Fragment, { children: [
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(EmbeddingCacheField, {
			values,
			setValue,
			entry
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
			className: "choice-grid two",
			role: "radiogroup",
			"aria-label": "Clustering technique",
			children: [{
				id: "kmeans",
				name: "K-means",
				sub: "Assign every sequence"
			}, {
				id: "hdbscan",
				name: "HDBSCAN",
				sub: "Discover groups + noise"
			}].map((item) => /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
				className: `choice-card ${values.clusterer === item.id ? "selected" : ""}`,
				onClick: () => setValue("clusterer", item.id),
				type: "button",
				"aria-pressed": values.clusterer === item.id,
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { className: "choice-dot" }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: item.name }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("small", { children: item.sub })
				]
			}, item.id))
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "field-grid two",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Clustering space",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
					...select("space"),
					children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: "pca",
						children: "PCA-reduced embeddings"
					}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: "full",
						children: "Full embedding dimensions"
					})]
				})
			}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Output tag",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
					...input("label"),
					placeholder: "Use embedding tag"
				})
			})]
		}),
		values.space === "pca" ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "subsection",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "segmented",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
					type: "button",
					className: values.pcaMode === "dims" ? "active" : "",
					onClick: () => setValue("pcaMode", "dims"),
					children: "Fixed dimensions"
				}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
					type: "button",
					className: values.pcaMode === "variance" ? "active" : "",
					onClick: () => setValue("pcaMode", "variance"),
					children: "Variance target"
				})]
			}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: values.pcaMode === "dims" ? "PCA dimensions" : "Variance fraction",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
					...input(values.pcaMode === "dims" ? "pcaDims" : "pcaVariance", "number"),
					step: values.pcaMode === "dims" ? 1 : .01
				})
			})]
		}) : null,
		values.clusterer === "kmeans" ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "subsection",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "segmented",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
					type: "button",
					className: values.kMode === "auto" ? "active" : "",
					onClick: () => setValue("kMode", "auto"),
					children: "Auto-select k"
				}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
					type: "button",
					className: values.kMode === "fixed" ? "active" : "",
					onClick: () => setValue("kMode", "fixed"),
					children: "Fixed k"
				})]
			}), values.kMode === "auto" ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "field-grid two",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
					label: "Minimum k",
					children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("kMin", "number") })
				}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
					label: "Maximum k",
					children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("kMax", "number") })
				})]
			}) : /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Number of clusters",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("k", "number") })
			})]
		}) : /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "field-grid two subsection",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Minimum cluster size",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("minClusterSize", "number") })
			}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Minimum samples",
				hint: "Blank follows cluster size",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
					...input("minSamples", "number"),
					placeholder: "Automatic"
				})
			})]
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
			checked: Boolean(values.analysis),
			onChange: (v) => setValue("analysis", v),
			label: "Generate Tier-2 analysis",
			description: "Profiles, enrichment, and representative sequences."
		}),
		values.analysis ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "field-grid two",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Representatives per cluster",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", { ...input("topN", "number") })
			}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "FDR threshold",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
					...input("fdr", "number"),
					step: "0.01"
				})
			})]
		}) : null,
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
			checked: Boolean(values.raw),
			onChange: (v) => setValue("raw", v),
			label: "Use raw embedding geometry",
			description: "Ignore the normalized sibling and cluster on raw vectors instead."
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
			checked: Boolean(values.force),
			onChange: (v) => setValue("force", v),
			label: "Replace an existing matching clustering"
		})
	] });
	if (tool === "boltz") {
		const selections = entry?.selections ?? [];
		const refOptions = ((values.selection ? selections.find((sel) => sel.name === values.selection) : selections[0]) ?? selections[0])?.ids ?? [];
		const datafileName = `${entry?.name ?? "<entry>"}.sse.tsv`;
		return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(import_jsx_runtime.Fragment, { children: [
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Selection to analyze",
				hint: "Exported from the visualizer's “Export selection for Boltz” button",
				children: selections.length ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
					...select("selection"),
					children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
						value: "",
						children: "Most recent selection"
					}), selections.map((sel) => /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("option", {
						value: sel.name,
						children: [sel.name, sel.count != null ? ` (${sel.count} seq)` : ""]
					}, sel.name))]
				}) : /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
					className: "notice",
					children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "No selections found" }), " Select points in the explorer and click “Export selection for Boltz”, then reload this page."]
				})
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "NVIDIA API key",
				hint: "Required; passed to this job without being written to job history",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
					...input("apiKey", "password"),
					placeholder: "nvapi-…"
				})
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Substrate SMILES",
				hint: "Optional; one per line. Present = holo prediction, empty = apo",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("textarea", {
					value: String(values.smiles ?? ""),
					onChange: (e) => setValue("smiles", e.target.value),
					placeholder: "One SMILES per line",
					rows: 2
				})
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
				label: "Ligand label",
				hint: "Optional; names the holo output columns/folders",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
					...input("smilesLabel"),
					placeholder: "e.g. UDP-Glc"
				})
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
				checked: Boolean(values.useMsa),
				onChange: (v) => setValue("useMsa", v),
				label: "Generate MSA (recommended)",
				description: "Uses ColabFold; more accurate but slower. Off runs single-sequence."
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("details", {
				className: "advanced",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("summary", { children: "Prediction parameters" }), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
					className: "field-grid two advanced-body",
					children: [
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
							label: "Recycling steps",
							children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
								...input("recyclingSteps", "number"),
								min: "1",
								max: "10"
							})
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
							label: "Sampling steps",
							children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
								...input("samplingSteps", "number"),
								min: "10",
								max: "500"
							})
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
							label: "Diffusion samples",
							children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
								...input("diffusionSamples", "number"),
								min: "1",
								max: "10"
							})
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
							label: "Step scale",
							children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
								...input("stepScale", "number"),
								step: "0.001",
								min: "0.1",
								max: "5"
							})
						})
					]
				})]
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
				checked: Boolean(values.rmsd),
				onChange: (v) => setValue("rmsd", v),
				label: "Measure structural RMSD",
				description: "After prediction, run the RMSD comparison (Kabsch superposition) over the predicted apo structures."
			}),
			values.rmsd ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "subsection",
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
						className: "subsection-title",
						children: "RMSD structural comparison"
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
						className: "notice",
						children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
							className: "notice-body",
							children: [
								"Every predicted apo structure is aligned to one ",
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: "reference structure" }),
								" and its RMSD is measured. Results are appended to ",
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("code", { children: datafileName }),
								" as ",
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("code", { children: "RMSD_vs_<reference>_r<rank>_<method>" }),
								" columns — reload the explorer to color by them."
							]
						})
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
						label: "Reference structure",
						hint: "The structure every other one is compared against — one of the analyzed sequences",
						children: refOptions.length ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
							...select("rmsdReference"),
							children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
								value: "",
								children: "Select a reference…"
							}), refOptions.map((id) => /* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
								value: id,
								children: id
							}, id))]
						}) : /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
							...input("rmsdReference"),
							placeholder: "e.g. OleD_S1"
						})
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
						className: "field-grid two",
						children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
							label: "Reference rank",
							hint: "Which predicted rank of the reference",
							children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
								...input("rmsdRefRank", "number"),
								min: "0"
							})
						}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
							label: "Alignment method",
							children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
								...select("rmsdMethod"),
								children: [
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
										value: "seq",
										children: "Sequence-guided"
									}),
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
										value: "ce",
										children: "Structure-based (CE)"
									}),
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
										value: "both",
										children: "Both"
									})
								]
							})
						})]
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
						label: "Compare against",
						hint: "Which query structures to measure against the reference",
						children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("select", {
							...select("rmsdScope"),
							children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
								value: "all",
								children: "All predicted apo structures"
							}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", {
								value: "selected",
								children: "Only this selection’s sequences"
							})]
						})
					})
				]
			}) : null,
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Toggle, {
				checked: Boolean(values.force),
				onChange: (v) => setValue("force", v),
				label: "Force re-run (ignore cached predictions)"
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
				className: "notice",
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
					className: "notice-body",
					children: [
						"Outputs: ranked ",
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("code", { children: ".cif" }),
						" structures under ",
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("code", { children: "structures/" }),
						", binding scores in ",
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("code", { children: "logs/boltz_log.csv" }),
						", and pTM/pLDDT (plus any RMSD) columns appended to ",
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("code", { children: datafileName }),
						". Use “Open structures folder” on the finished job to browse them."
					]
				})
			})
		] });
	}
	return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(import_jsx_runtime.Fragment, { children: [
		/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "notice success",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Entry ready" }), " 22 coordinate columns and 2 clustering systems detected."]
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsx)(Field, {
			label: "Dash server port",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("input", {
				...input("port", "number"),
				min: "1024",
				max: "65535"
			})
		}),
		/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "explorer-preview",
			children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
				className: "mini-plot",
				"aria-label": "Illustrative sequence-space scatter plot",
				children: Array.from({ length: 34 }).map((_, index) => /* @__PURE__ */ (0, import_jsx_runtime.jsx)("i", { style: {
					"--x": `${8 + index * 29 % 82}%`,
					"--y": `${10 + index * 47 % 75}%`,
					"--d": `${index % 5 * .04}s`
				} }, index))
			}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: "Sequence Space Explorer" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Coordinates, filters, layers, Boltz-2, RMSD, and exports" })] })]
		})
	] });
}
function LandingPage({ onEnter }) {
	return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
		className: "landing-page",
		children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
			className: "landing-copy",
			children: [
				/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
					className: "landing-kicker",
					children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {}), " Protein sequence intelligence"]
				}),
				/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("h1", { children: [
					"Sequence Space",
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("br", {}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("em", { children: "Explorer" })
				] }),
				/* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", { children: "Build, enrich, map, and analyze protein sequence spaces through one guided visual workflow." }),
				/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
					className: "landing-actions",
					children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
						type: "button",
						className: "landing-primary",
						onClick: onEnter,
						children: ["Enter pipeline ", /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "→" })]
					}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("span", {
						className: "landing-note",
						children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("b", { children: tools.length }), " integrated tools · terminal optional"]
					})]
				}),
				/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
					className: "landing-capabilities",
					"aria-label": "Pipeline capabilities",
					children: [
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Embed" }),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("i", {}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Reduce" }),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("i", {}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Cluster" }),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("i", {}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Explore" })
					]
				})
			]
		}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
			className: "landing-visual",
			"aria-label": "Stylized sequence-space visualization example",
			children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "hero-plot",
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { className: "density-field density-slate" }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { className: "density-field density-magenta" }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { className: "density-field density-gray" }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { className: "density-field density-teal" }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { className: "density-field density-purple" }),
					landingPoints.map((point, index) => /* @__PURE__ */ (0, import_jsx_runtime.jsx)("i", {
						className: `data-point population-${point.cluster}`,
						style: {
							"--x": `${point.x}%`,
							"--y": `${point.y}%`,
							"--s": `${point.size}px`,
							"--c": point.color,
							"--o": point.opacity
						}
					}, index))
				]
			})
		})]
	});
}
var terminalStatuses = new Set([
	"succeeded",
	"failed",
	"cancelled",
	"interrupted"
]);
function jobStatusLabel(job) {
	if (job.ready) return "Explorer ready";
	return {
		queued: "Queued",
		running: "Running",
		cancelling: "Stopping",
		succeeded: "Complete",
		failed: "Failed",
		cancelled: "Cancelled",
		interrupted: "Interrupted"
	}[job.status];
}
function formatCount(value) {
	return typeof value === "number" ? value.toLocaleString() : "—";
}
function Home() {
	const [tool, setTool] = (0, import_react.useState)("coordinates");
	const [entry, setEntry] = (0, import_react.useState)("EnzymeMiner_Selection_Table_ri4plk");
	const [entries, setEntries] = (0, import_react.useState)([]);
	const [values, setValues] = (0, import_react.useState)(initialValues);
	const [job, setJob] = (0, import_react.useState)(null);
	const [runnerState, setRunnerState] = (0, import_react.useState)("connecting");
	const [environmentReady, setEnvironmentReady] = (0, import_react.useState)(false);
	const [environmentMissing, setEnvironmentMissing] = (0, import_react.useState)([]);
	const [runnerError, setRunnerError] = (0, import_react.useState)("");
	const [submitting, setSubmitting] = (0, import_react.useState)(false);
	const [copied, setCopied] = (0, import_react.useState)(false);
	const [showTools, setShowTools] = (0, import_react.useState)(true);
	const [theme, setTheme] = (0, import_react.useState)(() => {
		if (typeof window === "undefined") return "modern";
		return window.localStorage.getItem("sse-theme") === "original" ? "original" : "modern";
	});
	const [showLanding, setShowLanding] = (0, import_react.useState)(true);
	const [entryFiles, setEntryFiles] = (0, import_react.useState)(null);
	const [showFullLog, setShowFullLog] = (0, import_react.useState)(false);
	const [shuttingDown, setShuttingDown] = (0, import_react.useState)(false);
	const explorerWindow = (0, import_react.useRef)(null);
	const active = tools.find((item) => item.id === tool);
	const activeEntry = entries.find((item) => item.name === entry);
	const previewCommand = (0, import_react.useMemo)(() => buildCommand(tool, values[tool], entry), [
		tool,
		values,
		entry
	]);
	const command = job?.tool === tool && job.entry === (tool === "initialize" ? String(values.initialize.name || entry) : entry) ? job.command : previewCommand;
	const running = job && [
		"queued",
		"running",
		"cancelling"
	].includes(job.status);
	const visualizerRunning = Boolean(job && job.tool === "visualizer" && [
		"queued",
		"running",
		"cancelling"
	].includes(job.status));
	const jobId = job?.id;
	const jobStatus = job?.status;
	const refreshEntries = (0, import_react.useCallback)(async () => {
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
	(0, import_react.useEffect)(() => {
		const timer = window.setTimeout(() => void refreshEntries(), 0);
		return () => window.clearTimeout(timer);
	}, [refreshEntries]);
	(0, import_react.useEffect)(() => {
		if (runnerState !== "online" || !entry || job) return;
		getJobs(entry).then((history) => {
			if (history[0]) setJob(history[0]);
		}).catch(() => void 0);
	}, [
		entry,
		job,
		runnerState
	]);
	(0, import_react.useEffect)(() => {
		if (!jobId || !jobStatus || terminalStatuses.has(jobStatus)) return;
		const interval = window.setInterval(async () => {
			try {
				const updated = await getJob(jobId);
				setJob(updated);
				if (terminalStatuses.has(updated.status)) refreshEntries();
			} catch (reason) {
				setRunnerError(reason instanceof Error ? reason.message : "Could not refresh the job");
			}
		}, 1e3);
		return () => window.clearInterval(interval);
	}, [
		jobId,
		jobStatus,
		refreshEntries
	]);
	(0, import_react.useEffect)(() => {
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
	const setValue = (key, value) => {
		setValues((current) => ({
			...current,
			[tool]: {
				...current[tool],
				[key]: value
			}
		}));
	};
	const chooseTheme = (nextTheme) => {
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
			const submitted = await submitJob({
				tool,
				entry,
				values: values[tool]
			});
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
		try {
			setJob(await cancelPipelineJob(job.id));
		} catch (reason) {
			setRunnerError(reason instanceof Error ? reason.message : "Could not stop the job");
		}
	};
	const shutDownApp = async () => {
		if (!window.confirm("Shut down Sequence Space Explorer? Any active pipeline job will be stopped.")) return;
		setRunnerError("");
		try {
			await requestPipelineShutdown();
			setShuttingDown(true);
		} catch (reason) {
			setRunnerError(reason instanceof Error ? reason.message : "Could not shut down the app");
		}
	};
	const copyCommand = async () => {
		try {
			await navigator.clipboard.writeText(command);
		} catch {}
		setCopied(true);
		window.setTimeout(() => setCopied(false), 1500);
	};
	const viewFiles = async () => {
		setRunnerError("");
		try {
			setEntryFiles(await getEntryFiles(entry));
		} catch (reason) {
			setRunnerError(reason instanceof Error ? reason.message : "Could not list entry files");
		}
	};
	const revealFolder = async (subpath) => {
		setRunnerError("");
		try {
			await revealEntryPath(entry, subpath);
		} catch (reason) {
			setRunnerError(reason instanceof Error ? reason.message : "Could not open the folder");
		}
	};
	const completedTool = (id) => Boolean(activeEntry && (id === "initialize" || id === "taxonomy" && activeEntry.taxonomyFields > 0 || id === "coordinates" && activeEntry.coordinates > 0 || id === "distance" && activeEntry.annotations > 0 || id === "cluster" && activeEntry.analyses > 0));
	return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("main", {
		className: `app-shell ${showLanding ? "landing-open" : ""}`,
		"data-theme": theme,
		children: [
			/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "ambient-space",
				"aria-hidden": "true",
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { className: "ambient-orbit orbit-one" }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { className: "ambient-orbit orbit-two" }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { className: "ambient-glow glow-one" }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { className: "ambient-glow glow-two" }),
					Array.from({ length: 76 }).map((_, index) => /* @__PURE__ */ (0, import_jsx_runtime.jsx)("i", { style: {
						"--x": `${2 + index * 37 % 92}%`,
						"--y": `${8 + index * 53 % 84}%`,
						"--s": `${2 + index % 4}px`,
						"--o": `${.16 + index % 6 * .09}`,
						"--delay": `${index % 11 * -.43}s`
					} }, index))
				]
			}),
			/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("header", {
				className: "topbar",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
					className: "brand",
					type: "button",
					onClick: () => setShowLanding(true),
					"aria-label": "Return to the Sequence Space Explorer landing page",
					children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
						className: "brand-mark",
						children: [
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {}),
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {}),
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {})
						]
					}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: "Sequence Space Explorer" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Pipeline control center" })] })]
				}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
					className: "top-actions",
					children: [
						/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
							className: `runner-badge ${runnerState === "online" && !environmentReady ? "offline" : runnerState}`,
							children: [
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {}),
								" ",
								runnerState === "online" ? environmentReady ? "Runner connected" : "Setup needed" : runnerState === "offline" ? "Runner offline" : "Connecting"
							]
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
							className: "theme-switch",
							role: "group",
							"aria-label": "Interface theme",
							children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
								type: "button",
								className: theme === "modern" ? "active" : "",
								"aria-pressed": theme === "modern",
								onClick: () => chooseTheme("modern"),
								children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "◐" }), " Modern"]
							}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
								type: "button",
								className: theme === "original" ? "active" : "",
								"aria-pressed": theme === "original",
								onClick: () => chooseTheme("original"),
								children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "○" }), " Original"]
							})]
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
							className: "icon-button",
							type: "button",
							"aria-label": "Reconnect to runner",
							onClick: () => void refreshEntries(),
							children: "↻"
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
							className: "shutdown-button",
							type: "button",
							onClick: () => void shutDownApp(),
							disabled: runnerState !== "online",
							children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
								"aria-hidden": "true",
								children: "⏻"
							}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("b", { children: "Shut down" })]
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
							className: "avatar",
							children: "WG"
						})
					]
				})]
			}),
			showLanding ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)(LandingPage, { onEnter: () => setShowLanding(false) }) : /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "workspace",
				children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("aside", {
					className: `sidebar ${showTools ? "" : "collapsed"}`,
					children: [
						/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
							className: "side-heading",
							children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Workspace" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
								type: "button",
								"aria-label": "Collapse tools",
								onClick: () => setShowTools(!showTools),
								children: showTools ? "‹" : "›"
							})]
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
							className: "entry-card",
							children: [
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
									className: "entry-icon",
									children: "E"
								}),
								/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("span", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("select", {
									"aria-label": "Selected entry",
									value: entry,
									onChange: (event) => {
										const embedding = entries.find((item) => item.name === event.target.value)?.embeddings[0]?.tag;
										setEntry(event.target.value);
										setJob(null);
										if (embedding) setValues((current) => ({
											...current,
											distance: {
												...current.distance,
												embedding
											},
											cluster: {
												...current.cluster,
												embedding
											}
										}));
									},
									disabled: !entries.length,
									children: entries.length ? entries.map((item) => /* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", { children: item.name }, item.name)) : /* @__PURE__ */ (0, import_jsx_runtime.jsx)("option", { children: "No entries found" })
								}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("small", { children: [formatCount(activeEntry?.sequences), " sequences"] })] }),
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("b", { children: "⌄" })
							]
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("nav", {
							"aria-label": "Pipeline tools",
							children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
								className: "nav-label",
								children: "Pipeline"
							}), tools.map((item) => /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
								type: "button",
								onClick: () => setTool(item.id),
								className: tool === item.id ? "active" : "",
								children: [
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
										className: "nav-mark",
										children: item.mark
									}),
									/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("span", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: item.label }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("small", { children: item.stage })] }),
									completedTool(item.id) ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("i", {
										className: "complete",
										children: "✓"
									}) : null
								]
							}, item.id))]
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
							className: "side-footer",
							children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
								className: `health ${runnerState === "online" && !environmentReady ? "offline" : runnerState}`,
								children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: runnerState === "online" ? environmentReady ? "Environment ready" : "Python setup needed" : "Runner unavailable" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("small", { children: runnerState === "online" ? environmentReady ? "Local pipeline access" : `${environmentMissing.length} packages missing` : "Start the local runner" })] })]
							}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
								type: "button",
								onClick: () => void refreshEntries(),
								children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "↻" }), " Refresh environment"]
							})]
						})
					]
				}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
					className: "content",
					children: [
						runnerError ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
							className: "runner-alert",
							role: "alert",
							children: [
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: runnerState === "offline" ? "Local runner not connected." : "Action needed." }),
								" ",
								runnerError,
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
									type: "button",
									onClick: () => void refreshEntries(),
									children: "Retry"
								})
							]
						}) : null,
						/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
							className: "entry-header",
							children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
									className: "eyebrow",
									children: "Selected entry"
								}),
								/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
									className: "entry-title-row",
									children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("h1", { children: entry || "No entry selected" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
										className: "ready-pill",
										children: activeEntry ? "Ready" : "Unavailable"
									})]
								}),
								/* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", { children: activeEntry ? `${activeEntry.source.toUpperCase()} source · Updated ${new Date(activeEntry.updatedAt).toLocaleString()}` : "Connect the local runner to discover entries." })
							] }), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
								className: "entry-actions",
								children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
									className: "secondary-button",
									type: "button",
									onClick: () => void viewFiles(),
									disabled: !activeEntry,
									children: "View files"
								}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
									className: "secondary-button",
									type: "button",
									onClick: () => setTool("visualizer"),
									disabled: !activeEntry,
									children: "Open explorer ↗"
								})]
							})]
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
							className: "stats-grid",
							children: [
								/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Sequences" }),
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: formatCount(activeEntry?.sequences) }),
									/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("small", { children: [formatCount(activeEntry?.queries), " marked queries"] })
								] }),
								/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Annotations" }),
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: formatCount(activeEntry?.annotations) }),
									/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("small", { children: [formatCount(activeEntry?.taxonomyFields), " taxonomy fields"] })
								] }),
								/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Coordinates" }),
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: formatCount(activeEntry?.coordinates) }),
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("small", { children: activeEntry?.embeddings.map((item) => item.tag).join(" · ") || "No embedding cache" })
								] }),
								/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Analyses" }),
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: formatCount(activeEntry?.analyses) }),
									/* @__PURE__ */ (0, import_jsx_runtime.jsx)("small", { children: activeEntry?.analyses ? "Cluster reports available" : "No cluster reports" })
								] })
							]
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
							className: "stepper",
							"aria-label": "Pipeline progress",
							children: [
								"Entry",
								"Enrich",
								"Embed",
								"Analyze",
								"Explore"
							].map((name, index) => {
								const done = Boolean(activeEntry && (index === 0 || index === 1 && activeEntry.annotations > 2 || index === 2 && activeEntry.coordinates > 0 || index === 3 && activeEntry.analyses > 0));
								return /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
									className: done ? "done" : index === 4 ? "current" : "",
									children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: done ? "✓" : index + 1 }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("b", { children: name })]
								}, name);
							})
						}),
						/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
							className: "panel-grid",
							children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
								className: "panel configure-panel",
								children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
									className: "panel-header",
									children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
										className: "tool-heading",
										children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: active.mark }), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [
											/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
												className: "eyebrow",
												children: active.stage
											}),
											/* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", { children: active.label }),
											/* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", { children: active.description })
										] })]
									}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
										className: "script-chip",
										children: active.script
									})]
								}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
									className: "panel-body",
									children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)(ConfigurationForm, {
										tool,
										values: values[tool],
										setValue,
										runnerOnline: runnerState === "online",
										entry: activeEntry
									})
								})]
							}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("aside", {
								className: "run-column",
								children: [
									/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
										className: "panel run-panel",
										children: [
											/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
												className: "panel-header compact",
												children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
													className: "eyebrow",
													children: "Review"
												}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", { children: "Run command" })] }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
													className: `simulation ${runnerState}`,
													children: runnerState === "online" ? "Validated locally" : "Runner required"
												})]
											}),
											/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
												className: "command-box",
												children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Generated command" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
													type: "button",
													onClick: copyCommand,
													children: copied ? "Copied" : "Copy"
												})] }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("code", { children: command })]
											}),
											/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
												className: "run-summary",
												children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Writes to" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: tool === "visualizer" ? "Local browser session" : tool === "initialize" ? `entries/${String(values.initialize.name || "<entry>")}/` : `${entry}.sse.tsv` })] }), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "Execution" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: tool === "visualizer" ? "Long-lived local service" : tool === "coordinates" ? "Model and cache dependent" : "Serialized for this entry" })] })]
											}),
											tool === "visualizer" && visualizerRunning ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
												className: "primary-button danger",
												type: "button",
												onClick: () => void cancelJob(),
												disabled: job?.status === "cancelling",
												children: ["⏻ ", job?.status === "cancelling" ? "Shutting down explorer…" : "Shut down explorer"]
											}) : /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("button", {
												className: "primary-button",
												type: "button",
												onClick: () => void startJob(),
												disabled: runnerState !== "online" || !environmentReady || submitting || Boolean(running) || tool !== "initialize" && !activeEntry,
												children: ["▶ ", submitting ? "Submitting…" : tool === "visualizer" ? "Start explorer" : `Run ${active.label.toLowerCase()}`]
											}),
											/* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", {
												className: "simulation-note",
												children: "Commands run without a shell in the configured Python environment. Writers for the same entry are queued to prevent data corruption."
											})
										]
									}),
									/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
										className: `panel job-panel ${job ? "visible" : ""}`,
										"aria-live": "polite",
										children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
											className: "panel-header compact",
											children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
												className: "eyebrow",
												children: "Latest job"
											}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", { children: job ? tools.find((item) => item.id === job.tool)?.label || job.tool : "No active run" })] }), job ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {
												className: `job-status ${job.status}`,
												children: jobStatusLabel(job)
											}) : null]
										}), job ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)(import_jsx_runtime.Fragment, { children: [
											/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
												className: "progress-meta",
												children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: job.ready ? "Explorer is accepting connections" : job.status === "succeeded" ? "Finished successfully" : job.status === "failed" ? `Exited with code ${job.exitCode ?? "?"}` : job.status === "queued" ? "Waiting for this entry" : job.status === "cancelled" ? "Stopped safely" : "Processing entry" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: job.progress === null ? "Live" : `${job.progress}%` })]
											}),
											/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
												className: `progress-track ${job.progress === null ? "indeterminate" : ""}`,
												children: /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { style: job.progress === null ? void 0 : { width: `${job.progress}%` } })
											}),
											/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
												className: "log-window",
												children: job.logs.slice(-8).map((line, index) => /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
													className: line.stream,
													children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: index === job.logs.slice(-8).length - 1 && ["running", "queued"].includes(job.status) ? "›" : line.stream === "stderr" ? "!" : "✓" }), line.text]
												}, `${line.at}-${index}`))
											}),
											/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
												className: "job-actions",
												children: [
													[
														"running",
														"queued",
														"cancelling"
													].includes(job.status) ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
														type: "button",
														onClick: () => void cancelJob(),
														disabled: job.status === "cancelling",
														children: "Cancel job"
													}) : /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
														type: "button",
														onClick: () => setJob(null),
														children: "Dismiss"
													}),
													job.ready && job.url ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
														type: "button",
														onClick: () => window.open(job.url, "_blank", "noopener,noreferrer"),
														children: "Open explorer ↗"
													}) : /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
														type: "button",
														onClick: () => setShowFullLog(true),
														disabled: !job.logs.length,
														children: "View full log"
													}),
													job.tool === "boltz" && job.status === "succeeded" ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
														type: "button",
														onClick: () => void revealFolder("structures"),
														children: "Open structures folder ↗"
													}) : null
												]
											})
										] }) : /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
											className: "empty-job",
											children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: "○" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", { children: "Configure a tool and run it to see its real progress, logs, and exit status here." })]
										})]
									}),
									/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
										className: "panel guidance-panel",
										children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
											className: "guidance-icon",
											children: "i"
										}), /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("strong", { children: "Workflow guidance" }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", { children: tool === "cluster" ? "Add taxonomy and external annotations before clustering if you want them included in enrichment tests. Clustering reads the embedding cache as it was written; it does not normalize." : tool === "coordinates" ? "This is the only tool that normalizes. Distances and clustering re-read the cache it writes, so the geometry chosen here is the geometry the whole entry gets." : tool === "distance" ? "Distances are Euclidean over the cache as it was written. On normalized vectors that is cosine geometry; on raw vectors it is not." : tool === "initialize" ? "Creation is bootstrap-only. Later changes should use the additive enrichment tools." : tool === "boltz" ? "Select points in the explorer and export them first. Prediction saves .cif structures under structures/ and appends pTM/pLDDT (and any RMSD) columns to the datafile; reload the explorer to see them." : "Every completed tool records its parameters in the entry manifest or run logs." })] })]
									})
								]
							})]
						})
					]
				})]
			}),
			entryFiles ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
				className: "modal-backdrop",
				role: "presentation",
				onMouseDown: () => setEntryFiles(null),
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
					className: "modal-panel",
					role: "dialog",
					"aria-modal": "true",
					"aria-label": `Files for ${entry}`,
					onMouseDown: (event) => event.stopPropagation(),
					children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
						className: "modal-header",
						children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
							className: "eyebrow",
							children: "Entry contents"
						}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", { children: entry })] }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
							type: "button",
							onClick: () => setEntryFiles(null),
							"aria-label": "Close",
							children: "×"
						})]
					}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
						className: "file-list",
						children: entryFiles.map((file) => /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("code", { children: file.path }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: file.size < 1024 ? `${file.size} B` : file.size < 1024 * 1024 ? `${(file.size / 1024).toFixed(1)} KB` : `${(file.size / 1024 / 1024).toFixed(1)} MB` })] }, file.path))
					})]
				})
			}) : null,
			showFullLog && job ? /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
				className: "modal-backdrop",
				role: "presentation",
				onMouseDown: () => setShowFullLog(false),
				children: /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("section", {
					className: "modal-panel log-modal",
					role: "dialog",
					"aria-modal": "true",
					"aria-label": "Full job log",
					onMouseDown: (event) => event.stopPropagation(),
					children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
						className: "modal-header",
						children: [/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", { children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
							className: "eyebrow",
							children: jobStatusLabel(job)
						}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("h2", { children: job.command })] }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("button", {
							type: "button",
							onClick: () => setShowFullLog(false),
							"aria-label": "Close",
							children: "×"
						})]
					}), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("div", {
						className: "full-log",
						children: job.logs.map((line, index) => /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
							className: line.stream,
							children: [/* @__PURE__ */ (0, import_jsx_runtime.jsx)("time", { children: new Date(line.at).toLocaleTimeString() }), /* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", { children: line.text })]
						}, `${line.at}-${index}`))
					})]
				})
			}) : null,
			shuttingDown ? /* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
				className: "shutdown-screen",
				role: "status",
				children: [
					/* @__PURE__ */ (0, import_jsx_runtime.jsxs)("div", {
						className: "brand-mark",
						"aria-hidden": "true",
						children: [
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {}),
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {}),
							/* @__PURE__ */ (0, import_jsx_runtime.jsx)("span", {})
						]
					}),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("h1", { children: "Shutting down…" }),
					/* @__PURE__ */ (0, import_jsx_runtime.jsx)("p", { children: "Pipeline jobs and local services are being stopped safely. You can close this browser tab." })
				]
			}) : null
		]
	});
}
//#endregion
export { Home as default };
