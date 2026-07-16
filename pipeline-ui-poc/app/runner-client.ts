export type EntrySummary = {
  name: string;
  source: string;
  sequences: number | null;
  queries: number | null;
  annotations: number;
  taxonomyFields: number;
  coordinates: number;
  analyses: number;
  embeddings: Array<{ tag: string; normalized: boolean }>;
  selections?: Array<{ name: string; count: number | null; created: string | null; ids?: string[] }>;
  updatedAt: string;
};

export type JobStatus = "queued" | "running" | "cancelling" | "succeeded" | "failed" | "cancelled" | "interrupted";

export type PipelineJob = {
  id: string;
  tool: string;
  entry: string;
  command: string;
  status: JobStatus;
  progress: number | null;
  logs: Array<{ at: string; stream: "stdout" | "stderr" | "system"; text: string }>;
  createdAt: string;
  startedAt?: string;
  finishedAt?: string;
  exitCode?: number | null;
  error?: string;
  ready?: boolean;
  url?: string;
};

export type EntryFile = { path: string; size: number; modifiedAt: string };

export const runnerUrl = (process.env.NEXT_PUBLIC_SSE_RUNNER_URL || "http://127.0.0.1:8788").replace(/\/$/, "");

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const response = await fetch(`${runnerUrl}${path}`, init);
  const payload = await response.json().catch(() => ({}));
  if (!response.ok) throw new Error(payload.error || `Runner request failed (${response.status})`);
  return payload as T;
}

export async function getRunnerState() {
  const [health, entries] = await Promise.all([
    request<{ ok: true; root: string; python: string; activeJobs: number; environment: { ready: boolean; missing: string[]; error?: string } }>("/api/health"),
    request<{ entries: EntrySummary[] }>("/api/entries"),
  ]);
  return { health, entries: entries.entries };
}

export async function getJobs(entry?: string) {
  const suffix = entry ? `?entry=${encodeURIComponent(entry)}` : "";
  return (await request<{ jobs: PipelineJob[] }>(`/api/jobs${suffix}`)).jobs;
}

export async function getJob(id: string) {
  return (await request<{ job: PipelineJob }>(`/api/jobs/${encodeURIComponent(id)}`)).job;
}

export async function submitJob(payload: { tool: string; entry: string; values: Record<string, string | boolean | number> }) {
  return (await request<{ job: PipelineJob }>("/api/jobs", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  })).job;
}

export async function cancelPipelineJob(id: string) {
  return (await request<{ job: PipelineJob }>(`/api/jobs/${encodeURIComponent(id)}`, { method: "DELETE" })).job;
}

export async function requestPipelineShutdown() {
  return await request<{ ok: true; managed: boolean }>("/api/shutdown", { method: "POST" });
}

export async function getEntryFiles(entry: string) {
  return (await request<{ files: EntryFile[] }>(`/api/entries/${encodeURIComponent(entry)}/files`)).files;
}

export async function revealEntryPath(entry: string, subpath: string) {
  return await request<{ ok: true; path: string }>(`/api/entries/${encodeURIComponent(entry)}/reveal`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ subpath }),
  });
}

export async function uploadRunnerFile(file: File, onProgress?: (progress: number) => void) {
  // XMLHttpRequest exposes upload progress; fetch still does not do so consistently.
  return await new Promise<{ path: string; name: string; size: number }>((resolve, reject) => {
    const xhr = new XMLHttpRequest();
    xhr.open("POST", `${runnerUrl}/api/uploads`);
    xhr.setRequestHeader("Content-Type", "application/octet-stream");
    xhr.setRequestHeader("X-File-Name", file.name);
    xhr.upload.onprogress = (event) => {
      if (event.lengthComputable) onProgress?.(Math.round((event.loaded / event.total) * 100));
    };
    xhr.onerror = () => reject(new Error("Could not reach the local runner"));
    xhr.onload = () => {
      let payload: { path?: string; name?: string; size?: number; error?: string } = {};
      try { payload = JSON.parse(xhr.responseText); } catch { /* handled below */ }
      if (xhr.status < 200 || xhr.status >= 300 || !payload.path) return reject(new Error(payload.error || "Upload failed"));
      resolve({ path: payload.path, name: payload.name || file.name, size: payload.size || file.size });
    };
    xhr.send(file);
  });
}
