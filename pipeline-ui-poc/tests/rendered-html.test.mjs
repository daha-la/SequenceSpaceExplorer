import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

async function render() {
  const workerUrl = new URL("../dist/server/index.js", import.meta.url);
  workerUrl.searchParams.set("test", `${process.pid}-${Date.now()}`);
  const { default: worker } = await import(workerUrl.href);
  return worker.fetch(
    new Request("http://localhost/", { headers: { accept: "text/html" } }),
    { ASSETS: { fetch: async () => new Response("Not found", { status: 404 }) } },
    { waitUntil() {}, passThroughOnException() {} },
  );
}

test("server-renders the Sequence Space Explorer application", async () => {
  const response = await render();
  assert.equal(response.status, 200);
  assert.match(response.headers.get("content-type") ?? "", /^text\/html\b/i);
  const html = await response.text();
  assert.match(html, /<title>Sequence Space Explorer · Pipeline Control Center<\/title>/i);
  assert.match(html, /Protein sequence intelligence/);
  assert.match(html, /Enter pipeline/);
  assert.doesNotMatch(html, /Your site is taking shape|Codex is working/);
});

test("wires real runner operations instead of simulated progress", async () => {
  const [page, runner, commands, launcher] = await Promise.all([
    readFile(new URL("../app/page.tsx", import.meta.url), "utf8"),
    readFile(new URL("../runner/server.mjs", import.meta.url), "utf8"),
    readFile(new URL("../runner/commands.mjs", import.meta.url), "utf8"),
    readFile(new URL("../Start Sequence Space Explorer.command", import.meta.url), "utf8"),
  ]);
  assert.match(page, /submitJob/);
  assert.match(page, /getJob/);
  assert.match(page, /cancelPipelineJob/);
  assert.match(page, /requestPipelineShutdown/);
  assert.match(page, /Shut down/);
  assert.match(page, /window\.open\("", "_blank"\)/);
  assert.match(page, /pendingWindow\.location\.replace/);
  assert.doesNotMatch(page, /jobLogs|simulates execution|production version would submit/i);
  assert.match(runner, /spawn\(pythonExecutable\(\), job\.args/);
  assert.match(runner, /writesEntry/);
  assert.match(runner, /\/api\/shutdown/);
  assert.match(commands, /WRITER_TOOLS/);
  assert.match(launcher, /pnpm run runner/);
  assert.match(launcher, /pnpm run dev/);
  assert.match(launcher, /SSE_SHUTDOWN_FILE/);
  assert.match(launcher, /APP_URL="http:\/\/localhost:3000\/"/);
  assert.match(launcher, /\/usr\/bin\/open "\$APP_URL"/);
});
