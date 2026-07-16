# Sequence Space Explorer pipeline control center

This app configures, launches, and monitors the Sequence Space Explorer scripts without requiring day-to-day terminal work.

It now provides:

- live discovery and summary statistics for entries in `../entries/`;
- validated execution of all eight command-line entry points, including the Boltz-2 structure-prediction + RMSD module;
- file uploads for initialization, external data, translators, and Foldseek JSON;
- persistent job history, live stdout/stderr, cancellation, and exit status;
- per-entry serialization for commands that write to SSE data;
- real explorer startup with a link to the local Dash application.

## Local development

### Double-click launch on macOS

Double-click **Start Sequence Space Explorer.command** in Finder. It starts the runner and interface, waits until both are healthy, opens `http://localhost:3000/` in the default browser, and keeps both services supervised in one Terminal window. Use **Shut down** in the app—or close that Terminal window—to stop both services and any active pipeline process.

Starting the explorer from the pipeline opens a waiting browser tab immediately. When the Dash visualizer is ready, that tab automatically navigates to its local address (normally `http://127.0.0.1:8051/`). If the browser blocks the tab, allow pop-ups for the local pipeline app and press **Start explorer** again.

The first launch may show a macOS security prompt. If macOS blocks the file, Control-click it in Finder, choose **Open**, and confirm once.

### Manual launch

The UI and the local pipeline runner are separate processes. Run both from this directory:

```bash
pnpm run runner
pnpm run dev
```

Then open the local URL printed by the development server. The runner listens only on `127.0.0.1:8788` by default.

The runner uses `SSE_PYTHON` when set, then checks `../.venv/bin/python` and `../venv/bin/python`, and finally falls back to `python3`. For example:

```bash
SSE_PYTHON=/path/to/your/sse/environment/bin/python pnpm run runner
```

Optional settings:

- `SSE_RUNNER_PORT` changes the runner port.
- `SSE_UI_ORIGINS` is a comma-separated allowlist of UI origins.
- `NEXT_PUBLIC_SSE_RUNNER_URL` points the UI at a non-default runner URL.

Job metadata is stored in the ignored `work/jobs.json` file. Uploaded inputs are stored in `work/uploads/`.

## Production build

```bash
pnpm run build
```

The hosted site remains a browser interface; pipeline execution intentionally stays on the machine that owns the datasets, Python environment, models, and GPU.
