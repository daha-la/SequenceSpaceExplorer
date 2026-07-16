#!/bin/zsh

set -u

APP_DIR="${0:A:h}"
ROOT_DIR="${APP_DIR:h}"
WORK_DIR="${APP_DIR}/work"
SHUTDOWN_FILE="${WORK_DIR}/shutdown.request"
APP_URL="http://localhost:3000/"
BUNDLED_NODE="/Users/wpgaard/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin"
BUNDLED_TOOLS="/Users/wpgaard/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/fallback"
RUNNER_PID=""
UI_PID=""

mkdir -p "$WORK_DIR"
rm -f "$SHUTDOWN_FILE"

if [[ -x "${BUNDLED_NODE}/node" && -x "${BUNDLED_TOOLS}/pnpm" ]]; then
  export PATH="${BUNDLED_NODE}:${BUNDLED_TOOLS}:${PATH}"
fi

pause_for_error() {
  echo
  echo "$1"
  echo
  read -k 1 "?Press any key to close this window..."
  echo
  exit 1
}

command -v node >/dev/null 2>&1 || pause_for_error "Node.js was not found. Open this project in Codex once so its bundled runtime is available, or install Node.js."
command -v pnpm >/dev/null 2>&1 || pause_for_error "pnpm was not found. Open this project in Codex once so its bundled runtime is available, or install pnpm."

if [[ -x "${ROOT_DIR}/.venv/bin/python" ]]; then
  export SSE_PYTHON="${ROOT_DIR}/.venv/bin/python"
elif [[ -x "${ROOT_DIR}/venv/bin/python" ]]; then
  export SSE_PYTHON="${ROOT_DIR}/venv/bin/python"
fi

export SSE_SHUTDOWN_FILE="$SHUTDOWN_FILE"
export PORT=3000

stop_process() {
  local pid="$1"
  [[ -n "$pid" ]] || return 0
  kill -0 "$pid" >/dev/null 2>&1 || return 0
  kill -TERM "$pid" >/dev/null 2>&1 || true
}

cleanup() {
  trap - EXIT INT TERM
  stop_process "$UI_PID"
  stop_process "$RUNNER_PID"
  [[ -n "$UI_PID" ]] && wait "$UI_PID" >/dev/null 2>&1 || true
  [[ -n "$RUNNER_PID" ]] && wait "$RUNNER_PID" >/dev/null 2>&1 || true
  rm -f "$SHUTDOWN_FILE"
}

trap cleanup EXIT INT TERM
cd "$APP_DIR" || exit 1

echo "Starting Sequence Space Explorer..."
echo "This window supervises the app and may be minimized."
echo

pnpm run runner >"${WORK_DIR}/runner.log" 2>&1 &
RUNNER_PID=$!
pnpm run dev >"${WORK_DIR}/ui.log" 2>&1 &
UI_PID=$!

ready=false
for _ in {1..60}; do
  if ! kill -0 "$RUNNER_PID" >/dev/null 2>&1; then
    tail -n 20 "${WORK_DIR}/runner.log" 2>/dev/null || true
    pause_for_error "The pipeline runner could not start."
  fi
  if ! kill -0 "$UI_PID" >/dev/null 2>&1; then
    tail -n 20 "${WORK_DIR}/ui.log" 2>/dev/null || true
    pause_for_error "The pipeline interface could not start."
  fi
  if curl -fsS "http://127.0.0.1:8788/api/health" >/dev/null 2>&1 && curl -fsS "$APP_URL" >/dev/null 2>&1; then
    ready=true
    break
  fi
  sleep 1
done

if [[ "$ready" != true ]]; then
  pause_for_error "The app did not become ready within one minute. Logs are available in pipeline-ui-poc/work/."
fi

if ! /usr/bin/open "$APP_URL"; then
  pause_for_error "The app started, but macOS could not open the default browser. Open ${APP_URL} manually."
fi
echo "Sequence Space Explorer is running."
echo "Use the Shut down button in the app, or close this Terminal window, to stop it."

while kill -0 "$RUNNER_PID" >/dev/null 2>&1 && kill -0 "$UI_PID" >/dev/null 2>&1; do
  [[ -f "$SHUTDOWN_FILE" ]] && break
  sleep 1
done

echo
echo "Shutting down Sequence Space Explorer..."
cleanup
echo "Shutdown complete."
