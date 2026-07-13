"""Pytest bootstrap for the SSE test suite.

SSE is currently run from the repo root rather than being pip-installed, so the
`sse_tools` package isn't on `sys.path` by default when pytest collects tests
from `tests/`. Adding the repo root here lets the tests do
`from sse_tools import common` exactly the way the scripts in `scripts/` do.

`scripts/` is added too so pipeline tests can `import sse_initialization` and
exercise `build_entry` (the creation contract) directly, the same way that
script imports itself when run.

(If SSE ever grows a real `pyproject.toml` and becomes `pip install -e .`, this
file can go away.)
"""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
for _p in (REPO_ROOT, REPO_ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))
