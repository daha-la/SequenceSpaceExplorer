#!/usr/bin/env python3
"""Development-sized coordinate run using the first 25 sequences only.

This invokes the normal sse_coordinates pipeline, but gives its embedding
cache, coordinate columns, figure, and log a ``_first25`` suffix so this small
cache cannot be confused with a complete entry embedding.
"""

import sys

from sse_coordinates import main


if __name__ == "__main__":
    sys.exit(main([*sys.argv[1:], "--limit", "25"]))
