#!/usr/bin/env bash
# Thin wrapper around convergence_study.py, run inside the pixi environment
# so numpy/matplotlib/gmsh are on PATH. See scripts/README.md for usage.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec pixi run python3 "$SCRIPT_DIR/convergence_study.py" "$@"
