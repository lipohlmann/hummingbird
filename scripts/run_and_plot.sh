#!/usr/bin/env bash
# Run a hummingbird case, then plot its scalar and angular flux with plot_csv.py.
# The results CSV is written next to the case JSON and the plots go in
# <case_dir>/plots/. See scripts/README.md for details.
set -euo pipefail

if [[ $# -ne 1 || ! -f $1 ]]; then
  echo "Usage: $0 <case.json>" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/../build/release/hummingbird"
if [[ ! -x $EXE ]]; then
  echo "Error: $EXE not found; build the project first (see README.md)." >&2
  exit 1
fi

CASE_DIR="$(cd "$(dirname "$1")" && pwd)"
CASE_JSON="$CASE_DIR/$(basename "$1")"

# The solver writes <problem.name>_results.csv to the cwd, so run from the case dir.
cd "$CASE_DIR"
NAME="$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["problem"]["name"])' "$CASE_JSON")"
CSV="${NAME}_results.csv"

"$EXE" "$CASE_JSON"

# plot_csv.py needs pandas/seaborn, which pixi's python lacks, so use system python3.
mkdir -p plots
python3 "$SCRIPT_DIR/plot_csv.py" "$CSV" --angular --output-dir plots
echo "Plots written to $CASE_DIR/plots"
