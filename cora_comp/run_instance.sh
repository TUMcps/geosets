#!/bin/bash
# CORA-COMP: run one instance (timed). Args: v1 <benchmark> <instance> <params-json> ... <results-file>
set -e
[ "$1" = v1 ] || { echo "Expected interface version 'v1', got '$1'"; exit 1; }
DIR="$(dirname "$0")"
PARAMS="$4"
RESULTS_FILE="${@: -1}"

case "$PARAMS" in
    *'"device": "gpu"'*) printf 'result\nunsupported\n' > "$RESULTS_FILE"; exit 0 ;;
esac
exec "$DIR/../.cora-venv/bin/python" "$DIR/run.py" "$PARAMS" "$RESULTS_FILE"
