#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_ROOT=$(cd -- "$HERE/.." && pwd)
RUNTIME_DIR="${RNA_RUNTIME_DIR:-$WORKFLOW_ROOT/wdl_runtime}"
PREFIX="${RNA_ENV_PREFIX:-$RUNTIME_DIR/env/resolveome-rna}"
CONDA_BASE="${CONDA_BASE:-$(cat "$RUNTIME_DIR/conda_base.txt" 2>/dev/null || true)}"
if [[ -z "$CONDA_BASE" || ! -s "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
  echo "[ERROR] Set CONDA_BASE or run setup/setup.sh first (runtime: $RUNTIME_DIR)" >&2
  exit 1
fi
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$PREFIX"
export R_LIBS_USER="$PREFIX/lib/R/library"
if [[ "${1:-}" == "--submit" ]]; then
  shift
  mkdir -p "$WORKFLOW_ROOT/logs"
  bsub -J resolveome_rna_driver -q "${LSF_QUEUE:-standard}" -n 2 -M 16000 \
    -R 'rusage[mem=16000] span[hosts=1]' -W 168:00 -cwd "$WORKFLOW_ROOT" \
    -oo "$WORKFLOW_ROOT/logs/driver.%J.out" -eo "$WORKFLOW_ROOT/logs/driver.%J.err" \
    bash "$HERE/run_all_resolveome_rna.sh" --run "$@"
else
  exec "$PREFIX/bin/python" "$HERE/run_all_resolveome_rna.py" "$@"
fi
