#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_ROOT=$(cd -- "$HERE/.." && pwd)
RUNTIME_DIR="${RNA_RUNTIME_DIR:-$WORKFLOW_ROOT/wdl_runtime}"
PYTHON="${RNA_PYTHON:-$RUNTIME_DIR/env/resolveome-rna/bin/python}"
[[ -x "$PYTHON" ]] || {
  echo "[ERROR] RNA Python is missing or not executable: $PYTHON" >&2
  exit 1
}
exec "$PYTHON" \
  "$HERE/audit_expanded_fastqs.py" \
  --manifest "${1:?RNA sample manifest required}" \
  --raw-plan "${2:?Raw input plan required}" \
  --output-dir "${3:?Per-sample audit directory required}" \
  --index "${LSB_JOBINDEX:?Must run as an LSF array task}"
