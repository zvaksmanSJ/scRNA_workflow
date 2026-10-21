#!/usr/bin/env bash
set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RNAW="$(cd "$HERE/.." && pwd)"
RUNTIME_DIR="${RNA_RUNTIME_DIR:-$RNAW/wdl_runtime}"
PREFIX="${RNA_ENV_PREFIX:-$RUNTIME_DIR/env/resolveome-rna}"
RNAR="${RESOLVEOME_RNA_RESULTS:?Set RESOLVEOME_RNA_RESULTS to a new, completed run directory}"
CONDA_BASE="${CONDA_BASE:-$(cat "$RUNTIME_DIR/conda_base.txt" 2>/dev/null || true)}"
STATUS="$RNAR/final_results/07_qc_provenance/downstream_run_status.txt"

record_failure() {
    local exit_code=$?
    printf 'FAILED %s exit_code=%s\n' \
        "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$exit_code" > "$STATUS"
    exit "$exit_code"
}
trap record_failure ERR

printf 'RUNNING %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$STATUS"
cd "$RNAW"
[[ -s "$CONDA_BASE/etc/profile.d/conda.sh" ]] || {
    echo "[ERROR] Conda activation script missing; set CONDA_BASE or run setup/setup.sh" >&2
    exit 1
}
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$PREFIX"
python extra_methods/quantify_transcripts.py --validate
python extra_methods/run_rna_downstream.py
printf 'PASS %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$STATUS"
