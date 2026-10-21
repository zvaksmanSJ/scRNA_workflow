#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" != "--submit" ]]; then
  echo "Usage: bash python_lsf_pipeline/submit_fastq_audit.sh --submit --manifest TSV --raw-plan JSON --audit-dir DIR [--queue NAME]" >&2
  exit 2
fi
shift
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
WORKFLOW_ROOT=$(cd -- "$HERE/.." && pwd)
RUNTIME_DIR="${RNA_RUNTIME_DIR:-$WORKFLOW_ROOT/wdl_runtime}"
PYTHON="${RNA_PYTHON:-$RUNTIME_DIR/env/resolveome-rna/bin/python}"
MANIFEST=""
RAW_PLAN=""
AUDIT_DIR=""
QUEUE="${LSF_QUEUE:-standard}"
while (($#)); do
  case "$1" in
    --manifest) MANIFEST="${2:?Missing manifest path}"; shift 2 ;;
    --raw-plan) RAW_PLAN="${2:?Missing raw plan path}"; shift 2 ;;
    --audit-dir) AUDIT_DIR="${2:?Missing audit output directory}"; shift 2 ;;
    --queue) QUEUE="${2:?Missing LSF queue}"; shift 2 ;;
    *) echo "[ERROR] Unknown option: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$AUDIT_DIR" && -s "$MANIFEST" && -s "$RAW_PLAN" && -x "$PYTHON" ]] || {
  echo "[ERROR] Manifest, raw plan, and RNA Python must exist: $MANIFEST $RAW_PLAN $PYTHON" >&2
  exit 2
}
command -v bsub >/dev/null || {
  echo "[ERROR] bsub is required for LSF array submission" >&2
  exit 1
}
AUDIT_DIR="$(mkdir -p "$AUDIT_DIR" && cd "$AUDIT_DIR" && pwd)"
LOG_DIR="$AUDIT_DIR/logs"

sample_count=$("$PYTHON" -c \
  'import csv,sys; print(sum(1 for _ in csv.DictReader(open(sys.argv[1]), delimiter="\t")))' \
  "$MANIFEST")
[[ "$sample_count" -gt 0 ]] || { echo "No samples in $MANIFEST" >&2; exit 2; }
mkdir -p "$AUDIT_DIR" "$LOG_DIR"
bsub -J "rna_expanded_fastq_audit[1-$sample_count]" \
  -q "$QUEUE" -n 1 -M 8000 -R 'rusage[mem=8000]' -W 24:00 \
  -cwd "$WORKFLOW_ROOT" \
  -oo "$LOG_DIR/audit.%J.%I.out" -eo "$LOG_DIR/audit.%J.%I.err" \
  bash "$HERE/audit_fastq_worker.sh" "$MANIFEST" "$RAW_PLAN" "$AUDIT_DIR" \
  | tee -a "$LOG_DIR/submission.log"
