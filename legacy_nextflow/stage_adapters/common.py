"""Shared import configuration for the stage scripts copied into a run tree."""
import os
from pathlib import Path
import sys

WORKFLOW = Path(os.environ.get(
    "RESOLVEOME_RNA_WORKFLOW",
    Path(__file__).resolve().parents[2],
)).resolve()
EXECUTOR = os.environ.get("RNA_NEXTFLOW_EXECUTOR", "lsf")
QUEUE = os.environ.get("RNA_LSF_QUEUE", "standard")
THREADS = int(os.environ.get("RNA_NEXTFLOW_THREADS", "8"))
MEMORY_MB = int(os.environ.get("RNA_NEXTFLOW_MEMORY_MB", "64000"))

for directory in ("", "python_lsf_pipeline", "setup", "legacy_nextflow", "extra_methods"):
    path = str(WORKFLOW / directory)
    if path not in sys.path:
        sys.path.insert(0, path)
