#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# Uses the prepared RNA manifest, GRCh38/Ensembl104, both depths, male-independent RNA
# quantification, unstranded counting, 8 CPUs and 64000 MB worker cap.
bash "$HERE/run_all_resolveome_rna.sh" --submit
