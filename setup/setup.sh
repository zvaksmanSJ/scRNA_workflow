#!/usr/bin/env bash
set -euo pipefail

SETUP_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd -- "$SETUP_DIR/.." && pwd)"
RUNTIME_DIR="${RNA_RUNTIME_DIR:-$WORKFLOW_ROOT/wdl_runtime}"
THREADS=8
MEMORY_MB=64000
SKIP_REFERENCES=0
INSTALL_RESEARCH_ADDONS=0

usage() {
  cat <<'EOF'
Set up the portable ResolveOME RNA runtime.

Usage: bash setup/setup.sh [options]

Options:
  --runtime-dir DIR   Install environments and references in DIR
                      (default: workflow_RNA/wdl_runtime)
  --threads N         Threads for index generation (default: 8)
  --memory-mb N       STAR index memory limit in MB (default: 64000)
  --skip-references   Install software only; do not download/build references
  --research-addons  Also install the optional RNA-only research environment
  -h, --help          Show this help

Prerequisites: Linux, conda or mamba, curl, gzip, tar, sha256sum, and
approximately 100 GB free disk plus at least 64 GB RAM when building STAR.
EOF
}

while (($#)); do
  case "$1" in
    --runtime-dir) RUNTIME_DIR="${2:?Missing directory after --runtime-dir}"; shift 2 ;;
    --threads) THREADS="${2:?Missing count after --threads}"; shift 2 ;;
    --memory-mb) MEMORY_MB="${2:?Missing amount after --memory-mb}"; shift 2 ;;
    --skip-references) SKIP_REFERENCES=1; shift ;;
    --research-addons) INSTALL_RESEARCH_ADDONS=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if ! [[ "$THREADS" =~ ^[1-9][0-9]*$ && "$MEMORY_MB" =~ ^[1-9][0-9]*$ ]]; then
  echo "[ERROR] --threads and --memory-mb must be positive integers" >&2
  exit 2
fi
if ((MEMORY_MB < 48000)); then
  echo "[ERROR] STAR index generation requires at least 48000 MB" >&2
  exit 2
fi

CONDA_BIN="${CONDA_EXE:-}"
if [[ -z "$CONDA_BIN" ]]; then
  CONDA_BIN="$(command -v conda || command -v mamba || true)"
fi
if [[ -z "$CONDA_BIN" ]]; then
  echo "[ERROR] Install Miniforge/Miniconda or Mamba and put conda/mamba on PATH" >&2
  exit 1
fi
for tool in curl gzip tar sha256sum; do
  command -v "$tool" >/dev/null || {
    echo "[ERROR] Required command not found: $tool" >&2
    exit 1
  }
done

RUNTIME_DIR="$(mkdir -p "$RUNTIME_DIR" && cd "$RUNTIME_DIR" && pwd)"
ENV_PREFIX="$RUNTIME_DIR/env/resolveome-rna"
SALMON_PREFIX="$RUNTIME_DIR/env/resolveome-rna-salmon"
RESEARCH_PREFIX="$RUNTIME_DIR/env/resolveome-rna-research"
mkdir -p "$RUNTIME_DIR/env" "$RUNTIME_DIR/bin" "$RUNTIME_DIR/packages"

create_env() {
  local prefix="$1" definition="$2"
  if [[ -x "$prefix/bin/python" || -x "$prefix/bin/salmon" ]]; then
    echo "[INFO] Environment already exists: $prefix"
  elif [[ -e "$prefix" ]]; then
    echo "[ERROR] Incomplete environment exists at $prefix; inspect it before removing it" >&2
    exit 1
  else
    "$CONDA_BIN" env create --yes --prefix "$prefix" --file "$definition"
  fi
}

create_env "$ENV_PREFIX" "$SETUP_DIR/environment.yml"
create_env "$SALMON_PREFIX" "$SETUP_DIR/analysis_salmon_environment.yml"
if ((INSTALL_RESEARCH_ADDONS)); then
  create_env "$RESEARCH_PREFIX" "$SETUP_DIR/research_addons_environment.yml"
fi

download_pinned() {
  local url="$1" destination="$2" expected_sha256="$3"
  if [[ ! -s "$destination" ]]; then
    curl --fail --location --retry 3 --output "$destination.part" "$url"
    mv "$destination.part" "$destination"
  fi
  local actual_sha256
  actual_sha256="$(sha256sum "$destination" | awk '{print $1}')"
  [[ "$expected_sha256" == "$actual_sha256" ]] || {
    echo "[ERROR] SHA-256 mismatch: $destination" >&2
    exit 1
  }
}

AMOR="$RUNTIME_DIR/packages/AMOR-0.1-2.tar.gz"
OHCHIBI="$RUNTIME_DIR/packages/ohchibi-9d3dac9.tar.gz"
HARRIETR="$RUNTIME_DIR/packages/harrietr-0.2.3.tar.gz"
download_pinned \
  "https://codeload.github.com/surh/AMOR/tar.gz/refs/tags/0.1-2" "$AMOR" \
  "dabd53a84524b246edf97ef9060e21658d437f53c11ff22631db921a4660a773"
download_pinned \
  "https://codeload.github.com/isaisg/ohchibi/tar.gz/9d3dac92081993508facc58877fdf0cfa8bba77c" "$OHCHIBI" \
  "61f8c7e007d519736e67cde949ce49054aee42a1d09bc893f3de615e1e93cc73"
download_pinned \
  "https://cloud.r-project.org/src/contrib/harrietr_0.2.3.tar.gz" "$HARRIETR" \
  "920179cba12072872eba6ec1de423889931f3a92155bf30d803212a93dd6db58"
REACTOME="$WORKFLOW_ROOT/analysis_resources/ReactomePathways.gmt.zip"
download_pinned \
  "https://reactome.org/download/current/ReactomePathways.gmt.zip" "$REACTOME" \
  "8c1dbc8578431da5d2d5118262718c60b553a9be3398e93658daa069e4a9afd4"

"$ENV_PREFIX/bin/R" CMD INSTALL "$AMOR" "$OHCHIBI" "$HARRIETR"
"$ENV_PREFIX/bin/Rscript" -e \
  'stopifnot(as.character(packageVersion("AMOR")) == "0.1.2", as.character(packageVersion("harrietr")) == "0.2.3", requireNamespace("ohchibi", quietly=TRUE))'
"$ENV_PREFIX/bin/STAR" --version
"$ENV_PREFIX/bin/fastp" --version
"$ENV_PREFIX/bin/seqkit" version
"$ENV_PREFIX/bin/samtools" --version | head -1
"$SALMON_PREFIX/bin/salmon" --version
if ((INSTALL_RESEARCH_ADDONS)); then
  for tool in STAR arriba suppa.py dexseq_prepare_annotation.py dexseq_count.py tin.py read_duplication.py; do
    if [[ ! -x "$RESEARCH_PREFIX/bin/$tool" ]]; then
      echo "[ERROR] Research-add-on tool missing after install: $tool" >&2
      exit 1
    fi
  done
fi

NXF_VER=26.04.6
export NXF_VER
if [[ ! -x "$RUNTIME_DIR/bin/nextflow" ]]; then
  (cd "$RUNTIME_DIR/bin" && curl --fail --location --retry 3 https://get.nextflow.io | bash)
fi
PATH="$ENV_PREFIX/bin:$PATH" "$RUNTIME_DIR/bin/nextflow" -version

if ((SKIP_REFERENCES == 0)); then
  RNA_RUNTIME_DIR="$RUNTIME_DIR" bash "$SETUP_DIR/prepare_references.sh" \
    --threads "$THREADS" --memory-mb "$MEMORY_MB"
fi

CONDA_BASE="$("$CONDA_BIN" info --base)"
printf '%s\n' "$CONDA_BASE" > "$RUNTIME_DIR/conda_base.txt"
"$CONDA_BIN" list --prefix "$ENV_PREFIX" --explicit > "$RUNTIME_DIR/environment.resolveome-rna.explicit.txt"
"$CONDA_BIN" list --prefix "$SALMON_PREFIX" --explicit > "$RUNTIME_DIR/environment.salmon.explicit.txt"
if ((INSTALL_RESEARCH_ADDONS)); then
  "$CONDA_BIN" list --prefix "$RESEARCH_PREFIX" --explicit \
    > "$RUNTIME_DIR/environment.research-addons.explicit.txt"
fi
python_bin="$ENV_PREFIX/bin/python"
"$python_bin" - "$RUNTIME_DIR/runtime_config.json" "$CONDA_BASE" "$ENV_PREFIX" \
  "$SALMON_PREFIX" "$RUNTIME_DIR/bin/nextflow" "$RUNTIME_DIR/references/Ensembl104" <<'PY'
import json
import sys
from pathlib import Path

output, conda_base, rna_env, salmon_env, nextflow, refs = sys.argv[1:]
refs = Path(refs)
config = {
    "conda_base": conda_base,
    "rna_environment": rna_env,
    "salmon_environment": salmon_env,
    "research_addons_environment": (
        str(Path(rna_env).parent / "resolveome-rna-research")
        if Path(Path(rna_env).parent / "resolveome-rna-research").is_dir()
        else None
    ),
    "nextflow": nextflow,
    "references": {
        "fasta": str(refs / "GRCh38.primary.fa"),
        "gtf": str(refs / "Homo_sapiens.GRCh38.104.chr.gtf"),
        "tx2gene": str(refs / "Ensembl104.tx2gene.tsv"),
        "genebody": str(refs / "Ensembl104.transcripts.bed"),
        "star_index": str(refs / "STAR"),
        "ref_flat": str(refs / "Ensembl104.refFlat.txt"),
        "salmon_index": str(refs / "salmon_index"),
    },
}
Path(output).write_text(json.dumps(config, indent=2) + "\n")
PY
printf '[PASS] RNA runtime installed under %s\n' "$RUNTIME_DIR"
