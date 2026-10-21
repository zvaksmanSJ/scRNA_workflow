#!/usr/bin/env bash
set -euo pipefail

SETUP_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd -- "$SETUP_DIR/.." && pwd)"
RUNTIME_DIR="${RNA_RUNTIME_DIR:-$WORKFLOW_ROOT/wdl_runtime}"
mkdir -p "$RUNTIME_DIR"
RUNTIME_DIR="$(cd -- "$RUNTIME_DIR" && pwd)"
THREADS=8
MEMORY_MB=64000

usage() {
  cat <<'EOF'
Download the public Ensembl 104 reference inputs and build RNA indices.

Usage: bash setup/prepare_references.sh [--threads N] [--memory-mb N]

This builds a STAR genome index and Salmon transcript index. STAR index
generation needs at least 48 GB RAM; increase --memory-mb to match available
memory. Files are placed under wdl_runtime/references/Ensembl104 by default.
EOF
}

while (($#)); do
  case "$1" in
    --threads) THREADS="${2:?Missing count after --threads}"; shift 2 ;;
    --memory-mb) MEMORY_MB="${2:?Missing amount after --memory-mb}"; shift 2 ;;
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

CONDA_BIN="${CONDA_EXE:-$(command -v conda || command -v mamba || true)}"
if [[ -z "$CONDA_BIN" ]]; then
  echo "[ERROR] conda or mamba is required; run setup/setup.sh first" >&2
  exit 1
fi
R_ENV="$RUNTIME_DIR/env/resolveome-rna"
SALMON_ENV="$RUNTIME_DIR/env/resolveome-rna-salmon"
[[ -x "$R_ENV/bin/STAR" && -x "$R_ENV/bin/samtools" ]] || {
  echo "[ERROR] RNA tool environment is missing; run setup/setup.sh first" >&2
  exit 1
}
[[ -x "$SALMON_ENV/bin/salmon" ]] || {
  echo "[ERROR] Salmon environment is missing; run setup/setup.sh first" >&2
  exit 1
}

REF_DIR="$RUNTIME_DIR/references/Ensembl104"
mkdir -p "$REF_DIR/downloads"

download_ensembl() {
  local group="$1" name="$2" destination="$3"
  local base="https://ftp.ensembl.org/pub/release-104"
  local url="$base/$group/$name"
  local checksums="$REF_DIR/downloads/${group##*/}.CHECKSUMS"
  if [[ ! -s "$checksums" ]]; then
    curl --fail --location --retry 3 \
      --output "$checksums.part" "$base/$group/CHECKSUMS"
    mv "$checksums.part" "$checksums"
  fi
  local expected expected_crc expected_size actual_crc actual_size actual_name
  expected="$(awk -v file="$name" '$3 == file {print $1, $2}' "$checksums")"
  if [[ -z "$expected" ]]; then
    echo "[ERROR] No Ensembl CHECKSUMS entry for $name" >&2
    return 1
  fi
  read -r expected_crc expected_size <<< "$expected"
  if [[ ! -s "$destination" ]]; then
    curl --fail --location --retry 3 \
      --output "$destination.part" "$url"
    read -r actual_crc actual_size actual_name < <(cksum "$destination.part")
    if [[ "$expected_crc" != "$actual_crc" || "$expected_size" != "$actual_size" ]]; then
      rm -f "$destination.part"
      echo "[ERROR] Ensembl checksum mismatch for $name" >&2
      return 1
    fi
    mv "$destination.part" "$destination"
  else
    read -r actual_crc actual_size actual_name < <(cksum "$destination")
    if [[ "$expected_crc" != "$actual_crc" || "$expected_size" != "$actual_size" ]]; then
      echo "[ERROR] Existing Ensembl file failed checksum validation: $destination" >&2
      return 1
    fi
  fi
  gzip --test "$destination"
  echo "[OK] $destination"
}

FASTA_GZ="$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_GZ="$REF_DIR/Homo_sapiens.GRCh38.104.chr.gtf.gz"
CDNA_GZ="$REF_DIR/Homo_sapiens.GRCh38.cdna.all.fa.gz"
NCRNA_GZ="$REF_DIR/Homo_sapiens.GRCh38.ncrna.fa.gz"
download_ensembl "fasta/homo_sapiens/dna" \
  "$(basename "$FASTA_GZ")" "$FASTA_GZ"
download_ensembl "gtf/homo_sapiens" \
  "$(basename "$GTF_GZ")" "$GTF_GZ"
download_ensembl "fasta/homo_sapiens/cdna" \
  "$(basename "$CDNA_GZ")" "$CDNA_GZ"
download_ensembl "fasta/homo_sapiens/ncrna" \
  "$(basename "$NCRNA_GZ")" "$NCRNA_GZ"

FASTA="$REF_DIR/GRCh38.primary.fa"
GTF="$REF_DIR/Homo_sapiens.GRCh38.104.chr.gtf"
TRANSCRIPTS="$REF_DIR/Ensembl104.transcripts.fa"
if [[ ! -s "$FASTA" ]]; then
  gzip --decompress --stdout "$FASTA_GZ" > "$FASTA.part"
  mv "$FASTA.part" "$FASTA"
fi
if [[ ! -s "$GTF" ]]; then
  gzip --decompress --stdout "$GTF_GZ" > "$GTF.part"
  mv "$GTF.part" "$GTF"
fi
"$R_ENV/bin/samtools" faidx "$FASTA"
RNA_RUNTIME_DIR="$RUNTIME_DIR" "$R_ENV/bin/python" \
  "$SETUP_DIR/step2_prepare_references.py" \
  --fasta "$FASTA" --gtf "$GTF" \
  --tx2gene "$REF_DIR/Ensembl104.tx2gene.tsv" \
  --genebody "$REF_DIR/Ensembl104.transcripts.bed" \
  --star-index "$REF_DIR/STAR" \
  --ref-flat "$REF_DIR/Ensembl104.refFlat.txt" \
  --threads "$THREADS" --memory-mb "$MEMORY_MB" \
  --derive-annotation --build

if [[ ! -s "$TRANSCRIPTS" ]]; then
  "$R_ENV/bin/python" "$SETUP_DIR/step2_prepare_references.py" \
    --combine-transcripts --cdna-fasta "$CDNA_GZ" --ncrna-fasta "$NCRNA_GZ" \
    --transcript-fasta "$TRANSCRIPTS"
fi

SALMON_INDEX="$REF_DIR/salmon_index"
if [[ -s "$SALMON_INDEX/versionInfo.json" ]]; then
  echo "[OK] Existing Salmon index retained: $SALMON_INDEX"
elif [[ -e "$SALMON_INDEX" ]]; then
  echo "[ERROR] Incomplete Salmon index at $SALMON_INDEX; inspect before removing" >&2
  exit 1
else
  "$SALMON_ENV/bin/salmon" index -t "$TRANSCRIPTS" -i "$SALMON_INDEX" \
    -k 31 -p "$THREADS"
fi

printf '[PASS] Public Ensembl 104 reference bundle ready at %s\n' "$REF_DIR"
