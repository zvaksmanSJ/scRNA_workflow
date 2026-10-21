version 1.0

# One independent full-depth transcript-abundance task per sample. The WDL
# scatters this task, pins the Salmon executable/index contract, and returns
# quant.sf plus its audit metadata; HTSeq remains the canonical gene-count path.
task SalmonQuant {
  input {
    String sample_id
    String results_dir
    File full_preprocessing_gate
    String salmon_executable
    String salmon_index
    Int threads
    Int memory_mb
    String queue
  }

  command <<<
    set -euo pipefail
    python - "~{full_preprocessing_gate}" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    report = json.load(stream)
if report.get("status") != "PASS" or report.get("depth") != "full":
    raise SystemExit("Full-depth preprocessing audit gate failed")
PY
    if ! [[ "~{sample_id}" =~ ^[A-Za-z0-9_.-]+$ ]]; then
      echo "[ERROR] Unsafe sample ID: ~{sample_id}" >&2
      exit 22
    fi
    r1="~{results_dir}/stage3_full_preprocessing/trimmed/~{sample_id}_R1_trim.fastq.gz"
    r2="~{results_dir}/stage3_full_preprocessing/trimmed/~{sample_id}_R2_trim.fastq.gz"
    out="~{results_dir}/final_results/02_transcript/quant_by_sample/~{sample_id}"
    test -s "$r1"
    test -s "$r2"
    if [ -e "$out/quant.sf" ]; then
      echo "[ERROR] Salmon output already exists; refusing to overwrite: $out/quant.sf" >&2
      exit 23
    fi
    mkdir -p "$out"
    "~{salmon_executable}" quant -i "~{salmon_index}" -l A \
      -1 "$r1" -2 "$r2" -p ~{threads} --validateMappings \
      --seqBias --gcBias -o "$out" > salmon.log 2>&1
    test -s "$out/quant.sf"
    cp "$out/quant.sf" quant.sf
    mkdir -p "~{results_dir}/logs_driver"
    cp salmon.log "~{results_dir}/logs_driver/wdl_salmon_~{sample_id}.log"
  >>>

  output {
    File quant_sf = "quant.sf"
    File log = "salmon.log"
  }

  runtime {
    cpu: threads
    memory: "64 GB"
    memory_mb: memory_mb
    lsf_queue: queue
    lsf_walltime: "72:00"
  }
}
