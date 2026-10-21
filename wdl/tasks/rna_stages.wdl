version 1.0

# Validates immutable inputs and the new run tree before any processing starts.
# Its PASS gate is shared by the initial QC-depth branch.
task RnaPreflight {
  input {
    String workflow_rna_dir
    String runtime_dir
    String results_dir
    String basej_source_dir
    File sample_sheet
    Array[String] sample_ids
    File ref_fasta
    File ref_gtf
    File tx2gene
    File genebody_bed
    String star_index_dir
    File ref_flat
    String conda_base
    String conda_prefix
    String queue
    Int memory_mb
  }

  command <<<
    set -euo pipefail
    workflow_dir="~{workflow_rna_dir}"
    results="~{results_dir}"
    export RESOLVEOME_RNA_WORKFLOW="$workflow_dir"
    export RNA_RUNTIME_DIR="~{runtime_dir}"
    export CONDA_BASE="~{conda_base}"
    export BASEJ_RNA_SOURCE="~{basej_source_dir}"
    for path in \
      "$results/stage3_qc_preprocessing" "$results/stage4_qc_alignment" \
      "$results/stage5_qc_counts" "$results/stage6_qc_matrix" \
      "$results/stage7_qc_reports_v2" "$results/stage3_full_preprocessing" \
      "$results/stage4_full_alignment" "$results/stage5_full_counts" \
      "$results/stage6_full_matrix" "$results/stage7_full_reports"; do
      if [ -e "$path" ]; then
        echo "[ERROR] Refusing to reuse existing stage output: $path" >&2
        echo "Use a fresh RNA project/results tree; this workflow never overwrites runs." >&2
        exit 20
      fi
    done
    test -s "$results/rna_cells.tsv"
    test -s "$results/common.py"
    test -s "$results/run_step3_preprocessing.py"
    test -s "$results/check_step7_outputs.py"
    source "~{conda_base}/etc/profile.d/conda.sh"
    conda activate "~{conda_prefix}"
    test -s "$basej_source_dir/main.nf"
    printf '%s\n' "$workflow_dir" > "$results/workflow_rna_dir.txt"
    printf '%s\n' "$runtime_dir" > "$results/runtime_dir.txt"
    printf '%s\n' "$basej_source_dir" > "$results/basej_source_dir.txt"
    python "$workflow_dir/legacy_nextflow/run_all_resolveome_rna.py" \
      --cells "~{sample_sheet}" --outdir "$results" \
      --env-prefix "~{conda_prefix}" --source-dir "$basej_source_dir" \
      --fasta "~{ref_fasta}" \
      --gtf "~{ref_gtf}" --tx2gene "~{tx2gene}" \
      --genebody "~{genebody_bed}" --star-index "~{star_index_dir}" \
      --ref-flat "~{ref_flat}" --executor lsf --queue "~{queue}" --preflight
    python "$results/check_step2_references.py"
    cp "$results/tool_preflight.json" preflight.json
    python - "~{sample_sheet}" '~{write_json(sample_ids)}' <<'PY'
import csv, json, sys
with open(sys.argv[1], newline="") as stream:
    rows = list(csv.DictReader(stream, delimiter="\t"))
expected = [row["sample_id"] for row in rows]
provided = json.loads(sys.argv[2])
if len(expected) < 5 or len(set(expected)) != len(expected) or provided != expected:
    raise SystemExit("WDL sample_ids must match at least five unique sample rows")
with open("gate.json", "w") as out:
    json.dump({"status": "PASS", "gate": "tools_and_references"}, out)
PY
  >>>

  output {
    File preflight_report = "preflight.json"
    File gate_report = "gate.json"
  }

  runtime {
    cpu: 1
    memory: "4 GB"
    memory_mb: memory_mb
    lsf_queue: queue
    lsf_walltime: "04:00"
  }
}

# Runs one existing BaseJumper stage through its adapter, then invokes the
# matching checker and writes the stage audit consumed by the next WDL call.
task RunRnaStage {
  input {
    String stage_name
    String results_dir
    String runner_name
    String checker_name
    String depth
    String runner_flags
    String checker_flags
    File upstream_gate
    String audit_path
    String conda_base
    String conda_prefix
    String queue
    Int memory_mb
    Int walltime_hours
    Int nextflow_threads
    Int nextflow_memory_mb
  }

  command <<<
    set -euo pipefail
    results="~{results_dir}"
    export RESOLVEOME_RNA_WORKFLOW="$(cat "$results/workflow_rna_dir.txt")"
    export RNA_RUNTIME_DIR="$(cat "$results/runtime_dir.txt")"
    export CONDA_BASE="~{conda_base}"
    export RNA_ENV_PREFIX="~{conda_prefix}"
    export BASEJ_RNA_SOURCE="$(cat "$results/basej_source_dir.txt")"
    export RNA_LSF_QUEUE="~{queue}"
    export RNA_NEXTFLOW_THREADS="~{nextflow_threads}"
    export RNA_NEXTFLOW_MEMORY_MB="~{nextflow_memory_mb}"
    python - "~{upstream_gate}" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    report = json.load(stream)
if report.get("status") != "PASS":
    raise SystemExit(f"Upstream audit gate is not PASS: {report.get('status')}")
PY
    source "~{conda_base}/etc/profile.d/conda.sh"
    conda activate "~{conda_prefix}"
    export RESOLVEOME_RNA_CELLS="$results/rna_cells.tsv"
    runner_flags=(~{runner_flags})
    checker_flags=(~{checker_flags})
    printf 'WDL stage started: %s depth=%s\\n' "~{stage_name}" "~{depth}"
    python "~{results_dir}/~{runner_name}" --depth "~{depth}" "${runner_flags[@]}" \
      > stage.log 2>&1
    python "~{results_dir}/~{checker_name}" --depth "~{depth}" "${checker_flags[@]}" \
      >> stage.log 2>&1
    python - "~{audit_path}" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    report = json.load(stream)
if report.get("status") != "PASS":
    raise SystemExit(f"Stage audit did not pass: {report.get('status')}")
PY
    cp "~{audit_path}" stage_audit.json
    mkdir -p "~{results_dir}/logs_driver"
    cp stage.log "~{results_dir}/logs_driver/wdl_~{stage_name}.log"
    printf 'WDL stage passed: %s depth=%s\\n' "~{stage_name}" "~{depth}" >> stage.log
  >>>

  output {
    File audit = "stage_audit.json"
    File log = "stage.log"
  }

  runtime {
    cpu: 1
    memory: "8 GB"
    memory_mb: memory_mb
    lsf_queue: queue
    lsf_walltime: "~{walltime_hours}:00"
  }
}

# Records successful per-sample Salmon task outputs in the run's provenance
# manifest; it does not re-quantify samples or replace canonical HTSeq counts.
task RecordWdlSalmonManifest {
  input {
    String results_dir
    String workflow_rna_dir
    File full_alignment_audit
    Array[File] quantification_files
    String conda_base
    String conda_prefix
    String queue
    Int memory_mb
  }

  command <<<
    set -euo pipefail
    python - "~{full_alignment_audit}" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    report = json.load(stream)
if report.get("status") != "PASS" or report.get("depth") != "full":
    raise SystemExit("Full-depth alignment audit gate failed")
PY
    for quant in ~{sep=" " quantification_files}; do
      test -s "$quant"
    done
    source "~{conda_base}/etc/profile.d/conda.sh"
    conda activate "~{conda_prefix}"
    export RESOLVEOME_RNA_RESULTS="~{results_dir}"
    python "~{workflow_rna_dir}/extra_methods/quantify_transcripts.py" --record-wdl-completions
    manifest="~{results_dir}/final_results/07_qc_provenance/salmon_jobs.tsv"
    cp "$manifest" salmon_jobs.tsv
  >>>

  output {
    File manifest = "salmon_jobs.tsv"
  }

  runtime {
    cpu: 1
    memory: "2 GB"
    memory_mb: memory_mb
    lsf_queue: queue
    lsf_walltime: "04:00"
  }
}

# Checks cohort membership, quant.sf structure, and transcript/reference
# consistency before allowing final downstream analysis.
task ValidateSalmon {
  input {
    String workflow_rna_dir
    String results_dir
    File full_alignment_audit
    File salmon_manifest
    String conda_base
    String conda_prefix
    String queue
    Int memory_mb
  }

  command <<<
    set -euo pipefail
    test -s "~{salmon_manifest}"
    results="~{results_dir}"
    export RESOLVEOME_RNA_WORKFLOW="~{workflow_rna_dir}"
    export RNA_RUNTIME_DIR="$(cat "$results/runtime_dir.txt")"
    export BASEJ_RNA_SOURCE="$(cat "$results/basej_source_dir.txt")"
    export RESOLVEOME_RNA_CELLS="$results/rna_cells.tsv"
    source "~{conda_base}/etc/profile.d/conda.sh"
    conda activate "~{conda_prefix}"
    export RESOLVEOME_RNA_RESULTS="~{results_dir}"
    python - "~{full_alignment_audit}" <<'PY'
import json, sys
with open(sys.argv[1]) as stream:
    report = json.load(stream)
if report.get("status") != "PASS" or report.get("depth") != "full":
    raise SystemExit("Full-depth alignment audit gate failed")
PY
    python "~{workflow_rna_dir}/extra_methods/quantify_transcripts.py" --validate
    cp "~{results_dir}/final_results/07_qc_provenance/salmon_output_audit.json" salmon_audit.json
  >>>

  output {
    File audit = "salmon_audit.json"
  }

  runtime {
    cpu: 1
    memory: "8 GB"
    memory_mb: memory_mb
    lsf_queue: queue
    lsf_walltime: "08:00"
  }
}

# Builds cohort-level descriptive expression summaries and the dashboard from
# audited full-depth counts, reports, and Salmon outputs.
task RunFinalAnalysis {
  input {
    String workflow_rna_dir
    String results_dir
    File full_reports_audit
    File salmon_audit
    String conda_base
    String conda_prefix
    String queue
    Int memory_mb
  }

  command <<<
    set -euo pipefail
    results="~{results_dir}"
    export RESOLVEOME_RNA_WORKFLOW="~{workflow_rna_dir}"
    export RNA_RUNTIME_DIR="$(cat "$results/runtime_dir.txt")"
    export BASEJ_RNA_SOURCE="$(cat "$results/basej_source_dir.txt")"
    export RESOLVEOME_RNA_CELLS="$results/rna_cells.tsv"
    python - "~{salmon_audit}" "~{full_reports_audit}" <<'PY'
import json, sys
for path in sys.argv[1:]:
    with open(path) as stream:
        report = json.load(stream)
    if report.get("status") != "PASS":
        raise SystemExit(f"Analysis input audit gate failed: {path}")
PY
    source "~{conda_base}/etc/profile.d/conda.sh"
    conda activate "~{conda_prefix}"
    export RESOLVEOME_RNA_RESULTS="~{results_dir}"
    python "~{workflow_rna_dir}/extra_methods/run_rna_downstream.py" > analysis.log 2>&1
    test -s "~{results_dir}/final_results/00_dashboard/index.html"
    cp "~{results_dir}/final_results/00_dashboard/index.html" dashboard.html
    mkdir -p "~{results_dir}/logs_driver"
    cp analysis.log "~{results_dir}/logs_driver/wdl_final_analysis.log"
  >>>

  output {
    File dashboard = "dashboard.html"
    File log = "analysis.log"
  }

  runtime {
    cpu: 2
    memory: "16 GB"
    memory_mb: memory_mb
    lsf_queue: queue
    lsf_walltime: "24:00"
  }
}

# Publishes the already-generated dashboard and results guide into the
# organized result tree using checked links/moves rather than duplicate copies.
task OrganizeRnaResults {
  input {
    String workflow_rna_dir
    String results_dir
    File final_dashboard
    String queue
    Int memory_mb
  }

  command <<<
    set -euo pipefail
    test -s "~{final_dashboard}"
    python "~{workflow_rna_dir}/extra_methods/organize_rna_results.py" \
      --root "~{results_dir}" --apply > organize.log 2>&1
    test -s "~{results_dir}/README.md"
    cp "~{results_dir}/README.md" results_readme.md
  >>>

  output {
    File readme = "results_readme.md"
    File log = "organize.log"
  }

  runtime {
    cpu: 1
    memory: "4 GB"
    memory_mb: memory_mb
    lsf_queue: queue
    lsf_walltime: "04:00"
  }
}
