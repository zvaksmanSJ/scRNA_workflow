## Top-level WDL DAG for the staged ResolveOME RNA workflow.
##
## Cromwell owns the stage order, audit gates, and Salmon scatter. The tested
## BaseJumper scientific tasks remain in a separately authorized Nextflow source and are
## invoked by the existing per-stage adapters. This preserves those methods
## while making WDL the user-facing workflow engine.
##
## Execution map:
## 1. RnaPreflight validates the run directory, sample set, tools, and references.
## 2. The QC-depth branch runs preprocessing, alignment, counts, matrix, and QC
##    reports; each stage audit gates the next stage.
## 3. The full-depth branch starts only after the QC report audit passes and
##    repeats preprocessing/alignment/counting/reporting on all reads.
## 4. SalmonQuant is scattered per sample; its manifest and output audit gate
##    final cohort analysis.
## 5. OrganizeRnaResults publishes the organized results/readme after analysis.
## Optional research add-ons are separate from this DAG and are not launched here.
version 1.0

import "wdl/tasks/rna_stages.wdl" as stages
import "wdl/tasks/rna_salmon.wdl" as salmon

workflow ResolveOmeRnaWorkflow {
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
    String salmon_executable
    String salmon_index
    String queue = "standard"
    Int stage_memory_mb = 8000
    Int nextflow_memory_mb = 64000
    Int nextflow_threads = 8
    Int salmon_memory_mb = 64000
    Int salmon_threads = 8
    Int run_serial_hours = 168
  }

  call stages.RnaPreflight as preflight {
    input:
      workflow_rna_dir = workflow_rna_dir,
      runtime_dir = runtime_dir,
      results_dir = results_dir,
      basej_source_dir = basej_source_dir,
      sample_sheet = sample_sheet,
      sample_ids = sample_ids,
      ref_fasta = ref_fasta,
      ref_gtf = ref_gtf,
      tx2gene = tx2gene,
      genebody_bed = genebody_bed,
      star_index_dir = star_index_dir,
      ref_flat = ref_flat,
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb
  }

  # The two-million-read QC branch is completed and audited before the
  # full-depth branch is allowed to start.
  call stages.RunRnaStage as qc_preprocessing {
    input:
      stage_name = "qc_step3_preprocessing",
      results_dir = results_dir,
      runner_name = "run_step3_preprocessing.py",
      checker_name = "check_step3_outputs.py",
      depth = "qc",
      runner_flags = "",
      checker_flags = "",
      upstream_gate = preflight.gate_report,
      audit_path = results_dir + "/step3_output_audit.json",
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb,
      walltime_hours = run_serial_hours,
      nextflow_threads = nextflow_threads,
      nextflow_memory_mb = nextflow_memory_mb
  }

  call stages.RunRnaStage as qc_alignment {
    input:
      stage_name = "qc_step4_alignment",
      results_dir = results_dir,
      runner_name = "run_step4_alignment.py",
      checker_name = "check_step4_outputs.py",
      depth = "qc",
      runner_flags = "",
      checker_flags = "",
      upstream_gate = qc_preprocessing.audit,
      audit_path = results_dir + "/step4_output_audit.json",
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb,
      walltime_hours = run_serial_hours,
      nextflow_threads = nextflow_threads,
      nextflow_memory_mb = nextflow_memory_mb
  }

  call stages.RunRnaStage as qc_counts {
    input:
      stage_name = "qc_step5_counts",
      results_dir = results_dir,
      runner_name = "run_step5_counts_qc.py",
      checker_name = "check_step5_outputs.py",
      depth = "qc",
      runner_flags = "",
      checker_flags = "",
      upstream_gate = qc_alignment.audit,
      audit_path = results_dir + "/step5_output_audit.json",
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb,
      walltime_hours = run_serial_hours,
      nextflow_threads = nextflow_threads,
      nextflow_memory_mb = nextflow_memory_mb
  }

  call stages.RunRnaStage as qc_matrix {
    input:
      stage_name = "qc_step6_matrix",
      results_dir = results_dir,
      runner_name = "run_step6_matrix.py",
      checker_name = "check_step6_outputs.py",
      depth = "qc",
      runner_flags = "",
      checker_flags = "",
      upstream_gate = qc_counts.audit,
      audit_path = results_dir + "/step6_output_audit.json",
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb,
      walltime_hours = run_serial_hours,
      nextflow_threads = nextflow_threads,
      nextflow_memory_mb = nextflow_memory_mb
  }

  call stages.RunRnaStage as qc_reports {
    input:
      stage_name = "qc_step7_reports",
      results_dir = results_dir,
      runner_name = "run_step7_reports.py",
      checker_name = "check_step7_outputs.py",
      depth = "qc",
      runner_flags = "--output-name stage7_qc_reports_v2",
      checker_flags = "",
      upstream_gate = qc_matrix.audit,
      audit_path = results_dir + "/step7_output_audit.json",
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb,
      walltime_hours = run_serial_hours,
      nextflow_threads = nextflow_threads,
      nextflow_memory_mb = nextflow_memory_mb
  }

  call stages.RunRnaStage as full_preprocessing {
    input:
      stage_name = "full_step3_preprocessing",
      results_dir = results_dir,
      runner_name = "run_step3_preprocessing.py",
      checker_name = "check_step3_outputs.py",
      depth = "full",
      runner_flags = "",
      checker_flags = "",
      upstream_gate = qc_reports.audit,
      audit_path = results_dir + "/step3_full_output_audit.json",
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb,
      walltime_hours = run_serial_hours,
      nextflow_threads = nextflow_threads,
      nextflow_memory_mb = nextflow_memory_mb
  }

  call stages.RunRnaStage as full_alignment {
    input:
      stage_name = "full_step4_alignment",
      results_dir = results_dir,
      runner_name = "run_step4_alignment.py",
      checker_name = "check_step4_outputs.py",
      depth = "full",
      runner_flags = "",
      checker_flags = "--quick",
      upstream_gate = full_preprocessing.audit,
      audit_path = results_dir + "/step4_full_output_audit_quick.json",
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb,
      walltime_hours = run_serial_hours,
      nextflow_threads = nextflow_threads,
      nextflow_memory_mb = nextflow_memory_mb
  }

  call stages.RunRnaStage as full_counts {
    input:
      stage_name = "full_step5_counts",
      results_dir = results_dir,
      runner_name = "run_step5_counts_qc.py",
      checker_name = "check_step5_outputs.py",
      depth = "full",
      runner_flags = "--accept-quick-audit",
      checker_flags = "",
      upstream_gate = full_alignment.audit,
      audit_path = results_dir + "/step5_full_output_audit.json",
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb,
      walltime_hours = run_serial_hours,
      nextflow_threads = nextflow_threads,
      nextflow_memory_mb = nextflow_memory_mb
  }

  call stages.RunRnaStage as full_matrix {
    input:
      stage_name = "full_step6_matrix",
      results_dir = results_dir,
      runner_name = "run_step6_matrix.py",
      checker_name = "check_step6_outputs.py",
      depth = "full",
      runner_flags = "",
      checker_flags = "",
      upstream_gate = full_counts.audit,
      audit_path = results_dir + "/step6_full_output_audit.json",
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb,
      walltime_hours = run_serial_hours,
      nextflow_threads = nextflow_threads,
      nextflow_memory_mb = nextflow_memory_mb
  }

  call stages.RunRnaStage as full_reports {
    input:
      stage_name = "full_step7_reports",
      results_dir = results_dir,
      runner_name = "run_step7_reports.py",
      checker_name = "check_step7_outputs.py",
      depth = "full",
      runner_flags = "",
      checker_flags = "",
      upstream_gate = full_matrix.audit,
      audit_path = results_dir + "/step7_full_output_audit.json",
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb,
      walltime_hours = run_serial_hours,
      nextflow_threads = nextflow_threads,
      nextflow_memory_mb = nextflow_memory_mb
  }

  scatter (i in range(length(sample_ids))) {
    call salmon.SalmonQuant as salmon_quant {
      input:
        sample_id = sample_ids[i],
        results_dir = results_dir,
        full_preprocessing_gate = full_preprocessing.audit,
        salmon_executable = salmon_executable,
        salmon_index = salmon_index,
        threads = salmon_threads,
        memory_mb = salmon_memory_mb,
        queue = queue
    }
  }

  call stages.RecordWdlSalmonManifest as salmon_manifest {
    input:
      results_dir = results_dir,
      workflow_rna_dir = workflow_rna_dir,
      full_alignment_audit = full_alignment.audit,
      quantification_files = salmon_quant.quant_sf,
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb
  }

  call stages.ValidateSalmon as validate_salmon {
    input:
      workflow_rna_dir = workflow_rna_dir,
      results_dir = results_dir,
      full_alignment_audit = full_alignment.audit,
      salmon_manifest = salmon_manifest.manifest,
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb
  }

  call stages.RunFinalAnalysis as final_analysis {
    input:
      workflow_rna_dir = workflow_rna_dir,
      results_dir = results_dir,
      full_reports_audit = full_reports.audit,
      salmon_audit = validate_salmon.audit,
      conda_base = conda_base,
      conda_prefix = conda_prefix,
      queue = queue,
      memory_mb = stage_memory_mb
  }

  call stages.OrganizeRnaResults as organize_results {
    input:
      workflow_rna_dir = workflow_rna_dir,
      results_dir = results_dir,
      final_dashboard = final_analysis.dashboard,
      queue = queue,
      memory_mb = stage_memory_mb
  }

  output {
    File qc_report_audit = qc_reports.audit
    File full_report_audit = full_reports.audit
    File salmon_audit = validate_salmon.audit
    File dashboard = final_analysis.dashboard
    File results_readme = organize_results.readme
  }
}
