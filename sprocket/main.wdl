version 1.2

import "tasks/rna_core.wdl" as core

workflow ResolveOmeRnaSprocketWorkflow {
  input {
    Array[String] sample_ids
    Array[File] r1_fastqs
    Array[File] r2_fastqs

    Directory star_index
    Directory salmon_index
    File ref_gtf
    File ref_flat
    File genebody_bed

    Int qc_target_individual_reads = 2000000
    Int full_target_individual_reads = 9223372036854775807
    Int sampling_seed = 12345
    Int preprocess_threads = 2
    Int alignment_threads = 8
    Int salmon_threads = 8
    String container = "docker://ghcr.io/replace-with-your-org/resolveome-rna:1.0.0"
  }

  call core.ValidateInputs as input_audit {
    input:
      sample_ids = sample_ids,
      r1_fastqs = r1_fastqs,
      r2_fastqs = r2_fastqs,
      image = container
  }

  scatter (i in range(length(sample_ids))) {
    call core.Preprocess as qc_preprocess {
      input:
        sample_id = sample_ids[i],
        r1 = r1_fastqs[i],
        r2 = r2_fastqs[i],
        target_individual_reads = qc_target_individual_reads,
        seed = sampling_seed,
        threads = preprocess_threads,
        image = container,
        upstream_audit = input_audit.audit
    }

    call core.Align as qc_align {
      input:
        sample_id = sample_ids[i],
        r1 = qc_preprocess.trimmed_r1,
        r2 = qc_preprocess.trimmed_r2,
        star_index = star_index,
        threads = alignment_threads,
        image = container
    }

    call core.CountAndQc as qc_count_qc {
      input:
        sample_id = sample_ids[i],
        bam = qc_align.bam,
        bai = qc_align.bai,
        gtf = ref_gtf,
        ref_flat = ref_flat,
        genebody_bed = genebody_bed,
        threads = preprocess_threads,
        image = container
    }
  }

  call core.BuildMatrix as qc_matrix {
    input:
      sample_ids = sample_ids,
      htseq_counts = qc_count_qc.htseq_counts,
      image = container
  }

  call core.MultiqcReport as qc_multiqc {
    input:
      fastp_jsons = qc_preprocess.fastp_json,
      fastp_htmls = qc_preprocess.fastp_html,
      fastqc_html = flatten([flatten(qc_preprocess.raw_fastqc_html), flatten(qc_preprocess.trimmed_fastqc_html)]),
      fastqc_zip = flatten([flatten(qc_preprocess.raw_fastqc_zip), flatten(qc_preprocess.trimmed_fastqc_zip)]),
      star_logs = qc_align.star_log_final,
      qualimap_dirs = qc_count_qc.qualimap_output,
      image = container
  }

  call core.BranchGate as qc_gate {
    input:
      audit = qc_matrix.audit,
      multiqc_report = qc_multiqc.report,
      image = container
  }

  scatter (i in range(length(sample_ids))) {
    call core.Preprocess as full_preprocess {
      input:
        sample_id = sample_ids[i],
        r1 = r1_fastqs[i],
        r2 = r2_fastqs[i],
        target_individual_reads = full_target_individual_reads,
        seed = sampling_seed,
        threads = preprocess_threads,
        image = container,
        upstream_audit = qc_gate.gate
    }

    call core.Align as full_align {
      input:
        sample_id = sample_ids[i],
        r1 = full_preprocess.trimmed_r1,
        r2 = full_preprocess.trimmed_r2,
        star_index = star_index,
        threads = alignment_threads,
        image = container
    }

    call core.CountAndQc as full_count_qc {
      input:
        sample_id = sample_ids[i],
        bam = full_align.bam,
        bai = full_align.bai,
        gtf = ref_gtf,
        ref_flat = ref_flat,
        genebody_bed = genebody_bed,
        threads = preprocess_threads,
        image = container
    }
  }

  call core.BuildMatrix as full_matrix {
    input:
      sample_ids = sample_ids,
      htseq_counts = full_count_qc.htseq_counts,
      image = container
  }

  call core.MultiqcReport as full_multiqc {
    input:
      fastp_jsons = full_preprocess.fastp_json,
      fastp_htmls = full_preprocess.fastp_html,
      fastqc_html = flatten([flatten(full_preprocess.raw_fastqc_html), flatten(full_preprocess.trimmed_fastqc_html)]),
      fastqc_zip = flatten([flatten(full_preprocess.raw_fastqc_zip), flatten(full_preprocess.trimmed_fastqc_zip)]),
      star_logs = full_align.star_log_final,
      qualimap_dirs = full_count_qc.qualimap_output,
      image = container
  }

  call core.BranchGate as full_gate {
    input:
      audit = full_matrix.audit,
      multiqc_report = full_multiqc.report,
      image = container
  }

  scatter (i in range(length(sample_ids))) {
    call core.SalmonQuant as salmon_quant {
      input:
        sample_id = sample_ids[i],
        r1 = full_preprocess.trimmed_r1[i],
        r2 = full_preprocess.trimmed_r2[i],
        salmon_index = salmon_index,
        threads = salmon_threads,
        image = container,
        upstream_gate = full_gate.gate
    }
  }

  output {
    File input_validation_report = input_audit.audit

    File qc_gene_count_matrix = qc_matrix.matrix
    File qc_gene_count_parquet = qc_matrix.matrix_parquet
    File qc_sample_expression_metrics = qc_matrix.sample_metrics
    File qc_branch_audit = qc_gate.gate
    File qc_multiqc_report = qc_multiqc.report
    Array[File] qc_bams = qc_align.bam
    Array[File] qc_htseq_counts = qc_count_qc.htseq_counts

    File full_gene_count_matrix = full_matrix.matrix
    File full_gene_count_parquet = full_matrix.matrix_parquet
    File full_sample_expression_metrics = full_matrix.sample_metrics
    File full_branch_audit = full_gate.gate
    File full_multiqc_report = full_multiqc.report
    Array[File] full_bams = full_align.bam
    Array[File] full_htseq_counts = full_count_qc.htseq_counts
    Array[File] salmon_quant_files = salmon_quant.quant
  }
}
