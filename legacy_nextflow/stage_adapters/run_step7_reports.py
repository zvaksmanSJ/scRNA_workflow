#!/usr/bin/env python3
"""Generate source PCA, consensus RNA QC, Parquet and MultiQC from audited products."""
import json
import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
from common import EXECUTOR, MEMORY_MB, QUEUE, THREADS, WORKFLOW
from run_all_resolveome_rna import PREFIX, NEXTFLOW, execute
from step1_prepare_samples import read_samples, write_source_csv
from step3_prepare_basejumper import SOURCE, prepare, write_config

STAGED_WORKFLOW = """
workflow {
    merged = channel.value([
        file("${params.matrix_root}/merged/df_gene_counts_starhtseq.tsv", checkIfExists: true),
        file("${params.matrix_root}/merged/df_gene_types_detected_summary_starhtseq.tsv", checkIfExists: true),
        file("${params.matrix_root}/merged/df_mt_gene_counts_starhtseq.tsv", checkIfExists: true)])
    PLOTTER_PCAHEATMAP(merged)
    CREATE_QC_REPORT(channel.fromPath("${params.qc_root}/tool_outputs/qualimap/qualimap_outdir_*",
                                   checkIfExists: true, type: 'dir').collect())
    logs = channel.fromPath("${params.align_root}/star/star_outdir_*/*Log.final.out",
                            checkIfExists: true).collect()
    counts = channel.fromPath("${params.qc_root}/tool_outputs/htseq/*.htseq_counts.tsv",
                              checkIfExists: true).collect()
    fastp = channel.fromPath("${params.preprocess_root}/trimmed/*_fastp.json",
                             checkIfExists: true).collect()
    reads = channel.fromPath("${params.preprocess_root}/sampled/*_read_counts.txt",
                             checkIfExists: true).collect()
    RNA_QC_PLOTS(logs, counts, fastp, reads, merged,
        file("${params.matrix_root}/matrix/matrix_gene_counts_starhtseq.txt", checkIfExists: true),
        CREATE_QC_REPORT.out.stats, file(params.input_csv, checkIfExists: true),
        params.dataset_id, params.workspace, params.workflow_id, params.staged_branch, params.pipeline_user)
    mqc = RNA_QC_PLOTS.out.mqc_metrics.flatten()
        .mix(fastp.flatten())
        .mix(PLOTTER_PCAHEATMAP.out.plots.flatten())
        .mix(RNA_QC_PLOTS.out.composition_jpg)
        .mix(channel.fromPath("${params.matrix_root}/matrix/HouseKeepingGenes_*_mqc.tsv",
                              checkIfExists: true))
        .mix(channel.fromPath("${params.matrix_root}/matrix/HKGenes_Expression__mqc.png",
                              checkIfExists: true))
        .mix(RNA_QC_PLOTS.out.summary_verdict).collect()
    MULTIQC(mqc, params.dataset_id, params.workspace, params.workflow_id, params.staged_branch,
            file(params.logo, checkIfExists: true))
}
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", choices=("qc", "full"), default="qc")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--output-name", choices=("stage7_qc_reports", "stage7_qc_reports_v2",
                                                "stage7_full_reports"))
    args = parser.parse_args()
    name = "step6_output_audit.json" if args.depth == "qc" else "step6_full_output_audit.json"
    audit = json.loads((HERE / name).read_text())
    rows = read_samples(HERE / "rna_cells.tsv")
    if audit["status"] != "PASS" or audit["samples"] != sorted(r["sample_id"] for r in rows):
        raise ValueError("Matrix audit must pass for the complete cohort")
    output_name = args.output_name or f"stage7_{args.depth}_reports"
    if not output_name.startswith(f"stage7_{args.depth}_"):
        parser.error("--output-name must match --depth")
    root = HERE / output_name
    root.mkdir(exist_ok=True)
    if (root / "submission_started").exists() and not args.resume:
        raise RuntimeError("Reports already started; inspect tasks before retrying")
    runtime = prepare(SOURCE, root / "runtime")
    main_file = root / "reports.nf"
    main_file.write_text((runtime / "main.nf").read_text().split("workflow {\n", 1)[0]
                         + STAGED_WORKFLOW)
    write_source_csv(rows, root / "rna_input.csv")
    config = write_config(runtime, root / "native.config", PREFIX, EXECUTOR, QUEUE,
                          THREADS, MEMORY_MB)
    with config.open("a") as stream:
        for name, directory in (("PLOTTER_PCAHEATMAP", "pca"),
                                ("CREATE_QC_REPORT", "qualimap_summary"),
                                ("RNA_QC_PLOTS", "consensus"), ("MULTIQC", "multiqc")):
            stream.write(f"""
process {{
    withName: '{name}' {{
        publishDir = [path: {{ "${{params.outputDir}}/{directory}" }}, mode: 'copy']
    }}
}}
""")
    (root / "submission_started").write_text("Inspect jobs and outputs before retrying.\n")
    execute([NEXTFLOW, "-log", root / "nextflow.log", "-c", config, "run", main_file,
             "--outputDir", root, "--input_csv", root / "rna_input.csv",
             "--matrix_root", HERE / f"stage6_{args.depth}_matrix",
             "--qc_root", HERE / f"stage5_{args.depth}_counts",
             "--align_root", HERE / f"stage4_{args.depth}_alignment",
             "--preprocess_root", HERE / f"stage3_{args.depth}_preprocessing",
             "--staged_branch", f"staged-{args.depth}",
             "--logo", runtime / "assets/bioskryb_logo-tagline.png",
             "-work-dir", root / "work", "-with-trace", root / "trace.tsv",
             "-with-report", root / "execution.html", *(["-resume"] if args.resume else [])])
    print("STEP 7 completed. Audit cohort metrics, QC labels, Parquet and reports "
          "before claiming this branch complete.", flush=True)


if __name__ == "__main__":
    main()
