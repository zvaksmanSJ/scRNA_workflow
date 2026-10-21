#!/usr/bin/env python3
"""Create source gene summaries and expression matrix from audited counts."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
from common import EXECUTOR, MEMORY_MB, QUEUE, THREADS, WORKFLOW
from run_all_resolveome_rna import PREFIX, NEXTFLOW, execute
from step1_prepare_samples import read_samples
from step2_prepare_references import DEFAULTS
from step3_prepare_basejumper import SOURCE, prepare, write_config

STAGED_WORKFLOW = """
workflow {
    counts = channel.fromPath(params.cells, checkIfExists: true)
        .splitCsv(header: true, sep: '\\t')
        .map { row -> [row.sample_id,
                      file("${params.count_dir}/${row.sample_id}.htseq_counts.tsv",
                           checkIfExists: true)] }
    CREATE_HTSEQ_SUMMARY(counts, file(params.tx2gene, checkIfExists: true))
    MERGE_HTSEQ_SUMMARY(CREATE_HTSEQ_SUMMARY.out.gene_counts.collect(),
                       CREATE_HTSEQ_SUMMARY.out.gene_types.collect(),
                       CREATE_HTSEQ_SUMMARY.out.mt_counts.collect())
    CREATE_HTSEQ_MATRIX(MERGE_HTSEQ_SUMMARY.out.merge_tsv,
                       params.workspace, params.workflow_id, params.dataset_id,
                       params.staged_branch, params.pipeline_user, params.genome)
}
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", choices=("qc", "full"), default="qc")
    args = parser.parse_args()
    name = "step5_output_audit.json" if args.depth == "qc" else "step5_full_output_audit.json"
    audit = json.loads((HERE / name).read_text())
    rows = read_samples(HERE / "rna_cells.tsv")
    if audit["status"] != "PASS" or sorted(r["sample"] for r in audit["samples"]) != sorted(
            r["sample_id"] for r in rows):
        raise ValueError("Complete count/QC audit must pass before matrix generation")
    root = HERE / f"stage6_{args.depth}_matrix"
    root.mkdir(exist_ok=True)
    if (root / "submission_started").exists():
        raise RuntimeError("Matrix stage already started; inspect existing jobs before retrying")
    runtime = prepare(SOURCE, root / "runtime")
    main_file = root / "matrix.nf"
    main_file.write_text((runtime / "main.nf").read_text().split("workflow {\n", 1)[0]
                         + STAGED_WORKFLOW)
    config = write_config(runtime, root / "native.config", PREFIX, EXECUTOR, QUEUE,
                          THREADS, MEMORY_MB)
    with config.open("a") as stream:
        stream.write("""
process {
    withName: 'CREATE_HTSEQ_SUMMARY' {
        publishDir = [path: { "${params.outputDir}/summaries" }, mode: 'copy']
    }
    withName: 'MERGE_HTSEQ_SUMMARY' {
        publishDir = [path: { "${params.outputDir}/merged" }, mode: 'copy']
    }
    withName: 'CREATE_HTSEQ_MATRIX' {
        publishDir = [path: { "${params.outputDir}/matrix" }, mode: 'copy']
    }
}
""")
    (root / "submission_started").write_text("Inspect tasks before another submission.\n")
    execute([NEXTFLOW, "-log", root / "nextflow.log", "-c", config, "run", main_file,
             "--cells", HERE / "rna_cells.tsv",
             "--count_dir", HERE / f"stage5_{args.depth}_counts/tool_outputs/htseq",
             "--staged_branch", f"staged-{args.depth}",
             "--outputDir", root, "--tx2gene", DEFAULTS["tx2gene"], "--genome", "GRCh38",
             "-work-dir", root / "work", "-with-trace", root / "trace.tsv",
             "-with-report", root / "execution.html"])
    print("STEP 6 processes completed. Inspect matrix sample membership, values and gene IDs "
          "before source PCA, scoring or reports.", flush=True)


if __name__ == "__main__":
    main()
