#!/usr/bin/env python3
"""Run source counts and per-cell RNA QC on verified BAMs; no cohort plotting."""
import json
import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
from common import EXECUTOR, MEMORY_MB, QUEUE, THREADS, WORKFLOW
from run_all_resolveome_rna import PREFIX, NEXTFLOW, execute
from step1_prepare_samples import read_samples
from step2_prepare_references import DEFAULTS, check_references
from step3_prepare_basejumper import SOURCE, prepare, write_config

STAGED_WORKFLOW = """
workflow {
    samples = channel.fromPath(params.cells, checkIfExists: true)
        .splitCsv(header: true, sep: '\\t')
        .map { row -> [row.sample_id,
                      file("${params.bam_dir}/${row.sample_id}.bam", checkIfExists: true),
                      file("${params.bam_dir}/${row.sample_id}.bam.bai", checkIfExists: true)] }
    HTSEQ_COUNTS(samples, file(params.gtf, checkIfExists: true))
    QUALIMAP_BAMRNA(samples, file(params.gtf, checkIfExists: true))
    PICARD_RNA(samples, file(params.ref_flat, checkIfExists: true))
    GENE_BODY_COVERAGE(samples, file(params.genebody_ref, checkIfExists: true))
}
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", choices=("qc", "full"), default="qc")
    parser.add_argument("--accept-quick-audit", action="store_true")
    args = parser.parse_args()
    name = "step4_output_audit" if args.depth == "qc" else "step4_full_output_audit"
    if args.accept_quick_audit:
        name += "_quick"
    audit = json.loads((HERE / f"{name}.json").read_text())
    if args.accept_quick_audit and (
            audit.get("audit_mode") != "quick" or audit.get("depth") != args.depth
            or len(audit.get("deep_checked_cells", [])) < 3):
        raise ValueError("Quick alignment audit lacks depth or deep-cell evidence")
    rows = read_samples(HERE / "rna_cells.tsv")
    if audit["status"] != "PASS" or sorted(r["sample"] for r in audit["samples"]) != sorted(
            r["sample_id"] for r in rows):
        raise ValueError("Complete alignment audit must pass before counting")
    check_references(DEFAULTS)
    root = HERE / f"stage5_{args.depth}_counts"
    root.mkdir(exist_ok=True)
    if (root / "submission_started").exists():
        raise RuntimeError("Counting/QC already started; inspect existing jobs before retrying")
    runtime = prepare(SOURCE, root / "runtime")
    main_file = root / "counts_qc.nf"
    main_file.write_text((runtime / "main.nf").read_text().split("workflow {\n", 1)[0]
                         + STAGED_WORKFLOW)
    config = write_config(runtime, root / "native.config", PREFIX, EXECUTOR, QUEUE,
                          THREADS, MEMORY_MB)
    (root / "submission_started").write_text("Inspect jobs and outputs before retrying.\n")
    execute([NEXTFLOW, "-log", root / "nextflow.log", "-c", config, "run", main_file,
             "--cells", HERE / "rna_cells.tsv", "--bam_dir", HERE / f"stage4_{args.depth}_alignment/bam",
             "--outputDir", root, "--gtf", DEFAULTS["gtf"],
             "--ref_flat", DEFAULTS["ref_flat"], "--genebody_ref", DEFAULTS["genebody"],
             "--strandedness", "no", "--picard_strand", "NONE",
             "-work-dir", root / "work", "-with-trace", root / "trace.tsv",
             "-with-report", root / "execution.html"])
    print("STEP 5 completed: inspect HTSeq counts, Qualimap, Picard and gene-body coverage "
          "before matrix generation and cohort plots.", flush=True)


if __name__ == "__main__":
    main()
