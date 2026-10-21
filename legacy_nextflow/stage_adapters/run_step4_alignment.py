#!/usr/bin/env python3
"""Align inspected reads using source STAR and samtools; no counting/report tasks."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
from common import EXECUTOR, MEMORY_MB, QUEUE, THREADS, WORKFLOW
from run_all_resolveome_rna import PREFIX, NEXTFLOW, execute
from step1_prepare_samples import read_samples, write_source_csv
from step2_prepare_references import DEFAULTS, check_references
from step3_prepare_basejumper import SOURCE, prepare, write_config

STAGED_WORKFLOW = """
workflow {
    reads = channel.fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> [row.biosampleName,
                      [file(row.read1, checkIfExists: true), file(row.read2, checkIfExists: true)]] }
    STAR_ALIGN(reads, file(params.star_index, checkIfExists: true))
    SAMTOOLS_INDEX_FILTER(
        STAR_ALIGN.out.bam.map { sid, bam, junction -> [sid, bam, junction] },
        file(params.ref_fasta, checkIfExists: true))
}
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", choices=("qc", "full"), default="qc")
    parser.add_argument("--accept-quick-audit", action="store_true")
    args = parser.parse_args()
    audit_name = "step3_output_audit" if args.depth == "qc" else "step3_full_output_audit"
    if args.accept_quick_audit:
        audit_name += "_quick"
    audit = json.loads((HERE / f"{audit_name}.json").read_text())
    if args.accept_quick_audit and (
            audit.get("audit_mode") != "quick" or audit.get("depth") != args.depth
            or len(audit.get("deep_checked_cells", [])) < 3):
        raise ValueError("Quick audit lacks the required depth and deep-cell evidence")
    rows = read_samples(HERE / "rna_cells.tsv")
    if audit["status"] != "PASS" or sorted(r["sample"] for r in audit["samples"]) != sorted(
            r["sample_id"] for r in rows):
        raise ValueError("Preprocessing audit has not passed for the complete cohort")
    check_references(DEFAULTS)
    root = HERE / f"stage4_{args.depth}_alignment"
    root.mkdir(exist_ok=True)
    if (root / "submission_started").exists():
        raise RuntimeError("Alignment was already started; inspect existing jobs before retrying")
    for row in rows:
        sid = row["sample_id"]
        row["r1"] = str(HERE / f"stage3_{args.depth}_preprocessing/trimmed" / f"{sid}_R1_trim.fastq.gz")
        row["r2"] = str(HERE / f"stage3_{args.depth}_preprocessing/trimmed" / f"{sid}_R2_trim.fastq.gz")
    (root / "preprocessing_audit_used.json").write_text(json.dumps({
        "path": str(HERE / f"{audit_name}.json"), "depth": args.depth,
        "audit_mode": audit.get("audit_mode", "comprehensive")}, indent=2) + "\n")
    runtime = prepare(SOURCE, root / "runtime")
    main_file = root / "alignment.nf"
    main_file.write_text((runtime / "main.nf").read_text().split("workflow {\n", 1)[0]
                         + STAGED_WORKFLOW)
    write_source_csv(rows, root / "rna_input.csv")
    config = write_config(runtime, root / "native.config", PREFIX, EXECUTOR, QUEUE,
                          THREADS, MEMORY_MB)
    with config.open("a") as stream:
        stream.write("""
process {
    withName: 'STAR_ALIGN' {
        publishDir = [path: { "${params.outputDir}/star" }, mode: 'copy']
    }
    withName: 'SAMTOOLS_INDEX_FILTER' {
        publishDir = [path: { "${params.outputDir}/bam" }, mode: 'copy',
                      pattern: '*.bam{,.bai}']
    }
}
""")
    (root / "submission_started").write_text("Inspect logs and trace before another submission.\n")
    execute([NEXTFLOW, "-log", root / "nextflow.log", "-c", config, "run", main_file,
             "--input_csv", root / "rna_input.csv", "--outputDir", root,
             "--star_index", DEFAULTS["star_index"], "--ref_fasta", DEFAULTS["fasta"],
             "-work-dir", root / "work", "-with-trace", root / "trace.tsv",
             "-with-report", root / "execution.html"])
    print("STEP 4 processes completed. Inspect STAR mapping and published BAM/index integrity "
          "before counting.", flush=True)


if __name__ == "__main__":
    main()
