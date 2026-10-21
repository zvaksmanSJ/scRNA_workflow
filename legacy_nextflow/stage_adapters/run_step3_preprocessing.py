#!/usr/bin/env python3
"""Run only source read preparation, trimming and FastQC; no alignment."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
from common import EXECUTOR, MEMORY_MB, QUEUE, THREADS, WORKFLOW
from run_all_resolveome_rna import PREFIX, NEXTFLOW, execute
from step1_prepare_samples import read_samples, write_source_csv
from step3_prepare_basejumper import SOURCE, prepare, write_config

STAGED_WORKFLOW = """
workflow {
    reads = channel.fromPath(params.input_csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> [row.biosampleName,
                      [file(row.read1, checkIfExists: true), file(row.read2, checkIfExists: true)]] }
    RAW_FASTQC(reads)
    SEQKIT_SAMPLE(reads.map { sid, mates -> [sid, mates, params.n_reads] },
                  params.seqkit_sample_seed)
    FASTP_TRIM(SEQKIT_SAMPLE.out.reads)
    TRIMMED_FASTQC(FASTP_TRIM.out.reads)
}
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", choices=("qc", "full"), default="qc")
    args = parser.parse_args()
    if args.depth == "full":
        report = json.loads((HERE / "step7_output_audit.json").read_text())
        if report["status"] != "PASS":
            raise ValueError("QC-depth final audit must pass before full-depth preparation")
    root = HERE / ("stage3_qc_preprocessing" if args.depth == "qc" else "stage3_full_preprocessing")
    root.mkdir(exist_ok=True)
    if (root / "submission_started").exists():
        raise RuntimeError("Preprocessing was already started; inspect its jobs before retrying")
    rows = read_samples(HERE / "rna_cells.tsv")
    runtime = prepare(SOURCE, root / "runtime")
    main_file = root / "preprocessing.nf"
    main_file.write_text((runtime / "main.nf").read_text().split("workflow {\n", 1)[0]
                         + STAGED_WORKFLOW)
    write_source_csv(rows, root / "rna_input.csv")
    config = write_config(runtime, root / "native.config", PREFIX, EXECUTOR, QUEUE,
                          THREADS, MEMORY_MB)
    with config.open("a") as stream:
        stream.write("""
process {
    withName: 'SEQKIT_SAMPLE' {
        publishDir = [path: { "${params.outputDir}/sampled" }, mode: 'copy']
    }
    withName: 'FASTP_TRIM' {
        publishDir = [path: { "${params.outputDir}/trimmed" }, mode: 'copy']
    }
}
""")
    (root / "submission_started").write_text("Inspect logs and trace before another submission.\n")
    execute([NEXTFLOW, "-log", root / "nextflow.log", "-c", config, "run", main_file,
             "--input_csv", root / "rna_input.csv", "--outputDir", root,
             "--n_reads", 2000000 if args.depth == "qc" else 9223372036854775807,
             "--seqkit_sample_seed", 12345,
             "--trim_poly_g", "false", "-work-dir", root / "work",
             "-with-trace", root / "trace.tsv",
             "-with-report", root / "execution.html"])
    print("STEP 3 processes completed. Inspect paired counts, fastp metrics and FastQC "
          "before any alignment.", flush=True)


if __name__ == "__main__":
    main()
