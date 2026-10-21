#!/usr/bin/env python3
"""Run matched-version Salmon quantification on the audited full-depth libraries."""
import argparse
import csv
import json
import math
import os
import re
import shlex
import subprocess
from pathlib import Path

RNAW = Path(__file__).resolve().parents[1]
RUNTIME = Path(os.environ.get("RNA_RUNTIME_DIR", RNAW / "wdl_runtime"))
RNAR = Path(os.environ.get("RESOLVEOME_RNA_RESULTS", RNAW / "wdl_runtime/results")).resolve()
PREFIX = Path(os.environ.get("RNA_SALMON_PREFIX", RUNTIME / "env/resolveome-rna-salmon"))
INDEX = Path(os.environ.get("RNA_SALMON_INDEX", RUNTIME / "references/Ensembl104/salmon_index"))
TX2GENE = Path(os.environ.get("RNA_TX2GENE", RUNTIME / "references/Ensembl104/Ensembl104.tx2gene.tsv"))
SAMPLES = Path(os.environ.get("RESOLVEOME_RNA_CELLS", RNAW / "data/rna_cells.local.tsv"))
OUT = RNAR / "final_results"
QUANTS = OUT / "02_transcript/quant_by_sample"
JOBS = OUT / "07_qc_provenance/salmon_jobs.tsv"


def read_samples(path=None):
    path = Path(path) if path is not None else SAMPLES
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or len({row["sample_id"] for row in rows}) != len(rows):
        raise ValueError("Sample manifest must contain unique sample IDs")
    return rows


def preflight():
    audit = json.loads(
        (RNAR / "step4_full_output_audit_quick.json").read_text())
    if audit.get("status") != "PASS" or audit.get("depth") != "full":
        raise RuntimeError("Full-depth alignment audit must pass before Salmon quantification")
    version_info = json.loads((INDEX / "versionInfo.json").read_text())
    if version_info.get("salmonVersion") != "1.6.0":
        raise RuntimeError(
            f"Shared Salmon index was built by {version_info.get('salmonVersion')}, "
            "not the pinned Salmon 1.6.0")
    if not TX2GENE.is_file():
        raise FileNotFoundError(TX2GENE)
    if not (PREFIX / "bin/salmon").is_file():
        raise FileNotFoundError(f"Salmon executable missing: {PREFIX / 'bin/salmon'}")
    rows = read_samples()
    inputs = []
    for row in rows:
        sid = row["sample_id"]
        trimmed = RNAR / "stage3_full_preprocessing/trimmed"
        r1 = trimmed / f"{sid}_R1_trim.fastq.gz"
        r2 = trimmed / f"{sid}_R2_trim.fastq.gz"
        if not r1.is_file() or not r2.is_file():
            raise FileNotFoundError(f"Missing full-depth trimmed mates for {sid}")
        inputs.append((sid, r1, r2))
    return inputs


def write_job_manifest(rows):
    temporary = JOBS.with_suffix(".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample_id", "job_id", "status", "output",
                                "stdout_log_pattern"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(JOBS)


def record_wdl_completions():
    inputs = preflight()
    if JOBS.exists():
        raise FileExistsError(f"Salmon submission manifest already exists: {JOBS}")
    rows = []
    for sid, _, _ in inputs:
        quant = QUANTS / sid / "quant.sf"
        if not quant.is_file() or not quant.stat().st_size:
            raise FileNotFoundError(f"WDL Salmon output is missing for {sid}: {quant}")
        rows.append({
            "sample_id": sid,
            "job_id": "Cromwell-WDL",
            "status": "WDL_COMPLETED",
            "output": str(QUANTS / sid),
            "stdout_log_pattern": "managed-by-cromwell",
        })
    write_job_manifest(rows)


def submit():
    inputs = preflight()
    if JOBS.exists():
        raise FileExistsError(f"Salmon submission manifest already exists: {JOBS}")
    QUANTS.mkdir(parents=True, exist_ok=True)
    JOBS.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for sid, r1, r2 in inputs:
        outdir = QUANTS / sid
        if outdir.exists():
            raise FileExistsError(f"Refusing to overwrite existing Salmon output: {outdir}")
        log = JOBS.parent / "salmon_logs" / f"{sid}.%J.out"
        log.parent.mkdir(parents=True, exist_ok=True)
        rows.append({"sample_id": sid, "job_id": "", "status": "NOT_SUBMITTED",
                     "output": str(outdir), "stdout_log_pattern": str(log)})
    write_job_manifest(rows)
    for index, (sid, r1, r2) in enumerate(inputs):
        outdir = QUANTS / sid
        log = JOBS.parent / "salmon_logs" / f"{sid}.%J.out"
        salmon = [
            str(PREFIX / "bin/salmon"), "quant", "-i", str(INDEX), "-l", "A",
            "-1", str(r1), "-2", str(r2), "-p", "8", "--validateMappings",
            "--seqBias", "--gcBias", "-o", str(outdir),
        ]
        shell_command = (
            f"export PATH={shlex.quote(str(PREFIX / 'bin'))}:\"$PATH\"; "
            f"exec {shlex.join(salmon)}"
        )
        command = [
            "bsub", "-q", "standard", "-J", f"rna_salmon_{sid}", "-n", "8",
            "-M", "64000", "-R", "rusage[mem=8000] span[hosts=1]",
            "-oo", str(log), "bash", "-lc", shell_command,
        ]
        result = subprocess.run(command, text=True, capture_output=True, check=True)
        match = re.search(r"<(\d+)>", result.stdout)
        if not match:
            raise RuntimeError(f"Could not read LSF job ID from bsub output: {result.stdout}")
        rows[index]["job_id"] = match.group(1)
        rows[index]["status"] = "SUBMITTED"
        write_job_manifest(rows)
        print(f"Submitted Salmon for {sid}: job {match.group(1)}", flush=True)
    print(f"Submission manifest: {JOBS}", flush=True)


def validate_outputs():
    inputs = preflight()
    if not JOBS.is_file():
        raise FileNotFoundError(f"Salmon submission manifest is missing: {JOBS}")
    with JOBS.open(newline="") as handle:
        jobs = list(csv.DictReader(handle, delimiter="\t"))
    valid_job = lambda row: (
        (row.get("status") == "SUBMITTED" and row.get("job_id", "").isdigit())
        or (row.get("status") == "WDL_COMPLETED"
            and row.get("job_id") == "Cromwell-WDL")
    )
    if ({row.get("sample_id") for row in jobs} != {sid for sid, _, _ in inputs}
            or any(not valid_job(row) for row in jobs)):
        raise ValueError("Salmon submission manifest is incomplete or inconsistent")
    summary = []
    expected_transcripts = None
    for sid, _, _ in inputs:
        quant = QUANTS / sid / "quant.sf"
        if not quant.is_file() or quant.stat().st_size == 0:
            raise FileNotFoundError(f"Salmon quant.sf not complete for {sid}: {quant}")
        with quant.open() as handle:
            header = handle.readline().rstrip("\n").split("\t")
            if header != ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]:
                raise ValueError(f"Unexpected Salmon quant.sf header for {sid}: {header}")
            total_rows = 0
            total_reads = 0.0
            total_tpm = 0.0
            transcript_ids = set()
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 5:
                    raise ValueError(f"Malformed quant.sf row for {sid}: {line[:200]}")
                name = fields[0]
                length, effective_length, tpm, num_reads = map(float, fields[1:])
                values = (length, effective_length, tpm, num_reads)
                if not name or name in transcript_ids or not all(
                        math.isfinite(value) for value in values):
                    raise ValueError(f"Invalid Salmon transcript row for {sid}")
                if length <= 0 or effective_length < 0 or tpm < 0 or num_reads < 0:
                    raise ValueError(f"Negative/invalid Salmon estimate for {sid}")
                transcript_ids.add(name)
                total_rows += 1
                total_reads += num_reads
                total_tpm += tpm
        if expected_transcripts is None:
            expected_transcripts = transcript_ids
        elif transcript_ids != expected_transcripts:
            raise ValueError(f"Salmon transcript reference differs across samples: {sid}")
        if total_rows < 1000 or total_reads <= 0 or not (900000 <= total_tpm <= 1100000):
            raise ValueError(
                f"Salmon output sanity check failed for {sid}: rows={total_rows}, "
                f"NumReads={total_reads:.1f}, summed TPM={total_tpm:.1f}")
        summary.append({"sample_id": sid, "transcripts": total_rows,
                        "estimated_reads": round(total_reads, 2),
                        "summed_tpm": round(total_tpm, 2), "status": "PASS"})
    result_path = OUT / "07_qc_provenance/salmon_output_audit.json"
    result_path.parent.mkdir(parents=True, exist_ok=True)
    result_path.write_text(json.dumps({
        "status": "PASS", "salmon_version": "1.6.0",
        "index_version": json.loads((INDEX / "versionInfo.json").read_text()),
        "samples": summary,
    }, indent=2) + "\n")
    print(f"PASS: Salmon quant.sf validated for {len(summary)} samples: {result_path}")


def main():
    global SAMPLES, RNAR, PREFIX, INDEX, TX2GENE, OUT, QUANTS, JOBS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=Path, default=SAMPLES)
    parser.add_argument("--results-dir", type=Path, default=RNAR)
    parser.add_argument("--salmon-prefix", type=Path, default=PREFIX)
    parser.add_argument("--salmon-index", type=Path, default=INDEX)
    parser.add_argument("--tx2gene", type=Path, default=TX2GENE)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--plan", action="store_true")
    action.add_argument("--submit", action="store_true")
    action.add_argument("--record-wdl-completions", action="store_true")
    action.add_argument("--validate", action="store_true")
    args = parser.parse_args()
    SAMPLES = args.cells.resolve()
    RNAR = args.results_dir.resolve()
    PREFIX = args.salmon_prefix.resolve()
    INDEX = args.salmon_index.resolve()
    TX2GENE = args.tx2gene.resolve()
    OUT = RNAR / "final_results"
    QUANTS = OUT / "02_transcript/quant_by_sample"
    JOBS = OUT / "07_qc_provenance/salmon_jobs.tsv"
    if args.plan:
        for sid, r1, r2 in preflight():
            cmd = [PREFIX / "bin/salmon", "quant", "-i", INDEX, "-l", "A",
                   "-1", r1, "-2", r2, "-p", "8", "--validateMappings",
                   "--seqBias", "--gcBias", "-o", QUANTS / sid]
            print(shlex.join(map(str, cmd)))
    elif args.submit:
        submit()
    elif args.record_wdl_completions:
        record_wdl_completions()
    else:
        validate_outputs()


if __name__ == "__main__":
    main()
