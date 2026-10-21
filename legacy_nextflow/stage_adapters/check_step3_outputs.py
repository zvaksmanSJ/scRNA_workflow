#!/usr/bin/env python3
"""Audit retained preprocessing reads, metrics and both FastQC reports per cell."""
import csv
import argparse
import gzip
import json
from pathlib import Path
import re
import sys
import zipfile

HERE = Path(__file__).resolve().parent
from common import WORKFLOW
from step1_prepare_samples import fastq_record, read_samples, validate_pairs


def spot_check_pairs(r1, r2, limit=1000):
    for path in (r1, r2):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Missing or empty FASTQ: {path}")
    with gzip.open(r1, "rt") as left, gzip.open(r2, "rt") as right:
        for index in range(limit):
            a, b = fastq_record(left), fastq_record(right)
            if a is None or b is None or a != b:
                raise ValueError(f"Invalid spot-check pair {index + 1}: {r1}")
    return limit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", choices=("qc", "full"), default="qc")
    parser.add_argument("--quick", action="store_true",
                        help="Spot-check all cells and reuse completed deep checks for at least 3 cells")
    parser.add_argument("--completed-deep-log", type=Path,
                        help="Output log from this depth's already-running comprehensive audit")
    args = parser.parse_args()
    if args.quick != bool(args.completed_deep_log):
        parser.error("--quick and --completed-deep-log must be supplied together")
    completed = {}
    if args.quick:
        for line in args.completed_deep_log.read_text().splitlines():
            match = re.fullmatch(
                r"(\S+): sampled=([\d,]+) trimmed=([\d,]+) retained=[\d.]+% Q30=[\d.]+%", line)
            if match:
                completed[match[1]] = tuple(int(value.replace(",", "")) for value in match.groups()[1:])
    root = HERE / ("stage3_qc_preprocessing" if args.depth == "qc" else "stage3_full_preprocessing")
    samples = read_samples(HERE / "rna_cells.tsv")
    trace = list(csv.DictReader((root / "trace.tsv").open(), delimiter="\t"))
    for process in ("RAW_FASTQC", "SEQKIT_SAMPLE", "FASTP_TRIM", "TRIMMED_FASTQC"):
        tasks = [row for row in trace if row["name"].startswith(process + " (")]
        names = [row["name"].split(" (", 1)[1].removesuffix(")") for row in tasks]
        if sorted(names) != sorted(row["sample_id"] for row in samples):
            raise ValueError(f"Task cohort mismatch: {process}")
        if any(row["status"] not in ("COMPLETED", "CACHED") or row["exit"] != "0" for row in tasks):
            raise ValueError(f"Failed preprocessing task: {process}")
    results = []
    for row in samples:
        sid = row["sample_id"]
        provenance = HERE / "fastqs" / f"{sid}.merge.json"
        counts = list(map(int, (root / "sampled" / f"{sid}_read_counts.txt").read_text().split()))
        if len(counts) != 2 or counts[0] % 2:
            raise ValueError(f"Invalid source/sample read-count record: {sid}")
        merged_pairs = counts[0] // 2
        if provenance.is_file():
            declared_pairs = int(json.loads(provenance.read_text())["merged_pairs"])
            if declared_pairs != merged_pairs:
                raise ValueError(f"Merge provenance pair count mismatch: {sid}")
        metrics = json.loads((root / "trimmed" / f"{sid}_fastp.json").read_text())
        before, after = metrics["summary"]["before_filtering"], metrics["summary"]["after_filtering"]
        if args.quick:
            for directory, suffix in (("sampled", "subsampled"), ("trimmed", "trim")):
                paths = [
                    root / directory / (f"{sid}_subsampled_{mate}.fastq.gz" if suffix == "subsampled"
                                        else f"{sid}_{mate}_trim.fastq.gz")
                    for mate in ("R1", "R2")]
                spot_check_pairs(*paths)
            if before["total_reads"] % 2 or after["total_reads"] % 2:
                raise ValueError(f"Odd paired-read totals: {sid}")
            sampled, trimmed = before["total_reads"] // 2, after["total_reads"] // 2
            if sid in completed and completed[sid] != (2 * sampled, 2 * trimmed):
                raise ValueError(f"Completed deep-check counts disagree: {sid}")
            if args.depth == "full":
                for mate in ("r1", "r2"):
                    copied = root / "sampled" / f"{sid}_subsampled_{mate.upper()}.fastq.gz"
                    if copied.stat().st_size != Path(row[mate]).stat().st_size:
                        raise ValueError(f"Full-depth copy size mismatch: {sid}/{mate}")
        else:
            sampled = validate_pairs(root / "sampled" / f"{sid}_subsampled_R1.fastq.gz",
                                     root / "sampled" / f"{sid}_subsampled_R2.fastq.gz")
            trimmed = validate_pairs(root / "trimmed" / f"{sid}_R1_trim.fastq.gz",
                                     root / "trimmed" / f"{sid}_R2_trim.fastq.gz")
        if counts != [2 * merged_pairs, 2 * sampled]:
            raise ValueError(f"Source read-count mismatch: {sid}")
        if args.depth == "full" and sampled != merged_pairs:
            raise ValueError(f"Full-depth branch unexpectedly subsampled reads: {sid}")
        if before["total_reads"] != 2 * sampled or after["total_reads"] != 2 * trimmed:
            raise ValueError(f"fastp read-count mismatch: {sid}")
        if not 0 < trimmed <= sampled:
            raise ValueError(f"Invalid retained counts: {sid}")
        if sum(metrics["filtering_result"].values()) != before["total_reads"]:
            raise ValueError(f"fastp filtering categories do not account for every read: {sid}")
        qc = {}
        for stage in ("raw_fastqc", "trimmed_fastqc"):
            directory = root / "extras" / stage / sid
            archives = list(directory.glob("*.zip"))
            if len(archives) != 2 or len(list(directory.glob("*.html"))) != 2:
                raise ValueError(f"Expected two FastQC HTML/ZIP reports: {stage}/{sid}")
            qc[stage] = {}
            for archive in archives:
                with zipfile.ZipFile(archive) as stream:
                    if stream.testzip() is not None:
                        raise ValueError(f"Corrupt FastQC archive: {archive}")
                    summary = next(name for name in stream.namelist() if name.endswith("/summary.txt"))
                    qc[stage][archive.name] = stream.read(summary).decode().splitlines()
        result = {"sample": sid, "sampled_reads": 2 * sampled, "trimmed_reads": 2 * trimmed,
                  "retained_fraction": trimmed / sampled, "Q30_before": before["q30_rate"],
                  "Q30_after": after["q30_rate"], "FastQC": qc}
        if args.quick:
            result["fastq_verification"] = (
                "completed_deep_check_reused" if sid in completed else "first_1000_pairs_only")
        results.append(result)
        print(f'{sid}: sampled={2 * sampled:,} trimmed={2 * trimmed:,} '
              f'retained={trimmed / sampled:.1%} Q30={after["q30_rate"]:.1%}', flush=True)
    report = {"status": "PASS", "samples": results}
    name = "step3_output_audit" if args.depth == "qc" else "step3_full_output_audit"
    if args.quick:
        deep_cells = [r["sample"] for r in results if r["sample"] in completed]
        if len(deep_cells) < 3:
            raise ValueError("Quick audit requires at least 3 completed deep checks")
        report.update(audit_mode="quick", depth=args.depth, deep_checked_cells=deep_cells,
                      deep_check_log=str(args.completed_deep_log.resolve()),
                      limitation="Remaining cells have only 1000-pair spot checks; full gzip integrity "
                                 "and independently counted totals are not verified for those cells.")
        name += "_quick"
    (HERE / f"{name}.json").write_text(json.dumps(report, indent=2) + "\n")
    if args.quick:
        print(f"STEP 3 QUICK PASS: {len(results)} cells checked; {len(deep_cells)} completed "
              "deep checks reused. Remaining FASTQs are spot-checked, not fully verified.")
        return
    print("STEP 3 PASS: all 64 tasks, 64 retained FASTQs, 16 count files, "
          "16 fastp JSONs and 64 FastQC HTML/ZIP pairs accounted for.")


if __name__ == "__main__":
    main()
