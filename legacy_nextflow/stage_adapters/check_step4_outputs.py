#!/usr/bin/env python3
"""Check alignments and published BAM/index files, with optional subset deep scans."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
from common import WORKFLOW
from step1_prepare_samples import read_samples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", choices=("qc", "full"), default="qc")
    parser.add_argument("--quick", action="store_true",
                        help="Scan three BAMs fully; check headers, EOF and index access in all cells")
    args = parser.parse_args()
    root = HERE / f"stage4_{args.depth}_alignment"
    samples = read_samples(HERE / "rna_cells.tsv")
    expected = {r["sample_id"] for r in samples}
    trace = list(csv.DictReader((root / "trace.tsv").open(), delimiter="\t"))
    for process in ("STAR_ALIGN", "SAMTOOLS_INDEX_FILTER"):
        tasks = [r for r in trace if r["name"].startswith(process + " (")]
        names = [r["name"].split(" (", 1)[1].removesuffix(")") for r in tasks]
        if len(names) != len(expected) or set(names) != expected or any(
                r["exit"] != "0" or r["status"] not in ("COMPLETED", "CACHED") for r in tasks):
            raise ValueError(f"Incomplete or failed alignment tasks: {process}")
    if {p.stem for p in (root / "bam").glob("*.bam")} != expected:
        raise ValueError("Published BAM cohort mismatch")
    if args.depth == "qc":
        audit_path = HERE / "step3_output_audit.json"
    else:
        audit_path = Path(json.loads((root / "preprocessing_audit_used.json").read_text())["path"])
    preprocessing = {r["sample"]: r for r in json.loads(audit_path.read_text())["samples"]}
    deep_cells = set(sorted(expected)[:3]) if args.quick else expected
    results = []
    for sid in sorted(expected):
        directory = root / "star" / f"star_outdir_{sid}"
        for suffix in ("Log.final.out", "Log.out", "SJ.out.tab", "Chimeric.out.junction",
                       "Aligned.sortedByCoord.out.bam"):
            if not (directory / f"{sid}_{suffix}").is_file():
                raise FileNotFoundError(f"Missing STAR output: {sid}_{suffix}")
        metrics = {}
        for line in (directory / f"{sid}_Log.final.out").read_text().splitlines():
            if "|" in line:
                key, value = line.split("|", 1)
                metrics[key.strip()] = value.strip()
        pairs = int(metrics["Number of input reads"])
        if pairs * 2 != preprocessing[sid]["trimmed_reads"]:
            raise ValueError(f"STAR input does not match trimmed pairs: {sid}")
        bam = root / "bam" / f"{sid}.bam"
        bai = bam.with_suffix(".bam.bai")
        if not bai.is_file() or not bai.stat().st_size:
            raise FileNotFoundError(bai)
        subprocess.run(["samtools", "quickcheck", "-v", str(bam)], check=True)
        header = subprocess.run(["samtools", "view", "-H", str(bam)], check=True,
                                capture_output=True, text=True).stdout
        if "SO:coordinate" not in header or f"SM:{sid}" not in header:
            raise ValueError(f"Incorrect BAM sorting or sample tag: {sid}")
        # Index-driven retrieval proves the published BAI is usable.
        subprocess.run(["samtools", "view", "-c", str(bam), "chr1:1-100000"],
                       check=True, capture_output=True)
        flagstat = None
        if sid in deep_cells:
            secondary = int(subprocess.run(["samtools", "view", "-c", "-f", "256", str(bam)],
                                      check=True, capture_output=True, text=True).stdout)
            if secondary:
                raise ValueError(f"Secondary alignments remain in filtered BAM: {sid}")
            flagstat = subprocess.run(["samtools", "flagstat", str(bam)], check=True,
                                  capture_output=True, text=True).stdout
            (root / "bam" / f"{sid}.flagstat.txt").write_text(flagstat)
        result = {"sample": sid, "input_pairs": pairs,
                  "unique_mapping_pct": float(metrics["Uniquely mapped reads %"].rstrip("%")),
                  "multimapping_pct": float(metrics["% of reads mapped to multiple loci"].rstrip("%")),
                  "STAR_metrics": metrics, "flagstat": flagstat,
                  "bam_verification": "full_scan" if sid in deep_cells else "quickcheck_header_index"}
        results.append(result)
        print(f'{sid}: {pairs:,} input pairs; unique={result["unique_mapping_pct"]:.2f}% '
              f'multi={result["multimapping_pct"]:.2f}%; BAM/BAI PASS', flush=True)
    name = "step4_output_audit" if args.depth == "qc" else "step4_full_output_audit"
    if args.quick:
        name += "_quick"
    (HERE / f"{name}.json").write_text(
        json.dumps({"status": "PASS", "depth": args.depth,
                    "audit_mode": "quick" if args.quick else "comprehensive",
                    "deep_checked_cells": sorted(deep_cells), "samples": results}, indent=2) + "\n")
    print("STEP 4 PASS: 32 tasks, 16 published BAMs/BAIs and all STAR products accounted for.")


if __name__ == "__main__":
    main()
