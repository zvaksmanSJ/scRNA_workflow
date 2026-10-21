#!/usr/bin/env python3
"""Prepare an isolated WDL results tree using symlinks to immutable RNA inputs."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import sys

RNAW = Path(__file__).resolve().parents[2]
STAGE_CODE = RNAW / "legacy_nextflow/stage_adapters"
RESULT_CODE = (
    "common.py",
    "run_step3_preprocessing.py",
    "run_step4_alignment.py",
    "run_step5_counts_qc.py",
    "run_step6_matrix.py",
    "run_step7_reports.py",
    "check_step2_references.py",
    "check_step3_outputs.py",
    "check_step4_outputs.py",
    "check_step5_outputs.py",
    "check_step6_outputs.py",
    "check_step7_outputs.py",
)


def prepare(destination: Path, sample_sheet: Path) -> None:
    destination = destination.resolve()
    if destination == RNAW or RNAW in destination.parents:
        raise ValueError(f"WDL results cannot be placed inside the workflow source tree: {destination}")
    if destination.exists():
        raise FileExistsError(
            f"Refusing to reuse an existing WDL results directory: {destination}")

    sys.path.insert(0, str(RNAW))
    from python_lsf_pipeline.step1_prepare_samples import read_samples

    sample_sheet = sample_sheet.resolve(strict=True)
    rows = read_samples(sample_sheet)
    if len(rows) < 5:
        raise ValueError("The staged cohort needs at least five samples for cohort analyses")
    source_code = {name: STAGE_CODE / name for name in RESULT_CODE}
    missing_code = [path for path in source_code.values() if not path.is_file()]
    if missing_code:
        raise FileNotFoundError(f"Required staged workflow code is missing: {missing_code[0]}")

    destination.mkdir(parents=True)
    fastq_dir = destination / "fastqs"
    fastq_dir.mkdir()
    with (destination / "rna_cells.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["sample_id", "r1", "r2", "group"], delimiter="\t")
        writer.writeheader()
        for row in rows:
            sid = row["sample_id"]
            original_r1 = Path(row["r1"])
            for mate in ("r1", "r2"):
                source = Path(row[mate]).resolve(strict=True)
                link = fastq_dir / f"{sid}_{mate.upper()}.fastq.gz"
                link.symlink_to(source)
                row[mate] = str(link)
            provenance = original_r1.parent / f"{sid}.merge.json"
            if provenance.is_file():
                (fastq_dir / f"{sid}.merge.json").symlink_to(provenance.resolve())
            writer.writerow(row)

    for name in RESULT_CODE:
        shutil.copy2(source_code[name], destination / name)

    (destination / "wdl_run_manifest.json").write_text(json.dumps({
        "status": "PREPARED",
        "workflow_rna_dir": str(RNAW),
        "source_sample_sheet": str(sample_sheet),
        "results_dir": str(destination),
        "sample_count": len(rows),
        "input_files_are_symlinks": True,
    }, indent=2) + "\n")
    print(f"Prepared isolated WDL results tree: {destination}")
    print(f"Linked validated merged FASTQs for {len(rows)} samples; no FASTQ data copied.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir", type=Path, required=True,
        help="New results directory; it must not exist and cannot be inside workflow_RNA")
    parser.add_argument("--sample-sheet", type=Path, required=True,
                        help="TSV with sample_id, r1, r2, and group columns")
    args = parser.parse_args()
    prepare(args.results_dir, args.sample_sheet)


if __name__ == "__main__":
    main()
