#!/usr/bin/env python3
"""Audit merged expanded-cohort FASTQs against their manifest and provenance."""
import argparse
import json
import os
from pathlib import Path
import sys

WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKFLOW_ROOT))
from python_lsf_pipeline.step1_prepare_samples import read_samples, validate_pairs


def audit_one(manifest_path, raw_plan_path, output_dir, index):
    rows = read_samples(manifest_path)
    if index < 1 or index > len(rows):
        raise ValueError(f"Sample index must be in 1..{len(rows)}: {index}")
    row = rows[index - 1]
    sample_id = row["sample_id"]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plan = json.loads(Path(raw_plan_path).read_text())
    planned = {cell["sample_id"]: cell for cell in plan["cells"]}
    if set(planned) != {sample["sample_id"] for sample in rows}:
        raise ValueError("Raw input plan and output manifest sample sets differ")
    expected = planned[sample_id]

    fastq_dir = Path(row["r1"]).parent
    provenance_path = fastq_dir / f"{sample_id}.merge.json"
    provenance = json.loads(provenance_path.read_text())
    if provenance.get("input") != expected:
        raise ValueError(f"Merge provenance input differs from raw plan for {sample_id}")

    lane_counts = provenance.get("source_pair_counts")
    lanes = expected.get("pairs", [])
    if (not isinstance(lane_counts, list) or len(lane_counts) != len(lanes)
            or any(not isinstance(count, int) or count <= 0 for count in lane_counts)):
        raise ValueError(f"Invalid source lane pair counts for {sample_id}")

    for pair in lanes:
        for mate in ("R1", "R2"):
            source = pair[mate]
            path = Path(source["path"])
            stat = path.stat()
            if (stat.st_size, stat.st_mtime_ns) != (source["size"], source["mtime_ns"]):
                raise ValueError(f"Raw input changed since planning: {path}")

    output_sizes = {}
    for key, path in (("r1", Path(row["r1"])), ("r2", Path(row["r2"]))):
        stat = path.stat()
        saved = provenance.get("outputs", {}).get(key)
        if saved != [stat.st_size, stat.st_mtime_ns]:
            raise ValueError(f"Output size/mtime differs from provenance for {path}")
        output_sizes[key] = stat.st_size

    merged_pairs = validate_pairs(row["r1"], row["r2"])
    expected_pairs = sum(lane_counts)
    if merged_pairs != expected_pairs or merged_pairs != provenance.get("merged_pairs"):
        raise ValueError(
            f"Pair-count mismatch for {sample_id}: merged={merged_pairs}, "
            f"source_sum={expected_pairs}, provenance={provenance.get('merged_pairs')}")

    partials = sorted(f"{path.name}" for path in fastq_dir.glob("*.partial"))
    if partials:
        raise ValueError(f"Partial output files exist: {partials}")

    result = {
        "status": "PASS",
        "sample_id": sample_id,
        "merged_pairs": merged_pairs,
        "source_lane_count": len(lanes),
        "source_pair_counts": lane_counts,
        "output_sizes": output_sizes,
        "provenance": str(provenance_path),
    }
    report_path = output_dir / f"{sample_id}.json"
    if report_path.exists():
        if json.loads(report_path.read_text()) != result:
            raise FileExistsError(f"Existing audit differs; inspect before replacing: {report_path}")
    else:
        temporary = report_path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(result, indent=2) + "\n")
        temporary.replace(report_path)
    print(f"[PASS] {sample_id}: {merged_pairs:,} read pairs")
    return result


def combine(manifest_path, raw_plan_path, output_dir, report_path):
    rows = read_samples(manifest_path)
    output_dir = Path(output_dir)
    sample_reports = []
    for index, row in enumerate(rows, start=1):
        report = output_dir / f"{row['sample_id']}.json"
        if not report.is_file():
            raise FileNotFoundError(f"Missing per-sample audit: {report}")
        result = json.loads(report.read_text())
        if result.get("status") != "PASS" or result.get("sample_id") != row["sample_id"]:
            raise ValueError(f"Invalid per-sample audit: {report}")
        for key, path in (("r1", row["r1"]), ("r2", row["r2"])):
            if result.get("output_sizes", {}).get(key) != Path(path).stat().st_size:
                raise ValueError(f"FASTQ changed since per-sample audit: {path}")
        sample_reports.append(result)

    fastq_dir = Path(rows[0]["r1"]).parent
    expected_fastqs = {row[key] for row in rows for key in ("r1", "r2")}
    observed_fastqs = {str(path.resolve()) for path in fastq_dir.glob("*.fastq.gz")}
    if expected_fastqs != observed_fastqs:
        raise ValueError(
            f"FASTQ inventory mismatch: expected {len(expected_fastqs)}, "
            f"observed {len(observed_fastqs)}")
    expected_provenance = {f"{row['sample_id']}.merge.json" for row in rows}
    observed_provenance = {path.name for path in fastq_dir.glob("*.merge.json")}
    if expected_provenance != observed_provenance:
        raise ValueError(
            f"Provenance inventory mismatch: expected {len(expected_provenance)}, "
            f"observed {len(observed_provenance)}")
    if list(fastq_dir.glob("*.partial")):
        raise ValueError("Partial FASTQ outputs remain in the output directory")

    raw_plan = json.loads(Path(raw_plan_path).read_text())
    if {cell["sample_id"] for cell in raw_plan.get("cells", [])} != {
            row["sample_id"] for row in rows}:
        raise ValueError("Raw input plan sample set differs from the output manifest")
    result = {
        "status": "PASS",
        "samples": len(rows),
        "merged_fastqs": len(observed_fastqs),
        "provenance_files": len(observed_provenance),
        "total_read_pairs": sum(item["merged_pairs"] for item in sample_reports),
        "samples_audited": sample_reports,
    }
    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = report_path.with_suffix(report_path.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2) + "\n")
    temporary.replace(report_path)
    print(f"[PASS] {len(rows)} samples, {len(observed_fastqs)} merged FASTQs; "
          f"audit: {report_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--raw-plan", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--index", type=int, default=int(os.environ.get("LSB_JOBINDEX", "0")))
    parser.add_argument("--combine-report", type=Path)
    args = parser.parse_args()
    raw_plan = args.raw_plan or args.manifest.with_name("raw_inputs.json")
    output_dir = args.output_dir or args.manifest.parent / "audit/fastq"
    if args.combine_report:
        combine(args.manifest, raw_plan, output_dir, args.combine_report)
    elif args.index > 0:
        audit_one(args.manifest, raw_plan, output_dir, args.index)
    else:
        parser.error("Provide --index (or run in an LSF array) or --combine-report")


if __name__ == "__main__":
    main()
