#!/usr/bin/env python3
"""Move completed RNA results into a navigable tree without copying data."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

RNAW = Path(__file__).resolve().parents[1]

DIRECTORY_MOVES = {
    "fastqs": "inputs/merged_fastqs",
    "stage3_qc_preprocessing": "steps/qc_depth/step3_preprocessing",
    "stage4_qc_alignment": "steps/qc_depth/step4_alignment",
    "stage5_qc_counts": "steps/qc_depth/step5_counts",
    "stage6_qc_matrix": "steps/qc_depth/step6_matrix",
    "stage7_qc_reports_v2": "steps/qc_depth/step7_reports",
    "stage7_qc_reports": "steps/qc_depth/superseded/step7_reports_initial",
    "stage3_full_preprocessing": "steps/full_depth/step3_preprocessing",
    "stage4_full_alignment": "steps/full_depth/step4_alignment",
    "stage5_full_counts": "steps/full_depth/step5_counts",
    "stage6_full_matrix": "steps/full_depth/step6_matrix",
    "stage7_full_reports": "steps/full_depth/step7_reports",
}

LOG_MOVES = {
    "logs_driver": "logs/driver",
    "logs_full_continuation": "logs/full_depth/continuation",
    "logs_full_step3": "logs/full_depth/step3_preprocessing",
    "logs_full_step3_audit": "logs/full_depth/step3_audit",
    "logs_full_step4": "logs/full_depth/step4_alignment",
    "logs_step1": "logs/input_merge",
    "logs_step3": "logs/qc_depth/step3_preprocessing",
    "logs_step4": "logs/qc_depth/step4_alignment",
    "logs_step5": "logs/qc_depth/step5_counts",
    "logs_step6": "logs/qc_depth/step6_matrix",
    "logs_step7": "logs/qc_depth/step7_reports",
    "logs_validation": "logs/validation",
}

FILE_MOVES = {
    "step1_input_inspection.tsv": "inputs/manifests/step1_input_inspection.tsv",
    "raw_inputs.json": "inputs/manifests/raw_inputs.json",
    "rna_cells.tsv": "inputs/manifests/rna_cells.tsv",
    "step2_reference_audit.json": "audits/reference/step2_reference_audit.json",
    "tool_preflight.json": "audits/preflight/tool_preflight.json",
    "step3_output_audit.json": "audits/qc_depth/step3_output_audit.json",
    "step4_output_audit.json": "audits/qc_depth/step4_output_audit.json",
    "step5_output_audit.json": "audits/qc_depth/step5_output_audit.json",
    "step6_output_audit.json": "audits/qc_depth/step6_output_audit.json",
    "step7_output_audit.json": "audits/qc_depth/step7_output_audit.json",
    "step3_full_output_audit.json": "audits/full_depth/step3_full_output_audit.json",
    "step3_full_output_audit_quick.json": "audits/full_depth/step3_full_output_audit_quick.json",
    "step4_full_output_audit_quick.json": "audits/full_depth/step4_full_output_audit_quick.json",
    "step5_full_output_audit.json": "audits/full_depth/step5_full_output_audit.json",
    "step6_full_output_audit.json": "audits/full_depth/step6_full_output_audit.json",
    "step7_full_output_audit.json": "audits/full_depth/step7_full_output_audit.json",
    "full_staged_status.txt": "logs/full_depth/staged_status.txt",
}

PACKAGE_LOG_MOVES = {
    "final_results/07_qc_provenance/logs_downstream": "logs/downstream/driver",
    "final_results/07_qc_provenance/salmon_logs": "logs/downstream/salmon",
}

PACKAGE_AUDIT_LINKS = {
    "step4_full_output_audit_quick.json": "audits/full_depth/step4_full_output_audit_quick.json",
    "step5_full_output_audit.json": "audits/full_depth/step5_full_output_audit.json",
    "step6_full_output_audit.json": "audits/full_depth/step6_full_output_audit.json",
    "step7_full_output_audit.json": "audits/full_depth/step7_full_output_audit.json",
    "salmon_output_audit.json": "final_results/07_qc_provenance/salmon_output_audit.json",
}

PACKAGE_REFERENCE_LINKS = {
    "final_results/09_references/ReactomePathways.gmt.zip":
        "analysis_resources/ReactomePathways.gmt.zip",
    "final_results/09_references/README.txt": "analysis_resources/README.txt",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def destination(root: Path, relative: str) -> Path:
    return root / relative


def source_for_move(root: Path, source: Path, target: Path) -> Path:
    if source.is_symlink():
        resolved = source.resolve(strict=True)
        if resolved != target.resolve(strict=True):
            raise RuntimeError(f"Compatibility link points to unexpected target: {source} -> {resolved}")
        return target
    if source.exists():
        if target.exists():
            raise FileExistsError(f"Both source and destination exist: {source} ; {target}")
        return source
    if target.exists():
        return target
    raise FileNotFoundError(f"Neither source nor destination exists: {source} ; {target}")


def check_moves(root: Path) -> list[tuple[Path, Path]]:
    moves = []
    for mapping in (DIRECTORY_MOVES, LOG_MOVES, FILE_MOVES):
        for old, new in mapping.items():
            src, dst = root / old, destination(root, new)
            if not src.exists() and not src.is_symlink() and not dst.exists():
                continue
            actual = source_for_move(root, src, dst)
            if actual != dst:
                if actual.stat().st_dev != root.stat().st_dev:
                    raise OSError(f"Move would cross filesystems: {actual}")
                moves.append((actual, dst))
    for old, new in PACKAGE_LOG_MOVES.items():
        src, dst = destination(root, old), destination(root, new)
        if not src.exists() and not src.is_symlink() and not dst.exists():
            continue
        actual = source_for_move(root, src, dst)
        if actual != dst:
            if actual.stat().st_dev != root.stat().st_dev:
                raise OSError(f"Move would cross filesystems: {actual}")
            moves.append((actual, dst))
    return moves


def check_duplicate_links(root: Path) -> list[tuple[Path, Path]]:
    links = []
    for mapping, canonical_root in (
        (PACKAGE_AUDIT_LINKS, root),
        (PACKAGE_REFERENCE_LINKS, RNAW),
    ):
        for old, new in mapping.items():
            duplicate = root / "final_results/07_qc_provenance/source_audits" / old \
                if mapping is PACKAGE_AUDIT_LINKS else root / old
            canonical = canonical_root / new
            if not canonical.exists():
                for move_source, move_target in FILE_MOVES.items():
                    if move_target == new:
                        canonical = root / move_source
                        break
            if not canonical.is_file():
                raise FileNotFoundError(f"Canonical reference is missing: {canonical}")
            if not duplicate.exists() and not duplicate.is_symlink():
                raise FileNotFoundError(f"Expected package reference is missing: {duplicate}")
            if duplicate.is_file() and sha256(duplicate) != sha256(canonical):
                raise ValueError(f"Refusing to replace non-identical package copy: {duplicate}")
            links.append((duplicate, canonical))
    for old, new in PACKAGE_LOG_MOVES.items():
        alias, canonical = root / old, destination(root, new)
        if not alias.exists() and not alias.is_symlink():
            continue
        if canonical.exists() and alias.is_dir() and not alias.is_symlink():
            raise FileExistsError(f"Both package log source and destination exist: {alias} ; {canonical}")
    return links


def make_link(link: Path, target: Path) -> None:
    target = target.resolve(strict=True)
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.is_symlink():
        if link.resolve(strict=True) == target:
            return
        link.unlink()
    elif link.exists():
        raise FileExistsError(f"Refusing to replace non-link path: {link}")
    link.symlink_to(os.path.relpath(target, link.parent))


def move_path(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        raise FileExistsError(f"Move destination appeared after preflight: {target}")
    os.rename(source, target)


def write_results_readme(root: Path) -> None:
    manifest = root / "rna_cells.tsv"
    if not manifest.is_file():
        raise FileNotFoundError(f"Organized sample manifest is missing: {manifest}")
    with manifest.open(newline="") as stream:
        sample_count = sum(1 for _ in csv.DictReader(stream, delimiter="\t"))
    if sample_count < 1:
        raise ValueError(f"Sample manifest has no rows: {manifest}")
    quick_audit_path = root / "step4_full_output_audit_quick.json"
    if quick_audit_path.is_file():
        quick_audit = json.loads(quick_audit_path.read_text())
        deep_count = len(quick_audit.get("deep_checked_cells", []))
        audit_note = (
            f"The full-depth alignment audit used {quick_audit.get('audit_mode', 'unknown')} "
            f"mode and deeply checked {deep_count} sample(s); consult the audit for its exact "
            "scope.")
    else:
        audit_note = "Consult the full-depth alignment audit for the exact validation scope."
    text = f"""# ResolveOME RNA results

The completed RNA workflow outputs are organized by input, depth branch, stage,
audit, and log type. `final_results/` is the user-facing cohort analysis package;
open `final_results/00_dashboard/index.html` for the figures and interpretation.

## Folder map

- `inputs/merged_fastqs/` — the validated merged paired-end FASTQs.
- `inputs/manifests/` — original-input inventory and merge inspection outputs.
- `steps/qc_depth/` — sampled-read QC branch, with one folder per pipeline step.
- `steps/full_depth/` — full-depth expression branch, with one folder per step.
- `steps/qc_depth/superseded/` — preserved first QC-report iteration; the accepted
  reports are in `steps/qc_depth/step7_reports/`.
- `audits/reference/`, `audits/qc_depth/`, `audits/full_depth/`, and
  `audits/preflight/` — canonical audit reports.
- `logs/` — driver, per-stage, validation, and downstream Salmon/analysis logs.
- `final_results/` — dashboard, figures, matrices, transcript quantification,
  junction/pathway summaries, QC, and provenance.

Each stage retains its Nextflow `work/` tree and published outputs for provenance
and resume support. These were moved in place, not copied or removed. Identical
audit and reference assets in the final package are symbolic links to their one
canonical copy.

For compatibility, former root paths such as `stage4_full_alignment`,
`logs_step4`, `step4_full_output_audit_quick.json`, `fastqs`, and `rna_cells.tsv`
remain as symbolic links to the organized locations. They do not create
additional data copies.

All {sample_count} samples are retained regardless of QC label. {audit_note}
This is an analysis output package, not an official St. Jude publication or
endorsement.
"""
    (root / "README.md").write_text(text)


def apply(root: Path, moves: list[tuple[Path, Path]],
          duplicate_links: list[tuple[Path, Path]]) -> None:
    for source, target in moves:
        move_path(source, target)

    for old, new in DIRECTORY_MOVES.items():
        if destination(root, new).exists():
            make_link(root / old, destination(root, new))
    for old, new in LOG_MOVES.items():
        if destination(root, new).exists():
            make_link(root / old, destination(root, new))
    for old, new in FILE_MOVES.items():
        if destination(root, new).exists():
            make_link(root / old, destination(root, new))

    for old, new in PACKAGE_LOG_MOVES.items():
        if destination(root, new).exists():
            make_link(root / old, destination(root, new))
    for duplicate, canonical in duplicate_links:
        if duplicate.is_symlink() and duplicate.resolve(strict=True) == canonical.resolve(strict=True):
            continue
        if duplicate.exists():
            duplicate.unlink()
        make_link(duplicate, canonical)

    write_results_readme(root)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--plan", action="store_true", help="validate and print the planned organization")
    mode.add_argument("--apply", action="store_true", help="move completed outputs and create compatibility links")
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    status = root / "final_results/07_qc_provenance/downstream_run_status.txt"
    if not status.is_file() or status.read_text().split(maxsplit=1)[0] != "PASS":
        raise RuntimeError(f"Downstream final-results status is not PASS: {status}")

    moves = check_moves(root)
    duplicate_links = check_duplicate_links(root)
    print(f"Validated result root: {root}")
    print(f"Same-filesystem moves (no copy): {len(moves)}")
    for source, target in moves:
        print(f"MOVE {source.relative_to(root)} -> {target.relative_to(root)}")
    print(f"Compatibility links: {len(DIRECTORY_MOVES) + len(LOG_MOVES) + len(FILE_MOVES)}")
    print(f"Verified duplicate references to replace with links: {len(duplicate_links)}")
    for duplicate, canonical in duplicate_links:
        print(f"LINK {duplicate.relative_to(root)} -> {canonical}")
    if args.apply:
        apply(root, moves, duplicate_links)
        print("Organization applied successfully.")
    else:
        print("Plan only; no filesystem changes made.")


if __name__ == "__main__":
    main()
