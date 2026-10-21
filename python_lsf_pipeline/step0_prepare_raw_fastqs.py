#!/usr/bin/env python3
"""Plan and validate per-cell merges of RNA lane/flowcell FASTQ pairs."""
import argparse
import csv
import fcntl
import gzip
import json
from pathlib import Path
import re
import shutil

from python_lsf_pipeline.step1_prepare_samples import fastq_record, validate_pairs


def plan(raw_dir, outdir, fastq_folder="fastqs"):
    raw_dir, outdir = Path(raw_dir).resolve(), Path(outdir).resolve()
    if fastq_folder not in {"fastq", "fastqs"}:
        raise ValueError("FASTQ output folder must be 'fastq' or 'fastqs'")
    groups = {}
    for path in sorted(raw_dir.glob("*RNA*.fastq.gz")):
        match = re.fullmatch(r"(.+)_S(\d+)_(L\d{3})_(R[12])_(\d+)\.fastq\.gz", path.name)
        if not match:
            raise ValueError(f"Unsupported RNA FASTQ name: {path.name}")
        label, sample, lane, mate, part = match.groups()
        label = label.replace("-", "_")
        identity = re.fullmatch(
            r"(SJNBL\d+)_(C\d+)_ResolOME_SY5Y_P[234]_RNA_C(\d+)", label)
        if not identity:
            raise ValueError(f"Unsupported RNA library identity: {label}")
        patient, library, cell = identity.groups()
        sid = f"SY5Y_Cell{int(cell)}_RNA_{patient}_{library}_RNA_C{int(cell)}"
        key = (sid, int(sample), lane, part)
        entry = groups.setdefault(key, {})
        if mate in entry:
            raise ValueError(f"Duplicate mate for {key}: {path}")
        entry[mate] = path
    if not groups:
        raise ValueError(f"No RNA FASTQs in {raw_dir}")
    cells, physical = {}, set()
    for (sid, sample, lane, part), mates in sorted(groups.items()):
        if set(mates) != {"R1", "R2"}:
            raise ValueError(f"Missing mate: {sid}, S{sample}, {lane}, {part}")
        with gzip.open(mates["R1"], "rt") as left, gzip.open(mates["R2"], "rt") as right:
            name = fastq_record(left)
            if name is None or name != fastq_record(right):
                raise ValueError(f"First read is not paired: {mates}")
        fields = name.split(":")
        if len(fields) < 7 or not fields[3].isdigit() or int(fields[3]) != int(lane[1:]):
            raise ValueError(f"Illumina read header and filename lane disagree: {mates['R1']}")
        run_lane = ":".join(fields[:4])
        key = (sid, run_lane, part)
        if key in physical:
            raise ValueError(f"Duplicate physical run/lane input: {key}")
        physical.add(key)
        cell = cells.setdefault(sid, {"sample_id": sid, "group": "SY5Y", "pairs": [],
                                     "r1": str(outdir / fastq_folder / f"{sid}_R1.merged.fastq.gz"),
                                     "r2": str(outdir / fastq_folder / f"{sid}_R2.merged.fastq.gz")})
        cell["pairs"].append({
            "sample_number": sample, "lane": lane, "run_lane": run_lane, "part": part,
            **{mate: {"path": str(path), "size": path.stat().st_size,
                      "mtime_ns": path.stat().st_mtime_ns} for mate, path in mates.items()}})
    return {"raw_directory": str(raw_dir), "cells": list(cells.values())}


def write_plan(manifest, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(manifest, indent=2) + "\n"
    path = outdir / "raw_inputs.json"
    if path.exists():
        if path.read_text() != encoded:
            raise ValueError("Raw input plan changed; use a new output directory")
    else:
        with path.open("x") as stream:
            stream.write(encoded)
    table = outdir / "rna_cells.tsv"
    if not table.exists():
        temporary = table.with_suffix(".tsv.tmp")
        with temporary.open("x") as stream:
            writer = csv.DictWriter(stream, delimiter="\t",
                                    fieldnames=["sample_id", "r1", "r2", "group"])
            writer.writeheader()
            writer.writerows({k: cell[k] for k in writer.fieldnames} for cell in manifest["cells"])
        temporary.replace(table)
    else:
        with table.open() as stream:
            existing = list(csv.DictReader(stream, delimiter="\t"))
        expected = [{k: cell[k] for k in ("sample_id", "r1", "r2", "group")}
                    for cell in manifest["cells"]]
        if existing != expected:
            raise ValueError("Existing RNA cell manifest differs from raw plan")


def merge_cell(manifest_path, index):
    manifest_path = Path(manifest_path)
    cells = json.loads(manifest_path.read_text())["cells"]
    if index < 1 or index > len(cells):
        raise ValueError(f"Cell index outside 1..{len(cells)}: {index}")
    cell = cells[index - 1]
    root = Path(cell["r1"]).parent
    root.mkdir(parents=True, exist_ok=True)
    marker = root / f'{cell["sample_id"]}.merge.json'
    with (root / f'{cell["sample_id"]}.lock').open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for pair in cell["pairs"]:
            for mate in ("R1", "R2"):
                item = pair[mate]
                stat = Path(item["path"]).stat()
                if (stat.st_size, stat.st_mtime_ns) != (item["size"], item["mtime_ns"]):
                    raise ValueError(f'Raw input changed: {item["path"]}')
        if marker.exists():
            saved = json.loads(marker.read_text())
            if saved["input"] != cell:
                raise ValueError(f"Existing merge provenance changed: {marker}")
            for mate in ("r1", "r2"):
                stat = Path(cell[mate]).stat()
                if [stat.st_size, stat.st_mtime_ns] != saved["outputs"][mate]:
                    raise ValueError(f"Existing merged input changed: {cell[mate]}")
            print(f'[OK] Validated merge retained: {cell["sample_id"]}', flush=True)
            return
        targets = [Path(cell[m]) for m in ("r1", "r2")]
        temporary = [p.with_suffix(p.suffix + ".partial") for p in targets]
        if any(p.exists() for p in targets + temporary):
            raise FileExistsError(f'Uncommitted merge exists: {cell["sample_id"]}; inspect before retry')
        counts = []
        try:
            with temporary[0].open("xb") as left, temporary[1].open("xb") as right:
                for pair in cell["pairs"]:
                    count = validate_pairs(pair["R1"]["path"], pair["R2"]["path"])
                    counts.append(count)
                    for mate, output in (("R1", left), ("R2", right)):
                        with Path(pair[mate]["path"]).open("rb") as stream:
                            shutil.copyfileobj(stream, output, length=8 * 1024 * 1024)
                    print(f'[OK] {cell["sample_id"]} {pair["run_lane"]}: {count} pairs', flush=True)
            merged_count = validate_pairs(*temporary)
            if merged_count != sum(counts):
                raise ValueError("Merged FASTQ count does not equal source lane counts")
            for source, target in zip(temporary, targets):
                source.replace(target)
            saved = {"input": cell, "source_pair_counts": counts, "merged_pairs": merged_count,
                     "outputs": {m: [Path(cell[m]).stat().st_size, Path(cell[m]).stat().st_mtime_ns]
                                 for m in ("r1", "r2")}}
            with marker.open("x") as stream:
                stream.write(json.dumps(saved, indent=2) + "\n")
        finally:
            for path in temporary:
                if path.exists():
                    path.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--raw-dir", type=Path)
    modes.add_argument("--merge-index", type=int)
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--fastq-folder", choices=("fastq", "fastqs"), default="fastqs")
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()
    if args.raw_dir:
        if not args.outdir:
            parser.error("--raw-dir requires --outdir")
        manifest = plan(args.raw_dir, args.outdir, args.fastq_folder)
        write_plan(manifest, args.outdir)
        print(f'[OK] Planned {len(manifest["cells"])} RNA cells; no raw input modified')
    else:
        if not args.manifest:
            parser.error("--merge-index requires --manifest")
        merge_cell(args.manifest, args.merge_index)


if __name__ == "__main__":
    main()
