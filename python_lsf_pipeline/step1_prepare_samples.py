#!/usr/bin/env python3
"""Prepare RNA manifests and validate paired FASTQ records without changing inputs."""
import argparse
import csv
import gzip
import re
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]


def read_samples(path):
    with Path(path).open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["sample_id", "r1", "r2", "group"]:
            raise ValueError("Expected TSV columns: sample_id, r1, r2, group")
        rows = list(reader)
    names = set()
    for row in rows:
        if None in row or any(not value for value in row.values()):
            raise ValueError("Empty field or malformed sample row")
        sid = row["sample_id"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sid) or sid in names:
            raise ValueError(f"Unsafe or duplicate sample ID: {sid}")
        names.add(sid)
        for mate in ("r1", "r2"):
            source = Path(row[mate])
            source = source if source.is_absolute() else Path(path).resolve().parent / source
            if any(c in str(source) for c in ("'", '"', "\n", "\r", " ", ";", "$", "`", "|")):
                raise ValueError(f"Unsupported shell characters in FASTQ path: {source}")
            if not source.is_file() or source.stat().st_size == 0:
                raise FileNotFoundError(source)
            row[mate] = str(source.resolve())
        if row["r1"] == row["r2"]:
            raise ValueError(f"Identical mate files: {sid}")
    if not rows:
        raise ValueError("No RNA samples")
    return rows


def discover(fastqs):
    rows = []
    for r1 in sorted(Path(fastqs).glob("*_RNA_*_R1.merged.fastq.gz")):
        sid = r1.name.removesuffix("_R1.merged.fastq.gz")
        rows.append({"sample_id": sid, "r1": str(r1.resolve()),
                     "r2": str(r1.with_name(sid + "_R2.merged.fastq.gz").resolve()),
                     "group": "SY5Y"})
    if not rows:
        raise ValueError(f"No merged RNA pairs found in {fastqs}")
    return rows


def write_manifest(rows, path):
    path = Path(path)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite a sample manifest: {path}")
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=["sample_id", "r1", "r2", "group"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    try:
        read_samples(temporary)
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def fastq_record(stream):
    header = stream.readline()
    if not header:
        return None
    sequence, plus, quality = (stream.readline().rstrip("\r\n") for _ in range(3))
    if not header.startswith("@") or not plus.startswith("+") or not sequence \
            or len(sequence) != len(quality):
        raise ValueError(f"Malformed FASTQ record: {header.strip()}")
    name = header.split()[0][1:]
    return re.sub(r"/[12]$", "", name)


def validate_pairs(r1, r2):
    count = 0
    with gzip.open(r1, "rt") as left, gzip.open(r2, "rt") as right:
        while True:
            a, b = fastq_record(left), fastq_record(right)
            if a is None and b is None:
                break
            if a is None or b is None or a != b:
                raise ValueError(f"Unpaired FASTQ at pair {count + 1}: {a!r}, {b!r}")
            count += 1
    if not count:
        raise ValueError("Empty FASTQ pair")
    return count


def write_source_csv(rows, path):
    with Path(path).open("w") as stream:
        writer = csv.writer(stream)
        writer.writerow(["biosampleName", "read1", "read2", "groups"])
        writer.writerows([row["sample_id"], row["r1"], row["r2"], row["group"]] for row in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastqs", type=Path, default=HERE.parent / "fastqs")
    parser.add_argument("--output", type=Path, default=HERE / "rna_cells.tsv")
    parser.add_argument("--validate-pair", nargs=2, metavar=("R1", "R2"))
    args = parser.parse_args()
    if args.validate_pair:
        print(validate_pairs(*args.validate_pair))
    else:
        write_manifest(discover(args.fastqs), args.output)
        print(f"Prepared {len(read_samples(args.output))} RNA samples: {args.output}")


if __name__ == "__main__":
    main()
