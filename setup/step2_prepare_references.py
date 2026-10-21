#!/usr/bin/env python3
"""Build a private complete STAR index and Picard refFlat from matching Ensembl 104 files."""
import argparse
import hashlib
import gzip
import json
import os
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
RNA_WORKFLOW = HERE.parent
RUNTIME = Path(os.environ.get("RNA_RUNTIME_DIR", RNA_WORKFLOW / "wdl_runtime"))
REFERENCE_DIR = RUNTIME / "references/Ensembl104"
DEFAULTS = {
    "fasta": str(REFERENCE_DIR / "GRCh38.primary.fa"),
    "gtf": str(REFERENCE_DIR / "Homo_sapiens.GRCh38.104.chr.gtf"),
    "tx2gene": str(REFERENCE_DIR / "Ensembl104.tx2gene.tsv"),
    "genebody": str(REFERENCE_DIR / "Ensembl104.transcripts.bed"),
    "star_index": str(REFERENCE_DIR / "STAR"),
    "ref_flat": str(REFERENCE_DIR / "Ensembl104.refFlat.txt"),
}


def require_file(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing/empty reference: {path}")


def check_references(refs, require_index=True):
    for name in ("fasta", "gtf", "tx2gene", "genebody"):
        require_file(refs[name])
    require_file(refs["fasta"] + ".fai")
    if require_index:
        require_file(refs["ref_flat"])
        for name in ("Genome", "SA", "SAindex", "genomeParameters.txt", "chrNameLength.txt"):
            require_file(Path(refs["star_index"]) / name)
        fai = {row.split("\t")[0]: int(row.split("\t")[1])
               for row in Path(refs["fasta"] + ".fai").read_text().splitlines()}
        index = dict(line.split("\t") for line in
                     (Path(refs["star_index"]) / "chrNameLength.txt").read_text().splitlines())
        if {key: int(value) for key, value in index.items()} != fai:
            raise ValueError("STAR index contigs/lengths do not match the FASTA")


def parse_attributes(text):
    return dict(re.findall(r'(\w+) "([^"]*)"', text))


def create_tx2gene(gtf, output):
    transcripts = {}
    with Path(gtf).open() as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in ("transcript", "exon"):
                continue
            attributes = parse_attributes(fields[8])
            transcript_id = attributes.get("transcript_id")
            gene_id = attributes.get("gene_id")
            if not transcript_id or not gene_id:
                raise ValueError("GTF transcript/exon lacks transcript_id or gene_id")
            transcripts[transcript_id] = (
                attributes.get("transcript_biotype", attributes.get("transcript_type", "unknown")),
                gene_id,
                attributes.get("gene_biotype", attributes.get("gene_type", "unknown")),
                attributes.get("gene_name", gene_id),
            )
    if not transcripts:
        raise ValueError(f"No transcript mappings found in GTF: {gtf}")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w") as stream:
        for transcript_id, values in sorted(transcripts.items()):
            stream.write("\t".join((transcript_id, *values)) + "\n")
    require_file(temporary)
    temporary.replace(output)
    return len(transcripts)


def create_genebody_bed(gtf, output):
    transcripts = {}
    with Path(gtf).open() as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attributes = parse_attributes(fields[8])
            transcript_id = attributes.get("transcript_id")
            gene_id = attributes.get("gene_id")
            if not transcript_id or not gene_id:
                raise ValueError("GTF exon lacks transcript_id or gene_id")
            start, end = int(fields[3]) - 1, int(fields[4])
            if start < 0 or end <= start:
                raise ValueError("Invalid GTF exon interval")
            key = (transcript_id, fields[0], fields[6], gene_id)
            transcripts.setdefault(key, []).append((start, end))
    if not transcripts:
        raise ValueError(f"No exons found in GTF: {gtf}")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w") as stream:
        for (transcript_id, chrom, strand, gene_id), exons in sorted(
                transcripts.items(), key=lambda item: (item[0][1], min(x[0] for x in item[1]))):
            exons = sorted(set(exons))
            start, end = exons[0][0], max(interval[1] for interval in exons)
            sizes = ",".join(str(stop - begin) for begin, stop in exons) + ","
            starts = ",".join(str(begin - start) for begin, _ in exons) + ","
            name = f"{gene_id}|{transcript_id}"
            stream.write("\t".join(map(str, (
                chrom, start, end, name, 0, strand, start, end, 0,
                len(exons), sizes, starts,
            ))) + "\n")
    require_file(temporary)
    temporary.replace(output)
    return len(transcripts)


def combine_transcript_fastas(inputs, output):
    seen = {}
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w") as out:
        for input_path in inputs:
            with gzip.open(input_path, "rt") as source:
                header = None
                sequence = []

                def emit_record():
                    if header is None:
                        return
                    transcript_id = header[1:].split(None, 1)[0]
                    if not transcript_id:
                        raise ValueError(f"Empty transcript ID in {input_path}")
                    sequence_text = "".join(sequence)
                    sequence_hash = hashlib.sha256(sequence_text.encode()).hexdigest()
                    prior = seen.get(transcript_id)
                    if prior is not None:
                        if prior != sequence_hash:
                            raise ValueError(
                                f"Conflicting duplicate transcript sequence: {transcript_id}")
                        return
                    seen[transcript_id] = sequence_hash
                    out.write(header + "\n")
                    for offset in range(0, len(sequence_text), 80):
                        out.write(sequence_text[offset:offset + 80] + "\n")

                for line in source:
                    if line.startswith(">"):
                        emit_record()
                        header, sequence = line.rstrip("\n"), []
                    elif header is not None:
                        sequence.append(line.strip())
                emit_record()
    require_file(temporary)
    temporary.replace(output)
    return len(seen)


def create_ref_flat(gtf, output):
    transcripts = {}
    with Path(gtf).open() as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in ("exon", "CDS"):
                continue
            attributes = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            tid = attributes.get("transcript_id")
            if not tid:
                raise ValueError("GTF exon/CDS lacks transcript_id")
            key = (attributes.get("gene_name", attributes["gene_id"]), tid,
                   fields[0], fields[6])
            record = transcripts.setdefault(key, {"exons": set(), "cds": []})
            start, end = int(fields[3]) - 1, int(fields[4])
            if start < 0 or end <= start:
                raise ValueError("Invalid GTF interval")
            if fields[2] == "exon":
                record["exons"].add((start, end))
            else:
                record["cds"].append((start, end))
    output = Path(output)
    temporary = output.with_suffix(".tmp")
    with temporary.open("w") as stream:
        for key, record in sorted(transcripts.items()):
            exons = sorted(record["exons"])
            if not exons:
                raise ValueError(f"Transcript has CDS but no exons: {key[1]}")
            start, end = exons[0][0], max(e[1] for e in exons)
            cds = record["cds"]
            cds_start = min(e[0] for e in cds) if cds else end
            cds_end = max(e[1] for e in cds) if cds else end
            stream.write("\t".join(map(str, [*key, start, end, cds_start, cds_end,
                                             len(exons), "".join(f"{e[0]}," for e in exons),
                                             "".join(f"{e[1]}," for e in exons)])) + "\n")
    require_file(temporary)
    temporary.replace(output)
    return len(transcripts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key, value in DEFAULTS.items():
        parser.add_argument("--" + key.replace("_", "-"), default=value)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--sjdb-overhang", type=int, default=100)
    parser.add_argument("--memory-mb", type=int, default=64000)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--derive-annotation", action="store_true",
                        help="Build tx2gene, transcript BED12 and refFlat from the GTF")
    parser.add_argument("--combine-transcripts", action="store_true")
    parser.add_argument("--cdna-fasta", type=Path)
    parser.add_argument("--ncrna-fasta", type=Path)
    parser.add_argument("--transcript-fasta", type=Path)
    args = parser.parse_args()
    if args.combine_transcripts:
        if not all((args.cdna_fasta, args.ncrna_fasta, args.transcript_fasta)):
            parser.error("--combine-transcripts requires all three transcript FASTA paths")
        count = combine_transcript_fastas(
            (args.cdna_fasta, args.ncrna_fasta), args.transcript_fasta)
        print(f"Combined unique transcript sequences: {count:,}")
        return
    refs = {key: str(Path(getattr(args, key)).resolve()) for key in DEFAULTS}
    if args.derive_annotation:
        print(f"tx2gene mappings: {create_tx2gene(refs['gtf'], refs['tx2gene']):,}")
        print(f"gene-body transcript models: {create_genebody_bed(refs['gtf'], refs['genebody']):,}")
    check_references(refs, require_index=False)
    if args.threads < 1 or args.memory_mb < 48000 or args.sjdb_overhang < 1:
        parser.error("Reference generation needs positive threads/overhang and at least 48000 MB")
    if args.build:
        index = Path(refs["star_index"])
        if index.exists() and any(index.iterdir()):
            check_references(refs)
            print("Existing complete reference retained; no overwrite")
            return
        temporary_index = index.with_name(index.name + ".building")
        if temporary_index.exists():
            raise FileExistsError(
                f"Interrupted STAR index build exists at {temporary_index}; inspect it before removing")
        temporary_index.mkdir(parents=True)
        subprocess.run(["STAR", "--runMode", "genomeGenerate", "--runThreadN", str(args.threads),
                        "--genomeDir", str(temporary_index), "--genomeFastaFiles", refs["fasta"],
                        "--sjdbGTFfile", refs["gtf"], "--sjdbOverhang", str(args.sjdb_overhang),
                        "--limitGenomeGenerateRAM", str(args.memory_mb * 1000000 - 4000000000)],
                       check=True, cwd=temporary_index)
        create_ref_flat(refs["gtf"], refs["ref_flat"])
        temporary_index.rename(index)
    check_references(refs)
    manifest = Path(refs["star_index"]).parent / "reference_manifest.json"
    manifest.write_text(json.dumps(refs, indent=2) + "\n")
    print(json.dumps(refs, indent=2))


if __name__ == "__main__":
    main()
