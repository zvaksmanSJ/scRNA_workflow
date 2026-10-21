#!/usr/bin/env python3
"""Opt-in synthetic RNA integration test; never reads production FASTQs."""
import argparse
import csv
import gzip
import random
import subprocess
import sys
from pathlib import Path

WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKFLOW_ROOT / "python_lsf_pipeline"))
sys.path.insert(0, str(WORKFLOW_ROOT / "setup"))

from step1_prepare_samples import write_manifest
from step2_prepare_references import create_ref_flat

HK = ["GAPDH", "ACTB", "RPL36", "HINT1", "TBP", "PPIA", "B2M", "HPRT1",
      "UBC", "RPL13A", "PGK1", "EIF3K", "RPLP0", "GUSB", "CLTC", "HMBS"]


def fixture(root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=False)
    rng = random.Random(1928)
    sequence = "".join(rng.choices("ACGT", k=2800000))
    fasta = root / "genome.fa"
    with fasta.open("w") as stream:
        stream.write(">chr1\n")
        for start in range(0, len(sequence), 60):
            stream.write(sequence[start:start + 60] + "\n")
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    genes = []
    gtf = root / "genes.gtf"
    mapping = root / "tx2gene.tsv"
    bed = root / "housekeeping.bed"
    with gtf.open("w") as annotation, mapping.open("w") as tx, bed.open("w") as body:
        for i in range(1400):
            start, end = i * 2000 + 200, i * 2000 + 1000
            gid, tid = f"ENSG{i:011d}", f"ENST{i:011d}"
            symbol = HK[i] if i < len(HK) else f"GENE{i}"
            attrs = (f'gene_id "{gid}"; transcript_id "{tid}"; gene_name "{symbol}"; '
                     'gene_biotype "protein_coding"; transcript_biotype "protein_coding";')
            for kind in ("transcript", "exon"):
                annotation.write(f"chr1\tsynthetic\t{kind}\t{start + 1}\t{end}\t.\t+\t.\t{attrs}\n")
            tx.write(f"{tid}\tprotein_coding\t{gid}\tprotein_coding\t{symbol}\n")
            body.write(f"chr1\t{start}\t{end}\t{symbol}\t0\t+\t{start}\t{end}\t0\t1\t800,\t0,\n")
            genes.append((start, end))
    flat = root / "refFlat.txt"
    create_ref_flat(gtf, flat)
    index = root / "star"
    index.mkdir()
    subprocess.run(["STAR", "--runMode", "genomeGenerate", "--runThreadN", "2",
                    "--genomeDir", str(index), "--genomeFastaFiles", str(fasta),
                    "--sjdbGTFfile", str(gtf), "--sjdbOverhang", "100",
                    "--genomeSAindexNbases", "9", "--genomeChrBinNbits", "12"],
                   check=True, cwd=index, stdout=subprocess.DEVNULL)
    rows = []
    complement = str.maketrans("ACGT", "TGCA")
    for cell in range(6):
        sid = f"SyntheticRNA{cell + 1}"
        left, right = root / f"{sid}_R1.gz", root / f"{sid}_R2.gz"
        with gzip.open(left, "wt") as r1, gzip.open(right, "wt") as r2:
            number = 0
            for start, end in genes:
                for _ in range(rng.randint(12, 28)):
                    number += 1
                    position = rng.randint(start, end - 400)
                    a = sequence[position:position + 150]
                    b = sequence[position + 250:position + 400].translate(complement)[::-1]
                    r1.write(f"@{sid}:{number}/1\n{a}\n+\n{'I' * 150}\n")
                    r2.write(f"@{sid}:{number}/2\n{b}\n+\n{'I' * 150}\n")
        rows.append({"sample_id": sid, "r1": str(left), "r2": str(right), "group": "synthetic"})
    manifest = root / "samples.tsv"
    write_manifest(rows, manifest)
    return ["--cells", str(manifest), "--fasta", str(fasta), "--gtf", str(gtf),
            "--tx2gene", str(mapping), "--genebody", str(bed),
            "--ref-flat", str(flat), "--star-index", str(index)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--mode", choices=("qc", "full", "both"), default="both")
    parser.add_argument("--reuse-fixture", action="store_true")
    parser.add_argument("--output-name", default="output")
    args = parser.parse_args()
    root = args.directory.resolve()
    inputs = (["--cells", str(root / "samples.tsv"), "--fasta", str(root / "genome.fa"),
               "--gtf", str(root / "genes.gtf"), "--tx2gene", str(root / "tx2gene.tsv"),
               "--genebody", str(root / "housekeeping.bed"), "--ref-flat", str(root / "refFlat.txt"),
               "--star-index", str(root / "star")] if args.reuse_fixture else fixture(root))
    output = root / args.output_name
    subprocess.run([sys.executable, str(WORKFLOW_ROOT / "legacy_nextflow/run_all_resolveome_rna.py"),
                    "--run", "--executor", "local", "--threads", "2", "--memory-mb", "48000",
                    "--mode", args.mode,
                    "--outdir", str(output), "--qc-reads", "45000", *inputs], check=True)
    if not (output / "workflow_dashboard.html").is_file():
        raise RuntimeError("Synthetic integration did not produce the final dashboard")
    print(f"[OK] Synthetic {args.mode} branches and audited dashboard: {output}")


if __name__ == "__main__":
    main()
