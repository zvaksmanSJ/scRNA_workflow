#!/usr/bin/env python3
"""Inspect matching RNA references before any real-data alignment."""
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
from common import WORKFLOW
from step2_prepare_references import DEFAULTS, check_references


def main():
    check_references(DEFAULTS)
    fai = {r.split("\t")[0]: int(r.split("\t")[1])
           for r in Path(DEFAULTS["fasta"] + ".fai").read_text().splitlines()}
    genes, transcripts, contigs = set(), set(), set()
    exon_genes = set()
    exons = 0
    with Path(DEFAULTS["gtf"]).open() as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                raise ValueError("Malformed GTF")
            chrom, start, end = fields[0], int(fields[3]), int(fields[4])
            if chrom not in fai or not 1 <= start <= end <= fai[chrom]:
                raise ValueError(f"GTF interval incompatible with FASTA: {fields[:5]}")
            contigs.add(chrom)
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            if "gene_id" in attrs:
                genes.add(attrs["gene_id"])
            if fields[2] == "exon":
                exons += 1
                transcripts.add(attrs["transcript_id"])
                exon_genes.add(attrs["gene_id"])
    flat_transcripts = set()
    with Path(DEFAULTS["ref_flat"]).open() as stream:
        for line in stream:
            fields = line.rstrip().split("\t")
            if len(fields) != 11:
                raise ValueError("Malformed refFlat")
            start, end, cds_start, cds_end, count = map(int, fields[4:9])
            starts = list(map(int, fields[9].rstrip(",").split(",")))
            ends = list(map(int, fields[10].rstrip(",").split(",")))
            if not (fields[2] in fai and 0 <= start < end <= fai[fields[2]]
                    and start <= cds_start <= cds_end <= end
                    and len(starts) == len(ends) == count
                    and all(start <= a < b <= end for a, b in zip(starts, ends))):
                raise ValueError(f"Invalid refFlat coordinates: {fields[:4]}")
            flat_transcripts.add(fields[1])
    if flat_transcripts != transcripts:
        raise ValueError("refFlat transcript set differs from GTF exons")
    mapping_genes = set()
    mapping_rows = 0
    with Path(DEFAULTS["tx2gene"]).open() as stream:
        for line in stream:
            fields = line.rstrip().split("\t")
            if len(fields) != 5:
                raise ValueError("Malformed five-column tx2gene")
            mapping_genes.add(fields[2].split(".")[0])
            mapping_rows += 1
    missing = exon_genes - mapping_genes
    if missing:
        raise ValueError(f"{len(missing)} countable GTF genes absent from normalized tx2gene: "
                         f"{sorted(missing)[:10]}")
    bed_rows = 0
    with Path(DEFAULTS["genebody"]).open() as stream:
        for line in stream:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if not (len(fields) >= 12 and fields[0] in fai
                    and 0 <= int(fields[1]) < int(fields[2]) <= fai[fields[0]]):
                raise ValueError(f"Invalid housekeeping BED interval: {fields[:3]}")
            bed_rows += 1
    report = {
        "status": "PASS", "FASTA_STAR_contigs": len(fai), "GTF_contigs": len(contigs),
        "GTF_genes": len(genes), "countable_GTF_genes": len(exon_genes), "GTF_exons": exons,
        "refFlat_transcripts": len(flat_transcripts), "tx2gene_rows": mapping_rows,
        "gene_ID_normalization": "Strip version suffix as in source htseq_to_gene_df.R",
        "housekeeping_BED_rows": bed_rows, "references": DEFAULTS,
    }
    (HERE / "step2_reference_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "references"}, indent=2))


if __name__ == "__main__":
    main()
