#!/usr/bin/env python3
"""Account for every per-cell count and RNA-QC product before matrix generation."""
import csv
import argparse
import json
import math
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
from common import WORKFLOW
from step1_prepare_samples import read_samples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", choices=("qc", "full"), default="qc")
    args = parser.parse_args()
    root = HERE / f"stage5_{args.depth}_counts"
    expected = {r["sample_id"] for r in read_samples(HERE / "rna_cells.tsv")}
    tasks = list(csv.DictReader((root / "trace.tsv").open(), delimiter="\t"))
    for process in ("HTSEQ_COUNTS", "QUALIMAP_BAMRNA", "PICARD_RNA", "GENE_BODY_COVERAGE"):
        selected = [r for r in tasks if r["name"].startswith(process + " (")]
        names = [r["name"].split(" (", 1)[1].removesuffix(")") for r in selected]
        if len(names) != len(expected) or set(names) != expected or any(
                r["exit"] != "0" or r["status"] not in ("COMPLETED", "CACHED") for r in selected):
            raise ValueError(f"Incomplete or failed tasks: {process}")
    results = []
    gene_set = None
    for sid in sorted(expected):
        counts, special = {}, {}
        with (root / "tool_outputs/htseq" / f"{sid}.htseq_counts.tsv").open() as stream:
            for line in stream:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 3:
                    raise ValueError(f"Incorrect HTSeq schema: {sid}")
                gene, symbol, count_text = fields
                count = int(count_text)
                target = special if gene.startswith("__") else counts
                if gene in target or count < 0:
                    raise ValueError(f"Duplicate/negative HTSeq count: {sid}/{gene}")
                target[gene] = count
        if set(special) != {"__no_feature", "__ambiguous", "__too_low_aQual",
                            "__not_aligned", "__alignment_not_unique"}:
            raise ValueError(f"HTSeq special categories missing: {sid}")
        if len(counts) != 60664 or (gene_set is not None and set(counts) != gene_set):
            raise ValueError(f"HTSeq gene cohort mismatch: {sid}")
        gene_set = set(counts)
        if not sum(counts.values()):
            raise ValueError(f"No assigned counts: {sid}")
        directory = root / "tool_outputs/qualimap" / f"qualimap_outdir_{sid}"
        for name in ("rnaseq_qc_results.txt", f"{sid}.report.pdf", f"{sid}.oc"):
            if not (directory / name).is_file() or not (directory / name).stat().st_size:
                raise FileNotFoundError(directory / name)
        text = (directory / "rnaseq_qc_results.txt").read_text()
        region = {}
        for name in ("exonic", "intronic", "intergenic"):
            match = re.search(rf"^\s*{name} =\s+([\d,]+) \(([\d.]+)%\)", text, re.M)
            if not match:
                raise ValueError(f"Missing Qualimap regional metric: {sid}/{name}")
            region[name] = float(match[2])
        if abs(sum(region.values()) - 100) > 0.05:
            raise ValueError(f"Qualimap regional percentages do not sum to 100: {sid}")
        metrics_path = root / "extras/picard" / sid / f"{sid}.rna_metrics"
        lines = metrics_path.read_text().splitlines()
        start = next(i for i, line in enumerate(lines) if line.startswith("PF_BASES\t"))
        picard = dict(zip(lines[start].split("\t"), lines[start + 1].split("\t")))
        fractions = [float(picard[name]) for name in
                     ("PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES", "PCT_INTERGENIC_BASES")]
        if any(not math.isfinite(x) or not 0 <= x <= 1 for x in fractions) or abs(sum(fractions) - 1) > 0.00001:
            raise ValueError(f"Invalid Picard regional fractions: {sid}")
        if int(picard["PF_ALIGNED_BASES"]) <= 0:
            raise ValueError(f"No Picard aligned bases: {sid}")
        flagstat = root / "extras/picard" / sid / f"{sid}.flagstat.txt"
        if not flagstat.is_file() or not flagstat.stat().st_size:
            raise FileNotFoundError(flagstat)
        body = root / "extras/gene_body" / sid
        with (body / "gene_body.geneBodyCoverage.txt").open() as stream:
            coverage = list(csv.reader(stream, delimiter="\t"))
        if len(coverage) != 2 or coverage[0][1:] != [str(n) for n in range(1, 101)]:
            raise ValueError(f"Incorrect gene-body percentile schema: {sid}")
        values = list(map(float, coverage[1][1:]))
        if coverage[1][0] != sid or len(values) != 100 or any(
                not math.isfinite(v) or v < 0 for v in values) or sum(values) <= 0:
            raise ValueError(f"Invalid gene-body coverage: {sid}")
        if not (body / "gene_body.geneBodyCoverage.curves.pdf").is_file():
            raise FileNotFoundError(body / "gene_body.geneBodyCoverage.curves.pdf")
        result = {"sample": sid, "genes": len(counts),
                  "detected_genes": sum(n > 0 for n in counts.values()),
                  "assigned_counts": sum(counts.values()), "HTSeq_special": special,
                  "Qualimap_region_pct": region, "Picard_mRNA_fraction": float(picard["PCT_MRNA_BASES"])}
        results.append(result)
        print(f'{sid}: detected={result["detected_genes"]:,} assigned={result["assigned_counts"]:,} '
              f'exonic={region["exonic"]:.2f}% intergenic={region["intergenic"]:.2f}%', flush=True)
    name = "step5_output_audit.json" if args.depth == "qc" else "step5_full_output_audit.json"
    (HERE / name).write_text(
        json.dumps({"status": "PASS", "samples": results}, indent=2) + "\n")
    print("STEP 5 PASS: all 64 tasks and required per-cell count/QC outputs accounted for. "
          "Picard ribosomal metrics are unavailable without rRNA intervals.")


if __name__ == "__main__":
    main()
