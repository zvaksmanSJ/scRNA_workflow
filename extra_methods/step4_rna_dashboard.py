#!/usr/bin/env python3
"""Audit BaseJumper RNA outputs and link the two analysis depths in an HTML dashboard."""
import argparse
import csv
import html
import json
import subprocess
from pathlib import Path


def audit_branch(root, expected):
    root = Path(root)
    index = root / "workflow_outputs/resolveome/RNA_WORKFLOW/index"
    metrics = root / "workflow_outputs/resolveome/RNA_WORKFLOW/metrics/rnaqc_metrics/rnaqc_all_metrics.tsv"
    required = [index / "bam.csv", index / "per_biosample_status.csv", metrics,
                root / "workflow_outputs/resolveome/RNA_WORKFLOW/reports/multiqc_report.html"]
    for path in required:
        if not path.is_file() or not path.stat().st_size:
            raise FileNotFoundError(f"Missing final output: {path}")
    with (index / "per_biosample_status.csv").open() as stream:
        statuses = list(csv.DictReader(stream))
    names = [row["biosampleName"] for row in statuses]
    if len(names) != len(expected) or set(names) != set(expected):
        raise ValueError(f"RNA status cohort mismatch: {names}")
    with metrics.open() as stream:
        metric_rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(metric_rows) != len(expected) or {row["biosample"] for row in metric_rows} != set(expected):
        raise ValueError("RNA metric cohort mismatch")
    bams = list(root.glob("bam/**/*.bam"))
    if len(bams) != len(expected) or {path.stem for path in bams} != set(expected):
        raise ValueError(f"Final BAM cohort mismatch: expected {sorted(expected)}, "
                         f"found {sorted(path.stem for path in bams)} under {root / 'bam'}")
    with (index / "bam.csv").open() as stream:
        bam_rows = list(csv.DictReader(stream))
    if len(bam_rows) != len(expected) or {row["biosampleName"] for row in bam_rows} != set(expected):
        raise ValueError("BAM index cohort mismatch")
    published = {bam.stem: bam for bam in bams}
    for row in bam_rows:
        sid = row["biosampleName"]
        for key, actual in (("bam", published[sid]), ("bai", published[sid].with_suffix(".bam.bai"))):
            recorded = Path(row[key])
            if not recorded.is_absolute():
                recorded = root / recorded
            if recorded.resolve() != actual.resolve():
                raise ValueError(f"BAM index does not reference published {key} for {sid}: {row[key]}")
    for bam in bams:
        if not bam.with_suffix(".bam.bai").is_file():
            raise FileNotFoundError(f"Missing BAM index: {bam}")
        subprocess.run(["samtools", "quickcheck", "-v", str(bam)], check=True)
    import pyarrow.parquet as pq
    count_cells = set()
    qc_cells = set()
    for path in (root / "tables").rglob("*.parquet"):
        schema = pq.read_schema(path)
        if "gene_id" in schema.names and "count" in schema.names:
            table = pq.ParquetFile(path).read(columns=["biosample", "count"]).to_pandas()
            if table["count"].isna().any() or (table["count"] < 0).any():
                raise ValueError(f"Invalid expression counts: {path}")
            count_cells.update(table["biosample"])
        elif "qc_status" in schema.names:
            qc_cells.update(pq.ParquetFile(path).read(columns=["biosample"])
                            .column("biosample").to_pylist())
    if count_cells != set(expected) or qc_cells != set(expected):
        raise ValueError("Expression/QC Parquet cohort mismatch")
    for sid in expected:
        for directory, pattern in (("raw_fastqc", "*.zip"), ("trimmed_fastqc", "*.zip"),
                                   ("picard", "*.rna_metrics"),
                                   ("gene_body", "*.geneBodyCoverage.txt")):
            if not list((root / "extras" / directory / sid).glob(pattern)):
                raise FileNotFoundError(f"Missing supplemental QC: {directory}/{sid}")
    return {"samples": len(names), "statuses": statuses}


def dashboard(outdir, branches, expected):
    outdir = Path(outdir)
    protocol = json.loads((outdir / "workflow_state.json").read_text())["fingerprint"]
    qc_reads = protocol["qc_reads"]
    if not isinstance(qc_reads, int) or qc_reads < 1:
        raise ValueError("Invalid QC read target in workflow state")
    report = {}
    sections = []
    for branch in branches:
        root = outdir / branch
        report[branch] = audit_branch(root, expected)
        subprocess.run(["multiqc", str(root / "extras"), str(root / "tool_outputs"),
                        "--force", "--outdir", str(root / "reports"),
                        "--filename", "multiqc_supplemental.html"], check=True)
        links = [
            ("BaseJumper consensus dashboard",
             f"{branch}/workflow_outputs/resolveome/RNA_WORKFLOW/reports/multiqc_report.html"),
            ("Source figures",
             f"{branch}/workflow_outputs/resolveome/RNA_WORKFLOW/qc_plots/"),
            ("All RNA metrics",
             f"{branch}/workflow_outputs/resolveome/RNA_WORKFLOW/metrics/rnaqc_metrics/rnaqc_all_metrics.tsv"),
            ("FastQC, Picard and tool QC", f"{branch}/reports/multiqc_supplemental.html"),
            ("Per-cell gene-body coverage", f"{branch}/extras/gene_body/"),
            ("Expression Parquet tables", f"{branch}/tables/"),
        ]
        listing = "".join(f'<li><a href="{html.escape(path)}">{html.escape(label)}</a></li>'
                          for label, path in links)
        sections.append(f"<h2>{html.escape(branch)}</h2><p>{len(expected)} samples audited.</p>"
                        f"<ul>{listing}</ul>")
    document = (
        '<!doctype html><html><head><meta charset="utf-8"><title>ResolveOME RNA</title>'
        '<style>body{font-family:sans-serif;max-width:1100px;margin:40px auto;padding:20px}'
        'li{margin:8px 0}h2{color:#12284c}</style></head><body>'
        '<h1>ResolveOME RNA workflow</h1>'
        f'<p>Configured QC target: approximately {qc_reads:,} individual reads per library; '
        'full-depth uses all '
        'reads. Counts are paired-fragment HTSeq counts, not UMIs. QC PASS does not '
        'exclude cells from the full-depth matrix. Source score cutoffs: mappability &gt;0.7, '
        'exonic &gt;0.7, intergenic &lt;0.1, mitochondrial &lt;0.1, '
        'protein-coding genes &gt;500. Score 4–5 PASS, 3 Borderline, 0–2 FAIL.</p>'
        + "".join(sections) + "</body></html>"
    )
    (outdir / "workflow_dashboard.html").write_text(document)
    (outdir / "output_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--branches", default="qc_2M,full_depth")
    args = parser.parse_args()
    from python_lsf_pipeline.step1_prepare_samples import read_samples
    dashboard(args.outdir, args.branches.split(","), [r["sample_id"] for r in read_samples(args.cells)])


if __name__ == "__main__":
    main()
