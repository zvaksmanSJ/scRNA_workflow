#!/usr/bin/env python3
"""Verify source QC metrics, score components and Parquet against audited inputs."""
import json
import argparse
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
from common import WORKFLOW


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", choices=("qc", "full"), default="qc")
    args = parser.parse_args()
    root = HERE / ("stage7_qc_reports_v2" if args.depth == "qc" else "stage7_full_reports")
    name = "step5_output_audit.json" if args.depth == "qc" else "step5_full_output_audit.json"
    previous = json.loads((HERE / name).read_text())
    counts = {row["sample"]: row for row in previous["samples"]}
    if args.depth == "qc":
        prep_path = HERE / "step3_output_audit.json"
    else:
        prep_path = Path(json.loads((HERE / "stage4_full_alignment/preprocessing_audit_used.json")
                                   .read_text())["path"])
    prep = {row["sample"]: row for row in json.loads(prep_path.read_text())["samples"]}
    expected = set(counts)
    trace = pd.read_csv(root / "trace.tsv", sep="\t")
    processes = {name.split(" (")[0] for name in trace["name"]}
    if len(trace) != 4 or processes != {"PLOTTER_PCAHEATMAP", "CREATE_QC_REPORT",
                                      "RNA_QC_PLOTS", "MULTIQC"} or not (
            (trace["exit"] == 0) & trace.status.isin(["COMPLETED", "CACHED"])).all():
        raise ValueError("Report tasks incomplete")
    metrics = pd.read_csv(root / "consensus/rnaqc_all_metrics.tsv", sep="\t")
    scores = pd.read_csv(root / "consensus/RNA-QC_ConsensusScores.txt", sep="\t").set_index("SampleId")
    statuses = pd.read_csv(root / "consensus/per_biosample_status.csv").set_index("biosampleName")
    for names in (metrics.biosample, scores.index, statuses.index):
        if len(names) != len(expected) or set(names) != expected:
            raise ValueError("Report cohort mismatch")
    result = []
    for _, row in metrics.iterrows():
        sid = row.biosample
        prior = counts[sid]
        if int(row.htseq_assigned) != prior["assigned_counts"]:
            raise ValueError(f"Reported assigned counts differ from HTSeq: {sid}")
        for field, category in (("htseq_no_feature", "__no_feature"),
                                ("htseq_ambiguous", "__ambiguous"),
                                ("htseq_alignment_not_unique", "__alignment_not_unique")):
            if int(row[field]) != prior["HTSeq_special"][category]:
                raise ValueError(f"Reported HTSeq category differs: {sid}/{field}")
        if int(row.final_reads) != prep[sid]["sampled_reads"] or int(
                row.fastp_reads_after_filter) != prep[sid]["trimmed_reads"]:
            raise ValueError(f"Report input/trimmed read counts differ: {sid}")
        verdicts = {
            "Verdict_PropMappability": int(row.prop_mappability > 0.7),
            "Verdict_PropExonic": int(row.prop_exonic > 0.7),
            "Verdict_PropIntergenic": int(row.prop_intergenic < 0.1),
            "Verdict_ProportionCountsMitochondrialGenes": int(row.mt_pct < 0.1),
            "Verdict_ProteinCodingGenesDetected": int(row.genes_protein_coding > 500),
        }
        score = sum(verdicts.values())
        label = "PASS" if score >= 4 else "Borderline" if score == 3 else "FAIL"
        if int(row.qc_score) != score or row.qc_status != label or statuses.loc[sid, "qc_status"] != label:
            raise ValueError(f"Incorrect QC score/status: {sid}")
        if int(scores.loc[sid, "CompositeScore"]) != score or any(
                int(scores.loc[sid, field]) != value for field, value in verdicts.items()):
            raise ValueError(f"Consensus component mismatch: {sid}")
        paths = list((root / "consensus").rglob(f"biosample={sid}/output.parquet"))
        if len(paths) != 1:
            raise ValueError(f"Missing/duplicate QC Parquet: {sid}")
        table = pq.ParquetFile(paths[0]).read().to_pandas()
        if len(table) != 1 or table.iloc[0].biosample != sid:
            raise ValueError(f"QC Parquet cohort mismatch: {sid}")
        for field in ("htseq_assigned", "qc_score", "qc_status", "genes_protein_coding"):
            if table.iloc[0][field] != row[field]:
                raise ValueError(f"QC Parquet metric differs: {sid}/{field}")
        result.append({"sample": sid, "score": score, "status": label,
                       "assigned_counts": int(row.htseq_assigned)})
    for relative in ("multiqc/multiqc_report.html", "pca/PCA_mqc.png", "pca/Heatmap__mqc.png",
                     "pca/MT_mqc.png", "pca/GenesDetected_mqc.png",
                     "consensus/composition_rnaqc.pdf", "consensus/RNAQC_composition_mqc.jpg",
                     "consensus/RNA-QC_ConsensusScores_SummaryTable_mqc.txt",
                     "consensus/summary_verdict_group.txt", "consensus/combined_selected_metrics.txt"):
        path = root / relative
        if not path.is_file() or not path.stat().st_size:
            raise FileNotFoundError(path)
    if len(list((root / "consensus").glob("*_selected_metrics_mqc.txt"))) != 16:
        raise ValueError("Missing per-cell MultiQC metric files")
    report = {"status": "PASS", "samples": result,
              "labels": metrics.qc_status.value_counts().to_dict(),
              "scope": f"{args.depth}-depth staged reports"}
    name = "step7_output_audit.json" if args.depth == "qc" else "step7_full_output_audit.json"
    (HERE / name).write_text(json.dumps(report, indent=2) + "\n")
    print("STEP 7 PASS:", report["labels"], "All 16 assigned counts, scores and QC Parquets verified.")


if __name__ == "__main__":
    main()
