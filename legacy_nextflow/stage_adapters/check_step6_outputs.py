#!/usr/bin/env python3
"""Verify every matrix value against original HTSeq counts and count Parquet."""
import csv
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
from common import WORKFLOW
from step1_prepare_samples import read_samples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--depth", choices=("qc", "full"), default="qc")
    args = parser.parse_args()
    root = HERE / f"stage6_{args.depth}_matrix"
    expected = {r["sample_id"] for r in read_samples(HERE / "rna_cells.tsv")}
    trace = list(csv.DictReader((root / "trace.tsv").open(), delimiter="\t"))
    for process, size in (("CREATE_HTSEQ_SUMMARY", 16), ("MERGE_HTSEQ_SUMMARY", 1),
                          ("CREATE_HTSEQ_MATRIX", 1)):
        tasks = [r for r in trace if r["name"].startswith(process + " (")]
        if len(tasks) != size or any(r["exit"] != "0" or r["status"] not in
                                     ("COMPLETED", "CACHED") for r in tasks):
            raise ValueError(f"Incomplete matrix tasks: {process}")
    for sid in expected:
        for prefix in ("df_gene_star_htseq_", "df_sum_detected_gene_", "df_mtcounts_star_htseq_"):
            path = root / "summaries" / f"{prefix}{sid}.tsv"
            if not path.is_file() or not path.stat().st_size:
                raise FileNotFoundError(path)
    matrix = pd.read_csv(root / "matrix/matrix_gene_counts_starhtseq.txt", sep="\t")
    if matrix.columns[0] != "gene_id_gene_symbol" or set(matrix.columns[1:]) != expected:
        raise ValueError("Matrix cohort/schema mismatch")
    if matrix.iloc[:, 0].duplicated().any():
        raise ValueError("Duplicate matrix gene labels")
    values = matrix.iloc[:, 1:].to_numpy()
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError("Invalid expression matrix counts")
    gene_ids = matrix.iloc[:, 0].str.rsplit("_", n=1).str[0]
    base_ids = gene_ids.str.split(".").str[0]
    merged = pd.read_csv(root / "merged/df_gene_counts_starhtseq.tsv", sep="\t")
    if merged.duplicated(["File", "gene_id"]).any():
        raise ValueError("Duplicate merged sample/gene counts")
    totals = {}
    for sid in sorted(expected):
        raw = pd.read_csv(HERE / f"stage5_{args.depth}_counts/tool_outputs/htseq" / f"{sid}.htseq_counts.tsv",
                          sep="\t", names=["gene", "symbol", "count"])
        raw = raw[~raw.gene.str.startswith("__")].set_index("gene")["count"]
        expected_values = raw.reindex(base_ids).to_numpy()
        if not np.array_equal(matrix[sid].to_numpy(), expected_values):
            raise ValueError(f"Matrix values differ from original HTSeq counts: {sid}")
        if int(matrix[sid].sum()) != int(raw.sum()):
            raise ValueError(f"Matrix loses counts: {sid}")
        summary = merged[merged.File == f"{sid}.htseq_counts.tsv"]
        if int(summary.countHTSeq.sum()) != int(raw.sum()):
            raise ValueError(f"Gene summary loses counts: {sid}")
        totals[sid] = int(raw.sum())
    path = root / "matrix/gene_count_summary/workspace=resolveome/workflow_id=RNA_WORKFLOW/output.parquet"
    parquet = pq.ParquetFile(path).read().to_pandas()
    if set(parquet.biosample) != expected or len(parquet) != len(matrix) * len(expected):
        raise ValueError("Count Parquet cohort/row count mismatch")
    if parquet.duplicated(["gene_id", "biosample"]).any():
        raise ValueError("Duplicate gene/sample rows in Parquet")
    expected_long = matrix.assign(gene_id=gene_ids).drop(columns=["gene_id_gene_symbol"]).melt(
        id_vars=["gene_id"], var_name="biosample", value_name="count")
    compared = expected_long.merge(parquet[["gene_id", "biosample", "count"]],
                                  on=["gene_id", "biosample"], validate="one_to_one", suffixes=("_m", "_p"))
    if len(compared) != len(expected_long) or not (compared.count_m == compared.count_p).all():
        raise ValueError("Count Parquet differs from matrix")
    for name in ("HouseKeepingGenes_Counts_mqc.tsv", "HouseKeepingGenes_CV_mqc.tsv",
                 "HKGenes_Expression__mqc.png"):
        product = root / "matrix" / name
        if not product.is_file() or not product.stat().st_size:
            raise FileNotFoundError(product)
    report = {"status": "PASS", "samples": sorted(expected), "matrix_genes": len(matrix),
              "Parquet_rows": len(parquet), "assigned_counts": totals,
              "all_matrix_values_match_original_HTSeq": True}
    name = "step6_output_audit.json" if args.depth == "qc" else "step6_full_output_audit.json"
    (HERE / name).write_text(json.dumps(report, indent=2) + "\n")
    print(f"STEP 6 PASS: {len(matrix):,} genes x {len(expected)} cells; "
          f"{len(parquet):,} Parquet rows; every count matches original HTSeq.")


if __name__ == "__main__":
    main()
