#!/usr/bin/env python3
"""Run isolated RNA-only QC, exon, splicing, and fusion research add-ons.

The command first gates on the full-depth report and Salmon audits, then emits
an execution plan. With --execute it runs RSeQC, DEXSeq, SUPPA2, and a
fusion-specific STAR/Arriba branch into a new results subdirectory. It never
changes canonical BAM/count outputs, performs no differential test, and uses no
matched DNA. The branch is optional and is not yet a Sprocket/WDL scatter stage.
"""
from __future__ import annotations

import argparse
import csv
from io import StringIO
import json
import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SAMPLE_PATTERN = re.compile(r"[A-Za-z0-9_.-]+")
OUTPUT_RELATIVE = Path("final_results/10_rna_research_addons")
BASEJ_RULES = (
    ("prop_mappability", ">", 0.70, "mappability must be >70%"),
    ("prop_exonic", ">", 0.70, "exonic fraction must be >70%"),
    ("prop_intergenic", "<", 0.10, "intergenic fraction must be <10%"),
    ("mt_pct", "<", 0.10, "mitochondrial assigned-count fraction must be <10%"),
    ("genes_protein_coding", ">", 500, "detected protein-coding genes must be >500"),
)


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"sample_id", "group"}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Sample sheet lacks required columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Sample sheet is empty: {path}")
    sample_ids = [row["sample_id"] for row in rows]
    if any(not SAMPLE_PATTERN.fullmatch(sample_id) for sample_id in sample_ids):
        raise ValueError("Sample IDs contain unsupported characters")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Sample IDs must be unique")
    if any(not row.get("group", "").strip() for row in rows):
        raise ValueError("Each sample must have a non-empty group label")
    return rows


def require_pass(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} is missing: {path}")
    with path.open() as stream:
        report = json.load(stream)
    if report.get("status") != "PASS":
        raise ValueError(f"{label} is not PASS: {report.get('status')}")


def basejumper_reasons(row: pd.Series) -> list[str]:
    reasons = []
    for column, operator, threshold, explanation in BASEJ_RULES:
        value = pd.to_numeric(pd.Series([row.get(column)]), errors="coerce").iloc[0]
        if pd.isna(value):
            reasons.append(f"{explanation} (metric missing)")
        elif (operator == ">" and value <= threshold) or (
                operator == "<" and value >= threshold):
            rendered = f"{value:.3g}" if column != "genes_protein_coding" else f"{int(value)}"
            reasons.append(f"{explanation} (observed {rendered})")
    return reasons


def source_qc_status(score: object) -> str:
    value = pd.to_numeric(pd.Series([score]), errors="coerce").iloc[0]
    if pd.isna(value):
        return "FAIL"
    if value >= 4:
        return "PASS"
    if value == 3:
        return "Borderline"
    return "FAIL"


def robust_low_outlier_scores(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    center = float(numeric.median())
    mad = float((numeric - center).abs().median())
    if mad == 0:
        return numeric.map(lambda value: 0.0 if value == center else (
            float("-inf") if value < center else float("inf")))
    return 0.6745 * (numeric - center) / mad


def improved_qc_status(
    base_status: str, base_reasons: Sequence[str], tin_z: float, sample_count: int,
) -> tuple[str, str]:
    reasons = list(base_reasons)
    if base_status == "FAIL":
        return "FAIL", "; ".join(["BaseJumper QC is FAIL", *reasons])
    if math.isnan(tin_z):
        return "FAIL", "Required transcript-integrity metric is missing or invalid"
    if sample_count < 8:
        reasons.append(
            f"Only {sample_count} samples; robust cohort TIN thresholds require at least 8")
        return "Borderline", "; ".join(reasons)
    if tin_z <= -3:
        reasons.append(f"TIN is a severe low cohort outlier (robust z={tin_z:.2f}, threshold <= -3)")
        return "FAIL", "; ".join(reasons)
    if tin_z <= -2:
        reasons.append(f"TIN is a low cohort outlier (robust z={tin_z:.2f}, threshold <= -2)")
        return "Borderline", "; ".join(reasons)
    if base_status == "Borderline":
        reasons.insert(0, "BaseJumper QC is Borderline")
        return "Borderline", "; ".join(reasons)
    return "PASS", "No BaseJumper failure or enhanced TIN outlier"


def load_tin_value(path: Path) -> float:
    if not path.is_file():
        raise FileNotFoundError(f"RSeQC TIN output missing: {path}")
    lines = [line for line in path.read_text().splitlines() if line.strip() and
             not line.startswith("##")]
    if len(lines) < 2:
        raise ValueError(f"RSeQC TIN output has no records: {path}")
    from io import StringIO

    frame = pd.read_csv(StringIO("\n".join(lines)), sep="\t")
    columns = [column for column in frame.columns
               if str(column).lstrip("#").casefold() in {"tin", "tin_score"}]
    if not columns:
        raise ValueError(f"RSeQC TIN output has no TIN column: {path}")
    values = pd.to_numeric(frame[columns[0]], errors="coerce").dropna()
    if values.empty or not np.isfinite(values.to_numpy()).all():
        raise ValueError(f"RSeQC TIN output has no finite scores: {path}")
    return float(values.median())


def run_qc_summaries(
    results_dir: Path, output_dir: Path, samples: Sequence[dict[str, str]],
) -> pd.DataFrame:
    basej_path = results_dir / "stage7_full_reports/consensus/rnaqc_all_metrics.tsv"
    basej = pd.read_csv(basej_path, sep="\t")
    if "biosample" not in basej or "qc_score" not in basej:
        raise ValueError(f"BaseJumper metrics lack required identity/score columns: {basej_path}")
    if basej.biosample.duplicated().any() or set(basej.biosample) != {
            row["sample_id"] for row in samples}:
        raise ValueError("BaseJumper metric sample membership differs from the sample sheet")
    basej = basej.set_index("biosample")
    rows = []
    for sample in samples:
        sample_id = sample["sample_id"]
        if sample_id not in basej.index:
            raise ValueError(f"BaseJumper metrics missing sample: {sample_id}")
        tin_path = output_dir / f"qc/tin/{sample_id}/{sample_id}.bam.TIN.xls"
        tin = load_tin_value(tin_path)
        source_row = basej.loc[sample_id]
        reasons = basejumper_reasons(source_row)
        status = source_qc_status(source_row.qc_score)
        if str(source_row.get("qc_status", status)) != status:
            raise ValueError(
                f"BaseJumper status/score disagree for {sample_id}: "
                f"{source_row.get('qc_status')} vs score {source_row.qc_score}")
        rows.append({
            "sample_id": sample_id,
            "group": sample["group"],
            "basejumper_status": status,
            "basejumper_score_0_to_5": int(source_row.qc_score),
            "basejumper_reasons": (
                "; ".join(reasons) if reasons else "All five BaseJumper criteria passed"),
            "median_transcript_integrity_number": tin,
        })
    metrics = pd.DataFrame(rows)
    metrics["tin_robust_z"] = robust_low_outlier_scores(
        metrics.median_transcript_integrity_number)
    if not np.isfinite(metrics.median_transcript_integrity_number.to_numpy()).all():
        raise ValueError("Per-sample median TIN values must all be finite")
    improved = metrics.apply(
        lambda row: improved_qc_status(
            row.basejumper_status,
            [] if row.basejumper_status == "PASS" else [row.basejumper_reasons],
            float(row.tin_robust_z), len(metrics),
        ),
        axis=1,
    )
    metrics["improved_status"] = [row[0] for row in improved]
    metrics["improved_reasons"] = [row[1] for row in improved]
    qc_dir = output_dir / "qc_enhanced"
    qc_dir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(qc_dir / "per_cell_qc_status.tsv", sep="\t", index=False)
    median_tin = float(metrics.median_transcript_integrity_number.median())
    mad_tin = float((metrics.median_transcript_integrity_number - median_tin).abs().median())
    (qc_dir / "thresholds.json").write_text(json.dumps({
        "basejumper": {
            "source": "BaseJumper 0-5 composite score",
            "status": "PASS for score >=4; Borderline for score ==3; FAIL for score <=2",
            "criteria": [
                {"metric": "prop_mappability", "operator": ">", "threshold": 0.70},
                {"metric": "prop_exonic", "operator": ">", "threshold": 0.70},
                {"metric": "prop_intergenic", "operator": "<", "threshold": 0.10},
                {"metric": "mt_pct", "operator": "<", "threshold": 0.10},
                {"metric": "genes_protein_coding", "operator": ">", "threshold": 500},
            ],
        },
        "enhanced": {
            "metric": "median RSeQC TIN per sample",
            "method": "modified robust z = 0.6745 * (sample median TIN - cohort median TIN) / MAD",
            "borderline": "robust z <= -2",
            "fail": "robust z <= -3",
            "minimum_samples_for_cohort_rule": 8,
            "cohort_median_tin": median_tin,
            "cohort_mad_tin": mad_tin,
            "automatic_exclusion": False,
            "interpretation": "cohort-relative research QC; no universal absolute TIN cutoff is asserted",
        },
    }, indent=2) + "\n")
    plot_tin_qc(metrics, qc_dir / "tin_by_sample.png")
    return metrics


def plot_tin_qc(metrics: pd.DataFrame, output: Path) -> None:
    ordered = metrics.sort_values("median_transcript_integrity_number")
    colors = {
        "PASS": "#197C83",
        "Borderline": "#D99A22",
        "FAIL": "#A23A3A",
    }
    fig, ax = plt.subplots(figsize=(max(10, 0.55 * len(ordered)), 5.4))
    ax.bar(
        ordered.sample_id,
        ordered.median_transcript_integrity_number,
        color=[colors[status] for status in ordered.improved_status],
        edgecolor="white",
        linewidth=0.7,
    )
    median = float(ordered.median_transcript_integrity_number.median())
    mad = float((ordered.median_transcript_integrity_number - median).abs().median())
    ax.axhline(median, color="#315A86", linestyle="-", label=f"Cohort median ({median:.1f})")
    if mad > 0:
        ax.axhline(median - 2 * mad / 0.6745, color="#D99A22", linestyle="--",
                   label="Borderline outlier threshold (robust z = -2)")
        ax.axhline(median - 3 * mad / 0.6745, color="#A23A3A", linestyle=":",
                   label="Fail outlier threshold (robust z = -3)")
    ax.set_ylabel("Median transcript integrity number (TIN)")
    ax.set_xlabel("Sample")
    ax.set_title("Transcript integrity by sample")
    ax.tick_params(axis="x", rotation=65)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def collect_exon_counts(
    output_dir: Path, samples: Sequence[dict[str, str]],
) -> pd.DataFrame:
    frames = []
    for row in samples:
        path = output_dir / f"exon_counts/{row['sample_id']}.counts.txt"
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"DEXSeq exon counts missing: {path}")
        with path.open() as stream:
            first_line = stream.readline().rstrip("\n").split("\t")
        if first_line[:3] == ["groupID", "featureID", "count"]:
            frame = pd.read_csv(path, sep="\t")
            frame = frame.rename(columns={
                "groupID": "gene_id", "featureID": "feature_id", "count": "count"})
        else:
            frame = pd.read_csv(
                path, sep="\t", header=None,
                names=["gene_id", "feature_id", "count"])
        required = {"gene_id", "feature_id", "count"}
        if not required.issubset(frame.columns):
            raise ValueError(f"Unexpected DEXSeq exon count schema: {path}")
        frame["sample_id"] = row["sample_id"]
        frame["count"] = pd.to_numeric(frame["count"], errors="raise")
        if frame["count"].lt(0).any() or not np.isfinite(frame["count"]).all():
            raise ValueError(f"Invalid DEXSeq exon count values: {path}")
        frames.append(frame[["gene_id", "feature_id", "count", "sample_id"]])
    counts = pd.concat(frames, ignore_index=True)
    if counts.duplicated(["gene_id", "feature_id", "sample_id"]).any():
        raise ValueError("Duplicate exon-bin/sample rows in DEXSeq counts")
    matrix = counts.pivot(
        index=["gene_id", "feature_id"], columns="sample_id", values="count")
    if matrix.isna().any().any():
        raise ValueError("DEXSeq exon-bin membership differs across samples")
    table = matrix.reset_index()
    table.columns.name = None
    table.to_csv(output_dir / "exon_counts/exon_bin_count_matrix.tsv",
                 sep="\t", index=False)
    fractions = matrix.div(matrix.groupby(level=0).transform("sum").replace(0, np.nan))
    fractions = fractions.fillna(0).reset_index()
    fractions.to_csv(output_dir / "exon_counts/exon_usage_fraction.tsv",
                     sep="\t", index=False)
    plot_exon_bins(matrix, output_dir / "exon_counts/top_variable_exon_bins.png")
    return table


def collect_duplication_summary(
    output_dir: Path, samples: Sequence[dict[str, str]],
) -> pd.DataFrame:
    rows = []
    for sample in samples:
        sample_id = sample["sample_id"]
        path = output_dir / f"qc/duplication/{sample_id}.dup.pos.DupRate.xls"
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"RSeQC positional duplication table missing: {path}")
        lines = [line for line in path.read_text().splitlines()
                 if line.strip() and not line.startswith("#")]
        if not lines:
            raise ValueError(f"RSeQC positional duplication table is empty: {path}")
        table = pd.read_csv(StringIO("\n".join(lines)), sep="\t", header=None)
        if table.shape[1] < 2:
            raise ValueError(f"Unexpected RSeQC positional duplication schema: {path}")
        occurrence = pd.to_numeric(table.iloc[:, 0], errors="coerce")
        unique_positions = pd.to_numeric(table.iloc[:, 1], errors="coerce")
        valid = occurrence.notna() & unique_positions.notna()
        occurrence = occurrence[valid].astype(float)
        unique_positions = unique_positions[valid].astype(float)
        if occurrence.empty or occurrence.le(0).any() or unique_positions.lt(0).any():
            raise ValueError(f"Invalid RSeQC positional duplication values: {path}")
        if not np.equal(occurrence, np.floor(occurrence)).all():
            raise ValueError(f"RSeQC occurrence values must be integers: {path}")
        distinct_positions = float(unique_positions.sum())
        mapped_reads = float((occurrence * unique_positions).sum())
        if mapped_reads <= 0 or distinct_positions > mapped_reads:
            raise ValueError(f"Invalid duplication totals in RSeQC output: {path}")
        rows.append({
            "sample_id": sample_id,
            "distinct_mapped_positions": int(distinct_positions),
            "reads_at_mapped_positions": int(mapped_reads),
            "estimated_positional_duplication_percent":
                100 * (1 - distinct_positions / mapped_reads),
            "method": "RSeQC positional occurrence distribution; not UMI molecule counting",
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(output_dir / "qc/duplication/positional_duplication_summary.tsv",
                   sep="\t", index=False)
    plot_duplication(summary, output_dir / "qc/duplication/positional_duplication.png")
    return summary


def plot_duplication(summary: pd.DataFrame, output: Path) -> None:
    ordered = summary.sort_values("estimated_positional_duplication_percent")
    fig, ax = plt.subplots(figsize=(max(10, 0.55 * len(ordered)), 5.4))
    ax.bar(ordered.sample_id, ordered.estimated_positional_duplication_percent,
           color="#315A86", edgecolor="white", linewidth=0.7)
    ax.set_ylabel("Estimated positional duplication (%)")
    ax.set_xlabel("Sample")
    ax.set_title("RSeQC positional read-duplication estimate")
    ax.tick_params(axis="x", rotation=65)
    ax.grid(axis="y", color="#E8ECEF", linewidth=0.6)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_exon_bins(matrix: pd.DataFrame, output: Path) -> None:
    if matrix.empty:
        raise ValueError("DEXSeq produced no exon bins")
    values = np.log2(matrix.astype(float) + 1)
    variances = values.var(axis=1).sort_values(ascending=False)
    selected = values.loc[variances.head(min(40, len(variances))).index]
    fig, ax = plt.subplots(figsize=(max(9, 0.55 * selected.shape[1]), 7))
    image = ax.imshow(selected.to_numpy(), aspect="auto", interpolation="nearest",
                      cmap="vlag" if "vlag" in plt.colormaps() else "coolwarm")
    ax.set_xticks(np.arange(len(selected.columns)), selected.columns, rotation=60, ha="right")
    ax.set_yticks(np.arange(len(selected.index)), [
        f"{gene}:{feature}" for gene, feature in selected.index], fontsize=7)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Gene: exon-bin feature")
    ax.set_title("Most variable exon-bin counts (log2[count + 1])")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.7)
    colorbar.set_label("log2 exon-bin count + 1")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def collect_psi(output_dir: Path, samples: Sequence[dict[str, str]]) -> pd.DataFrame:
    paths = sorted((output_dir / "splicing").glob("*.psi"))
    if not paths:
        raise FileNotFoundError("SUPPA2 produced no event PSI outputs")
    frames = []
    sample_ids = [row["sample_id"] for row in samples]
    for path in paths:
        frame = pd.read_csv(path, sep="\t")
        if frame.empty:
            continue
        event_column = frame.columns[0]
        missing = set(sample_ids).difference(frame.columns)
        if missing:
            raise ValueError(
                f"SUPPA2 event PSI file lacks sample columns {sorted(missing)}: {path}")
        values = frame[sample_ids].apply(pd.to_numeric, errors="coerce")
        if values.notna().sum().sum() == 0:
            raise ValueError(f"SUPPA2 event PSI file has no numeric values: {path}")
        table = pd.DataFrame({
            "event_id": frame[event_column].astype(str),
            "event_type": (
                re.search(r"_(SE|SS|MX|RI|FL)(?:_|$)", path.stem).group(1)
                if re.search(r"_(SE|SS|MX|RI|FL)(?:_|$)", path.stem)
                else "unspecified"),
        })
        table = pd.concat([table, values], axis=1)
        frames.append(table)
    if not frames:
        raise ValueError("SUPPA2 event PSI outputs contain no events")
    events = pd.concat(frames, ignore_index=True)
    if events.event_id.duplicated().any():
        raise ValueError("Duplicate SUPPA2 event IDs across event types")
    events.to_csv(output_dir / "splicing/event_psi_by_sample.tsv",
                  sep="\t", index=False, na_rep="NA")
    values = events.set_index("event_id")[sample_ids]
    variances = values.var(axis=1, skipna=True).sort_values(ascending=False)
    top = values.loc[variances.head(min(40, len(variances))).index]
    top.reset_index().to_csv(
        output_dir / "splicing/variable_event_psi.tsv",
        sep="\t", index=False, na_rep="NA")
    plot_psi_heatmap(top, output_dir / "splicing/variable_event_psi.png")
    return events


def plot_psi_heatmap(values: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(max(9, 0.55 * values.shape[1]), 7))
    image = ax.imshow(values.to_numpy(dtype=float), aspect="auto",
                      interpolation="nearest", vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(np.arange(len(values.columns)), values.columns, rotation=60, ha="right")
    ax.set_yticks(np.arange(len(values.index)), values.index, fontsize=7)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Splicing event ID")
    ax.set_title("Variable event PSI across samples")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.7)
    colorbar.set_label("Percent-spliced-in estimate (PSI)")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def parse_breakpoint(value: object) -> tuple[str, int] | None:
    match = re.fullmatch(r"(?:chr)?([^:]+):([\d,]+)", str(value))
    if not match:
        return None
    try:
        return match.group(1), int(match.group(2).replace(",", ""))
    except ValueError:
        return None


def summarize_fusions(
    output_dir: Path, samples: Sequence[dict[str, str]], fasta: Path,
) -> pd.DataFrame:
    rows = []
    for sample in samples:
        path = output_dir / f"fusions/{sample['sample_id']}.fusions.tsv"
        if not path.is_file():
            raise FileNotFoundError(f"Arriba fusion output missing: {path}")
        lines = [line for line in path.read_text().splitlines() if line.strip()]
        header_index = next((
            index for index, line in enumerate(lines)
            if {"gene1", "gene2"}.issubset(
                {field.lstrip("#") for field in line.split("\t")})
        ), None)
        if header_index is None:
            raise ValueError(f"Arriba fusion output has no recognized header: {path}")
        header = lines[header_index].lstrip("#")
        frame = pd.read_csv(
            StringIO("\n".join([header, *lines[header_index + 1:]])), sep="\t")
        if frame.empty:
            continue
        required = {
            "gene1", "gene2", "breakpoint1", "breakpoint2", "confidence",
            "filters", "split_reads1", "split_reads2", "discordant_mates",
        }
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"Arriba output lacks evidence columns {sorted(missing)}: {path}")
        for record in frame.to_dict("records"):
            confidence = str(record["confidence"]).lower()
            filters = str(record["filters"])
            split1 = int(record["split_reads1"])
            split2 = int(record["split_reads2"])
            discordant = int(record["discordant_mates"])
            evidence_reads = split1 + split2 + discordant
            if confidence == "high" and filters in {"", ".", "no_filter"} and evidence_reads > 0:
                tier = "HIGH"
            elif confidence == "medium" and filters in {"", ".", "no_filter"} and evidence_reads > 0:
                tier = "REVIEW"
            else:
                tier = "LOW_OR_FILTERED"
            rows.append({
                "sample_id": sample["sample_id"],
                "gene1": str(record["gene1"]),
                "gene2": str(record["gene2"]),
                "breakpoint1": str(record["breakpoint1"]),
                "breakpoint2": str(record["breakpoint2"]),
                "confidence": confidence,
                "filters": filters,
                "split_reads1": split1,
                "split_reads2": split2,
                "discordant_mates": discordant,
                "total_supporting_reads": evidence_reads,
                "evidence_tier": tier,
                "threshold_basis": (
                    "Arriba default filters/confidence; passed file and >=1 reported supporting read"
                ),
            })
    columns = [
        "sample_id", "gene1", "gene2", "breakpoint1", "breakpoint2",
        "confidence", "filters", "split_reads1", "split_reads2",
        "discordant_mates", "total_supporting_reads", "evidence_tier",
        "threshold_basis",
    ]
    result = pd.DataFrame(rows, columns=columns)
    fusion_dir = output_dir / "fusions"
    result.to_csv(fusion_dir / "fusion_evidence_summary.tsv",
                  sep="\t", index=False)
    plot_fusion_map(result, fasta, fusion_dir / "fusion_chromosome_map.png")
    return result


def plot_fusion_map(fusions: pd.DataFrame, fasta: Path, output: Path) -> None:
    fai = Path(str(fasta) + ".fai")
    if not fai.is_file():
        raise FileNotFoundError(f"FASTA index required for fusion map: {fai}")
    lengths = {}
    for line in fai.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) >= 2:
            chrom = fields[0].removeprefix("chr")
            if re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)", chrom):
                lengths[chrom] = int(fields[1])
    chromosomes = [str(i) for i in range(1, 23)] + ["X", "Y", "M", "MT"]
    chromosomes = [chrom for chrom in chromosomes if chrom in lengths]
    if not chromosomes:
        raise ValueError("FASTA index has no standard human chromosomes")
    if fusions.empty:
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.axis("off")
        ax.text(0.5, 0.58, "No fusion candidates reported by Arriba",
                ha="center", va="center", fontsize=16, color="#24313B")
        ax.text(0.5, 0.4, "The passed-candidate output was empty; discarded calls remain available.",
                ha="center", va="center", fontsize=10, color="#64717B")
        fig.tight_layout()
        fig.savefig(output, dpi=180)
        plt.close(fig)
        return

    chrom_y = {chrom: len(chromosomes) - i - 1 for i, chrom in enumerate(chromosomes)}
    fig, ax = plt.subplots(figsize=(12, max(6, len(chromosomes) * 0.34)))
    for chrom in chromosomes:
        length_mb = lengths[chrom] / 1_000_000
        y = chrom_y[chrom]
        ax.hlines(y, 0, length_mb, color="#B9C3CB", linewidth=5, zorder=1)
        ax.text(-4, y, f"chr{chrom}", ha="right", va="center", fontsize=8)
    colors = {"HIGH": "#6B1E3D", "REVIEW": "#D99A22", "LOW_OR_FILTERED": "#9AA5AE"}
    plotted = 0
    for fusion in fusions.to_dict("records"):
        first = parse_breakpoint(fusion["breakpoint1"])
        second = parse_breakpoint(fusion["breakpoint2"])
        if first is None or second is None or first[0] not in chrom_y or second[0] not in chrom_y:
            continue
        x1, x2 = first[1] / 1_000_000, second[1] / 1_000_000
        y1, y2 = chrom_y[first[0]], chrom_y[second[0]]
        color = colors.get(fusion["evidence_tier"], colors["LOW_OR_FILTERED"])
        ax.annotate(
            "", xy=(x2, y2), xytext=(x1, y1),
            arrowprops={"arrowstyle": "-", "connectionstyle": "arc3,rad=0.2",
                        "color": color, "alpha": 0.75, "linewidth": 1.6},
            zorder=2,
        )
        ax.scatter([x1, x2], [y1, y2], color=color, s=19, zorder=3)
        ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.18,
                f"{fusion['gene1']}–{fusion['gene2']}",
                fontsize=7, ha="center", color=color)
        plotted += 1
    handles = [
        plt.Line2D([0], [0], color=colors[tier], lw=2, label=label)
        for tier, label in (("HIGH", "Arriba high confidence"),
                            ("REVIEW", "Arriba medium confidence; review"),
                            ("LOW_OR_FILTERED", "Low/filtered call"))
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right", fontsize=8)
    max_length = max(lengths[chrom] for chrom in chromosomes) / 1_000_000
    ax.set_xlim(-8, max_length * 1.04)
    ax.set_ylim(-0.7, len(chromosomes) - 0.3)
    ax.set_yticks([])
    ax.set_xlabel("Genomic coordinate within each chromosome (Mb)")
    ax.set_title(f"RNA fusion breakpoints across chromosomes ({plotted} plotted)")
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.grid(axis="x", color="#E8ECEF", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def tool_path(name: str) -> str:
    located = shutil.which(name)
    if located is None:
        raise FileNotFoundError(f"Required research-add-on tool not found on PATH: {name}")
    return located


def build_plan(
    results_dir: Path,
    gtf: Path,
    fasta: Path,
    star_index: Path,
    genebody_bed: Path,
    arriba_blacklist: Path,
    arriba_known_fusions: Path,
    arriba_protein_domains: Path,
    samples: Sequence[dict[str, str]],
    output_dir: Path,
    threads: int,
) -> list[dict[str, object]]:
    star = tool_path("STAR")
    arriba = tool_path("arriba")
    tin = tool_path("tin.py")
    duplication = tool_path("read_duplication.py")
    prepare_dexseq = tool_path("dexseq_prepare_annotation.py")
    count_dexseq = tool_path("dexseq_count.py")
    suppa = tool_path("suppa.py")

    plan: list[dict[str, object]] = []
    annotation = output_dir / "exon_counts/Ensembl104.flattened.gff"
    plan.append({
        "stage": "prepare_exon_annotation",
        "command": [prepare_dexseq, str(gtf), str(annotation)],
    })
    for row in samples:
        sample_id = row["sample_id"]
        bam = results_dir / f"stage4_full_alignment/bam/{sample_id}.bam"
        r1 = results_dir / (
            f"stage3_full_preprocessing/trimmed/{sample_id}_R1_trim.fastq.gz")
        r2 = results_dir / (
            f"stage3_full_preprocessing/trimmed/{sample_id}_R2_trim.fastq.gz")
        plan.extend([
            {
                "stage": "transcript_integrity",
                "sample_id": sample_id,
                "command": [tin, "-i", str(bam), "-r", str(genebody_bed)],
                "cwd": str(output_dir / f"qc/tin/{sample_id}"),
            },
            {
                "stage": "positional_duplication_estimate",
                "sample_id": sample_id,
                "command": [
                    duplication, "-i", str(bam), "-o",
                    str(output_dir / f"qc/duplication/{sample_id}"),
                ],
            },
            {
                "stage": "exon_bin_counts",
                "sample_id": sample_id,
                "command": [
                    count_dexseq, "-p", "yes", "-f", "bam", "-r", "pos",
                    "-s", "no", str(annotation), str(bam),
                    str(output_dir / f"exon_counts/{sample_id}.counts.txt"),
                ],
            },
        ])

        fusion_prefix = output_dir / f"fusions/star/{sample_id}_"
        fusion_bam = Path(f"{fusion_prefix}Aligned.out.bam")
        plan.append({
            "stage": "fusion_alignment",
            "sample_id": sample_id,
            "command": [
                star, "--genomeDir", str(star_index), "--runThreadN", str(threads),
                "--readFilesCommand", "zcat", "--readFilesIn", str(r1), str(r2),
                "--outFileNamePrefix", str(fusion_prefix),
                "--outSAMtype", "BAM", "Unsorted",
                "--outSAMunmapped", "Within", "--outSAMattributes", "All",
                "--outFilterMultimapNmax", "50", "--peOverlapNbasesMin", "10",
                "--alignSplicedMateMapLminOverLmate", "0.5",
                "--alignSJstitchMismatchNmax", "5", "-1", "5", "5",
                "--chimSegmentMin", "10", "--chimOutType", "WithinBAM", "HardClip",
                "--chimJunctionOverhangMin", "10", "--chimScoreDropMax", "30",
                "--chimScoreJunctionNonGTAG", "0", "--chimScoreSeparation", "1",
                "--chimSegmentReadGapMax", "3", "--chimMultimapNmax", "50",
            ],
        })
        plan.append({
            "stage": "fusion_call",
            "sample_id": sample_id,
            "command": [
                arriba, "-x", str(fusion_bam), "-g", str(gtf), "-a", str(fasta),
                "-b", str(arriba_blacklist), "-k", str(arriba_known_fusions),
                "-t", str(arriba_known_fusions), "-p", str(arriba_protein_domains),
                "-s", "no", "-o",
                str(output_dir / f"fusions/{sample_id}.fusions.tsv"),
                "-O", str(output_dir / f"fusions/{sample_id}.discarded.tsv"),
            ],
        })

    expression = output_dir / "splicing/transcript_tpm_matrix.tsv"
    plan.append({
        "stage": "transcript_splicing_events",
        "command": [
            suppa, "generateEvents", "-i", str(gtf),
            "-o", str(output_dir / "splicing/Ensembl104"),
            "-f", "ioe", "-e", "SE", "SS", "MX", "RI", "FL",
        ],
        "expression_matrix": str(expression),
    })
    return plan


def run_command(command: Sequence[str], cwd: Path, log_path: Path) -> None:
    cwd.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        subprocess.run(
            list(command),
            cwd=cwd,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
            text=True,
        )


def load_transcript_tpm(results_dir: Path, samples: Sequence[dict[str, str]]) -> dict[str, dict[str, float]]:
    matrix: dict[str, dict[str, float]] = {}
    all_ids: set[str] | None = None
    for row in samples:
        sample_id = row["sample_id"]
        quant_path = results_dir / f"final_results/02_transcript/quant_by_sample/{sample_id}/quant.sf"
        if not quant_path.is_file():
            raise FileNotFoundError(f"Salmon quantification missing: {quant_path}")
        values: dict[str, float] = {}
        with quant_path.open(newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if not {"Name", "TPM"}.issubset(reader.fieldnames or ()):
                raise ValueError(f"Unexpected Salmon quant.sf schema: {quant_path}")
            for item in reader:
                transcript = item["Name"]
                abundance = float(item["TPM"])
                if not transcript or not math.isfinite(abundance) or abundance < 0.0:
                    raise ValueError(f"Invalid Salmon TPM row in {quant_path}")
                if transcript in values:
                    raise ValueError(f"Duplicate transcript ID in {quant_path}: {transcript}")
                values[transcript] = abundance
        ids = set(values)
        if all_ids is None:
            all_ids = ids
        elif ids != all_ids:
            raise ValueError(f"Salmon transcript sets differ for sample {sample_id}")
        for transcript, abundance in values.items():
            matrix.setdefault(transcript, {})[sample_id] = abundance
    return matrix


def write_tpm_matrix(
    results_dir: Path, samples: Sequence[dict[str, str]], output: Path,
) -> None:
    matrix = load_transcript_tpm(results_dir, samples)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["transcript_id", *[row["sample_id"] for row in samples]])
        for transcript in sorted(matrix):
            writer.writerow([
                transcript,
                *[format(matrix[transcript][row["sample_id"]], ".10g")
                  for row in samples],
            ])


def execute_plan(
    plan: Sequence[dict[str, object]],
    results_dir: Path,
    samples: Sequence[dict[str, str]],
    output_dir: Path,
    fasta: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=False)
    plan_path = output_dir / "execution_plan.json"
    plan_path.write_text(json.dumps(list(plan), indent=2) + "\n")
    sample_by_id = {row["sample_id"]: row for row in samples}
    for index, item in enumerate(plan):
        stage = str(item["stage"])
        sample_id = str(item.get("sample_id", "cohort"))
        command = item.get("command")
        if not isinstance(command, list) or not all(isinstance(arg, str) for arg in command):
            raise ValueError(f"Invalid command in generated plan: {item}")
        if stage == "prepare_exon_annotation":
            cwd = output_dir / "exon_counts"
        elif stage == "transcript_splicing_events":
            write_tpm_matrix(
                results_dir, samples,
                output_dir / "splicing/transcript_tpm_matrix.tsv",
            )
            cwd = output_dir / "splicing"
            log = output_dir / f"logs/{index:03d}_{stage}_cohort.log"
            run_command(command, cwd, log)
            event_files = sorted(cwd.glob("*.ioe"))
            if not event_files:
                raise FileNotFoundError(
                    f"SUPPA2 did not produce event files in {cwd}")
            expression = cwd / "transcript_tpm_matrix.tsv"
            for event_path in event_files:
                psi_path = event_path.with_suffix(".psi")
                psi_command = [
                    tool_path("suppa.py"), "psiPerEvent",
                    "-i", str(event_path), "-e", str(expression),
                    "-o", str(psi_path),
                ]
                psi_log = output_dir / f"logs/{index:03d}_psi_{event_path.stem}.log"
                run_command(psi_command, cwd, psi_log)
                if not psi_path.is_file() or psi_path.stat().st_size == 0:
                    raise FileNotFoundError(f"SUPPA2 PSI output missing: {psi_path}")
            continue
        elif stage == "transcript_integrity":
            cwd = output_dir / f"qc/tin/{sample_id}"
        elif stage == "positional_duplication_estimate":
            cwd = output_dir / "qc/duplication"
        elif stage == "exon_bin_counts":
            cwd = output_dir / "exon_counts"
        elif stage == "fusion_alignment":
            cwd = output_dir / f"fusions/star/{sample_id}"
        elif stage == "fusion_call":
            cwd = output_dir / "fusions"
        else:
            raise ValueError(f"Unexpected research-add-on stage: {stage}")
        log = output_dir / f"logs/{index:03d}_{stage}_{sample_id}.log"
        run_command(command, cwd, log)
        if stage == "positional_duplication_estimate":
            expected = output_dir / f"qc/duplication/{sample_id}.dup.pos.DupRate.xls"
            if not expected.is_file() or expected.stat().st_size == 0:
                raise FileNotFoundError(
                    f"RSeQC positional duplication output missing or empty: {expected}")
        elif stage == "fusion_call":
            expected = output_dir / f"fusions/{sample_id}.fusions.tsv"
            if not expected.is_file():
                raise FileNotFoundError(f"Arriba fusion output missing: {expected}")

    run_qc_summaries(results_dir, output_dir, samples)
    collect_duplication_summary(output_dir, samples)
    collect_exon_counts(output_dir, samples)
    collect_psi(output_dir, samples)
    summarize_fusions(output_dir, samples, fasta)
    groups: dict[str, int] = {}
    for row in sample_by_id.values():
        groups[row["group"]] = groups.get(row["group"], 0) + 1
    readme = [
        "# RNA-only research add-ons",
        "",
        "This package contains supplementary analyses of full-depth RNA data.",
        "It does not modify standard workflow outputs and does not use matched DNA.",
        "",
        "## Contents",
        "",
        "- `qc/tin/`: transcript integrity profiles from RSeQC.",
        "- `qc/duplication/`: positional duplication estimates; not UMI molecule counts.",
        "- `qc_enhanced/`: per-sample BaseJumper and enhanced QC status, reasons, thresholds and TIN figure.",
        "- `exon_counts/`: DEXSeq flattened exon-bin counts.",
        "- `splicing/`: Ensembl 104 event definitions and per-sample SUPPA2 PSI outputs.",
        "- `fusions/`: fusion-specific STAR alignments, Arriba results, read evidence and chromosome map.",
        "- `logs/`: command logs; `execution_plan.json` records exact arguments.",
        "",
        "Per-condition differential testing is intentionally not run by this module.",
        f"Sample groups present: {json.dumps(groups, sort_keys=True)}.",
        "Fusion and splice outputs are research candidates, not clinical findings.",
        "No RNA variant calling or RNA-DNA identity comparison is performed.",
        "",
    ]
    (output_dir / "README.md").write_text("\n".join(readme))
    (output_dir / "run_status.json").write_text(json.dumps({
        "status": "PASS",
        "analysis_scope": "RNA-only descriptive research add-ons",
        "samples": len(samples),
        "differential_tests_run": False,
        "dna_used": False,
        "clinical_interpretation": False,
    }, indent=2) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--star-index", required=True, type=Path)
    parser.add_argument("--genebody-bed", required=True, type=Path)
    parser.add_argument("--arriba-blacklist", required=True, type=Path)
    parser.add_argument("--arriba-known-fusions", required=True, type=Path)
    parser.add_argument("--arriba-protein-domains", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--execute", action="store_true",
        help="Run the RNA-only add-ons; omit to print the audited execution plan only",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_dir = args.results_dir.resolve(strict=True)
    files = {
        "GTF": args.gtf,
        "FASTA": args.fasta,
        "STAR index": args.star_index,
        "gene-body BED": args.genebody_bed,
        "Arriba blacklist": args.arriba_blacklist,
        "Arriba known fusions": args.arriba_known_fusions,
        "Arriba protein domains": args.arriba_protein_domains,
    }
    for label, path in files.items():
        if not path.exists() or (path.is_file() and path.stat().st_size == 0):
            raise FileNotFoundError(f"{label} is missing or empty: {path}")
    if args.threads < 1:
        raise ValueError("--threads must be a positive integer")
    sample_sheet = results_dir / "rna_cells.tsv"
    samples = read_samples(sample_sheet)
    require_pass(results_dir / "step7_full_output_audit.json", "Full-depth report audit")
    require_pass(
        results_dir / "final_results/07_qc_provenance/salmon_output_audit.json",
        "Salmon audit",
    )
    output_dir = results_dir / OUTPUT_RELATIVE
    if output_dir.exists():
        raise FileExistsError(
            f"Refusing to overwrite research add-on results: {output_dir}")
    for row in samples:
        sample_id = row["sample_id"]
        for path in (
            results_dir / f"stage4_full_alignment/bam/{sample_id}.bam",
            results_dir / (
                f"stage3_full_preprocessing/trimmed/{sample_id}_R1_trim.fastq.gz"),
            results_dir / (
                f"stage3_full_preprocessing/trimmed/{sample_id}_R2_trim.fastq.gz"),
        ):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"Required full-depth RNA input is missing: {path}")

    plan = build_plan(
        results_dir, args.gtf.resolve(strict=True), args.fasta.resolve(strict=True),
        args.star_index.resolve(strict=True), args.genebody_bed.resolve(strict=True),
        args.arriba_blacklist.resolve(strict=True),
        args.arriba_known_fusions.resolve(strict=True),
        args.arriba_protein_domains.resolve(strict=True), samples,         output_dir,
        args.threads,
    )
    if args.execute:
        execute_plan(
            plan, results_dir, samples, output_dir, args.fasta.resolve(strict=True))
        print(f"PASS: RNA-only add-on results: {output_dir}")
    else:
        print(json.dumps({
            "status": "PLAN_ONLY",
            "samples": len(samples),
            "output_dir": str(output_dir),
            "commands": plan,
        }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
