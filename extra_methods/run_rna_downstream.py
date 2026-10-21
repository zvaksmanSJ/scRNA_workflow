#!/usr/bin/env python3
"""Create descriptive, single-condition RNA analyses and a self-contained dashboard.

The full run validates audited HTSeq, BaseJumper, Salmon, and alignment inputs;
builds normalized gene/transcript summaries, pathway and junction summaries,
and labeled figures; then writes the dashboard and results guide. The refresh
mode reuses completed package tables to rebuild only the dashboard, including
the optional research-add-on tabs when their run-status audit is PASS. This
module does not fit condition contrasts or exclude samples based on QC labels.
"""
from __future__ import annotations

import csv
import hashlib
import html
import json
import os
import re
import zipfile
from datetime import date, datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import hypergeom, spearmanr

RNAW = Path(__file__).resolve().parents[1]
RUNTIME = Path(os.environ.get("RNA_RUNTIME_DIR", RNAW / "wdl_runtime"))
RNAR = Path(os.environ.get("RESOLVEOME_RNA_RESULTS", RNAW / "wdl_runtime/results")).resolve()
OUT = RNAR / "final_results"
MATRIX = RNAR / "stage6_full_matrix/matrix/matrix_gene_counts_starhtseq.txt"
SAMPLES = Path(os.environ.get("RESOLVEOME_RNA_CELLS", RNAW / "data/rna_cells.local.tsv"))
STAR = RNAR / "stage4_full_alignment/star"
BASEJUMPER_METRICS = RNAR / "stage7_full_reports/consensus/rnaqc_all_metrics.tsv"
BASEJUMPER_REPORTS = RNAR / "stage7_full_reports"
TX2GENE = Path(os.environ.get(
    "RNA_TX2GENE", RUNTIME / "references/Ensembl104/Ensembl104.tx2gene.tsv"))
RESOURCES = Path(os.environ.get("RNA_ANALYSIS_RESOURCES", RNAW / "analysis_resources"))
REACTOME_SHA256 = "8c1dbc8578431da5d2d5118262718c60b553a9be3398e93658daa069e4a9afd4"
BASEJUMPER_SOURCE = Path(os.environ.get("BASEJ_RNA_SOURCE", RNAW / "vendor/basej-rnaqc"))
BASEJUMPER_COMMIT_FILE = Path(os.environ.get(
    "BASEJ_RNA_COMMIT_FILE", RNAW / "vendor/UPSTREAM_COMMIT"))

DIRS = (
    "00_dashboard", "01_expression", "02_transcript/quant_by_sample",
    "03_splicing", "04_pathways", "05_cell_line_profiles", "06_figures",
    "07_qc_provenance", "08_large_files", "09_references",
)

COLORS = {
    "maroon": "#6B1E3D",
    "gold": "#D99A22",
    "teal": "#197C83",
    "blue": "#315A86",
    "light": "#F4F6F8",
    "ink": "#24313B",
    "muted": "#64717B",
    "borderline": "#D99A22",
    "pass": "#197C83",
}

FIGURES = {
    "pca_samples.png": (
        "Sample PCA",
        "Principal component analysis of the 500 most variable genes using log2 "
        "size-factor-normalized counts. Each point is one independent biological "
        "replicate. Point color denotes the existing source QC label; these labels "
        "do not define biological groups."
    ),
    "sample_correlation.png": (
        "Sample-to-sample expression correlation",
        "Hierarchically clustered Spearman correlations across all count-matrix genes "
        "after log2 size-factor normalization. Higher correlation indicates more "
        "similar rank ordering of gene expression; it is descriptive, not a replicate "
        "exclusion test."
    ),
    "expression_distributions.png": (
        "Normalized gene-count distributions",
        "Per-replicate distributions of log2 size-factor-normalized gene counts. "
        "The plot helps compare library-wide expression distributions; it is not a "
        "differential-expression analysis."
    ),
    "variable_gene_heatmap.png": (
        "Most variable genes across replicates",
        "Row-scaled expression of the 40 genes with the greatest variance among "
        "replicates after filtering for mean raw count of at least 10. Red indicates "
        "above-row-average expression and blue below-row-average expression."
    ),
    "mean_expression_top_genes.png": (
        "Highest mean-abundance genes",
        "The 20 genes with the highest mean size-factor-normalized counts across "
        "the cohort. This is an abundance summary, not a cell-identity marker test."
    ),
    "reactome_variable_gene_enrichment.png": (
        "Reactome pathways enriched among variable genes",
        "Over-representation of Reactome gene sets among the 500 most variable "
        "expressed genes, tested against all expressed genes (mean raw count >=10). "
        "This identifies pathways represented among variable genes; it does not "
        "indicate pathway activation or condition-associated change."
    ),
    "reactome_sample_scores.png": (
        "Descriptive Reactome expression scores",
        "For each pathway, the plotted value is the mean gene-wise z-score of "
        "log2-normalized counts within each replicate. Scores show relative "
        "within-cohort expression patterns and are not formal ssGSEA or activity "
        "estimates."
    ),
    "transcript_tpm_heatmap.png": (
        "Variable transcript abundance",
        "Row-scaled Salmon TPM for the 40 transcripts with the greatest variance "
        "across replicates. TPM supports abundance comparisons within a sample; "
        "differential transcript usage was not tested."
    ),
    "junction_support_heatmap.png": (
        "Variable splice-junction read support",
        "Log10-transformed unique STAR splice-junction read counts for recurrent "
        "junctions with at least 20 reads overall and support in at least three "
        "replicates. These descriptive counts are not PSI values or differential "
        "splicing results."
    ),
}


def sample_short(sample: str) -> str:
    match = re.search(r"_Cell(\d+)_", sample)
    return f"Cell {match.group(1)}" if match else sample


def ensure_output_tree() -> None:
    for directory in DIRS:
        (OUT / directory).mkdir(parents=True, exist_ok=True)


def ensure_reference_link(source: Path, target: Path) -> None:
    source = source.resolve(strict=True)
    target.parent.mkdir(parents=True, exist_ok=True)
    relative_target = Path(os.path.relpath(source, target.parent))
    if target.is_symlink():
        if target.resolve(strict=True) != source:
            raise FileExistsError(f"Existing reference link points elsewhere: {target}")
        return
    if target.exists():
        if not target.is_file() or hashlib.sha256(target.read_bytes()).digest() != hashlib.sha256(
                source.read_bytes()).digest():
            raise FileExistsError(f"Refusing to replace non-identical reference file: {target}")
        target.unlink()
    target.symlink_to(relative_target)


def load_basejumper_metrics(samples: list[str]) -> pd.DataFrame:
    if not BASEJUMPER_METRICS.is_file():
        raise FileNotFoundError(f"BaseJumper per-sample metrics missing: {BASEJUMPER_METRICS}")
    frame = pd.read_csv(BASEJUMPER_METRICS, sep="\t")
    required = {
        "biosample", "qc_status", "qc_score", "genes_mitochondrial",
        "mt_gene_counts", "total_gene_counts", "mt_pct",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"BaseJumper metrics lack required fields: {sorted(missing)}")
    if frame.biosample.duplicated().any() or set(frame.biosample) != set(samples):
        raise ValueError("BaseJumper per-sample metrics do not match the sample manifest")
    if frame[list(required - {"biosample", "qc_status"})].isna().any().any():
        raise ValueError("BaseJumper QC metrics contain missing numeric values")
    return frame.set_index("biosample").loc[samples].reset_index()


def format_basejumper_value(column: str, value: object) -> str:
    if pd.isna(value):
        return "NA"
    fraction_fields = {
        "fastp_q20_rate", "fastp_q30_rate", "fastp_gc_content",
        "fastp_q30_rate_after", "fastp_duplication_rate",
        "fastp_pct_filtered", "fastp_pct_low_quality", "fastp_pct_too_short",
        "htseq_assigned_pct", "htseq_no_feature_pct", "htseq_ambiguous_pct",
        "htseq_alignment_not_unique_pct", "mt_pct", "rrna_pct",
        "prop_mappability", "prop_exonic", "prop_intergenic", "prop_intronic",
    }
    star_percent_fields = {
        "star_uniquely_mapped_pct", "star_multimapped_pct",
        "star_unmapped_mismatches_pct", "star_unmapped_tooshort_pct",
        "star_unmapped_other_pct", "star_chimeric_pct",
        "star_mismatch_rate", "star_deletion_rate", "star_insertion_rate",
    }
    if column in fraction_fields:
        return f"{float(value) * 100:.2f}%"
    if column in star_percent_fields:
        return f"{float(value):.2f}%"
    if isinstance(value, (float, np.floating)):
        return f"{float(value):,.3f}".rstrip("0").rstrip(".")
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"
    return str(value)


def render_basejumper_table(frame: pd.DataFrame, columns: list[str]) -> str:
    heading_html = "".join(
        f'<th title="{html.escape(column)}">{html.escape(column.replace("_", " "))}</th>'
        for column in columns)
    row_html = []
    for row in frame[columns].itertuples(index=False, name=None):
        cells = "".join(
            f"<td>{html.escape(format_basejumper_value(column, value))}</td>"
            for column, value in zip(columns, row))
        row_html.append(f"<tr>{cells}</tr>")
    return (
        f'<div class="table-scroll"><table class="metrics-table">'
        f"<thead><tr>{heading_html}</tr></thead><tbody>"
        f"{''.join(row_html)}</tbody></table></div>")


def render_status_table(frame: pd.DataFrame, columns: list[str]) -> str:
    headings = "".join(
        f'<th title="{html.escape(column)}">{html.escape(column.replace("_", " "))}</th>'
        for column in columns)
    rows = []
    for values in frame[columns].itertuples(index=False, name=None):
        cells = []
        for column, value in zip(columns, values):
            text = html.escape(format_basejumper_value(column, value))
            if column.endswith("status"):
                css_status = str(value).casefold()
                css_status = "borderline" if css_status == "borderline" else css_status
                text = f'<span class="qc-badge qc-{html.escape(css_status)}">{text}</span>'
            cells.append(f"<td>{text}</td>")
        rows.append(f"<tr>{''.join(cells)}</tr>")
    return (
        '<div class="table-scroll"><table class="metrics-table">'
        f"<thead><tr>{headings}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>")


def audit_sample_ids(audit: dict, name: str) -> set[str]:
    rows = audit.get("samples")
    if not isinstance(rows, list):
        raise ValueError(f"{name} audit has no sample list")
    if all(isinstance(row, str) for row in rows):
        sample_ids = rows
    elif all(isinstance(row, dict) for row in rows):
        sample_ids = [row.get("sample", row.get("sample_id")) for row in rows]
    else:
        raise ValueError(f"{name} audit has an unsupported sample-list schema")
    if any(not isinstance(sample, str) or not sample for sample in sample_ids):
        raise ValueError(f"{name} audit has missing sample identifiers")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError(f"{name} audit has duplicate sample identifiers")
    return set(sample_ids)


def validate_prerequisites() -> tuple[pd.DataFrame, pd.DataFrame, dict, dict]:
    required_audits = {
        "step4": RNAR / "step4_full_output_audit_quick.json",
        "step5": RNAR / "step5_full_output_audit.json",
        "step6": RNAR / "step6_full_output_audit.json",
        "step7": RNAR / "step7_full_output_audit.json",
        "salmon": OUT / "07_qc_provenance/salmon_output_audit.json",
    }
    audits = {name: json.loads(path.read_text()) for name, path in required_audits.items()}
    for name, audit in audits.items():
        if audit.get("status") != "PASS":
            raise RuntimeError(f"Required full-depth {name} audit is not PASS")
    alignment = audits["step4"]
    if alignment.get("depth") != "full" or alignment.get("audit_mode") != "quick":
        raise ValueError("Expected the documented full-depth quick alignment audit")
    metadata = pd.read_csv(SAMPLES, sep="\t")
    if list(metadata.columns) != ["sample_id", "r1", "r2", "group"]:
        raise ValueError("Unexpected sample sheet schema")
    if len(metadata) < 5 or metadata.sample_id.duplicated().any():
        raise ValueError("At least five unique biological replicates are required")
    expected_samples = set(metadata.sample_id)
    for name, audit in audits.items():
        if audit_sample_ids(audit, name) != expected_samples:
            raise ValueError(f"Full-depth {name} audit and sample sheet differ")
    if metadata.group.nunique() != 1:
        raise ValueError(
            "This report is configured for one biological condition; use a separate "
            "design-aware differential-expression workflow for multiple groups")
    counts = pd.read_csv(MATRIX, sep="\t")
    if "gene_id_gene_symbol" not in counts:
        raise ValueError("Unexpected full-depth matrix gene identifier column")
    samples = [c for c in counts.columns if c != "gene_id_gene_symbol"]
    expected = set(metadata.sample_id)
    if set(samples) != expected:
        raise ValueError("Count matrix and sample manifest have different membership")
    if counts.gene_id_gene_symbol.duplicated().any():
        raise ValueError("Duplicate gene rows in full-depth matrix")
    numeric = counts[samples].to_numpy()
    if not np.isfinite(numeric).all() or (numeric < 0).any():
        raise ValueError("Counts must be finite and nonnegative")
    if not np.equal(numeric, np.floor(numeric)).all():
        raise ValueError("The authoritative HTSeq matrix must contain integer counts")
    return counts, metadata, audits["step7"], audits


def basejumper_settings_audit(audits: dict) -> dict:
    source_main = BASEJUMPER_SOURCE / "main.nf"
    runtime_dir = RNAR / "stage4_full_alignment/runtime"
    runtime_main = runtime_dir / "main.nf"
    manifest_path = runtime_dir / "runtime_manifest.json"
    source_text = source_main.read_text()
    runtime_text = runtime_main.read_text()
    manifest = json.loads(manifest_path.read_text())
    source_hash = hashlib.sha256(source_main.read_bytes()).hexdigest()
    if manifest.get("sha256", {}).get("main.nf") != source_hash:
        raise ValueError("Full-depth runtime does not record the current upstream BaseJumper source")
    current_hashes = {
        str(path.relative_to(BASEJUMPER_SOURCE)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(BASEJUMPER_SOURCE.rglob("*"))
        if path.is_file() and path.suffix in (".nf", ".config", ".R", ".yaml")
    }
    if manifest.get("sha256") != current_hashes:
        raise ValueError("Vendored BaseJumper source no longer matches the full-depth runtime manifest")
    if audits["step4"].get("depth") != "full" or not audits["step6"].get(
            "all_matrix_values_match_original_HTSeq"):
        raise ValueError("Parity audit requires passing full-depth alignment and matrix audits")

    def normalized(text: str) -> str:
        return re.sub(r"\s+", " ", text)

    source = normalized(source_text)
    runtime = normalized(runtime_text)
    checks = {
        "fastp paired-end adapter detection": (
            "--detect_adapter_for_pe", "--detect_adapter_for_pe"),
        "STAR two-pass Basic": ("--twopassMode Basic", "--twopassMode Basic"),
        "STAR sorted BAM": (
            "--outSAMtype BAM SortedByCoordinate",
            "--outSAMtype BAM SortedByCoordinate"),
        "STAR multimapping limit": (
            "--outFilterMultimapNmax 10", "--outFilterMultimapNmax 10"),
        "STAR primary alignment preference": (
            "--outSAMprimaryFlag OneBestScore",
            "--outSAMprimaryFlag OneBestScore"),
        "STAR chimeric segment threshold": (
            "--chimSegmentMin 12", "--chimSegmentMin 12"),
        "STAR chimeric overhang threshold": (
            "--chimJunctionOverhangMin 8", "--chimJunctionOverhangMin 8"),
        "STAR chimeric junction format": (
            "--chimOutJunctionFormat 1", "--chimOutJunctionFormat 1"),
        "samtools secondary-record filter": (
            "samtools view -F 256 -h", "samtools view -F 256 -h"),
        "HTSeq unstranded gene-level union counting": (
            "-f bam -r pos -s no -t exon -i gene_id --additional-attr=gene_name -m union",
            "-f bam -r pos -s ${params.strandedness} -t exon -i gene_id "
            "--additional-attr=gene_name -m union"),
        "HTSeq ignores secondary alignments": (
            "--secondary-alignments=ignore", "--secondary-alignments=ignore"),
    }
    results = [
        {"check": name, "upstream_source": expected in source,
         "adapted_runtime": adapted in runtime,
         "status": "PASS" if expected in source and adapted in runtime else "FAIL"}
        for name, (expected, adapted) in checks.items()
    ]
    if any(row["status"] != "PASS" for row in results):
        raise RuntimeError("Upstream/adapted BaseJumper scientific-settings parity check failed")
    source_tree_hashes = manifest.get("sha256", {})
    return {
        "status": "PASS_SETTINGS_LEVEL_WITH_LIMITATIONS",
        "upstream_commit": BASEJUMPER_COMMIT_FILE.read_text().strip(),
        "source_main_sha256": source_hash,
        "runtime_main_sha256": hashlib.sha256(runtime_main.read_bytes()).hexdigest(),
        "runtime_adapter_sha256": manifest.get("adapter_sha256"),
        "source_manifest_file_count": len(source_tree_hashes),
        "scientific_setting_checks": results,
        "scope": (
            "Compared the vendored upstream source with the actual staged full-depth "
            "runtime and verified the accepted full-depth output audits."
        ),
        "known_adaptations": [
            "Native conda/LSF execution and bounded resource configuration.",
            "Full-depth branch copies all reads instead of applying the QC subsample.",
            "Runtime parser reads the final HTSeq column and raises on malformed counts.",
            "STAR BAM sort memory is explicitly bounded.",
            "Output publication, input validation, and source QC parser handling are adapted.",
        ],
        "not_claimed": (
            "No separate untouched stock-BaseJumper production rerun or bitwise "
            "output equivalence is claimed; tool/container versions differ."
        ),
    }


def split_gene_labels(values: pd.Series) -> pd.DataFrame:
    pieces = values.astype(str).str.split("_", n=1, expand=True)
    if pieces.shape[1] != 2:
        raise ValueError("Expected matrix labels formatted as gene_id_gene_symbol")
    return pd.DataFrame({"gene_id": pieces[0], "gene_symbol": pieces[1]})


def median_ratio_size_factors(counts: pd.DataFrame) -> pd.Series:
    """DESeq-style median ratios using nonzero counts (the poscounts variant)."""
    values = counts.to_numpy(dtype=np.float64)
    positive = values > 0
    log_values = np.zeros_like(values)
    np.log(values, out=log_values, where=positive)
    n_positive = positive.sum(axis=1)
    usable = n_positive > 0
    geometric_means = np.zeros(values.shape[0], dtype=np.float64)
    geometric_means[usable] = np.exp(
        log_values[usable].sum(axis=1) / n_positive[usable])
    ratios = np.full_like(values, np.nan)
    np.divide(values, geometric_means[:, None], out=ratios,
              where=positive & usable[:, None])
    factors = np.nanmedian(ratios, axis=0)
    if not np.isfinite(factors).all() or (factors <= 0).any():
        raise ValueError("Could not estimate positive size factors for every sample")
    factors /= np.exp(np.mean(np.log(factors)))
    return pd.Series(factors, index=counts.columns, name="size_factor")


def benjamini_hochberg(p_values) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError("P-values must be a finite one-dimensional vector in [0, 1]")
    n = len(p)
    if n == 0:
        return p.copy()
    order = np.argsort(p)
    adjusted = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(n, dtype=float)
    result[order] = np.minimum(adjusted, 1)
    return result


def load_reactome_sets(path: Path) -> dict[str, set[str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Reactome gene-set archive missing: {path}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != REACTOME_SHA256:
        raise ValueError(f"Reactome archive checksum mismatch: {digest}")
    with zipfile.ZipFile(path) as archive:
        names = [name for name in archive.namelist() if name.endswith(".gmt")]
        if len(names) != 1:
            raise ValueError(f"Expected one Reactome GMT file, found {names}")
        result = {}
        with archive.open(names[0]) as source:
            for raw in source:
                fields = raw.decode("utf-8").rstrip("\r\n").split("\t")
                if len(fields) >= 4:
                    label = f"{fields[0]} [{fields[1]}]"
                    result[label] = {gene.upper() for gene in fields[2:] if gene}
        if not result:
            raise ValueError("Reactome GMT contained no usable gene sets")
        return result


def run_pathway_enrichment(
        metrics: pd.DataFrame, gene_sets: dict[str, set[str]], samples: list[str],
        log_counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    expressed = metrics.loc[metrics.mean_raw_count >= 10]
    background = set(expressed.gene_symbol.str.upper()) - {""}
    variable = set(expressed.nlargest(min(500, len(expressed)), "variance_log2").gene_symbol.str.upper())
    variable &= background
    if len(background) < 1000 or len(variable) < 100:
        raise ValueError("Too few expressed/variable genes for a useful pathway summary")
    rows = []
    usable_sets = {}
    for name, genes in gene_sets.items():
        overlap_set = genes & background
        if not 10 <= len(overlap_set) <= 500:
            continue
        overlap = len(overlap_set & variable)
        p_value = float(hypergeom.sf(overlap - 1, len(background), len(overlap_set),
                                     len(variable)))
        rows.append({"pathway": name, "pathway_genes_in_background": len(overlap_set),
                     "variable_gene_overlap": overlap,
                     "variable_gene_ratio": overlap / len(variable),
                     "p_value": p_value})
        usable_sets[name] = overlap_set
    enrichment = pd.DataFrame(rows)
    if enrichment.empty:
        raise ValueError("No Reactome gene sets met the expressed-background size limits")
    enrichment["fdr_bh"] = benjamini_hochberg(enrichment.p_value)
    enrichment = enrichment.sort_values(["fdr_bh", "p_value", "pathway"]).reset_index(drop=True)

    symbols = (metrics.set_index("gene_id").reindex(log_counts.index)
               .gene_symbol.fillna("").str.upper())
    z = log_counts.copy()
    z = z.sub(z.mean(axis=1), axis=0)
    standard_deviation = log_counts.std(axis=1, ddof=0).replace(0, np.nan)
    z = z.div(standard_deviation, axis=0).replace([np.inf, -np.inf], np.nan)
    symbol_z = z.assign(_symbol=symbols.values).groupby("_symbol").mean()
    score_rows = []
    for name, members in usable_sets.items():
        available = [symbol for symbol in members if symbol in symbol_z.index]
        if len(available) >= 10:
            values = symbol_z.loc[available, samples].mean(axis=0)
            score_rows.append({"pathway": name, **values.to_dict(),
                               "genes_scored": len(available)})
    scores = pd.DataFrame(score_rows).set_index("pathway")
    if scores.empty:
        raise ValueError("Could not calculate sample-level Reactome descriptive scores")
    return enrichment, scores


def load_star_junctions(samples: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_tables = []
    summary = []
    for sample in samples:
        paths = list(STAR.rglob(f"{sample}_SJ.out.tab"))
        if len(paths) != 1:
            raise FileNotFoundError(f"Expected one STAR splice table for {sample}; found {len(paths)}")
        table = pd.read_csv(paths[0], sep="\t", header=None, dtype={0: str})
        if table.shape[1] != 9:
            raise ValueError(f"STAR splice table must have 9 columns: {paths[0]}")
        table.columns = ["chrom", "intron_start", "intron_end", "strand",
                         "motif", "annotated", "unique_reads", "multi_reads",
                         "max_overhang"]
        table["unique_reads"] = pd.to_numeric(table.unique_reads, errors="raise").astype("int64")
        table = table.groupby(["chrom", "intron_start", "intron_end", "strand"],
                              as_index=False).agg(
            unique_reads=("unique_reads", "sum"),
            multi_reads=("multi_reads", "sum"),
            annotated=("annotated", "max"),
            motif=("motif", "first"),
            max_overhang=("max_overhang", "max"))
        table["sample_id"] = sample
        sample_tables.append(table)
        summary.append({
            "sample_id": sample,
            "junctions_reported": len(table),
            "unique_junction_reads": int(table.unique_reads.sum()),
            "annotated_junctions": int(table.annotated.eq(1).sum()),
        })
    long = pd.concat(sample_tables, ignore_index=True)
    counts = long.pivot_table(index=["chrom", "intron_start", "intron_end", "strand"],
                              columns="sample_id", values="unique_reads",
                              fill_value=0, aggfunc="sum")
    counts = counts.reindex(columns=samples, fill_value=0).astype("int64")
    strand = {0: "?", 1: "+", 2: "-"}
    counts.index = counts.index.map(
        lambda x: f"{x[0]}:{x[1]}-{x[2]}:{strand[int(x[3])]}")
    counts.index.name = "junction"
    summary_df = pd.DataFrame(summary)
    return counts, summary_df


def load_salmon_matrices(
        samples: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    quant_root = OUT / "02_transcript/quant_by_sample"
    tpm = {}
    num_reads = {}
    expected_transcripts = None
    for sample in samples:
        path = quant_root / sample / "quant.sf"
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing Salmon quantification for {sample}; submit/validate with "
                "quantify_transcripts.py first")
        quant = pd.read_csv(path, sep="\t")
        expected = ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]
        if list(quant.columns) != expected or quant.Name.duplicated().any():
            raise ValueError(f"Invalid Salmon quant.sf: {path}")
        transcript_ids = set(quant.Name.astype(str))
        if expected_transcripts is None:
            expected_transcripts = transcript_ids
        elif transcript_ids != expected_transcripts:
            raise ValueError(f"Salmon transcript reference differs across samples: {path}")
        tpm[sample] = quant.set_index("Name").TPM
        num_reads[sample] = quant.set_index("Name").NumReads
    tpm_df = pd.DataFrame(tpm).fillna(0)
    reads_df = pd.DataFrame(num_reads).fillna(0)
    tx2gene = pd.read_csv(TX2GENE, sep="\t", header=None,
                          names=["transcript_id", "transcript_biotype", "gene_id",
                                 "gene_biotype", "gene_symbol"], dtype=str)
    tx2gene["transcript_id"] = tx2gene.transcript_id.str.replace(r"\.\d+$", "", regex=True)
    tx2gene = tx2gene.drop_duplicates("transcript_id").set_index("transcript_id")
    quant_ids = pd.Index(tpm_df.index.astype(str))
    map_ids = quant_ids.str.replace(r"\.\d+$", "", regex=True)
    mapped = tx2gene.reindex(map_ids)
    mapped.index = quant_ids
    gene_ids = mapped.gene_id
    if gene_ids.nunique() < 1000:
        raise ValueError("Salmon transcript-to-gene mapping matched too few annotated genes")
    keep = gene_ids.notna()
    mapped_tpm = tpm_df.loc[keep].copy()
    mapped_reads = reads_df.loc[keep].copy()
    mapped_tpm["gene_id"] = gene_ids.loc[keep].values
    mapped_reads["gene_id"] = gene_ids.loc[keep].values
    gene_tpm = mapped_tpm.groupby("gene_id", sort=True)[samples].sum()
    gene_reads = mapped_reads.groupby("gene_id", sort=True)[samples].sum()
    mapped_total = gene_reads.sum(axis=0)
    all_total = reads_df.sum(axis=0)
    mapping_summary = pd.DataFrame({
        "sample_id": samples,
        "estimated_reads_mapped_to_gene_annotation": mapped_total.reindex(samples).values,
        "all_salmon_estimated_reads_including_decoys": all_total.reindex(samples).values,
        "fraction_estimated_reads_gene_mapped":
            (mapped_total / all_total).reindex(samples).values,
        "transcripts_quantified": len(tpm_df.index),
        "annotated_transcripts": int(keep.sum()),
    })
    if (mapping_summary.fraction_estimated_reads_gene_mapped <= 0).any():
        raise ValueError("Salmon estimated reads did not map to annotated transcripts")
    return tpm_df.loc[keep], reads_df.loc[keep], gene_tpm, mapping_summary


def set_plot_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 13,
        "axes.titleweight": "semibold",
        "axes.labelsize": 10,
        "axes.edgecolor": "#CBD2D9",
        "axes.labelcolor": COLORS["ink"],
        "text.color": COLORS["ink"],
        "xtick.color": COLORS["muted"],
        "ytick.color": COLORS["muted"],
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    })


def save_figure(fig, name: str, captions: list[dict]) -> None:
    path = OUT / "06_figures" / name
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    captions.append({"figure": name, "title": FIGURES[name][0],
                     "explanation": FIGURES[name][1]})


def plot_expression_figures(
        metrics: pd.DataFrame, log_counts: pd.DataFrame, samples: list[str],
        qc: pd.DataFrame, captions: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    short = [sample_short(s) for s in samples]
    variable = metrics.loc[metrics.mean_raw_count >= 10].nlargest(
        min(500, int((metrics.mean_raw_count >= 10).sum())), "variance_log2")
    variable_ids = variable.gene_id.tolist()
    data = log_counts.loc[variable_ids, samples]
    centered = data.sub(data.mean(axis=1), axis=0)
    u, singular, _ = np.linalg.svd(centered.T.to_numpy(), full_matrices=False)
    scores = u * singular
    pct = singular ** 2 / np.sum(singular ** 2) * 100
    score_df = pd.DataFrame({"sample_id": samples, "PC1": scores[:, 0],
                             "PC2": scores[:, 1],
                             "PC1_variance_percent": pct[0],
                             "PC2_variance_percent": pct[1]})
    score_df = score_df.merge(qc[["sample_id", "score", "status"]], on="sample_id", validate="one_to_one")
    fig, ax = plt.subplots(figsize=(9.2, 6.1))
    styles = [("PASS", COLORS["teal"], "o"), ("Borderline", COLORS["borderline"], "D"),
              ("FAIL", "#A23A3A", "X")]
    for status, color, marker in styles:
        subset = score_df.loc[score_df.status.eq(status)]
        if subset.empty:
            continue
        ax.scatter(subset.PC1, subset.PC2, s=68, c=color, marker=marker,
                   edgecolors="white", linewidths=0.8, label=f"{status} QC label", zorder=3)
        for _, row in subset.iterrows():
            ax.annotate(sample_short(row.sample_id), (row.PC1, row.PC2),
                        xytext=(5, 5), textcoords="offset points", fontsize=8)
    ax.axhline(0, color="#DDE2E7", lw=0.8, zorder=0)
    ax.axvline(0, color="#DDE2E7", lw=0.8, zorder=0)
    ax.set_xlabel(f"PC1 ({pct[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pct[1]:.1f}% variance)")
    ax.set_title("Replicate expression structure")
    ax.text(0.01, -0.17, "PCA on the 500 most variable genes; QC color is annotation only.",
            transform=ax.transAxes, fontsize=8, color=COLORS["muted"])
    ax.legend(frameon=False, loc="best")
    save_figure(fig, "pca_samples.png", captions)

    corr = spearmanr(log_counts[samples].to_numpy(), axis=0).statistic
    if not np.isfinite(corr).all():
        raise ValueError("Non-finite sample correlation")
    distance = np.clip(1 - corr, 0, 2)
    np.fill_diagonal(distance, 0)
    order = leaves_list(linkage(squareform(distance, checks=False), method="average"))
    ordered_corr = corr[np.ix_(order, order)]
    ordered_labels = [short[i] for i in order]
    fig, ax = plt.subplots(figsize=(9.5, 8.4))
    image = ax.imshow(ordered_corr, cmap="viridis", vmin=0.7, vmax=1, aspect="auto")
    ax.set_xticks(range(len(samples)), ordered_labels, rotation=55, ha="right", fontsize=8)
    ax.set_yticks(range(len(samples)), ordered_labels, fontsize=8)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Sample")
    ax.set_title("Spearman correlation of normalized expression")
    cbar = fig.colorbar(image, ax=ax, shrink=0.78, pad=0.02)
    cbar.set_label("Spearman correlation coefficient")
    save_figure(fig, "sample_correlation.png", captions)
    corr_df = pd.DataFrame(corr, index=samples, columns=samples)
    corr_df.index.name = "sample_id"

    fig, ax = plt.subplots(figsize=(11.0, 5.8))
    distributions = [log_counts[s].to_numpy() for s in samples]
    parts = ax.violinplot(distributions, showmeans=False, showmedians=True,
                          showextrema=False, widths=0.82)
    for body in parts["bodies"]:
        body.set_facecolor(COLORS["teal"])
        body.set_edgecolor(COLORS["teal"])
        body.set_alpha(0.62)
    parts["cmedians"].set_color(COLORS["maroon"])
    parts["cmedians"].set_linewidth(1.5)
    ax.set_xticks(range(1, len(samples) + 1), short, rotation=45, ha="right")
    ax.set_xlabel("Biological replicate")
    ax.set_ylabel("log2(size-factor-normalized count + 1)")
    ax.set_title("Normalized expression distributions")
    save_figure(fig, "expression_distributions.png", captions)

    top40 = variable.head(40).copy()
    z = log_counts.loc[top40.gene_id, samples]
    z = z.sub(z.mean(axis=1), axis=0).div(
        log_counts.loc[top40.gene_id, samples].std(axis=1, ddof=0).replace(0, np.nan),
        axis=0).fillna(0).clip(-2.5, 2.5)
    row_labels = [f"{r.gene_symbol} ({r.gene_id})" for r in top40.itertuples()]
    fig, ax = plt.subplots(figsize=(11.0, 9.0))
    image = ax.imshow(z.to_numpy(), cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_xticks(range(len(samples)), short, rotation=45, ha="right")
    ax.set_yticks(range(len(row_labels)), row_labels, fontsize=7)
    ax.set_xlabel("Biological replicate")
    ax.set_ylabel("Gene (symbol; Ensembl ID)")
    ax.set_title("Top 40 most variable expressed genes")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.74)
    colorbar.set_label("Gene-wise z-score of log2 normalized counts")
    save_figure(fig, "variable_gene_heatmap.png", captions)
    return score_df, corr_df


def plot_abundance(metrics: pd.DataFrame, captions: list[dict]) -> None:
    top = metrics.nlargest(20, "mean_normalized_count").sort_values(
        "mean_normalized_count", ascending=True)
    labels = [f"{row.gene_symbol} ({row.gene_id})" for row in top.itertuples()]
    fig, ax = plt.subplots(figsize=(10.5, 7.0))
    ax.barh(labels, top.mean_normalized_count, color=COLORS["maroon"])
    ax.set_xlabel("Mean size-factor-normalized gene count")
    ax.set_ylabel("Gene (symbol; Ensembl ID)")
    ax.set_title("Twenty highest-abundance genes in the cohort")
    save_figure(fig, "mean_expression_top_genes.png", captions)


def plot_pathways(
        enrichment: pd.DataFrame, scores: pd.DataFrame,
        samples: list[str], captions: list[dict]) -> None:
    top = enrichment.head(20).copy().sort_values("variable_gene_ratio")
    labels = [s.split(" [")[0] for s in top.pathway]
    fig, ax = plt.subplots(figsize=(11.0, 8.0))
    marker_sizes = 40 + top.variable_gene_overlap.to_numpy() * 4
    point_colors = -np.log10(top.fdr_bh.clip(lower=1e-300).to_numpy())
    scatter = ax.scatter(top.variable_gene_ratio, labels, s=marker_sizes,
                         c=point_colors, cmap="magma", edgecolors="white",
                         linewidths=0.6)
    ax.set_xlabel("Variable genes overlapping pathway / all selected variable genes")
    ax.set_ylabel("Reactome pathway")
    ax.set_title("Reactome over-representation among variable genes")
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label("-log10 Benjamini–Hochberg FDR")
    ax.text(0.01, -0.13, "Point size = number of overlapping variable genes.",
            transform=ax.transAxes, fontsize=8, color=COLORS["muted"])
    save_figure(fig, "reactome_variable_gene_enrichment.png", captions)

    sample_scores = scores.loc[:, samples]
    top_names = sample_scores.var(axis=1).nlargest(
        min(35, len(sample_scores))).index
    top_scores = sample_scores.loc[top_names]
    display = top_scores.copy()
    display.index = [name.split(" [")[0] for name in display.index]
    fig, ax = plt.subplots(figsize=(11.0, 9.0))
    image = ax.imshow(display.to_numpy(), cmap="RdBu_r", vmin=-1.5, vmax=1.5,
                      aspect="auto")
    samples = list(display.columns)
    ax.set_xticks(range(len(samples)), [sample_short(s) for s in samples],
                  rotation=45, ha="right")
    ax.set_yticks(range(len(display)), display.index, fontsize=7)
    ax.set_xlabel("Biological replicate")
    ax.set_ylabel("Reactome pathway")
    ax.set_title("Relative pathway expression scores")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.76)
    colorbar.set_label("Mean gene-wise z-score (descriptive)")
    save_figure(fig, "reactome_sample_scores.png", captions)


def plot_transcripts(tpm: pd.DataFrame, captions: list[dict]) -> pd.DataFrame:
    sample_mean = tpm.mean(axis=1)
    variance = tpm.var(axis=1, ddof=0)
    summary = pd.DataFrame({"transcript_id": tpm.index,
                            "mean_TPM": sample_mean.values,
                            "variance_TPM": variance.values})
    keep = summary.mean_TPM >= 0.1
    top = summary.loc[keep].nlargest(min(40, int(keep.sum())), "variance_TPM")
    if top.empty:
        raise ValueError("No expressed transcripts available for the TPM heatmap")
    data = tpm.loc[top.transcript_id]
    z = np.log2(data + 0.1)
    z = z.sub(z.mean(axis=1), axis=0).div(
        z.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0).clip(-2.5, 2.5)
    fig, ax = plt.subplots(figsize=(11, 8.8))
    image = ax.imshow(z.to_numpy(), cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_xticks(range(len(z.columns)), [sample_short(s) for s in z.columns],
                  rotation=45, ha="right")
    ax.set_yticks(range(len(z.index)), list(z.index), fontsize=6)
    ax.set_xlabel("Biological replicate")
    ax.set_ylabel("Salmon transcript ID")
    ax.set_title("Most variable expressed transcripts")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.75)
    colorbar.set_label("Transcript-wise z-score of log2(TPM + 0.1)")
    save_figure(fig, "transcript_tpm_heatmap.png", captions)
    return summary


def plot_junctions(counts: pd.DataFrame, captions: list[dict]) -> pd.DataFrame:
    prevalence = counts.gt(0).sum(axis=1)
    eligible = counts.loc[counts.sum(axis=1).ge(20) & prevalence.ge(3)]
    if eligible.empty:
        raise ValueError("No recurrent splice junctions met the plotting thresholds")
    metrics = pd.DataFrame({
        "junction": eligible.index,
        "total_unique_reads": eligible.sum(axis=1).values,
        "replicates_with_support": eligible.gt(0).sum(axis=1).values,
        "mean_unique_reads": eligible.mean(axis=1).values,
        "variance_unique_reads": eligible.var(axis=1, ddof=0).values,
    })
    top = metrics.nlargest(min(35, len(metrics)), "variance_unique_reads")
    data = eligible.loc[top.junction]
    transformed = np.log10(data + 1)
    fig, ax = plt.subplots(figsize=(11.0, 9.2))
    image = ax.imshow(transformed.to_numpy(), cmap="YlGnBu", aspect="auto")
    ax.set_xticks(range(len(data.columns)), [sample_short(s) for s in data.columns],
                  rotation=45, ha="right")
    ax.set_yticks(range(len(data.index)), list(data.index), fontsize=6)
    ax.set_xlabel("Biological replicate")
    ax.set_ylabel("Genomic splice junction (chromosome:start-end:strand)")
    ax.set_title("Most variable recurrent splice-junction support")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.75)
    colorbar.set_label("log10(unique STAR junction reads + 1)")
    save_figure(fig, "junction_support_heatmap.png", captions)
    return metrics.sort_values("variance_unique_reads", ascending=False)


def save_tsv(frame: pd.DataFrame, path: Path, index: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=index, na_rep="NA")


def build_dashboard(
        metrics: pd.DataFrame, qc: pd.DataFrame, size_factors: pd.Series,
        enrichment: pd.DataFrame, salmon_summary: pd.DataFrame,
        junction_summary: pd.DataFrame, audits: dict, captions: list[dict],
        basejumper_metrics: pd.DataFrame) -> None:
    dashboard_dir = OUT / "00_dashboard"
    figure_dir = OUT / "06_figures"
    logo_path = OUT / "09_references/stjude_logo.svg"
    source_logo = Path(os.environ.get("RNA_LOGO", RESOURCES / "stjude_logo.svg"))
    reactome_archive = RESOURCES / "ReactomePathways.gmt.zip"
    if not reactome_archive.is_file():
        raise FileNotFoundError(reactome_archive)
    if source_logo.is_file():
        ensure_reference_link(source_logo, logo_path)
        logo = logo_path.read_text()
        if "<script" in logo.lower() or "<!doctype" in logo.lower():
            raise ValueError("Logo SVG contains unexpected active markup")
        logo_markup = (
            '<img src="../09_references/stjude_logo.svg" '
            'alt="St. Jude Children&#39;s Research Hospital">')
    else:
        logo_markup = '<div class="brand-mark" aria-label="ResolveOME RNA">ResolveOME<br>RNA</div>'
    ensure_reference_link(reactome_archive, OUT / "09_references/ReactomePathways.gmt.zip")
    ensure_reference_link(RESOURCES / "README.txt", OUT / "09_references/README.txt")
    qc_counts = qc.status.value_counts().to_dict()
    sample_count = len(qc)
    condition = str(qc.condition.iloc[0])
    borderline = ", ".join(qc.loc[qc.status.eq("Borderline"), "sample_short"]) or "None"
    enriched_significant = int((enrichment.fdr_bh < 0.05).sum())
    basejumper_by_sample = basejumper_metrics.set_index("biosample")
    mt_summary = basejumper_metrics[
        ["biosample", "mt_gene_counts", "total_gene_counts", "mt_pct",
         "genes_mitochondrial", "qc_status", "qc_score"]].copy()
    mt_summary = mt_summary.rename(columns={
        "biosample": "Sample",
        "mt_gene_counts": "MT gene-assigned reads",
        "total_gene_counts": "Total gene-assigned reads",
        "mt_pct": "MT reads (% of gene-assigned reads)",
        "genes_mitochondrial": "MT genes detected",
        "qc_status": "BaseJumper QC",
        "qc_score": "QC score (0-5)",
    })
    mt_summary_html = render_basejumper_table(
        mt_summary, list(mt_summary.columns))
    all_metric_columns = list(basejumper_metrics.columns)
    all_metrics_html = render_basejumper_table(basejumper_metrics, all_metric_columns)
    basejumper_figure_specs = [
        ("PCA of top 500 genes", "BaseJumper PCA projection of the top 500 highly variable genes.",
         "pca/PCA_mqc.png"),
        ("Expression heatmap", "BaseJumper clustered expression heatmap for its variable-gene set.",
         "pca/Heatmap__mqc.png"),
        ("Mitochondrial count fraction", "BaseJumper plot of mitochondrial gene counts divided by total assigned gene counts; the red line marks 10%.",
         "pca/MT_mqc.png"),
        ("Protein-coding genes detected", "Number of protein-coding genes with detected counts in each sample.",
         "pca/GenesDetected_mqc.png"),
        ("BaseJumper sample-name view", "Original BaseJumper sample-label view used with its cohort QC plots.",
         "pca/samplenames_mqc.png"),
        ("RNA composition QC", "BaseJumper consensus composition plot; see the per-sample values below for exact percentages.",
         "consensus/RNAQC_composition_mqc.jpg"),
    ]
    basejumper_figures_html = []
    for title, explanation, relative_path in basejumper_figure_specs:
        asset = BASEJUMPER_REPORTS / relative_path
        if not asset.is_file():
            raise FileNotFoundError(f"BaseJumper report figure missing: {asset}")
        href = os.path.relpath(asset, dashboard_dir)
        basejumper_figures_html.append(
            f'<article class="figure-card"><h3>{html.escape(title)}</h3>'
            f'<a href="{html.escape(href)}"><img src="{html.escape(href)}" '
            f'alt="{html.escape(title)}"></a><p>{html.escape(explanation)}</p></article>')
    figures_html = []
    for caption in captions:
        name = caption["figure"]
        figures_html.append(
            f'<article class="figure-card"><h3>{html.escape(caption["title"])}</h3>'
            f'<a href="../06_figures/{html.escape(name)}">'
            f'<img src="../06_figures/{html.escape(name)}" alt="{html.escape(caption["title"])}">'
            f'</a><p>{html.escape(caption["explanation"])}</p></article>')
    result_links = [
        ("Raw gene-count matrix",
         "../01_expression/raw_gene_counts.tsv",
         "Authoritative integer HTSeq gene counts; use for count-based models."),
        ("Size-factor normalized counts",
         "../01_expression/normalized_gene_counts.tsv",
         "Descriptive normalized counts; not a DESeq2 model result."),
        ("Log2 normalized expression",
         "../01_expression/log2_normalized_gene_counts.tsv",
         "Exploratory expression values used in cohort plots."),
        ("Salmon transcript TPM matrix",
         "../02_transcript/transcript_tpm.tsv",
         "Transcript abundance estimates; no differential transcript usage test."),
        ("Salmon transcript NumReads matrix",
         "../02_transcript/transcript_numreads.tsv",
         "Estimated transcript-level counts; not interchangeable with HTSeq gene counts."),
        ("Salmon transcript-to-gene mapping summary",
         "../02_transcript/transcript_mapping_summary.tsv",
         "Estimated-read annotation coverage and transcript mapping counts by replicate."),
        ("STAR splice-junction count matrix",
         "../03_splicing/junction_unique_read_counts.tsv",
         "Descriptive unique read support; not PSI or differential-splicing results."),
        ("Per-replicate STAR junction summary",
         "../03_splicing/sample_junction_summary.tsv",
         "Number of reported junctions and unique junction reads in each replicate."),
        ("Reactome variable-gene enrichment",
         "../04_pathways/reactome_variable_gene_enrichment.tsv",
         "Exploratory over-representation among variable genes, with FDR."),
        ("Per-replicate pathway expression scores",
         "../04_pathways/reactome_descriptive_sample_scores.tsv",
         "Mean gene-wise z-score by set and replicate; descriptive only."),
        ("Top replicate-enriched genes",
         "../05_cell_line_profiles/sample_enriched_genes.tsv",
         "Within-cohort descriptive outliers; not validated cell-line markers."),
        ("Full-depth MultiQC",
         "../../stage7_full_reports/multiqc/multiqc_report.html",
         "Upstream QC and tool summary."),
        ("BaseJumper complete per-sample metrics",
         "../../stage7_full_reports/consensus/rnaqc_all_metrics.tsv",
         f"All {len(basejumper_metrics.columns)} BaseJumper QC fields for all "
         f"{len(basejumper_metrics)} samples, including mitochondrial and rRNA counts/fractions, "
         "Fastp, STAR, HTSeq, Qualimap and QC scores."),
        ("BaseJumper selected QC metrics",
         "../../stage7_full_reports/consensus/combined_selected_metrics.txt",
         "BaseJumper's selected summary table with QC score, MT fraction, protein-coding genes, read totals, Fastp quality/filtering and intronic fraction."),
        ("Qualimap raw RNA-seq metrics",
         "../../stage7_full_reports/qualimap_summary/qualimap_stats_mqc.csv",
         "Per-sample raw Qualimap alignment and feature-region counts, including exonic/intronic/intergenic, rRNA and junction read counts."),
        ("BaseJumper consensus verdict",
         "../../stage7_full_reports/consensus/summary_verdict_group.txt",
         "Original BaseJumper group-level consensus verdict."),
        ("BaseJumper consensus scores",
         "../../stage7_full_reports/consensus/RNA-QC_ConsensusScores.txt",
         "Original BaseJumper per-sample consensus scoring output."),
        ("BaseJumper score distribution",
         "../../stage7_full_reports/consensus/RNA-QC_ConsensusScores_SummaryTable_mqc.txt",
         "Sample counts and proportions in each BaseJumper QC score category."),
        ("Full-depth output audit",
         "../07_qc_provenance/analysis_summary.json",
         "Workflow audits, analysis methods and cohort status."),
        ("BaseJumper settings parity audit",
         "../07_qc_provenance/basejumper_settings_parity.json",
         "Checks adapted runtime scientific settings against the vendored source; not a separate stock run."),
        ("Alignment audit (quick)",
         "../07_qc_provenance/source_audits/step4_full_output_audit_quick.json",
         "Full cohort quick checks plus the recorded deep-check subset."),
        ("Counts/QC audit",
         "../07_qc_provenance/source_audits/step5_full_output_audit.json",
         "Per-sample count and QC output verification."),
        ("Matrix audit",
         "../07_qc_provenance/source_audits/step6_full_output_audit.json",
         "Full-depth cohort matrix membership and count-agreement verification."),
        ("Reports/QC audit",
         "../07_qc_provenance/source_audits/step7_full_output_audit.json",
         "Final source QC scores and report output verification."),
        ("Salmon output audit",
         "../07_qc_provenance/source_audits/salmon_output_audit.json",
         "All-sample matched-index Salmon quantification validation."),
    ]
    figure_captions = pd.DataFrame(captions)
    save_tsv(figure_captions, OUT / "07_qc_provenance/figure_captions.tsv")
    fdr_wording = ("No Reactome sets passed 5% FDR." if enriched_significant == 0
                   else f"{enriched_significant} Reactome sets passed 5% FDR.")
    alignment_mode = audits["step4"].get("audit_mode", "not recorded")
    addon_dir = RNAR / "final_results/10_rna_research_addons"
    addon_status_path = addon_dir / "run_status.json"
    addon_status = None
    if addon_status_path.is_file():
        addon_status = json.loads(addon_status_path.read_text())
    if addon_status and addon_status.get("status") == "PASS":
        enhanced_path = addon_dir / "qc_enhanced/per_cell_qc_status.tsv"
        enhanced = pd.read_csv(enhanced_path, sep="\t")
        required_enhanced = {
            "sample_id", "basejumper_status", "basejumper_score_0_to_5",
            "basejumper_reasons", "median_transcript_integrity_number",
            "tin_robust_z", "improved_status", "improved_reasons",
        }
        if not required_enhanced.issubset(enhanced.columns):
            raise ValueError(f"Enhanced QC table has an unexpected schema: {enhanced_path}")
        if enhanced.sample_id.duplicated().any() or set(enhanced.sample_id) != set(qc.sample_id):
            raise ValueError("Enhanced QC sample membership differs from the dashboard cohort")
        enhanced_columns = [
            "sample_id", "basejumper_status", "basejumper_score_0_to_5",
            "basejumper_reasons", "median_transcript_integrity_number",
            "tin_robust_z", "improved_status", "improved_reasons",
        ]
        enhanced_table_html = render_status_table(enhanced, enhanced_columns)
        tin_figure_html = (
            '<article class="figure-card"><h3>Enhanced transcript-integrity QC</h3>'
            '<a href="../10_rna_research_addons/qc_enhanced/tin_by_sample.png">'
            '<img src="../10_rna_research_addons/qc_enhanced/tin_by_sample.png" '
            'alt="Median transcript integrity number by sample"></a>'
            '<p>Each bar is the sample median RSeQC TIN. Colors show the combined enhanced-QC '
            'status; dashed thresholds are cohort-relative robust-z flags, not universal '
            'absolute cutoffs. No sample is automatically excluded.</p></article>')
        duplication_figure = addon_dir / "qc/duplication/positional_duplication.png"
        duplication_figure_html = ""
        if duplication_figure.is_file():
            duplication_figure_html = (
                '<article class="figure-card"><h3>Positional duplication estimate</h3>'
                '<a href="../10_rna_research_addons/qc/duplication/positional_duplication.png">'
                '<img src="../10_rna_research_addons/qc/duplication/positional_duplication.png" '
                'alt="RSeQC positional duplication estimates by sample"></a>'
                '<p>Estimated from RSeQC positional occurrence distributions. Identical mapping '
                'positions do not prove PCR duplication in non-UMI data; these values are not '
                'molecule counts and are not used as an exclusion threshold.</p></article>')
        enhanced_qc_html = (
            '<section><h2>Per-cell BaseJumper and enhanced QC</h2>'
            '<p>BaseJumper status is retained separately from enhanced status. Enhanced status '
            'adds low-TIN cohort outlier flags (robust z ≤−2 borderline; ≤−3 fail; minimum '
            '8 samples). The table gives the reason for every borderline/fail call. These '
            'are research QC annotations, not clinical thresholds.</p>'
            f'<div class="figure-grid">{tin_figure_html}{duplication_figure_html}</div>'
            f'<h3>Statuses and call reasons</h3>{enhanced_table_html}'
            '<p class="small">BaseJumper failures remain failures in the enhanced label; '
            'enhanced checks do not remove cells.</p></section>')
        exon_figure = addon_dir / "exon_counts/top_variable_exon_bins.png"
        psi_figure = addon_dir / "splicing/variable_event_psi.png"
        splicing_cards = []
        if exon_figure.is_file():
            splicing_cards.append(
                '<article class="figure-card"><h3>Variable exon-bin counts</h3>'
                '<a href="../10_rna_research_addons/exon_counts/top_variable_exon_bins.png">'
                '<img src="../10_rna_research_addons/exon_counts/top_variable_exon_bins.png" '
                'alt="Most variable DEXSeq exon-bin counts"></a>'
                '<p>Top variable flattened exon bins from DEXSeq counts; values are log2(count + 1). '
                'This is descriptive exon-level quantification, not differential exon usage.</p>'
                '</article>')
        if psi_figure.is_file():
            splicing_cards.append(
                '<article class="figure-card"><h3>Variable event PSI</h3>'
                '<a href="../10_rna_research_addons/splicing/variable_event_psi.png">'
                '<img src="../10_rna_research_addons/splicing/variable_event_psi.png" '
                'alt="Variable SUPPA2 event PSI across samples"></a>'
                '<p>SUPPA2 percent-spliced-in estimates for the most variable events across '
                'the cohort (scale 0–1). No differential-splicing test was performed.</p>'
                '</article>')
        splicing_html = (
            '<section><h2>Exon and splicing add-ons</h2>'
            '<p>DEXSeq exon-bin counts and SUPPA2 event PSI are descriptive summaries across '
            'the cohort. With one condition, no differential exon usage or differential '
            'splicing inference is performed.</p>'
            f'<div class="figure-grid">{"".join(splicing_cards)}</div></section>')
        fusion_path = addon_dir / "fusions/fusion_evidence_summary.tsv"
        fusion_figure = addon_dir / "fusions/fusion_chromosome_map.png"
        fusion_rows_html = (
            '<p>Fusion summary is unavailable in the completed add-on output.</p>')
        fusion_count = 0
        if fusion_path.is_file():
            fusions = pd.read_csv(fusion_path, sep="\t")
            fusion_count = len(fusions)
            if fusion_count:
                tier_order = {"HIGH": 0, "REVIEW": 1, "LOW_OR_FILTERED": 2}
                fusions["_tier_order"] = fusions.evidence_tier.map(tier_order).fillna(3)
                fusions = fusions.sort_values(
                    ["_tier_order", "total_supporting_reads"], ascending=[True, False])
                columns = [
                    "sample_id", "gene1", "gene2", "breakpoint1", "breakpoint2",
                    "confidence", "filters", "split_reads1", "split_reads2",
                    "discordant_mates", "total_supporting_reads", "evidence_tier",
                ]
                fusion_rows_html = render_basejumper_table(
                    fusions.head(200), columns)
                if fusion_count > 200:
                    fusion_rows_html += (
                        f'<p class="small">Showing top 200 of {fusion_count} candidates; '
                        'download the complete evidence table below.</p>')
            else:
                fusion_rows_html = '<p>No candidates passed Arriba filters.</p>'
        fusion_image_html = ""
        if fusion_figure.is_file():
            fusion_image_html = (
                '<article class="figure-card"><h3>Fusion chromosome map</h3>'
                '<a href="../10_rna_research_addons/fusions/fusion_chromosome_map.png">'
                '<img src="../10_rna_research_addons/fusions/fusion_chromosome_map.png" '
                'alt="RNA fusion candidate breakpoints mapped to chromosomes"></a>'
                '<p>Lines connect candidate breakpoints on their respective chromosomes; '
                'colors represent Arriba confidence/filter tiers. The table reports Arriba '
                'split-read and discordant-mate evidence. This is a research screen, not a '
                'clinical fusion call.</p></article>')
        fusion_html = (
            f'<section><h2>Fusion candidates ({fusion_count})</h2>'
            '<p>Evidence tiers preserve Arriba confidence and filter results. No independent '
            'minimum read threshold was invented: reported supporting reads and native Arriba '
            'filters are displayed transparently.</p>'
            f'<div class="figure-grid">{fusion_image_html}</div>'
            '<h3>Candidate evidence</h3>' + fusion_rows_html + '</section>')
        addon_downloads = [
            ("Enhanced per-cell QC statuses and reasons",
             "../10_rna_research_addons/qc_enhanced/per_cell_qc_status.tsv",
             "BaseJumper and enhanced status/reason calls with TIN robust-z values."),
            ("Positional duplication estimates",
             "../10_rna_research_addons/qc/duplication/positional_duplication_summary.tsv",
             "RSeQC mapped-position duplication estimate; not UMI molecule counting."),
            ("DEXSeq exon-bin count matrix",
             "../10_rna_research_addons/exon_counts/exon_bin_count_matrix.tsv",
             "Flattened exon-bin counts by sample; descriptive, no differential test."),
            ("SUPPA2 per-event PSI matrix",
             "../10_rna_research_addons/splicing/event_psi_by_sample.tsv",
             "Per-sample event PSI estimates; no differential-splicing test."),
            ("Arriba fusion evidence",
             "../10_rna_research_addons/fusions/fusion_evidence_summary.tsv",
             "Fusion breakpoints, native confidence/filter labels, and read support."),
        ]
        result_links.extend(
            (label, href, description) for label, href, description in addon_downloads
            if (OUT / href.removeprefix("../")).is_file())
    else:
        if addon_status:
            status_text = (
                f"Optional RNA research add-ons are not complete (status: "
                f"{addon_status.get('status', 'unknown')}). No enhanced QC, exon/splicing, "
                "or fusion conclusions are displayed.")
        elif addon_dir.exists():
            status_text = (
                "An add-on output directory exists without a successful run-status record. "
                "Treat it as an incomplete run; no enhanced QC, exon/splicing, or fusion "
                "conclusions are displayed.")
        else:
            status_text = (
                "Optional RNA research add-ons have not been run. Enhanced QC, exon/splicing, "
                "and fusion tabs will populate only after a separately authorized completed run.")
        enhanced_qc_html = f'<section><h2>Enhanced QC</h2><p>{html.escape(status_text)}</p></section>'
        splicing_html = f'<section><h2>Exon and splicing add-ons</h2><p>{html.escape(status_text)}</p></section>'
        fusion_html = f'<section><h2>Fusion candidates</h2><p>{html.escape(status_text)}</p></section>'
    link_rows = "".join(
        f'<tr><td><a href="{html.escape(href)}">{html.escape(label)}</a></td>'
        f'<td>{html.escape(description)}</td></tr>'
        for label, href, description in result_links)
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Full-depth RNA-seq cohort analysis</title>
<style>
:root{{--maroon:{COLORS['maroon']};--teal:{COLORS['teal']};--gold:{COLORS['gold']};--ink:{COLORS['ink']};--muted:{COLORS['muted']};--paper:#f3f5f7;--line:#dde3e8}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:15px/1.58 "Segoe UI",Arial,sans-serif}}
header{{background:#fff;border-bottom:4px solid var(--maroon);padding:22px max(24px,calc((100% - 1280px)/2));display:flex;align-items:center;gap:28px}}
header img{{width:235px;max-height:78px}}.brand-mark{{font-size:22px;font-weight:700;line-height:1.05;color:var(--maroon);border:2px solid var(--maroon);border-radius:8px;padding:12px 17px}}header .title{{border-left:1px solid var(--line);padding-left:24px}}h1{{font-size:27px;line-height:1.2;margin:0 0 5px;color:var(--maroon)}}header p{{margin:0;color:var(--muted)}}
main{{max-width:1280px;margin:26px auto;padding:0 22px}}section{{background:#fff;border:1px solid var(--line);border-radius:12px;padding:24px;margin:0 0 22px;box-shadow:0 2px 8px #13212c0a}}
h2{{font-size:21px;color:var(--maroon);margin:0 0 14px}}h3{{font-size:16px;margin:0 0 12px;color:var(--ink)}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(165px,1fr));gap:12px}}.stat{{background:#fafbfc;border:1px solid var(--line);border-radius:9px;padding:15px}}.stat b{{display:block;color:var(--maroon);font-size:25px}}.stat span{{color:var(--muted);font-size:13px}}
.notice{{border-left:4px solid var(--gold);background:#fff8e9;padding:14px 16px;border-radius:4px;margin:10px 0}}
.figure-grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,490px),1fr));gap:18px}}.figure-card{{border:1px solid var(--line);border-radius:9px;padding:16px;background:white;min-width:0}}.figure-card img{{width:100%;height:auto;border:1px solid #edf0f3;border-radius:5px}}.figure-card p{{font-size:13px;color:#4e5b65;margin:11px 0 0}}
table{{border-collapse:collapse;width:100%;font-size:14px}}th,td{{padding:10px 12px;text-align:left;border-bottom:1px solid var(--line);vertical-align:top}}th{{color:var(--maroon);background:#fafbfc}}a{{color:#145b75;text-decoration:none}}a:hover{{text-decoration:underline}}
.table-scroll{{overflow:auto;max-height:72vh;border:1px solid var(--line);border-radius:8px}}.metrics-table{{font-size:12px;white-space:nowrap}}.metrics-table th{{position:sticky;top:0;z-index:1}}details{{margin-top:16px}}summary{{cursor:pointer;color:#145b75;font-weight:600;padding:8px 0}}
.qc-badge{{display:inline-block;border-radius:12px;padding:2px 9px;font-weight:700}}.qc-pass{{background:#e4f3ef;color:#14695e}}.qc-borderline{{background:#fff2d4;color:#8c5b00}}.qc-fail{{background:#f8e5e5;color:#982c2c}}
.tab-nav{{display:flex;flex-wrap:wrap;gap:8px;margin:0 0 20px}}.tab-button{{border:1px solid var(--line);border-radius:20px;background:#fff;color:var(--ink);padding:9px 15px;font:600 13px "Segoe UI",Arial,sans-serif;cursor:pointer}}.tab-button[aria-selected="true"]{{background:var(--maroon);border-color:var(--maroon);color:#fff}}.tab-panel{{display:none}}.tab-panel.active{{display:block}}
.small{{font-size:13px;color:var(--muted)}}.attribution{{color:#4b5862;font-size:13px;margin-top:12px}}footer{{text-align:center;padding:8px 20px 32px;color:var(--muted);font-size:12px}}
@media(max-width:650px){{header{{align-items:flex-start;flex-direction:column}}header .title{{border:0;border-top:1px solid var(--line);padding:14px 0 0}}main{{padding:0 12px}}section{{padding:17px}}}}
</style></head><body>
<header>{logo_markup}
<div class="title"><h1>Full-depth RNA-seq cohort analysis</h1>
<p>{html.escape(condition)} cohort &nbsp;|&nbsp; {sample_count} biological replicates &nbsp;|&nbsp; GRCh38 / Ensembl 104</p>
<p class="attribution">Method/design by Zalman Vaksman PhD and Department of Computational Biology</p></div></header>
<main>
<nav class="tab-nav" role="tablist" aria-label="RNA analysis output types">
<button class="tab-button" id="tab-overview-button" role="tab" aria-controls="tab-overview" aria-selected="true">Overview</button>
<button class="tab-button" id="tab-expression-button" role="tab" aria-controls="tab-expression" aria-selected="false">Expression</button>
<button class="tab-button" id="tab-baseqc-button" role="tab" aria-controls="tab-baseqc" aria-selected="false">BaseJumper QC</button>
<button class="tab-button" id="tab-enhanced-button" role="tab" aria-controls="tab-enhanced" aria-selected="false">Enhanced QC</button>
<button class="tab-button" id="tab-splicing-button" role="tab" aria-controls="tab-splicing" aria-selected="false">Exon &amp; splicing</button>
<button class="tab-button" id="tab-fusions-button" role="tab" aria-controls="tab-fusions" aria-selected="false">Fusions</button>
<button class="tab-button" id="tab-data-button" role="tab" aria-controls="tab-data" aria-selected="false">Data &amp; methods</button>
</nav>
<div class="tab-panel active" id="tab-overview" role="tabpanel" aria-labelledby="tab-overview-button">
<section><h2>Study overview</h2><div class="cards">
<div class="stat"><b>{sample_count}</b><span>biological replicates</span></div>
<div class="stat"><b>{len(metrics):,}</b><span>genes in count matrix</span></div>
<div class="stat"><b>{len(metrics.loc[metrics.mean_raw_count.gt(0)]):,}</b><span>genes detected in ≥1 replicate</span></div>
<div class="stat"><b>{qc_counts.get('PASS',0)} / {qc_counts.get('Borderline',0)}</b><span>PASS / Borderline source QC labels</span></div>
<div class="stat"><b>{float(size_factors.min()):.2f}–{float(size_factors.max()):.2f}</b><span>median-ratio size-factor range</span></div>
<div class="stat"><b>{len(salmon_summary):,}</b><span>Salmon transcript quantifications audited</span></div>
</div><p>All samples are recorded as one condition ({html.escape(str(qc.condition.iloc[0]))}) and treated as independent biological replicates. The workflow therefore performs cohort characterization, not a condition contrast.</p>
<div class="notice"><b>Interpretation:</b> {fdr_wording} The pathway analysis is an exploratory enrichment of highly variable expressed genes across replicates; it is not evidence of activation or differential regulation. Differential expression and differential splicing were not run because there is only one condition.</div>
<p><b>Borderline QC labels:</b> {html.escape(borderline)}. All {sample_count} samples are retained; labels are reported, not used as exclusions.</p></section>
</div>
<div class="tab-panel" id="tab-expression" role="tabpanel" aria-labelledby="tab-expression-button">
<section><h2>Figures and interpretation</h2><div class="figure-grid">{''.join(figures_html)}</div></section>
</div>
<div class="tab-panel" id="tab-baseqc" role="tabpanel" aria-labelledby="tab-baseqc-button">
<section><h2>BaseJumper QC report</h2>
<p>These are the original full-depth BaseJumper QC plots and per-sample summary outputs. The dashboard includes every field from the BaseJumper full-metrics table, with the source report values preserved.</p>
<div class="figure-grid">{''.join(basejumper_figures_html)}</div>
<h3>Mitochondrial reads and QC by sample</h3>
<p>MT reads are counts assigned by HTSeq to genes whose symbol begins <code>MT-</code>. The percentage is MT-assigned gene counts divided by all assigned gene counts, not the fraction of all sequenced reads or all STAR-aligned reads.</p>
{mt_summary_html}
<details><summary>Show all {len(all_metric_columns)} BaseJumper metrics for all {len(basejumper_metrics)} samples</summary>
<p class="small">The complete source table is also available as <a href="../../stage7_full_reports/consensus/rnaqc_all_metrics.tsv">rnaqc_all_metrics.tsv</a>. Fractions are displayed as percentages; STAR fields already expressed in percent are not rescaled.</p>
{all_metrics_html}</details>
<p class="small">BaseJumper's original <a href="../../stage7_full_reports/multiqc/multiqc_report.html">MultiQC report</a>, <a href="../../stage7_full_reports/consensus/composition_rnaqc.pdf">composition PDF</a>, selected-metrics table, and Qualimap raw metric table are linked under Downloadable results.</p>
</section>
</div>
<div class="tab-panel" id="tab-enhanced" role="tabpanel" aria-labelledby="tab-enhanced-button">{enhanced_qc_html}</div>
<div class="tab-panel" id="tab-splicing" role="tabpanel" aria-labelledby="tab-splicing-button">{splicing_html}</div>
<div class="tab-panel" id="tab-fusions" role="tabpanel" aria-labelledby="tab-fusions-button">{fusion_html}</div>
<div class="tab-panel" id="tab-data" role="tabpanel" aria-labelledby="tab-data-button">
<section><h2>Downloadable results</h2><table><thead><tr><th>Output</th><th>Contents and intended use</th></tr></thead><tbody>{link_rows}</tbody></table>
<p class="small">Large final alignments are linked under <a href="../08_large_files/bam/">08_large_files/bam</a> to avoid duplicating hundreds of gigabytes. Their source files remain in the audited stage directory.</p></section>
<section><h2>Methods and limitations</h2>
<ul><li>Exploratory gene-expression analyses use the full-depth HTSeq raw-count matrix. Size factors use a median-ratio, positive-count method; log2 normalized counts are for visualization, not a fitted count model.</li>
<li>PCA uses 500 most variable genes after mean raw-count filtering. Sample correlation is Spearman correlation across count-matrix genes.</li>
<li>Reactome over-representation compares the 500 most variable expressed genes with a background of genes having mean raw count ≥10; Benjamini–Hochberg FDR is reported. Per-sample pathway scores are mean gene-wise z-scores.</li>
<li>Salmon transcript abundance and STAR splice-junction support are descriptive. Optional SUPPA2 PSI estimates and DEXSeq exon-bin counts are presented as research summaries only; no differential transcript-usage, exon-usage, or splicing test was fit.</li>
<li>Enhanced QC combines the source BaseJumper status with cohort-relative low-TIN outlier flags. Positional duplication is a mapping-position estimate, not UMI molecule counting; neither method automatically excludes a sample.</li>
<li>Optional fusion candidates use a separate STAR chimeric alignment and Arriba's confidence/filter labels and supporting-read fields. Candidate calls require independent review and are not clinical findings.</li>
<li>The cohort is bulk paired-end RNA-seq, not droplet/UMI single-cell RNA-seq. This workflow records fastp's estimated sequence-duplication rate but does not remove PCR duplicates; no UMI correction or coordinate-based deduplication was performed. In UMI-based single-cell workflows, PCR copies are typically collapsed as one molecule using cell barcode plus UMI (and gene/feature assignment); identical read sequences alone are not sufficient evidence of a duplicate.</li>
<li>The full-depth alignment audit mode and sample-level scope are recorded in the workflow audit reports.</li>
{("<li>The St. Jude logo is used for identification only. This is a research-workflow report, not an official St. Jude publication or endorsement.</li>" if source_logo.is_file() else "")}</ul>
<p class="small">Analysis generated {date.today().isoformat()} · Reactome release and archive checksum are recorded in 09_references/README.txt · Alignment audit mode: {html.escape(str(alignment_mode))} · BaseJumper per-sample QC metrics: {len(basejumper_by_sample)}.</p>
</section></div></main><footer>ResolveOME RNA workflow · Full-depth outputs and provenance are organized alongside this dashboard.</footer>
<script>
document.querySelectorAll(".tab-button").forEach((button) => {{
  button.addEventListener("click", () => {{
    document.querySelectorAll(".tab-button").forEach((item) => item.setAttribute("aria-selected", "false"));
    document.querySelectorAll(".tab-panel").forEach((panel) => panel.classList.remove("active"));
    button.setAttribute("aria-selected", "true");
    document.getElementById(button.getAttribute("aria-controls")).classList.add("active");
  }});
}});
</script></body></html>"""
    (dashboard_dir / "index.html").write_text(page, encoding="utf-8")


def run() -> Path:
    ensure_output_tree()
    counts_with_names, metadata, step7, audits = validate_prerequisites()
    parity = basejumper_settings_audit(audits)
    provenance = OUT / "07_qc_provenance/source_audits"
    provenance.mkdir(parents=True, exist_ok=True)
    audit_sources = {
        "step4_full_output_audit_quick.json": RNAR / "step4_full_output_audit_quick.json",
        "step5_full_output_audit.json": RNAR / "step5_full_output_audit.json",
        "step6_full_output_audit.json": RNAR / "step6_full_output_audit.json",
        "step7_full_output_audit.json": RNAR / "step7_full_output_audit.json",
        "salmon_output_audit.json": OUT / "07_qc_provenance/salmon_output_audit.json",
    }
    for name, source in audit_sources.items():
        ensure_reference_link(source, provenance / name)
    ensure_reference_link(
        BASEJUMPER_METRICS, provenance / "basejumper_full_metrics.tsv")
    (OUT / "07_qc_provenance/basejumper_settings_parity.json").write_text(
        json.dumps(parity, indent=2) + "\n")
    gene_labels = split_gene_labels(counts_with_names.gene_id_gene_symbol)
    samples = metadata.sample_id.tolist()
    raw = counts_with_names[samples].copy()
    raw.index = gene_labels.gene_id
    raw.index.name = "gene_id"
    if raw.index.duplicated().any():
        raise ValueError("Duplicate Ensembl gene IDs after parsing count matrix")
    size_factors = median_ratio_size_factors(raw)
    normalized = raw.div(size_factors, axis=1)
    log_counts = np.log2(normalized + 1)
    metrics = gene_labels.copy()
    metrics["mean_raw_count"] = raw.mean(axis=1).values
    metrics["mean_normalized_count"] = normalized.mean(axis=1).values
    metrics["variance_log2"] = log_counts.var(axis=1, ddof=0).values
    metrics["sd_log2"] = log_counts.std(axis=1, ddof=0).values
    metrics["detection_fraction"] = raw.gt(0).mean(axis=1).values
    metrics["cv_normalized_count"] = (
        normalized.std(axis=1, ddof=0) /
        normalized.mean(axis=1).replace(0, np.nan)).values
    metrics = metrics.sort_values("mean_normalized_count", ascending=False)

    sample_qc = pd.DataFrame(step7["samples"]).rename(columns={"sample": "sample_id"})
    if set(sample_qc.sample_id) != set(samples):
        raise ValueError("Final source QC score table and sample sheet differ")
    sample_qc = sample_qc.merge(metadata[["sample_id", "group"]], on="sample_id",
                                validate="one_to_one")
    sample_qc["sample_short"] = sample_qc.sample_id.map(sample_short)
    sample_qc = sample_qc.rename(columns={"group": "condition"})
    basejumper_metrics = load_basejumper_metrics(samples)

    save_tsv(raw.reset_index().merge(gene_labels.assign(gene_id=gene_labels.gene_id),
                                    on="gene_id", how="left").drop_duplicates("gene_id"),
             OUT / "01_expression/raw_gene_counts.tsv")
    save_tsv(normalized.rename_axis("gene_id").reset_index(),
             OUT / "01_expression/normalized_gene_counts.tsv")
    save_tsv(log_counts.rename_axis("gene_id").reset_index(),
             OUT / "01_expression/log2_normalized_gene_counts.tsv")
    save_tsv(metrics, OUT / "01_expression/gene_expression_summary.tsv")
    save_tsv(pd.DataFrame({"sample_id": samples, "sample_short": [sample_short(s) for s in samples],
                           "size_factor": size_factors.reindex(samples).values}),
             OUT / "01_expression/size_factors.tsv")
    save_tsv(sample_qc, OUT / "07_qc_provenance/sample_qc_status.tsv")
    save_tsv(metadata, OUT / "07_qc_provenance/sample_metadata.tsv")

    captions = []
    set_plot_style()
    score_df, corr_df = plot_expression_figures(metrics, log_counts, samples,
                                                 sample_qc, captions)
    save_tsv(score_df, OUT / "01_expression/pca_coordinates.tsv")
    save_tsv(corr_df.reset_index(), OUT / "01_expression/spearman_correlation.tsv")
    plot_abundance(metrics, captions)

    gene_sets = load_reactome_sets(RESOURCES / "ReactomePathways.gmt.zip")
    enrichment, pathway_scores = run_pathway_enrichment(
        metrics, gene_sets, samples, log_counts)
    save_tsv(enrichment, OUT / "04_pathways/reactome_variable_gene_enrichment.tsv")
    save_tsv(pathway_scores.reset_index(),
             OUT / "04_pathways/reactome_descriptive_sample_scores.tsv")
    plot_pathways(enrichment, pathway_scores, samples, captions)

    tpm, transcript_reads_df, gene_tpm, salmon_mapping = load_salmon_matrices(samples)
    transcript_summary = plot_transcripts(tpm, captions)
    save_tsv(tpm.rename_axis("transcript_id").reset_index(),
             OUT / "02_transcript/transcript_tpm.tsv")
    save_tsv(transcript_reads_df.rename_axis("transcript_id").reset_index(),
             OUT / "02_transcript/transcript_numreads.tsv")
    save_tsv(gene_tpm.rename_axis("gene_id").reset_index(),
             OUT / "02_transcript/gene_aggregated_tpm.tsv")
    save_tsv(gene_tpm.reset_index().merge(
        gene_labels.drop_duplicates("gene_id"), on="gene_id", how="left"),
        OUT / "02_transcript/gene_aggregated_tpm_with_symbols.tsv")
    save_tsv(transcript_summary.sort_values("variance_TPM", ascending=False),
             OUT / "02_transcript/transcript_variability.tsv")
    save_tsv(salmon_mapping, OUT / "02_transcript/transcript_mapping_summary.tsv")

    junction_counts, junction_summary = load_star_junctions(samples)
    save_tsv(junction_counts.reset_index(),
             OUT / "03_splicing/junction_unique_read_counts.tsv")
    junction_metrics = plot_junctions(junction_counts, captions)
    save_tsv(junction_metrics, OUT / "03_splicing/recurrent_junction_variability.tsv")
    save_tsv(junction_summary, OUT / "03_splicing/sample_junction_summary.tsv")

    gene_z = log_counts.sub(log_counts.mean(axis=1), axis=0).div(
        log_counts.std(axis=1, ddof=0).replace(0, np.nan), axis=0).fillna(0)
    enriched_rows = []
    expressed_gene_ids = metrics.loc[metrics.mean_raw_count.ge(10), "gene_id"]
    for col in samples:
        local_z = gene_z.loc[expressed_gene_ids, col]
        top = local_z.nlargest(10)
        for rank, (gene_id, z_value) in enumerate(top.items(), start=1):
            row = metrics.loc[metrics.gene_id.eq(gene_id)].iloc[0]
            enriched_rows.append({
                "sample_id": col, "sample_short": sample_short(col), "rank": rank,
                "gene_id": gene_id, "gene_symbol": row.gene_symbol,
                "log2_normalized_count": log_counts.loc[gene_id, col],
                "within_cohort_gene_z_score": z_value,
                "mean_raw_count": row.mean_raw_count,
                "interpretation": "descriptive replicate-enriched expression; not a validated marker",
            })
    save_tsv(pd.DataFrame(enriched_rows),
             OUT / "05_cell_line_profiles/sample_enriched_genes.tsv")
    save_tsv(metrics.nlargest(100, "mean_normalized_count"),
             OUT / "05_cell_line_profiles/top_cohort_abundance.tsv")

    analysis_summary = {
        "status": "PASS",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "samples": len(samples),
        "design": {
            "condition_count": int(metadata.group.nunique()),
            "condition": metadata.group.iloc[0],
            "replicate_interpretation": f"{len(samples)} independent biological replicates of one condition",
            "differential_expression": "NOT RUN: no condition contrast exists",
            "differential_splicing": "NOT RUN: no condition contrast exists",
        },
        "gene_matrix": {
            "genes": len(raw),
            "detected_in_at_least_one_sample": int(raw.gt(0).any(axis=1).sum()),
            "all_values_match_audited_HTSeq_matrix": bool(np.array_equal(
                raw[samples].to_numpy(), counts_with_names[samples].to_numpy())),
        },
        "qc_labels": sample_qc.status.value_counts().to_dict(),
        "borderline_samples": sample_qc.loc[
            sample_qc.status.eq("Borderline"), "sample_id"].tolist(),
        "alignment_audit": {
            "status": audits["step4"]["status"],
            "mode": audits["step4"].get("audit_mode"),
            "deep_checked_cells": audits["step4"].get("deep_checked_cells", []),
        },
        "basejumper_settings_parity": parity["status"],
        "pathway_analysis": {
            "database": "ReactomePathways.gmt",
            "source_url": "https://reactome.org/download/current/ReactomePathways.gmt.zip",
            "archive_sha256": REACTOME_SHA256,
            "method": "hypergeometric over-representation of 500 most variable expressed genes",
            "background": "genes with mean raw HTSeq count >=10",
            "fdr_bh_lt_0_05": int((enrichment.fdr_bh < 0.05).sum()),
            "interpretation": "descriptive variability; not differential pathway activation",
        },
        "transcript_quantification": {
            "tool": "Salmon 1.6.0",
            "shared_index_version": "1.6.0",
            "interpretation": "descriptive transcript abundance; no differential usage test",
        },
        "splicing": {
            "source": "STAR SJ.out.tab unique read support",
            "junctions": int(len(junction_counts)),
            "interpretation": "descriptive junction support; not PSI or differential splicing",
        },
        "figure_files": [row["figure"] for row in captions],
    }
    bam_link = OUT / "08_large_files/bam"
    source_bam = RNAR / "stage4_full_alignment/bam"
    target = Path("../../stage4_full_alignment/bam")
    if bam_link.is_symlink():
        if bam_link.resolve(strict=False) != (bam_link.parent / target).resolve(strict=False):
            raise FileExistsError(f"Existing BAM package link points elsewhere: {bam_link}")
    elif bam_link.exists():
        raise FileExistsError(f"Refusing to replace existing alignment package path: {bam_link}")
    else:
        bam_link.symlink_to(target)
    try:
        if not bam_link.resolve().is_dir() or not source_bam.is_dir():
            raise RuntimeError("Full-depth BAM package link does not resolve")
    except OSError as exc:
        raise RuntimeError(f"Could not verify BAM package link: {exc}") from exc
    build_dashboard(metrics, sample_qc, size_factors, enrichment, salmon_mapping,
                    junction_summary, audits, captions, basejumper_metrics)
    (OUT / "07_qc_provenance/analysis_summary.json").write_text(
        json.dumps(analysis_summary, indent=2) + "\n")
    (OUT / "07_qc_provenance/analysis_summary.txt").write_text(
        "Full-depth cohort characterization completed.\n"
        f"All {len(samples)} samples are one condition with biological replication. "
        "Differential expression and differential splicing were not run.\n"
        f"Final QC labels: {sample_qc.status.value_counts().to_dict()}\n"
        f"Reactome variability-enrichment FDR<0.05 terms: "
        f"{analysis_summary['pathway_analysis']['fdr_bh_lt_0_05']}\n"
        f"STAR junction rows summarized: {len(junction_counts)}\n"
        f"Full-depth quick alignment audit: {alignment_mode_for(audits)}\n")
    write_package_readme(analysis_summary)
    return OUT / "00_dashboard/index.html"


def alignment_mode_for(audits: dict) -> str:
    return str(audits["step4"].get("audit_mode", "not recorded"))


def refresh_dashboard() -> Path:
    qc = pd.read_csv(OUT / "07_qc_provenance/sample_qc_status.tsv", sep="\t")
    samples = qc.sample_id.tolist()
    metrics = pd.read_csv(OUT / "01_expression/gene_expression_summary.tsv", sep="\t")
    size_factor_table = pd.read_csv(OUT / "01_expression/size_factors.tsv", sep="\t")
    size_factors = pd.Series(
        size_factor_table.size_factor.to_numpy(),
        index=size_factor_table.sample_id.to_numpy())
    enrichment = pd.read_csv(
        OUT / "04_pathways/reactome_variable_gene_enrichment.tsv", sep="\t")
    salmon_summary = pd.read_csv(
        OUT / "02_transcript/transcript_mapping_summary.tsv", sep="\t")
    junction_summary = pd.read_csv(
        OUT / "03_splicing/sample_junction_summary.tsv", sep="\t")
    audit_paths = {
        "step4": RNAR / "step4_full_output_audit_quick.json",
        "step5": RNAR / "step5_full_output_audit.json",
        "step6": RNAR / "step6_full_output_audit.json",
        "step7": RNAR / "step7_full_output_audit.json",
        "salmon": OUT / "07_qc_provenance/salmon_output_audit.json",
    }
    audits = {name: json.loads(path.read_text()) for name, path in audit_paths.items()}
    captions = pd.read_csv(
        OUT / "07_qc_provenance/figure_captions.tsv", sep="\t").to_dict("records")
    basejumper_metrics = load_basejumper_metrics(samples)
    ensure_reference_link(
        BASEJUMPER_METRICS, OUT / "07_qc_provenance/basejumper_full_metrics.tsv")
    build_dashboard(metrics, qc, size_factors, enrichment, salmon_summary,
                    junction_summary, audits, captions, basejumper_metrics)
    return OUT / "00_dashboard/index.html"


def write_package_readme(summary: dict) -> None:
    content = f"""# Full-depth RNA-seq cohort results

This package contains descriptive analyses of {summary['samples']} {summary['design']['condition']}
biological replicates from the audited full-depth ResolveOME RNA workflow. All samples are one condition;
there is no differential-expression, differential-splicing, or differential
transcript-usage contrast.

Open `00_dashboard/index.html` for the figure-led report and download links.

## Folder guide

- `00_dashboard/` — local HTML dashboard with labeled figures and interpretations.
- `01_expression/` — audited raw HTSeq counts, normalized matrices, PCA and gene summaries.
- `02_transcript/` — Salmon 1.6.0 transcript TPM/NumReads and gene-aggregated TPM.
- `03_splicing/` — STAR unique junction support and descriptive variability.
- `04_pathways/` — Reactome variability enrichment and within-cohort expression scores.
- `05_cell_line_profiles/` — descriptive replicate-enriched and cohort-abundance tables.
- `06_figures/` — presentation-quality PNG figures used by the dashboard.
- `07_qc_provenance/` — source/runtime BaseJumper settings parity and audit references,
  methods, sample metadata/QC and analysis summary.
- `08_large_files/bam/` — symbolic link to the existing full-depth BAM folder; no BAMs copied.
- `09_references/` — link to the Reactome release archive; institution marks are optional and supplied separately.

## Interpretation limits

All samples are retained, including any existing Borderline QC labels.
Pathway enrichment is exploratory over-representation among variable genes;
it is not pathway activation. Per-sample scores are mean gene-wise z-scores,
not GSVA/ssGSEA. Transcript abundance and splice-junction support are
descriptive, not formal differential tests. The input is bulk RNA-seq, not
droplet/UMI single-cell RNA-seq.

The full-depth alignment audit is quick mode. Its exact scope and deep-checked
cells are linked from `07_qc_provenance/source_audits/` to the central
`../audits/full_depth/` files.
`07_qc_provenance/basejumper_settings_parity.json` compares the vendored source
with the adapted runtime's scientific flags. It is settings-level validation,
not a separate stock BaseJumper execution or a bitwise output-equivalence claim.

The report credits “Method/design by Zalman Vaksman PhD and Department of
Computational Biology.” Any institution logo is used for identification only;
this is not an official institution publication or endorsement.

Analysis status: {summary['status']}<br>
Genes in audited matrix: {summary['gene_matrix']['genes']:,}<br>
Samples: {summary['samples']}
"""
    (OUT / "README.md").write_text(content, encoding="utf-8")


def main() -> None:
    import argparse

    global BASEJUMPER_COMMIT_FILE, BASEJUMPER_METRICS, BASEJUMPER_REPORTS
    global BASEJUMPER_SOURCE, MATRIX, OUT, RESOURCES, RNAR, SAMPLES, STAR, TX2GENE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=RNAR)
    parser.add_argument("--cells", type=Path, default=SAMPLES)
    parser.add_argument("--tx2gene", type=Path, default=TX2GENE)
    parser.add_argument("--analysis-resources", type=Path, default=RESOURCES)
    parser.add_argument("--basej-source", type=Path, default=BASEJUMPER_SOURCE)
    parser.add_argument("--basej-commit-file", type=Path, default=BASEJUMPER_COMMIT_FILE)
    parser.add_argument(
        "--refresh-dashboard", action="store_true",
        help="refresh only dashboard content from completed results and BaseJumper QC outputs")
    args = parser.parse_args()
    RNAR = args.results_dir.resolve()
    OUT = RNAR / "final_results"
    MATRIX = RNAR / "stage6_full_matrix/matrix/matrix_gene_counts_starhtseq.txt"
    SAMPLES = args.cells.resolve()
    STAR = RNAR / "stage4_full_alignment/star"
    BASEJUMPER_METRICS = RNAR / "stage7_full_reports/consensus/rnaqc_all_metrics.tsv"
    BASEJUMPER_REPORTS = RNAR / "stage7_full_reports"
    TX2GENE = args.tx2gene.resolve()
    RESOURCES = args.analysis_resources.resolve()
    BASEJUMPER_SOURCE = args.basej_source.resolve()
    BASEJUMPER_COMMIT_FILE = args.basej_commit_file.resolve()
    if args.refresh_dashboard:
        dashboard = refresh_dashboard()
        print(f"Dashboard refreshed: {dashboard}")
        return

    status_path = OUT / "07_qc_provenance/analysis_run_status.json"
    summary_path = OUT / "07_qc_provenance/analysis_summary.json"
    if status_path.exists() or summary_path.exists():
        raise FileExistsError(
            "Analysis run state already exists; inspect the existing package before rerunning")
    ensure_output_tree()
    status_path.write_text(json.dumps({
        "status": "RUNNING",
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }, indent=2) + "\n")
    try:
        dashboard = run()
    except Exception as exc:
        status_path.write_text(json.dumps({
            "status": "FAILED",
            "failed_utc": datetime.now(timezone.utc).isoformat(),
            "error_type": type(exc).__name__,
            "error": str(exc),
        }, indent=2) + "\n")
        raise
    status_path.write_text(json.dumps({
        "status": "PASS",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
    }, indent=2) + "\n")
    print(f"Full-depth analysis package complete: {OUT}")
    print(f"Dashboard: {dashboard}")


if __name__ == "__main__":
    main()
