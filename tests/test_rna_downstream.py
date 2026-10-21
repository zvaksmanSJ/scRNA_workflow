import importlib
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "extra_methods"))
import numpy as np
import pandas as pd

import run_rna_downstream as analysis


class RNADownstreamTests(unittest.TestCase):
    def test_results_root_can_be_overridden_for_isolated_wdl_runs(self):
        with tempfile.TemporaryDirectory() as temporary:
            try:
                with patch.dict(os.environ, {"RESOLVEOME_RNA_RESULTS": temporary}):
                    importlib.reload(analysis)
                    self.assertEqual(analysis.RNAR, Path(temporary).resolve())
                    self.assertEqual(analysis.OUT, Path(temporary).resolve() / "final_results")
            finally:
                importlib.reload(analysis)

    def test_size_factors_are_positive_and_geometric_mean_one(self):
        counts = pd.DataFrame({
            "a": [100, 50, 0, 80],
            "b": [200, 25, 10, 40],
            "c": [50, 100, 20, 160],
        })
        factors = analysis.median_ratio_size_factors(counts)
        self.assertTrue((factors > 0).all())
        self.assertAlmostEqual(float(np.exp(np.log(factors).mean())), 1.0)

    def test_benjamini_hochberg_adjustment(self):
        adjusted = analysis.benjamini_hochberg([0.01, 0.04, 0.03, 0.002])
        np.testing.assert_allclose(adjusted, [0.02, 0.04, 0.04, 0.008])

    def test_audit_sample_ids_supports_rows_and_string_lists(self):
        expected = {"cell_a", "cell_b"}
        self.assertEqual(
            analysis.audit_sample_ids({"samples": [{"sample": "cell_a"},
                                                   {"sample_id": "cell_b"}]}, "step"),
            expected)
        self.assertEqual(
            analysis.audit_sample_ids({"samples": ["cell_a", "cell_b"]}, "step"),
            expected)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            analysis.audit_sample_ids({"samples": ["cell_a", "cell_a"]}, "step")

    def test_pathway_scores_align_gene_symbols_to_gene_ids(self):
        gene_ids = [f"ENSG{i:06d}" for i in range(1200)]
        metrics = pd.DataFrame({
            "gene_id": gene_ids,
            "gene_symbol": [f"G{i}" for i in range(1200)],
            "mean_raw_count": 20,
            "variance_log2": np.arange(1200, 0, -1),
        })
        values = np.column_stack((np.arange(1200), np.arange(1200)[::-1]))
        log_counts = pd.DataFrame(values, index=gene_ids, columns=["cell_a", "cell_b"])
        gene_sets = {"Test pathway [R-HSA-1]": {f"G{i}" for i in range(30)}}
        enrichment, scores = analysis.run_pathway_enrichment(
            metrics, gene_sets, ["cell_a", "cell_b"], log_counts)
        self.assertEqual(len(enrichment), 1)
        self.assertEqual(scores.loc["Test pathway [R-HSA-1]", "genes_scored"], 30)
        self.assertTrue(np.isfinite(scores[["cell_a", "cell_b"]].to_numpy()).all())

    def test_pathway_heatmap_uses_only_sample_columns(self):
        enrichment = pd.DataFrame({
            "pathway": ["A [R-HSA-1]"],
            "variable_gene_ratio": [0.1],
            "variable_gene_overlap": [20],
            "fdr_bh": [0.03],
        })
        scores = pd.DataFrame({
            "cell_a": [0.5],
            "cell_b": [-0.5],
            "genes_scored": [100],
        }, index=["A [R-HSA-1]"])
        figures = []

        def capture(fig, _name, _captions):
            figures.append(fig)

        with patch.object(analysis, "save_figure", side_effect=capture):
            analysis.plot_pathways(enrichment, scores, ["cell_a", "cell_b"], [])
        try:
            self.assertEqual(figures[1].axes[0].images[0].get_array().shape[1], 2)
        finally:
            for figure in figures:
                analysis.plt.close(figure)

    def test_salmon_mapping_is_vectorized_and_excludes_unmapped_targets(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            quant_root = root / "02_transcript/quant_by_sample"
            samples = ["cell_a", "cell_b"]
            quant_rows = []
            mapping_rows = []
            for index in range(1000):
                transcript = f"tx{index}.1"
                mapping_rows.append(
                    [transcript, "protein_coding", f"ENSG{index:06d}",
                     "protein_coding", f"G{index}"])
                quant_rows.append([transcript, 100, 90, 1.0, 10.0])
            quant_rows.append(["unmapped_decoy", 100, 90, 0.0, 100.0])
            tx2gene = root / "tx2gene.tsv"
            pd.DataFrame(mapping_rows).to_csv(
                tx2gene, sep="\t", header=False, index=False)
            for sample in samples:
                sample_dir = quant_root / sample
                sample_dir.mkdir(parents=True)
                pd.DataFrame(
                    quant_rows,
                    columns=["Name", "Length", "EffectiveLength", "TPM", "NumReads"],
                ).to_csv(sample_dir / "quant.sf", sep="\t", index=False)
            with patch.object(analysis, "OUT", root), \
                    patch.object(analysis, "TX2GENE", tx2gene):
                tpm, reads, gene_tpm, summary = analysis.load_salmon_matrices(samples)
            self.assertNotIn("unmapped_decoy", tpm.index)
            self.assertEqual(len(tpm), 1000)
            self.assertEqual(len(reads), 1000)
            self.assertEqual(gene_tpm.loc["ENSG000000", "cell_a"], 1.0)
            self.assertAlmostEqual(summary.fraction_estimated_reads_gene_mapped.iloc[0],
                                   10000 / 10100)

    def test_star_junction_strand_zero_is_not_labeled_minus(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samples = ["cell_a", "cell_b"]
            rows = [
                ["chr1", 10, 20, 0, 1, 0, 5, 0, 10],
                ["chr1", 30, 40, 1, 1, 1, 8, 0, 12],
                ["chr1", 50, 60, 2, 1, 1, 9, 0, 15],
            ]
            for sample in samples:
                pd.DataFrame(rows).to_csv(
                    root / f"{sample}_SJ.out.tab", sep="\t", header=False, index=False)
            with patch.object(analysis, "STAR", root):
                counts, _summary = analysis.load_star_junctions(samples)
            self.assertIn("chr1:10-20:?", counts.index)
            self.assertIn("chr1:30-40:+", counts.index)
            self.assertIn("chr1:50-60:-", counts.index)

    def test_final_results_readme_documents_one_condition_limit(self):
        with tempfile.TemporaryDirectory() as temporary:
            with patch.object(analysis, "OUT", Path(temporary)):
                analysis.write_package_readme({
                    "status": "PASS",
                    "gene_matrix": {"genes": 40_019},
                    "samples": 16,
                    "design": {"condition": "SY5Y"},
                })
                content = (Path(temporary) / "README.md").read_text()
            self.assertIn("no differential-expression", content)
            self.assertIn("08_large_files/bam", content)
            self.assertIn("not an official institution publication or endorsement", content)

    def test_status_table_color_codes_each_qc_call_and_escapes_reasons(self):
        table = analysis.render_status_table(pd.DataFrame({
            "sample_id": ["cell_a", "cell_b", "cell_c"],
            "basejumper_status": ["PASS", "Borderline", "FAIL"],
            "basejumper_reasons": ["", "low <70%", "missing & invalid"],
        }), ["sample_id", "basejumper_status", "basejumper_reasons"])
        self.assertIn("qc-pass", table)
        self.assertIn("qc-borderline", table)
        self.assertIn("qc-fail", table)
        self.assertIn("low &lt;70%", table)
        self.assertIn("missing &amp; invalid", table)

    def test_basejumper_percent_metrics_preserve_source_units(self):
        self.assertEqual(
            analysis.format_basejumper_value("mt_pct", 0.0617467), "6.17%")
        self.assertEqual(
            analysis.format_basejumper_value("star_uniquely_mapped_pct", 77.35),
            "77.35%")
        self.assertEqual(
            analysis.format_basejumper_value("fastp_duplication_rate", 0.557441),
            "55.74%")
        frame = pd.DataFrame({
            "biosample": ["cell_a"],
            "mt_gene_counts": [1234],
            "mt_pct": [0.0617],
        })
        table = analysis.render_basejumper_table(frame, list(frame.columns))
        self.assertIn("1,234", table)
        self.assertIn("6.17%", table)
        self.assertIn('title="mt_gene_counts"', table)

    def test_basejumper_settings_parity_audit_checks_runtime_copy(self):
        source_text = """--detect_adapter_for_pe
--twopassMode Basic
--outSAMtype BAM SortedByCoordinate
--outFilterMultimapNmax 10
--outSAMprimaryFlag OneBestScore
--chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1
samtools view -F 256 -h
-f bam -r pos -s no -t exon -i gene_id --additional-attr=gene_name -m union
--secondary-alignments=ignore
"""
        runtime_text = source_text.replace(
            "-s no", "-s ${params.strandedness}")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "workflow_RNA/vendor/basej-rnaqc"
            source.mkdir(parents=True)
            source_main = source / "main.nf"
            source_main.write_text(source_text)
            (root / "workflow_RNA/vendor/UPSTREAM_COMMIT").write_text("test-commit\n")
            runtime = root / "RNA_results/stage4_full_alignment/runtime"
            runtime.mkdir(parents=True)
            runtime_main = runtime / "main.nf"
            runtime_main.write_text(runtime_text)
            (runtime / "runtime_manifest.json").write_text(json.dumps({
                "sha256": {"main.nf": hashlib.sha256(source_main.read_bytes()).hexdigest()},
                "adapter_sha256": "test-adapter-hash",
            }))
            with patch.object(analysis, "RNAW", root / "workflow_RNA"), \
                    patch.object(analysis, "BASEJUMPER_SOURCE", source), \
                    patch.object(analysis, "BASEJUMPER_COMMIT_FILE",
                                 root / "workflow_RNA/vendor/UPSTREAM_COMMIT"), \
                    patch.object(analysis, "RNAR", root / "RNA_results"):
                result = analysis.basejumper_settings_audit({
                    "step4": {"depth": "full"},
                    "step6": {"all_matrix_values_match_original_HTSeq": True},
                })
        self.assertEqual(result["status"], "PASS_SETTINGS_LEVEL_WITH_LIMITATIONS")
        self.assertTrue(all(row["status"] == "PASS"
                            for row in result["scientific_setting_checks"]))
        self.assertIn("No separate untouched stock", result["not_claimed"])


if __name__ == "__main__":
    unittest.main()
