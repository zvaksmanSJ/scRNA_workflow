import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "extra_methods"))
import run_rna_research_addons as addons


class RNAResearchAddonTests(unittest.TestCase):
    def test_enhanced_status_uses_basejumper_and_cohort_tin_flags(self):
        self.assertEqual(
            addons.improved_qc_status("PASS", [], -2.2, 12)[0], "Borderline")
        self.assertEqual(
            addons.improved_qc_status("PASS", [], -3.2, 12)[0], "FAIL")
        self.assertEqual(
            addons.improved_qc_status("FAIL", ["bad MT fraction"], 0, 12)[0],
            "FAIL")
        self.assertEqual(
            addons.improved_qc_status("PASS", [], float("inf"), 12)[0], "PASS")
        self.assertEqual(
            addons.improved_qc_status("PASS", [], 0, 4)[0], "Borderline")

    def test_positional_duplication_summary_and_plot(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "qc/duplication"
            output.mkdir(parents=True)
            (output / "cell_a.dup.pos.DupRate.xls").write_text(
                "occurrence\tuniquely_mapped_reads\n1\t80\n2\t10\n")
            summary = addons.collect_duplication_summary(
                root, [{"sample_id": "cell_a", "group": "test"}])
            self.assertEqual(summary.loc[0, "distinct_mapped_positions"], 90)
            self.assertEqual(summary.loc[0, "reads_at_mapped_positions"], 100)
            self.assertAlmostEqual(
                summary.loc[0, "estimated_positional_duplication_percent"], 10)
            self.assertTrue((output / "positional_duplication.png").is_file())

    def test_exon_count_matrix_and_within_gene_fractions(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            count_dir = root / "exon_counts"
            count_dir.mkdir()
            (count_dir / "cell_a.counts.txt").write_text(
                "groupID\tfeatureID\tcount\nG1\tE001\t3\nG1\tE002\t7\n")
            (count_dir / "cell_b.counts.txt").write_text(
                "G1\tE001\t5\nG1\tE002\t5\n")
            matrix = addons.collect_exon_counts(root, [
                {"sample_id": "cell_a", "group": "test"},
                {"sample_id": "cell_b", "group": "test"},
            ])
            self.assertEqual(matrix.shape, (2, 4))
            fractions = pd.read_csv(count_dir / "exon_usage_fraction.tsv", sep="\t")
            self.assertAlmostEqual(fractions.loc[0, "cell_a"], 0.3)
            self.assertAlmostEqual(fractions.loc[1, "cell_b"], 0.5)

    def test_arriba_header_and_evidence_are_parsed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fusion_dir = root / "fusions"
            fusion_dir.mkdir()
            (fusion_dir / "cell_a.fusions.tsv").write_text(
                "#gene1\tgene2\tbreakpoint1\tbreakpoint2\tconfidence\tfilters\t"
                "split_reads1\tsplit_reads2\tdiscordant_mates\n"
                "ALK\tEML4\tchr2:29,200,000\tchr2:42,500,000\thigh\tno_filter\t"
                "4\t3\t8\n")
            fasta = root / "reference.fa"
            Path(f"{fasta}.fai").write_text(
                "chr1\t1000000\t0\t0\t0\nchr2\t243000000\t0\t0\t0\n")
            result = addons.summarize_fusions(
                root, [{"sample_id": "cell_a", "group": "test"}], fasta)
            self.assertEqual(result.loc[0, "evidence_tier"], "HIGH")
            self.assertEqual(result.loc[0, "total_supporting_reads"], 15)
            self.assertTrue((fusion_dir / "fusion_chromosome_map.png").is_file())


if __name__ == "__main__":
    unittest.main()
