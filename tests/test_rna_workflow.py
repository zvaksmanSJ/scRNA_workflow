import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import textwrap
import unittest
import sys
from unittest.mock import patch

WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
for directory in (
        "legacy_nextflow", "python_lsf_pipeline", "setup", "extra_methods"):
    sys.path.insert(0, str(WORKFLOW_ROOT / directory))

import run_all_resolveome_rna as workflow
import step1_prepare_samples as samples
import step2_prepare_references as refs
import step3_prepare_basejumper as runtime
import step0_prepare_raw_fastqs as raw
import step4_rna_dashboard as reports


class RNAWorkflowTests(unittest.TestCase):
    def test_staged_alignment_depth_routing_and_quick_gate(self):
        path = WORKFLOW_ROOT / "legacy_nextflow/stage_adapters/run_step4_alignment.py"
        sys.path.insert(0, str(path.parent))
        spec = importlib.util.spec_from_file_location("staged_alignment", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        for depth in ("qc", "full"):
            with self.subTest(depth=depth), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                prepared = root / "runtime"
                prepared.mkdir()
                (prepared / "main.nf").write_text("process EXAMPLE {}\nworkflow {\n}\n")
                rows = [{"sample_id": f"cell{i}", "r1": "old_r1", "r2": "old_r2",
                         "group": "test"} for i in range(3)]
                name = "step3_output_audit" if depth == "qc" else "step3_full_output_audit_quick"
                audit = {"status": "PASS", "samples": [{"sample": r["sample_id"]} for r in rows]}
                if depth == "full":
                    audit.update(audit_mode="quick", depth="full",
                                 deep_checked_cells=[r["sample_id"] for r in rows])
                (root / f"{name}.json").write_text(json.dumps(audit))
                argv = ["driver"] if depth == "qc" else [
                    "driver", "--depth", "full", "--accept-quick-audit"]
                with patch.object(module, "HERE", root), \
                        patch.object(module, "read_samples", return_value=rows), \
                        patch.object(module, "check_references"), \
                        patch.object(module, "prepare", return_value=prepared), \
                        patch.object(module, "write_config", return_value=root / "native.config"), \
                        patch.object(module, "execute") as execute, patch.object(sys, "argv", argv):
                    module.main()
                    command = execute.call_args.args[0]
                    self.assertEqual(command[command.index("--outputDir") + 1],
                                     root / f"stage4_{depth}_alignment")
                    manifest = (root / f"stage4_{depth}_alignment/rna_input.csv").read_text()
                    self.assertIn(f"stage3_{depth}_preprocessing/trimmed", manifest)
                    self.assertNotIn("old_r1", manifest)
                    with self.assertRaisesRegex(RuntimeError, "already started"):
                        module.main()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "step3_full_output_audit_quick.json").write_text(json.dumps({
                "status": "PASS", "audit_mode": "quick", "depth": "qc",
                "deep_checked_cells": ["a", "b", "c"]}))
            with patch.object(module, "HERE", root), patch.object(sys, "argv", [
                    "driver", "--depth", "full", "--accept-quick-audit"]):
                with self.assertRaisesRegex(ValueError, "depth and deep-cell"):
                    module.main()

    @unittest.skipUnless(runtime.SOURCE.is_dir(),
                         "requires separately authorized BaseJumper source")
    def test_report_parser_reads_htseq_three_column_counts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prepared = runtime.prepare(runtime.SOURCE, root / "runtime")
            source = (prepared / "main.nf").read_text()
            block = source[source.index('    htseq_files = glob.glob'):
                           source.index('    # ========== Parse Gene Detection Metrics')]
            # Unescape the embedded Nextflow string to obtain the executed Python.
            block = textwrap.dedent(block).replace('\\\\t', '\\t')
            counts = root / "sample.htseq_counts.tsv"
            counts.write_text("ENSG1\tGENE1\t12\nENSG2\tGENE2\t3\n"
                              "__no_feature\t\t5\n__ambiguous\t\t2\n"
                              "__alignment_not_unique\t\t1\n")
            import glob
            summary = {}
            with patch.object(glob, "glob", return_value=[str(counts)]):
                exec(block, {"glob": glob, "sample": "sample", "summary": summary})
            self.assertEqual(summary["htseq_assigned"], 15)
            self.assertEqual(summary["htseq_total_reads"], 23)
            counts.write_text("ENSG1\tGENE1\tnot_a_count\n")
            with patch.object(glob, "glob", return_value=[str(counts)]):
                with self.assertRaisesRegex(ValueError, "Could not parse HTSeq"):
                    exec(block, {"glob": glob, "sample": "sample", "summary": {}})

    def write_fastq(self, path, names):
        with gzip.open(path, "wt") as stream:
            for name in names:
                stream.write(f"@{name}\nACGT\n+\nIIII\n")

    def test_pair_validation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_fastq(root / "r1.gz", ["a/1", "b/1"])
            self.write_fastq(root / "r2.gz", ["a/2", "b/2"])
            self.assertEqual(samples.validate_pairs(root / "r1.gz", root / "r2.gz"), 2)
            self.write_fastq(root / "r2.gz", ["a/2"])
            with self.assertRaisesRegex(ValueError, "Unpaired"):
                samples.validate_pairs(root / "r1.gz", root / "r2.gz")

    def test_corrupt_gzip_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in ("r1.gz", "r2.gz"):
                self.write_fastq(root / name, ["a"])
            data = bytearray((root / "r2.gz").read_bytes())
            data[-8] ^= 1
            (root / "r2.gz").write_bytes(data)
            with self.assertRaises(gzip.BadGzipFile):
                samples.validate_pairs(root / "r1.gz", root / "r2.gz")

    def test_plan_never_executes(self):
        with patch.object(workflow, "execute") as execute:
            workflow.main(["--plan", "--cells", "/tmp/rna_samples.tsv",
                           "--outdir", "/tmp/rna_workflow_plan"])
            execute.assert_not_called()

    def test_depths_and_resource_flags(self):
        args = workflow.parse_args([
            "--plan", "--cells", "/tmp/rna_samples.tsv",
            "--outdir", "/tmp/rna_workflow_plan"])
        qc = workflow.branch_command(args, "qc_2M")
        full = workflow.branch_command(args, "full_depth")
        self.assertEqual(qc[qc.index("--n_reads") + 1], 2000000)
        self.assertEqual(full[full.index("--n_reads") + 1], 9223372036854775807)
        self.assertIn("--ref_flat", qc)
        self.assertIn("--max_memory", qc)
        self.assertNotIn("--gender", qc)

    def test_unsafe_output_rejected(self):
        for output in (workflow.HERE, workflow.HERE.parent / "DNA_results",
                       workflow.HERE.parent / "results/old"):
            with self.assertRaises(SystemExit):
                workflow.parse_args([
                    "--cells", "/tmp/rna_samples.tsv", "--outdir", str(output)])

    def test_reference_flat_coordinates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "genes.gtf").write_text(
                'chr1\tx\texon\t11\t20\t.\t+\t.\tgene_id "g"; transcript_id "t"; gene_name "G";\n'
                'chr1\tx\texon\t31\t40\t.\t+\t.\tgene_id "g"; transcript_id "t"; gene_name "G";\n'
                'chr1\tx\tCDS\t15\t35\t.\t+\t0\tgene_id "g"; transcript_id "t"; gene_name "G";\n')
            refs.create_ref_flat(root / "genes.gtf", root / "refFlat")
            self.assertEqual((root / "refFlat").read_text(),
                             "G\tt\tchr1\t+\t10\t40\t14\t35\t2\t10,30,\t20,40,\n")

    def test_public_reference_annotations_derive_from_one_gtf(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gtf = root / "genes.gtf"
            gtf.write_text(
                'chr1\tx\ttranscript\t11\t40\t.\t+\t.\tgene_id "g"; transcript_id "t"; '
                'gene_name "GENE"; gene_biotype "protein_coding"; '
                'transcript_biotype "protein_coding";\n'
                'chr1\tx\texon\t11\t20\t.\t+\t.\tgene_id "g"; transcript_id "t"; '
                'gene_name "GENE"; gene_biotype "protein_coding"; '
                'transcript_biotype "protein_coding";\n'
                'chr1\tx\texon\t31\t40\t.\t+\t.\tgene_id "g"; transcript_id "t"; '
                'gene_name "GENE"; gene_biotype "protein_coding"; '
                'transcript_biotype "protein_coding";\n')
            tx2gene = root / "tx2gene.tsv"
            bed = root / "transcripts.bed"
            self.assertEqual(refs.create_tx2gene(gtf, tx2gene), 1)
            self.assertEqual(refs.create_genebody_bed(gtf, bed), 1)
            self.assertEqual(
                tx2gene.read_text(),
                "t\tprotein_coding\tg\tprotein_coding\tGENE\n")
            self.assertEqual(
                bed.read_text(),
                "chr1\t10\t40\tg|t\t0\t+\t10\t40\t0\t2\t10,10,\t0,20,\n")

    @unittest.skipUnless(runtime.SOURCE.is_dir(),
                         "requires separately authorized BaseJumper source")
    def test_runtime_adaptation_retains_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            before = runtime.source_hashes(runtime.SOURCE)
            prepared = runtime.prepare(runtime.SOURCE, Path(temporary) / "runtime")
            text = (prepared / "main.nf").read_text()
            self.assertNotIn("/usr/local/bin/", text)
            self.assertIn("--twopassMode Basic", text)
            self.assertIn("--secondary-alignments=ignore", text)
            self.assertIn("RAW_FASTQC(ch_fastq_input)", text)
            self.assertIn("PICARD_RNA(", text)
            self.assertIn("count = int(parts[-1])", text)
            self.assertNotIn('count = int(parts[1]) if parts[1].isdigit() else 0', text)
            self.assertIn('sample.bam >> "bam/${params.workspace}/rna/tool=star"', text)
            self.assertIn('sample.bai >> "bam/${params.workspace}/rna/tool=star"', text)
            self.assertIn("--validate-pair", text)
            self.assertIn("export FINAL_READS=\\$((2 *", text)
            self.assertEqual(runtime.source_hashes(runtime.SOURCE), before)
            self.assertEqual(runtime.prepare(runtime.SOURCE, prepared), prepared)
            config = Path(temporary) / "native.config"
            conda_base = Path(temporary) / "conda"
            activation = conda_base / "etc/profile.d/conda.sh"
            activation.parent.mkdir(parents=True)
            activation.write_text("")
            runtime.write_config(prepared, config, workflow.PREFIX, "lsf", "standard",
                                 8, 64000, conda_base)
            config_text = config.read_text()
            self.assertIn(str(activation), config_text)
            self.assertIn("errorStrategy = 'terminate'", config_text)
            self.assertIn("resourceLimits = [cpus: 8, memory: '64000 MB']", config_text)
            self.assertIn("docker.enabled = false", config_text)
            self.assertIn("executor.perTaskReserve = true", config_text)
            self.assertNotIn("clusterOptions", config_text)
            self.assertNotIn("maxForks", config_text)
            self.assertNotIn("queueSize", config_text)
            self.assertIn("workflow.output.mode = 'copy'", config_text)

    def test_upstream_change_is_not_silently_accepted(self):
        with self.assertRaisesRegex(ValueError, "Upstream compatibility"):
            runtime.replace_once("different source", "required marker", "replacement")

    def test_raw_flowcell_merges_preserve_all_pairs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = root / "raw"
            inputs.mkdir()
            for sample, run in ((1, 102), (17, 101)):
                for lane in (1, 2):
                    for mate in (1, 2):
                        path = inputs / (
                            f"SJNBL999999_C9999-ResolOME-SY5Y_P2-RNA_C1_S{sample}"
                            f"_L{lane:03d}_R{mate}_001.fastq.gz")
                        self.write_fastq(path, [f"LH:{run}:FC{run}:{lane}:1101:1:1/{mate}"])
            manifest = raw.plan(inputs, root / "output")
            self.assertEqual(len(manifest["cells"]), 1)
            self.assertEqual(len(manifest["cells"][0]["pairs"]), 4)
            raw.write_plan(manifest, root / "output")
            raw.merge_cell(root / "output/raw_inputs.json", 1)
            cell = manifest["cells"][0]
            self.assertEqual(samples.validate_pairs(cell["r1"], cell["r2"]), 4)
            raw.merge_cell(root / "output/raw_inputs.json", 1)
            self.assertEqual(len(list(inputs.glob("*.gz"))), 8)
            first = next(inputs.glob("*.gz"))
            first.write_bytes(first.read_bytes() + b"x")
            with self.assertRaisesRegex(ValueError, "Raw input changed"):
                raw.merge_cell(root / "output/raw_inputs.json", 1)

    def test_extended_batch_raw_plan_uses_existing_names_and_fastq_folder(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = root / "raw"
            inputs.mkdir()
            libraries = (
                ("SJNBL999999_C0017", "P3", 17, 1, 101, (1, 2)),
                ("SJNBL999999_C0033", "P4", 33, 2, 102, (1,)),
            )
            for accession, plate, cell, sample, run, lanes in libraries:
                for lane in lanes:
                    for mate in (1, 2):
                        path = inputs / (
                            f"{accession}-ResolOME-SY5Y_{plate}-RNA_C{cell}_S{sample}"
                            f"_L{lane:03d}_R{mate}_001.fastq.gz")
                        self.write_fastq(
                            path, [f"LH:{run}:FC{run}:{lane}:1101:1:1/{mate}"])

            manifest = raw.plan(inputs, root / "expanded", "fastq")
            cells = {cell["sample_id"]: cell for cell in manifest["cells"]}
            self.assertEqual(set(cells), {
                "SY5Y_Cell17_RNA_SJNBL999999_C0017_RNA_C17",
                "SY5Y_Cell33_RNA_SJNBL999999_C0033_RNA_C33",
            })
            self.assertEqual(len(cells["SY5Y_Cell17_RNA_SJNBL999999_C0017_RNA_C17"]["pairs"]), 2)
            self.assertEqual(len(cells["SY5Y_Cell33_RNA_SJNBL999999_C0033_RNA_C33"]["pairs"]), 1)
            self.assertTrue(cells["SY5Y_Cell17_RNA_SJNBL999999_C0017_RNA_C17"]["r1"]
                            .endswith("/fastq/SY5Y_Cell17_RNA_SJNBL999999_C0017_RNA_C17_R1.merged.fastq.gz"))

    def test_duplicate_physical_flowcell_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for sample in (1, 17):
                for mate in (1, 2):
                    self.write_fastq(root / (
                        f"SJNBL999999_C9999-ResolOME-SY5Y_P2-RNA_C1_S{sample}"
                        f"_L001_R{mate}_001.fastq.gz"), [f"LH:101:FC:1:1101:1:1/{mate}"])
            with self.assertRaisesRegex(ValueError, "Duplicate physical"):
                raw.plan(root, root / "output")

    def test_dashboard_uses_actual_qc_read_target(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "workflow_state.json").write_text(json.dumps({"fingerprint": {"qc_reads": 45000}}))
            with patch.object(reports, "audit_branch", return_value={"samples": 1}), \
                    patch.object(reports.subprocess, "run"):
                reports.dashboard(root, ["qc_2M"], ["cell"])
            self.assertIn("45,000 individual reads", (root / "workflow_dashboard.html").read_text())


if __name__ == "__main__":
    unittest.main()
