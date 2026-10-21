import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "wdl/scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "extra_methods"))
import prepare_wdl_run
import build_inputs_json
import build_run_inputs
import organize_rna_results


class WdlPreparationTests(unittest.TestCase):
    def test_results_readme_uses_run_cohort_and_audit_scope(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (root / "rna_cells.tsv").open("w", newline="") as stream:
                writer = csv.DictWriter(
                    stream, fieldnames=["sample_id", "r1", "r2", "group"],
                    delimiter="\t")
                writer.writeheader()
                writer.writerows({
                    "sample_id": f"cell_{index}", "r1": "r1.fastq.gz",
                    "r2": "r2.fastq.gz", "group": "test",
                } for index in range(6))
            (root / "step4_full_output_audit_quick.json").write_text(json.dumps({
                "status": "PASS", "audit_mode": "quick",
                "deep_checked_cells": ["cell_0", "cell_2"],
            }))
            organize_rna_results.write_results_readme(root)
            text = (root / "README.md").read_text()
            self.assertIn("All 6 samples are retained", text)
            self.assertIn("deeply checked 2 sample(s)", text)

    def test_organizer_links_reactome_from_source_without_requiring_a_logo(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "results"
            source = Path(temporary) / "workflow"
            assets = source / "analysis_resources"
            assets.mkdir(parents=True)
            references = root / "final_results/09_references"
            references.mkdir(parents=True)
            for result_path, asset_path in organize_rna_results.PACKAGE_REFERENCE_LINKS.items():
                canonical = assets / asset_path.removeprefix("analysis_resources/")
                canonical.write_text(asset_path)
                (root / result_path).write_text(asset_path)
            with patch.object(organize_rna_results, "RNAW", source), \
                    patch.object(organize_rna_results, "PACKAGE_AUDIT_LINKS", {}), \
                    patch.object(organize_rna_results, "PACKAGE_LOG_MOVES", {}):
                links = organize_rna_results.check_duplicate_links(root)
            self.assertEqual(len(links), 2)
            self.assertTrue(all(canonical.parent == assets for _, canonical in links))

    def test_prepare_copies_stage_code_and_links_immutable_fastqs(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            fastqs = project / "inputs/merged_fastqs"
            fastqs.mkdir(parents=True)
            sample_sheet = project / "rna_cells.tsv"
            columns = ["sample_id", "r1", "r2", "group"]
            samples = [f"cell_{index}" for index in range(5)]
            with sample_sheet.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
                writer.writeheader()
                rows = []
                for sample in samples:
                    mates = []
                    for mate in ("R1", "R2"):
                        path = fastqs / f"{sample}_{mate}.fastq.gz"
                        with gzip.open(path, "wt") as output:
                            output.write(f"@{sample}/1\nACGT\n+\nIIII\n")
                        mates.append(path)
                    rows.append({"sample_id": sample, "r1": str(mates[0]),
                                 "r2": str(mates[1]), "group": "test"})
                writer.writerows(rows)
            for sample in samples:
                (fastqs / f"{sample}.merge.json").write_text("{}\n")

            destination = project / "RNA_results_wdl"
            prepare_wdl_run.prepare(destination, sample_sheet)

            self.assertEqual(
                (destination / "fastqs/cell_0_R1.fastq.gz").resolve(),
                (fastqs / "cell_0_R1.fastq.gz").resolve())
            self.assertTrue((destination / "fastqs/cell_0.merge.json").is_symlink())
            copied_runner = destination / "run_step3_preprocessing.py"
            self.assertFalse(copied_runner.is_symlink())
            self.assertEqual(copied_runner.resolve().parent, destination.resolve())
            state = json.loads((destination / "wdl_run_manifest.json").read_text())
            self.assertEqual(state["status"], "PREPARED")
            self.assertEqual(state["sample_count"], 5)
            self.assertTrue(state["input_files_are_symlinks"])

            references = project / "references"
            references.mkdir()
            for name in ("fasta", "gtf", "tx2gene", "genebody", "ref_flat"):
                (references / name).write_text("fixture\n")
            star_index = references / "STAR"
            star_index.mkdir()
            (star_index / "Genome").write_text("fixture\n")
            salmon_index = references / "salmon"
            salmon_index.mkdir()
            (salmon_index / "versionInfo.json").write_text("{}\n")
            basej_source = project / "basej"
            basej_source.mkdir()
            (basej_source / "main.nf").write_text("workflow {}\n")
            runtime_config = project / "runtime_config.json"
            runtime_config.write_text(json.dumps({
                "conda_base": str(project / "conda"),
                "rna_environment": str(project / "env/rna"),
                "salmon_environment": str(project / "env/salmon"),
                "references": {
                    "fasta": str(references / "fasta"),
                    "gtf": str(references / "gtf"),
                    "tx2gene": str(references / "tx2gene"),
                    "genebody": str(references / "genebody"),
                    "star_index": str(star_index),
                    "ref_flat": str(references / "ref_flat"),
                    "salmon_index": str(salmon_index),
                },
            }))
            inputs = build_inputs_json.build_inputs(
                destination, basej_source, runtime_config)
            self.assertEqual(inputs["ResolveOmeRnaWorkflow.sample_ids"], samples)
            self.assertEqual(
                inputs["ResolveOmeRnaWorkflow.basej_source_dir"], str(basej_source))
            self.assertEqual(
                inputs["ResolveOmeRnaWorkflow.runtime_dir"], str(project))

    def test_build_sprocket_inputs_uses_setup_references_and_sample_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastqs = root / "fastqs"
            fastqs.mkdir()
            rows = []
            for index in range(2):
                r1 = fastqs / f"cell_{index}_R1.fastq.gz"
                r2 = fastqs / f"cell_{index}_R2.fastq.gz"
                r1.write_bytes(b"fixture")
                r2.write_bytes(b"fixture")
                rows.append({
                    "sample_id": f"cell_{index}", "r1": str(r1), "r2": str(r2),
                    "group": "test",
                })
            sample_sheet = root / "samples.tsv"
            with sample_sheet.open("w", newline="") as stream:
                writer = csv.DictWriter(
                    stream, fieldnames=["sample_id", "r1", "r2", "group"],
                    delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

            refs = root / "references"
            refs.mkdir()
            for name in ("fasta", "gtf", "tx2gene", "genebody", "ref_flat"):
                (refs / name).write_text("fixture\n")
            star_index = refs / "STAR"
            star_index.mkdir()
            (star_index / "Genome").write_text("fixture\n")
            salmon_index = refs / "salmon"
            salmon_index.mkdir()
            (salmon_index / "versionInfo.json").write_text("{}\n")
            runtime_config = root / "runtime_config.json"
            runtime_config.write_text(json.dumps({
                "references": {
                    "fasta": str(refs / "fasta"),
                    "gtf": str(refs / "gtf"),
                    "tx2gene": str(refs / "tx2gene"),
                    "genebody": str(refs / "genebody"),
                    "star_index": str(star_index),
                    "ref_flat": str(refs / "ref_flat"),
                    "salmon_index": str(salmon_index),
                },
            }))

            inputs = build_run_inputs.build_sprocket_inputs(
                sample_sheet, runtime_config, "docker://example/rna:1.0")
            prefix = "ResolveOmeRnaSprocketWorkflow."
            self.assertEqual(inputs[prefix + "sample_ids"], ["cell_0", "cell_1"])
            self.assertEqual(
                inputs[prefix + "r1_fastqs"], [row["r1"] for row in rows])
            self.assertEqual(inputs[prefix + "r2_fastqs"], [row["r2"] for row in rows])
            self.assertEqual(inputs[prefix + "star_index"], str(star_index))
            self.assertEqual(inputs[prefix + "salmon_index"], str(salmon_index))
            self.assertEqual(inputs[prefix + "container"], "docker://example/rna:1.0")

    def test_input_json_writer_refuses_overwrite_and_workflow_tree(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "inputs.json"
            build_run_inputs.write_json(output, {"sample_ids": ["cell_1"]})
            with self.assertRaises(FileExistsError):
                build_run_inputs.write_json(output, {"sample_ids": ["cell_2"]})
            with patch.object(build_run_inputs, "RNAW", Path(temporary)):
                with self.assertRaisesRegex(ValueError, "inside workflow source"):
                    build_run_inputs.write_json(Path(temporary) / "generated.json", {})

    def test_build_both_input_files_and_prepare_new_results_tree(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastqs = root / "fastqs"
            fastqs.mkdir()
            sheet = root / "samples.tsv"
            rows = []
            for index in range(5):
                mates = []
                for mate in ("R1", "R2"):
                    path = fastqs / f"cell_{index}_{mate}.fastq.gz"
                    path.write_bytes(b"fixture")
                    mates.append(path)
                rows.append({
                    "sample_id": f"cell_{index}", "r1": str(mates[0]),
                    "r2": str(mates[1]), "group": "test",
                })
            with sheet.open("w", newline="") as stream:
                writer = csv.DictWriter(
                    stream, fieldnames=["sample_id", "r1", "r2", "group"],
                    delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)

            refs = root / "references"
            refs.mkdir()
            for name in ("fasta", "gtf", "tx2gene", "genebody", "ref_flat"):
                (refs / name).write_text("fixture\n")
            star_index = refs / "STAR"
            star_index.mkdir()
            (star_index / "Genome").write_text("fixture\n")
            salmon_index = refs / "salmon"
            salmon_index.mkdir()
            (salmon_index / "versionInfo.json").write_text("{}\n")
            runtime_config = root / "runtime_config.json"
            runtime_config.write_text(json.dumps({
                "conda_base": str(root / "conda"),
                "rna_environment": str(root / "env/rna"),
                "salmon_environment": str(root / "env/salmon"),
                "references": {
                    "fasta": str(refs / "fasta"),
                    "gtf": str(refs / "gtf"),
                    "tx2gene": str(refs / "tx2gene"),
                    "genebody": str(refs / "genebody"),
                    "star_index": str(star_index),
                    "ref_flat": str(refs / "ref_flat"),
                    "salmon_index": str(salmon_index),
                },
            }))
            basej_source = root / "basej"
            basej_source.mkdir()
            (basej_source / "main.nf").write_text("workflow {}\n")
            results_dir = root / "new_run"
            args = type("Args", (), {
                "workflow": "both",
                "sample_sheet": sheet,
                "runtime_config": runtime_config,
                "results_dir": results_dir,
                "basej_source_dir": basej_source,
                "sprocket_container": "docker://example/rna:1.0",
                "sprocket_output": None,
            })()

            outputs = build_run_inputs.build_run(args)
            self.assertEqual(
                {path.name for path in outputs},
                {"wdl_inputs.json", "sprocket_inputs.json"})
            self.assertTrue((results_dir / "fastqs/cell_0_R1.fastq.gz").is_symlink())
            cromwell = json.loads((results_dir / "wdl_inputs.json").read_text())
            sprocket = json.loads((results_dir / "sprocket_inputs.json").read_text())
            self.assertEqual(
                cromwell["ResolveOmeRnaWorkflow.sample_ids"],
                [row["sample_id"] for row in rows])
            self.assertEqual(
                sprocket["ResolveOmeRnaSprocketWorkflow.sample_ids"],
                [row["sample_id"] for row in rows])
            with self.assertRaises(FileExistsError):
                build_run_inputs.build_run(args)


if __name__ == "__main__":
    unittest.main()
