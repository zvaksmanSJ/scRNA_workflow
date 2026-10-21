import importlib
import json
import os
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "extra_methods"))
import quantify_transcripts as salmon


class SalmonWorkflowTests(unittest.TestCase):
    def test_results_root_can_be_overridden_for_isolated_wdl_runs(self):
        with tempfile.TemporaryDirectory() as temporary:
            try:
                with patch.dict(os.environ, {"RESOLVEOME_RNA_RESULTS": temporary}):
                    importlib.reload(salmon)
                    self.assertEqual(salmon.RNAR, Path(temporary).resolve())
            finally:
                importlib.reload(salmon)

    def test_partial_submission_manifest_is_persisted(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            samples = [("cell_a", root / "a_R1.gz", root / "a_R2.gz"),
                       ("cell_b", root / "b_R1.gz", root / "b_R2.gz")]
            salmon_rows = [
                {"sample_id": sid, "job_id": "", "status": "NOT_SUBMITTED",
                 "output": str(root / sid), "stdout_log_pattern": str(root / f"{sid}.%J.out")}
                for sid, _, _ in samples
            ]
            success = subprocess.CompletedProcess(
                args=["bsub"], returncode=0,
                stdout="Job <1234> is submitted to queue <standard>.\n", stderr="")
            failure = subprocess.CalledProcessError(1, ["bsub"], stderr="scheduler error")
            with patch.object(salmon, "preflight", return_value=samples), \
                    patch.object(salmon, "OUT", root), \
                    patch.object(salmon, "QUANTS", root / "quant"), \
                    patch.object(salmon, "JOBS", root / "jobs.tsv"), \
                    patch.object(salmon, "PREFIX", root / "env"), \
                    patch("subprocess.run", side_effect=[success, failure]):
                with self.assertRaises(subprocess.CalledProcessError):
                    salmon.submit()
            with (root / "jobs.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["status"], "SUBMITTED")
            self.assertEqual(rows[0]["job_id"], "1234")
            self.assertEqual(rows[1]["status"], "NOT_SUBMITTED")

    def test_validate_outputs_checks_all_samples_and_reference_consistency(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rnar = root / "RNA_results"
            rnar.mkdir()
            out = rnar / "final_results"
            quants = out / "02_transcript/quant_by_sample"
            index = root / "salmon_index"
            prefix = root / "env"
            tx2gene = root / "tx2gene.tsv"
            index.mkdir()
            (index / "versionInfo.json").write_text(json.dumps({"salmonVersion": "1.6.0"}))
            (prefix / "bin").mkdir(parents=True)
            (prefix / "bin/salmon").touch()
            tx2gene.touch()
            (rnar / "step4_full_output_audit_quick.json").write_text(json.dumps({
                "status": "PASS", "depth": "full",
            }))
            valid_quant = (
                "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
                + "".join(f"tx{number}\t100\t90\t1000\t5\n"
                          for number in range(1000))
            )
            samples = []
            for number in range(16):
                sample = f"cell_{number}"
                samples.append({"sample_id": sample})
                trimmed = rnar / f"stage3_full_preprocessing/trimmed/{sample}"
                trimmed.parent.mkdir(parents=True, exist_ok=True)
                trimmed.with_name(f"{sample}_R1_trim.fastq.gz").touch()
                trimmed.with_name(f"{sample}_R2_trim.fastq.gz").touch()
                quant_dir = quants / sample
                quant_dir.mkdir(parents=True)
                (quant_dir / "quant.sf").write_text(valid_quant)
            jobs = out / "07_qc_provenance/salmon_jobs.tsv"
            jobs.parent.mkdir(parents=True)
            with jobs.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=["sample_id", "job_id", "status", "output",
                                        "stdout_log_pattern"], delimiter="\t")
                writer.writeheader()
                writer.writerows({
                    "sample_id": row["sample_id"], "job_id": str(index + 100),
                    "status": "SUBMITTED", "output": str(quants / row["sample_id"]),
                    "stdout_log_pattern": "salmon.log",
                } for index, row in enumerate(samples))
            with patch.object(salmon, "read_samples", return_value=samples), \
                    patch.object(salmon, "RNAR", rnar), \
                    patch.object(salmon, "OUT", out), \
                    patch.object(salmon, "QUANTS", quants), \
                    patch.object(salmon, "JOBS", jobs), \
                    patch.object(salmon, "INDEX", index), \
                    patch.object(salmon, "PREFIX", prefix), \
                    patch.object(salmon, "TX2GENE", tx2gene):
                salmon.validate_outputs()
                audit = json.loads(
                    (out / "07_qc_provenance/salmon_output_audit.json").read_text())
                self.assertEqual(audit["status"], "PASS")
                self.assertEqual(len(audit["samples"]), 16)
                duplicate_quant = (
                    "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
                    + "".join(f"tx{0 if number == 1 else number}\t100\t90\t1000\t5\n"
                              for number in range(1000))
                )
                (quants / "cell_1/quant.sf").write_text(duplicate_quant)
                with self.assertRaisesRegex(ValueError, "Invalid Salmon transcript row"):
                    salmon.validate_outputs()


if __name__ == "__main__":
    unittest.main()
