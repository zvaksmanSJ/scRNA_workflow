import gzip
import json
from pathlib import Path
import tempfile
import unittest

from python_lsf_pipeline.audit_expanded_fastqs import audit_one, combine


class ExpandedFastqAuditTests(unittest.TestCase):
    @staticmethod
    def write_fastq(path, mate):
        with gzip.open(path, "wt") as stream:
            stream.write(f"@INST:1:FC:1:1101:1:1/{mate}\nACGT\n+\n!!!!\n")

    def test_audits_one_cell_and_builds_cohort_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq_dir = root / "fastq"
            fastq_dir.mkdir()
            sid = "SY5Y_Cell17_RNA_SJNBL999999_C0017_RNA_C17"
            r1 = fastq_dir / f"{sid}_R1.merged.fastq.gz"
            r2 = fastq_dir / f"{sid}_R2.merged.fastq.gz"
            self.write_fastq(r1, 1)
            self.write_fastq(r2, 2)
            source_r1 = root / "source_R1.fastq.gz"
            source_r2 = root / "source_R2.fastq.gz"
            self.write_fastq(source_r1, 1)
            self.write_fastq(source_r2, 2)
            source_pairs = []
            for lane, mate_paths in enumerate(((source_r1, source_r2),), start=1):
                source_pairs.append({
                    "sample_number": 1,
                    "lane": f"L{lane:03d}",
                    "run_lane": f"RUN:1:FC:{lane}",
                    "part": "001",
                    "R1": {"path": str(mate_paths[0]), "size": mate_paths[0].stat().st_size,
                           "mtime_ns": mate_paths[0].stat().st_mtime_ns},
                    "R2": {"path": str(mate_paths[1]), "size": mate_paths[1].stat().st_size,
                           "mtime_ns": mate_paths[1].stat().st_mtime_ns},
                })
            cell = {"sample_id": sid, "group": "SY5Y", "pairs": source_pairs,
                    "r1": str(r1), "r2": str(r2)}
            manifest = root / "rna_cells.tsv"
            manifest.write_text(
                "sample_id\tr1\tr2\tgroup\n"
                f"{sid}\t{r1}\t{r2}\tSY5Y\n")
            raw_plan = root / "raw_inputs.json"
            raw_plan.write_text(json.dumps({"cells": [cell]}))
            provenance = {
                "input": cell,
                "source_pair_counts": [1],
                "merged_pairs": 1,
                "outputs": {
                    "r1": [r1.stat().st_size, r1.stat().st_mtime_ns],
                    "r2": [r2.stat().st_size, r2.stat().st_mtime_ns],
                },
            }
            (fastq_dir / f"{sid}.merge.json").write_text(json.dumps(provenance))
            reports = root / "audit"

            result = audit_one(manifest, raw_plan, reports, 1)
            report_path = root / "fastq_audit.json"
            combine(manifest, raw_plan, reports, report_path)

            self.assertEqual(result["status"], "PASS")
            cohort = json.loads(report_path.read_text())
            self.assertEqual(cohort["status"], "PASS")
            self.assertEqual(cohort["samples"], 1)
            self.assertEqual(cohort["merged_fastqs"], 2)
            self.assertEqual(cohort["total_read_pairs"], 1)


if __name__ == "__main__":
    unittest.main()
