#!/usr/bin/env python3
"""Build the single Cromwell inputs JSON from the validated cohort and references."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

RNAW = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RNAW))
from python_lsf_pipeline.step1_prepare_samples import read_samples

RUNTIME = Path(os.environ.get("RNA_RUNTIME_DIR", RNAW / "wdl_runtime"))
RUNTIME_CONFIG = RUNTIME / "runtime_config.json"


def load_runtime_config(runtime_config: Path) -> dict:
    config_path = runtime_config.resolve(strict=True)
    config = json.loads(config_path.read_text())
    refs = config["references"]
    for key in ("fasta", "gtf", "tx2gene", "genebody", "ref_flat"):
        path = Path(refs[key]).resolve(strict=True)
        if not path.is_file():
            raise FileNotFoundError(f"Missing configured RNA reference: {path}")
        refs[key] = str(path)
    for key, marker in (
        ("star_index", "Genome"),
        ("salmon_index", "versionInfo.json"),
    ):
        index = Path(refs[key]).resolve(strict=True)
        if not index.is_dir() or not (index / marker).is_file():
            raise FileNotFoundError(
                f"Missing configured RNA index or marker {marker}: {index}")
        refs[key] = str(index)
    return config


def build_inputs(results_dir: Path, basej_source_dir: Path,
                 runtime_config: Path = RUNTIME_CONFIG) -> dict:
    results_dir = results_dir.resolve(strict=True)
    manifest = results_dir / "wdl_run_manifest.json"
    if not manifest.is_file():
        raise FileNotFoundError(
            f"Prepare an isolated WDL results tree first: {results_dir}")
    state = json.loads(manifest.read_text())
    if state.get("status") != "PREPARED" or state.get("workflow_rna_dir") != str(RNAW):
        raise ValueError(f"Invalid WDL run manifest: {manifest}")
    basej_source_dir = basej_source_dir.resolve(strict=True)
    if not (basej_source_dir / "main.nf").is_file():
        raise FileNotFoundError(f"BaseJumper Nextflow source is incomplete: {basej_source_dir}")
    config = load_runtime_config(runtime_config)
    samples = read_samples(results_dir / "rna_cells.tsv")
    if len(samples) < 5:
        raise ValueError("The cohort workflow requires at least five samples")
    refs = config["references"]
    inputs = {
        "ResolveOmeRnaWorkflow.workflow_rna_dir": str(RNAW),
        "ResolveOmeRnaWorkflow.runtime_dir": str(runtime_config.resolve().parent),
        "ResolveOmeRnaWorkflow.results_dir": str(results_dir),
        "ResolveOmeRnaWorkflow.basej_source_dir": str(basej_source_dir),
        "ResolveOmeRnaWorkflow.sample_sheet": str(results_dir / "rna_cells.tsv"),
        "ResolveOmeRnaWorkflow.sample_ids": [row["sample_id"] for row in samples],
        "ResolveOmeRnaWorkflow.ref_fasta": refs["fasta"],
        "ResolveOmeRnaWorkflow.ref_gtf": refs["gtf"],
        "ResolveOmeRnaWorkflow.tx2gene": refs["tx2gene"],
        "ResolveOmeRnaWorkflow.genebody_bed": refs["genebody"],
        "ResolveOmeRnaWorkflow.star_index_dir": refs["star_index"],
        "ResolveOmeRnaWorkflow.ref_flat": refs["ref_flat"],
        "ResolveOmeRnaWorkflow.conda_base": config["conda_base"],
        "ResolveOmeRnaWorkflow.conda_prefix": config["rna_environment"],
        "ResolveOmeRnaWorkflow.salmon_executable": str(
            Path(config["salmon_environment"]) / "bin/salmon"),
        "ResolveOmeRnaWorkflow.salmon_index": refs["salmon_index"],
    }
    return inputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument(
        "--basej-source-dir", type=Path,
        default=Path(os.environ["BASEJ_RNA_SOURCE"]) if os.environ.get("BASEJ_RNA_SOURCE") else None,
        required=not bool(os.environ.get("BASEJ_RNA_SOURCE")),
    )
    parser.add_argument("--runtime-config", type=Path, default=RUNTIME_CONFIG)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build_inputs(args.results_dir, args.basej_source_dir, args.runtime_config)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {len(result)} WDL workflow inputs: {args.output}")


if __name__ == "__main__":
    main()
