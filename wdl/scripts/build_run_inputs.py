#!/usr/bin/env python3
"""Build RNA workflow input JSON files from a sample sheet and setup config."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

RNAW = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(RNAW))
sys.path.insert(0, str(SCRIPTS))

from python_lsf_pipeline.step1_prepare_samples import read_samples
from build_inputs_json import build_inputs, load_runtime_config
from prepare_wdl_run import prepare

RUNTIME = Path(os.environ.get("RNA_RUNTIME_DIR", RNAW / "wdl_runtime"))
RUNTIME_CONFIG = RUNTIME / "runtime_config.json"
SPROCKET_WORKFLOW = "ResolveOmeRnaSprocketWorkflow"


def build_sprocket_inputs(sample_sheet: Path, runtime_config: Path,
                          container: str) -> dict:
    if not container.strip():
        raise ValueError("A non-empty Sprocket container URI is required")
    samples = read_samples(sample_sheet.resolve(strict=True))
    config = load_runtime_config(runtime_config)
    refs = config["references"]
    return {
        f"{SPROCKET_WORKFLOW}.sample_ids": [row["sample_id"] for row in samples],
        f"{SPROCKET_WORKFLOW}.r1_fastqs": [row["r1"] for row in samples],
        f"{SPROCKET_WORKFLOW}.r2_fastqs": [row["r2"] for row in samples],
        f"{SPROCKET_WORKFLOW}.star_index": refs["star_index"],
        f"{SPROCKET_WORKFLOW}.salmon_index": refs["salmon_index"],
        f"{SPROCKET_WORKFLOW}.ref_gtf": refs["gtf"],
        f"{SPROCKET_WORKFLOW}.ref_flat": refs["ref_flat"],
        f"{SPROCKET_WORKFLOW}.genebody_bed": refs["genebody"],
        f"{SPROCKET_WORKFLOW}.container": container.strip(),
    }


def write_json(path: Path, inputs: dict) -> None:
    path = path.resolve()
    if path == RNAW or RNAW in path.parents:
        raise ValueError(f"Refusing to write generated inputs inside workflow source: {path}")
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite an inputs JSON: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(inputs, stream, indent=2)
        stream.write("\n")


def build_run(args: argparse.Namespace) -> list[Path]:
    sample_sheet = args.sample_sheet.resolve(strict=True)
    config_path = args.runtime_config.resolve(strict=True)
    load_runtime_config(config_path)
    workflow = args.workflow
    outputs: list[Path] = []
    sprocket_inputs = None
    sprocket_output = args.sprocket_output

    if workflow in ("sprocket", "both"):
        if not args.sprocket_container:
            raise ValueError("--sprocket-container is required for Sprocket inputs")
        sprocket_inputs = build_sprocket_inputs(
            sample_sheet, config_path, args.sprocket_container)
        if workflow == "sprocket":
            if sprocket_output is None:
                raise ValueError("--sprocket-output is required for --workflow sprocket")
        elif sprocket_output is None:
            sprocket_output = args.results_dir / "sprocket_inputs.json"

    cromwell_inputs = None
    if workflow in ("cromwell", "both"):
        if args.results_dir is None:
            raise ValueError("--results-dir is required for Cromwell inputs")
        if args.basej_source_dir is None:
            raise ValueError(
                "Provide --basej-source-dir or set BASEJ_RNA_SOURCE for Cromwell inputs")
        basej_source = args.basej_source_dir.resolve(strict=True)
        if not (basej_source / "main.nf").is_file():
            raise FileNotFoundError(
                f"BaseJumper Nextflow source is incomplete: {basej_source}")
        results_dir = args.results_dir.resolve()
        cromwell_output = results_dir / "wdl_inputs.json"
        if cromwell_output.exists():
            raise FileExistsError(f"Refusing to overwrite an inputs JSON: {cromwell_output}")
        if sprocket_output is not None and sprocket_output.exists():
            raise FileExistsError(f"Refusing to overwrite an inputs JSON: {sprocket_output}")
        if sprocket_output is not None:
            sprocket_resolved = sprocket_output.resolve()
            if sprocket_resolved == cromwell_output.resolve():
                raise ValueError("Cromwell and Sprocket JSON outputs must be different files")
            if sprocket_resolved == RNAW or RNAW in sprocket_resolved.parents:
                raise ValueError(
                    f"Refusing to write generated inputs inside workflow source: {sprocket_resolved}")
        prepare(results_dir, sample_sheet)
        cromwell_inputs = build_inputs(results_dir, basej_source, config_path)
        write_json(cromwell_output, cromwell_inputs)
        outputs.append(cromwell_output)

    if sprocket_inputs is not None and sprocket_output is not None:
        write_json(sprocket_output, sprocket_inputs)
        outputs.append(sprocket_output)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workflow", choices=("cromwell", "sprocket", "both"), required=True,
        help="Which workflow input schema(s) to generate")
    parser.add_argument(
        "--sample-sheet", type=Path, required=True,
        help="TSV with sample_id, r1, r2, and group columns")
    parser.add_argument("--runtime-config", type=Path, default=RUNTIME_CONFIG)
    parser.add_argument(
        "--results-dir", type=Path,
        help="New results directory for Cromwell; must not already exist")
    parser.add_argument(
        "--basej-source-dir", type=Path,
        default=Path(os.environ["BASEJ_RNA_SOURCE"])
        if os.environ.get("BASEJ_RNA_SOURCE") else None)
    parser.add_argument(
        "--sprocket-container",
        help="Container URI accessible to Sprocket, e.g. docker://registry/image:tag")
    parser.add_argument(
        "--sprocket-output", type=Path,
        help="Output path for Sprocket JSON; required for --workflow sprocket")
    args = parser.parse_args()

    if args.workflow in ("cromwell", "both") and args.results_dir is None:
        parser.error("--results-dir is required for --workflow cromwell or both")
    if (args.workflow in ("cromwell", "both")
            and args.basej_source_dir is None):
        parser.error(
            "Provide --basej-source-dir or set BASEJ_RNA_SOURCE for Cromwell inputs")
    if args.workflow == "sprocket" and args.sprocket_output is None:
        parser.error("--sprocket-output is required for --workflow sprocket")
    if args.workflow == "sprocket" and not args.sprocket_container:
        parser.error("--sprocket-container is required for --workflow sprocket")
    if args.workflow == "both" and not args.sprocket_container:
        parser.error("--sprocket-container is required for --workflow both")
    if args.workflow not in ("sprocket", "both") and args.sprocket_output is not None:
        parser.error("--sprocket-output is only used with --workflow sprocket or both")

    paths = build_run(args)
    print(f"Prepared inputs for {len(read_samples(args.sample_sheet))} sample(s):")
    for path in paths:
        print(f"  {path}")


if __name__ == "__main__":
    main()
