#!/usr/bin/env python3
"""Run BaseJumper RNA-QC and a separate full-depth expression branch on LSF."""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
for package_dir in ("python_lsf_pipeline", "setup", "legacy_nextflow", "extra_methods"):
    sys.path.insert(0, str(HERE / package_dir))

from python_lsf_pipeline.step1_prepare_samples import read_samples, write_source_csv
from setup.step2_prepare_references import DEFAULTS, check_references
from legacy_nextflow.step3_prepare_basejumper import SOURCE, prepare, source_hashes, write_config
from extra_methods.step4_rna_dashboard import dashboard

RUNTIME = Path(os.environ.get("RNA_RUNTIME_DIR", HERE / "wdl_runtime"))
PREFIX = Path(os.environ.get("RNA_ENV_PREFIX", RUNTIME / "env/resolveome-rna"))
NEXTFLOW = Path(os.environ.get("RNA_NEXTFLOW", RUNTIME / "bin/nextflow"))
TOOLS = ("STAR", "samtools", "fastp", "seqkit", "htseq-count", "qualimap", "fastqc",
         "picard", "geneBody_coverage.py", "multiqc", "Rscript", "java")
PACKAGES = ("ohchibi", "AMOR", "dplyr", "tidyr", "reshape2", "pheatmap",
            "magrittr", "tibble", "ggtree", "optparse", "pals", "ggrepel", "ggpmisc",
            "gt", "egg", "paletteer")


def execute(argv, **kwargs):
    print("[CMD]", shlex.join(map(str, argv)), flush=True)
    try:
        return subprocess.run(list(map(str, argv)), text=True, check=True, **kwargs)
    except subprocess.CalledProcessError as exc:
        for output in (exc.stdout, exc.stderr):
            if output:
                print(output, file=sys.stderr, flush=True)
        raise


def references(args):
    return {name: str(getattr(args, name).resolve()) for name in DEFAULTS}


def preflight(args):
    rows = read_samples(args.cells)
    if len(rows) < 5:
        raise ValueError("Source RNA cohort plotting requires at least 5 samples; "
                         "use a separate targeted tool smoke test for smaller fixtures")
    if Path(sys.prefix).resolve() != args.env_prefix.resolve():
        raise RuntimeError(f"Activate {args.env_prefix}; selected runtime is {sys.prefix}")
    tools = {name: shutil.which(name) for name in TOOLS}
    if args.executor == "lsf":
        tools.update({name: shutil.which(name) for name in ("bsub", "bjobs")})
    if any(path is None for path in tools.values()):
        raise RuntimeError(f"Missing RNA tools: {[k for k, v in tools.items() if v is None]}")
    report = {"created_utc": datetime.now(timezone.utc).isoformat(), "tools": tools,
              "samples": len(rows), "python": sys.executable, "references": references(args)}
    report["versions"] = {}
    for name, options in (("STAR", ["--version"]), ("fastp", ["--version"]),
                          ("samtools", ["--version"]), ("seqkit", ["version"]),
                          ("htseq-count", ["--version"]), ("fastqc", ["--version"]),
                          ("multiqc", ["--version"])):
        result = execute([tools[name], *options], capture_output=True)
        report["versions"][name] = (result.stdout + result.stderr).strip()
    if "0.20.1" not in report["versions"]["fastp"] or "2.7.6a" not in report["versions"]["STAR"]:
        raise RuntimeError("Pinned fastp 0.20.1 and STAR 2.7.6a are required")
    execute([sys.executable, "-c", "import pandas,pyarrow,numpy,scipy"], capture_output=True)
    code = "for(p in c(" + ",".join(json.dumps(p) for p in PACKAGES) + (
        ")) {if(!requireNamespace(p,quietly=TRUE)) stop(paste('Missing package',p));"
        "cat(p,as.character(packageVersion(p)),'\\n')};stopifnot(capabilities('jpeg'))"
    )
    report["R_packages"] = execute(["Rscript", "-e", code], capture_output=True).stdout
    execute([args.nextflow, "-version"], capture_output=True)
    check_references(references(args), require_index=not args.check_tools)
    report["complete_reference_checked"] = not args.check_tools
    return report


def branch_command(args, branch):
    root = args.outdir
    target = args.qc_reads
    # Full depth is represented by an unreachable target, so the source copies all reads.
    if branch == "full_depth":
        target = 9223372036854775807
    refs = references(args)
    return [
        args.nextflow, "-log", root / "logs_driver" / f"{branch}.nextflow.log",
        "-c", root / "scripts" / "rna_native.config",
        "run", root / "runtime/basej-rnaqc/main.nf",
        "--input_csv", root / "rna_input.csv", "--outputDir", root / branch,
        "--platform", "Illumina", "--architecture", "x86", "--genome", "GRCh38",
        "--genomes_base", HERE.parent / "pipeline_resources",
        "--star_index", refs["star_index"], "--gtf", refs["gtf"],
        "--tx2gene", refs["tx2gene"], "--genebody_ref", refs["genebody"],
        "--ref_fasta", refs["fasta"], "--ref_flat", refs["ref_flat"],
        "--n_reads", target, "--seqkit_sample_seed", args.seed,
        "--strandedness", args.strandedness,
        "--picard_strand", {"no": "NONE", "yes": "FIRST_READ_TRANSCRIPTION_STRAND",
                            "reverse": "SECOND_READ_TRANSCRIPTION_STRAND"}[args.strandedness],
        "--trim_poly_g", str(args.trim_poly_g).lower(),
        "--max_cpus", args.threads, "--max_memory", f"{args.memory_mb}.MB",
        "-work-dir", root / f"work_{branch}",
        "-with-trace", root / "logs_driver" / f"{branch}.trace.tsv",
        "-with-report", root / "logs_driver" / f"{branch}.execution.html",
        "-with-timeline", root / "logs_driver" / f"{branch}.timeline.html",
        "-resume",
    ]


def branches(args):
    return ["qc_2M", "full_depth"] if args.mode == "both" else [
        "qc_2M" if args.mode == "qc" else "full_depth"]


def perform(args):
    root = args.outdir
    os.environ.setdefault("RNA_RUNTIME_DIR", str(RUNTIME))
    (root / "scripts").mkdir(exist_ok=True)
    runtime = prepare(args.source_dir, root / "runtime/basej-rnaqc")
    write_config(runtime, root / "scripts/rna_native.config", args.env_prefix,
                 args.executor, args.queue, args.threads, args.memory_mb)
    fingerprint = {
        "cells_sha256": hashlib.sha256(args.cells.read_bytes()).hexdigest(),
        "references": references(args), "mode": args.mode, "qc_reads": args.qc_reads,
        "seed": args.seed, "strandedness": args.strandedness, "trim_poly_g": args.trim_poly_g,
        "source": source_hashes(args.source_dir),
        "custom_code": {
            str(path.relative_to(HERE)): hashlib.sha256(path.read_bytes()).hexdigest()
            for directory in ("python_lsf_pipeline", "setup", "legacy_nextflow",
                              "extra_methods")
            for path in sorted((HERE / directory).rglob("*.py"))
        },
        "env_manifest": hashlib.sha256((HERE / "setup/environment.yml").read_bytes()).hexdigest(),
    }
    state_path = root / "workflow_state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {"fingerprint": fingerprint}
    if state["fingerprint"] != fingerprint:
        raise ValueError("Inputs/protocol/code changed: use a new output root")

    def persist():
        temporary = state_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(state, indent=2) + "\n")
        temporary.replace(state_path)

    write_source_csv(read_samples(args.cells), root / "rna_input.csv")
    for branch in branches(args):
        status = state.get(branch, {}).get("status")
        if status == "DONE":
            continue
        if status in ("STARTED", "FAILED") and not args.resume:
            raise RuntimeError(f"{branch} requires inspection of existing LSF jobs and --resume")
        state[branch] = {"status": "STARTED", "started_utc": datetime.now(timezone.utc).isoformat()}
        persist()
        try:
            with (root / "logs_driver" / f"{branch}.driver.log").open("a") as log:
                execute(branch_command(args, branch), stdout=log, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            state[branch]["status"] = "FAILED"
            persist()
            raise
        state[branch]["status"] = "DONE"
        persist()
    dashboard(root, branches(args), [row["sample_id"] for row in read_samples(args.cells)])
    state["output_audit"] = {"status": "DONE"}
    persist()


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--env-prefix", type=Path, default=PREFIX)
    parser.add_argument("--source-dir", type=Path, default=SOURCE)
    parser.add_argument("--nextflow", type=Path, default=NEXTFLOW)
    for key, value in DEFAULTS.items():
        parser.add_argument("--" + key.replace("_", "-"), type=Path, default=Path(value))
    parser.add_argument("--executor", choices=("lsf", "local"), default="lsf")
    parser.add_argument("--queue", default="standard")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--memory-mb", type=int, default=64000)
    parser.add_argument("--mode", choices=("both", "qc", "full"), default="both")
    parser.add_argument("--qc-reads", type=int, default=2000000,
                        help="Individual reads, R1+R2, not read pairs")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--strandedness", choices=("no", "yes", "reverse"), default="no")
    parser.add_argument("--trim-poly-g", action="store_true",
                        help="Force poly-G trimming; default preserves fastp auto-detection")
    parser.add_argument("--resume", action="store_true")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--plan", action="store_true")
    mode.add_argument("--preflight", action="store_true", help="Default: check full prerequisites")
    mode.add_argument("--check-tools", action="store_true", help="Check before index preparation")
    mode.add_argument("--run", action="store_true")
    args = parser.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", args.queue):
        parser.error("Unsafe queue name")
    if args.threads < 2 or args.memory_mb < 48000 or args.qc_reads < 2:
        parser.error("Need >=2 CPUs, >=48000 MB and >=2 QC reads")
    for name in ("env_prefix", "source_dir", "nextflow", "outdir", "cells"):
        value = getattr(args, name).resolve()
        if any(c in str(value) for c in ("'", '"', "\n", "\r", " ", ";", "$", "`")):
            parser.error(f"Unsupported shell characters in --{name.replace('_', '-')}")
        setattr(args, name, value)
    for name in DEFAULTS:
        value = getattr(args, name).resolve()
        if any(c in str(value) for c in ("'", '"', "\n", "\r", " ", ";", "$", "`")):
            parser.error(f"Unsupported shell characters in --{name.replace('_', '-')}")
        setattr(args, name, value)
    protected = [
        HERE,
        HERE.parent / "RNA_results",
        HERE.parent / "RNA_results_Expended",
        HERE.parent / "DNA_results",
        HERE.parent / "results",
    ]
    if args.outdir in HERE.parents or any(
            args.outdir == path or path in args.outdir.parents for path in protected):
        parser.error("RNA output must not overwrite code or DNA/previous result trees")
    return args


def main(argv=None):
    args = parse_args(argv)
    if args.plan:
        for branch in branches(args):
            print(shlex.join(map(str, branch_command(args, branch))))
        return
    report = preflight(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "logs_driver").mkdir(exist_ok=True)
    (args.outdir / "tool_preflight.json").write_text(json.dumps(report, indent=2) + "\n")
    print("[OK] Tool check" if args.check_tools else "[OK] Full reference/tool preflight")
    if args.run:
        with (args.outdir / "workflow.lock").open("w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            perform(args)


if __name__ == "__main__":
    main()
