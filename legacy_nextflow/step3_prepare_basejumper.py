#!/usr/bin/env python3
"""Prepare a versioned local BaseJumper runtime, preserving the upstream checkout."""
import argparse
import hashlib
import json
import os
import shlex
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
SOURCE = Path(os.environ.get("BASEJ_RNA_SOURCE", HERE / "vendor/basej-rnaqc")).resolve()

EXTRAS = r'''
process RAW_FASTQC {
    tag "${sample_name}"
    publishDir { "${params.outputDir}/extras/raw_fastqc/${sample_name}" }, mode: 'copy'
    input:
    tuple val(sample_name), path(reads)
    output:
    path("*_fastqc.{html,zip}")
    script:
    """
    fastqc -t ${task.cpus} -o . ${reads[0]} ${reads[1]}
    """
}
process TRIMMED_FASTQC {
    tag "${sample_name}"
    publishDir { "${params.outputDir}/extras/trimmed_fastqc/${sample_name}" }, mode: 'copy'
    input:
    tuple val(sample_name), path(reads)
    output:
    path("*_fastqc.{html,zip}")
    script:
    """
    fastqc -t ${task.cpus} -o . ${reads[0]} ${reads[1]}
    """
}
process PICARD_RNA {
    tag "${sample_name}"
    publishDir { "${params.outputDir}/extras/picard/${sample_name}" }, mode: 'copy'
    input:
    tuple val(sample_name), path(bam), path(bai)
    path(ref_flat)
    output:
    path("*.rna_metrics")
    path("*.flagstat.txt")
    script:
    """
    samtools quickcheck -v ${bam}
    samtools flagstat ${bam} > ${sample_name}.flagstat.txt
    picard -Xmx8g CollectRnaSeqMetrics I=${bam} O=${sample_name}.rna_metrics \\
      REF_FLAT=${ref_flat} STRAND_SPECIFICITY=${params.picard_strand}
    """
}
process GENE_BODY_COVERAGE {
    tag "${sample_name}"
    publishDir { "${params.outputDir}/extras/gene_body/${sample_name}" }, mode: 'copy'
    input:
    tuple val(sample_name), path(bam), path(bai)
    path(bed)
    output:
    path("gene_body.geneBodyCoverage.txt")
    path("gene_body.geneBodyCoverage.curves.pdf")
    script:
    """
    geneBody_coverage.py -r ${bed} -i ${bam} -o gene_body
    """
}
'''


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"Upstream compatibility check failed: expected one {old[:80]!r}")
    return text.replace(old, new, 1)


def source_hashes(source):
    return {str(path.relative_to(source)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(Path(source).rglob("*"))
            if path.is_file() and path.suffix in (".nf", ".config", ".R", ".yaml")}


def prepare(source, runtime):
    source, runtime = Path(source).resolve(), Path(runtime).resolve()
    if source == runtime or source in runtime.parents or runtime in source.parents:
        raise ValueError("Runtime and source must be separate directory trees")
    manifest = {"upstream": str(source), "sha256": source_hashes(source),
                "adapter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if runtime.exists():
        recorded = json.loads((runtime / "runtime_manifest.json").read_text())
        if recorded != manifest:
            raise ValueError("Runtime source/adapter changed: use a new run directory")
        return runtime
    main = (source / "main.nf").read_text()
    scripts = runtime / "scripts"
    main = main.replace("/usr/local/bin/", str(scripts) + "/")
    main = replace_once(
        main, '                        count = int(parts[1]) if parts[1].isdigit() else 0',
        '                        count = int(parts[-1])\n'
        '                        if count < 0:\n'
        '                            raise ValueError(f"Negative HTSeq count for {gene_id}")')
    main = replace_once(
        main, '            print(f"Warning: Could not parse HTSeq for {sample}: {e}")',
        '            raise ValueError(f"Could not parse HTSeq for {sample}: {e}") from e')
    main = replace_once(main, "--detect_adapter_for_pe \\\\",
                        "--detect_adapter_for_pe ${params.trim_poly_g ? '--trim_poly_g' : ''} \\\\")
    main = replace_once(main, "-f bam -r pos -s no -t exon",
                        "-f bam -r pos -s ${params.strandedness} -t exon")
    main = replace_once(main, "--runThreadN ${task.cpus} \\\\",
                        "--runThreadN ${task.cpus} \\\\\n"
                        "        --limitBAMsortRAM ${Math.min(8000000000L, task.memory.toBytes().intdiv(4))} \\\\")
    main = replace_once(
        main, "    # Count reads per mate in parallel\n"
        "    ( zcat '${r1}' | wc -l | awk '{print int(\\$1/4)}' > read1.txt ) &\n"
        "    ${count_r2}\n    wait",
        "    # Validate gzip CRC, record structure, mate counts and names.\n"
        f"    pair_count=\\$(python {HERE / 'python_lsf_pipeline/step1_prepare_samples.py'} "
        "--validate-pair '${r1}' '${r2}')\n"
        "    echo \\$pair_count > read1.txt\n    echo \\$pair_count > read2.txt")
    main = replace_once(
        main,
        "export FINAL_READS=\\$(zcat '${sample_name}_subsampled_R1.fastq.gz' | "
        "wc -l | awk '{print int(\\$1/4)}')",
        f"export FINAL_READS=\\$((2 * \\$(python {HERE / 'python_lsf_pipeline/step1_prepare_samples.py'} "
        "--validate-pair '${sample_name}_subsampled_R1.fastq.gz' "
        "'${sample_name}_subsampled_R2.fastq.gz')))")
    main = replace_once(
        main, "    // FASTQ preprocessing pipeline\n",
        "    RAW_FASTQC(ch_fastq_input)\n\n    // FASTQ preprocessing pipeline\n")
    main = replace_once(
        main, "        FASTP_TRIM(SEQKIT_SAMPLE.out.reads)\n",
        "        FASTP_TRIM(SEQKIT_SAMPLE.out.reads)\n"
        "        TRIMMED_FASTQC(FASTP_TRIM.out.reads)\n")
    main = replace_once(
        main, "    // Quantification and QC\n",
        "    PICARD_RNA(SAMTOOLS_INDEX_FILTER.out.bam_bai, "
        "file(params.ref_flat, checkIfExists: true))\n"
        "    GENE_BODY_COVERAGE(SAMTOOLS_INDEX_FILTER.out.bam_bai, "
        "file(params.genebody_ref, checkIfExists: true))\n\n    // Quantification and QC\n")
    main = replace_once(
        main, "        --java-mem-size=${task.memory.toGiga()}G",
        "        --java-mem-size=${Math.max(1, task.memory.toGiga() - 8)}G")
    main = replace_once(
        main, '        path "bam/${params.workspace}/rna/tool=star"',
        '        path { sample ->\n'
        '            sample.bam >> "bam/${params.workspace}/rna/tool=star"\n'
        '            sample.bai >> "bam/${params.workspace}/rna/tool=star"\n'
        '        }')
    offset = main.index("workflow {\n")
    main = main[:offset] + EXTRAS + "\n" + main[offset:]
    shutil.copytree(source, runtime, ignore=shutil.ignore_patterns(".git", "work", ".nextflow"))
    scripts.mkdir()
    for directory in ("custom_r_nf_rnaseq", "custom_r_qcplots"):
        for path in (source / "container" / directory / "scripts").glob("*.R"):
            contents = path.read_text().replace("/usr/local/bin/", str(scripts) + "/")
            (scripts / path.name).write_text(contents)
    (runtime / "main.nf").write_text(main)
    (runtime / "runtime_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return runtime


def write_config(runtime, output, prefix, executor, queue, threads, memory_mb,
                 conda_base=None):
    def quote(value):
        return "'" + str(value).replace("\\", "\\\\").replace("'", "\\'") + "'"
    if threads < 2 or memory_mb < 48000:
        raise ValueError("At least 2 CPUs and 48000 MB required")
    if conda_base is None:
        conda_base = os.environ.get("CONDA_BASE")
    runtime_dir = Path(os.environ.get("RNA_RUNTIME_DIR", HERE / "wdl_runtime"))
    if conda_base is None and (runtime_dir / "conda_base.txt").is_file():
        conda_base = (runtime_dir / "conda_base.txt").read_text().strip()
    if conda_base is None:
        conda_base = Path(prefix).resolve().parent.parent
    conda_script = Path(conda_base).resolve() / "etc/profile.d/conda.sh"
    if not conda_script.is_file():
        raise FileNotFoundError(f"Conda activation script is missing: {conda_script}")
    activate = (
        f"source {shlex.quote(str(conda_script))}; "
        f"conda activate {shlex.quote(str(Path(prefix).resolve()))}; "
        f"export R_LIBS_USER={shlex.quote(str(Path(prefix).resolve() / 'lib/R/library'))}; "
        "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1"
    )
    text = f"""
includeConfig {quote(Path(runtime) / 'nextflow.config')}
docker.enabled = false
singularity.enabled = false
conda.enabled = false
workflow.output.mode = 'copy'
process {{
    withName: '.*' {{
        executor = {quote(executor)}
        queue = {quote(queue)}
        container = null
        beforeScript = {quote(activate)}
        errorStrategy = 'terminate'
        maxRetries = 0
        cpus = 2
        memory = '16 GB'
        time = '48 h'
        resourceLimits = [cpus: {threads}, memory: '{memory_mb} MB']
    }}
    withName: 'STAR_ALIGN|QUALIMAP_BAMRNA' {{
        cpus = {threads}
        memory = '{memory_mb} MB'
    }}
    withName: 'RAW_FASTQC|TRIMMED_FASTQC' {{
        memory = '4 GB'
    }}
    withName: 'FASTP_TRIM' {{
        publishDir = [path: {{ "${{params.outputDir}}/tool_outputs/fastp" }},
                      mode: 'copy', pattern: '*fastp*']
    }}
    withName: 'STAR_ALIGN' {{
        publishDir = [path: {{ "${{params.outputDir}}/tool_outputs/star" }},
                      mode: 'copy', pattern: '*Log.final.out']
    }}
    withName: 'HTSEQ_COUNTS' {{
        publishDir = [path: {{ "${{params.outputDir}}/tool_outputs/htseq" }},
                      mode: 'copy']
    }}
    withName: 'QUALIMAP_BAMRNA' {{
        publishDir = [path: {{ "${{params.outputDir}}/tool_outputs/qualimap" }},
                      mode: 'copy']
    }}
}}
executor.perTaskReserve = true
params.architecture = 'x86'
params.genome = 'GRCh38'
params.workspace = 'resolveome'
params.dataset_id = 'SY5Y_RNA'
params.workflow_id = 'RNA_WORKFLOW'
"""
    # Exact selectors also override upstream exact selectors that ignore failures.
    import re
    for name in re.findall(r"^process (\w+) \{", (Path(runtime) / "main.nf").read_text(), re.M):
        compute = name in ("STAR_ALIGN", "QUALIMAP_BAMRNA")
        cpus = threads if compute else 2
        memory = memory_mb if compute else (4000 if "FASTQC" in name else 16000)
        text += f"""
process {{
    withName: {quote(name)} {{
        executor = {quote(executor)}
        queue = {quote(queue)}
        container = null
        beforeScript = {quote(activate)}
        errorStrategy = 'terminate'
        maxRetries = 0
        cpus = {cpus}
        memory = '{memory} MB'
        time = '48 h'
        resourceLimits = [cpus: {threads}, memory: '{memory_mb} MB']
    }}
}}
"""
    Path(output).write_text(text)
    return Path(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=SOURCE)
    parser.add_argument("--runtime-dir", type=Path, required=True)
    args = parser.parse_args()
    print(prepare(args.source_dir, args.runtime_dir))


if __name__ == "__main__":
    main()
