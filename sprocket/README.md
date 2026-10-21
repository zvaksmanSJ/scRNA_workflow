# RNA Sprocket port

This directory contains a WDL 1.2 workflow for Sprocket's LSF+Apptainer
backend. The workflow calls task-level tools directly inside a task container;
it does not launch Nextflow or Cromwell. It is separate from the existing
BaseJumper-adapted production runs and does not reuse or overwrite their
results.

**Status: static WDL validation only.** `miniwdl check sprocket/main.wdl`
passes. Sprocket is not installed on this cluster, Docker is unavailable,
and the custom image has not been built. The WDL therefore has not been
executed, its container has not been validated, and its outputs have not been
certified equivalent to the accepted BaseJumper reports. Apptainer is
available, but requires a built/published OCI image or a supported local
image source.

## Layout

```text
sprocket/
  main.wdl                         # QC-depth + full-depth + Salmon workflow
  sprocket.toml                    # LSF/Apptainer backend settings
  example_inputs.json              # one-cell example; not a production manifest
  tasks/rna_core.wdl               # direct SeqKit/fastp/STAR/HTSeq/QC/matrix tasks
  containers/rna-workflow.Dockerfile
```

## Methods in the port

For each sample, the workflow runs these task-level steps:

1. Validate cohort array lengths and sample IDs.
2. FastQC on the original reads, SeqKit sampling, fastp adapter detection,
   and FastQC on the trimmed reads.
3. STAR two-pass alignment with the BaseJumper STAR options, then samtools
   filtering with `-F 256`. Secondary alignments are removed; supplementary
   records are retained. No duplicate removal is performed.
4. HTSeq gene counting with `-f bam -r pos -s no -t exon -i gene_id`,
   `gene_name`, union mode, and secondary-alignment ignore behavior.
5. Qualimap RNA-seq QC, Picard RNA metrics, and RSeQC gene-body coverage.
6. Cohort gene-count matrix, long-form Parquet counts, sample expression
   summaries, and MultiQC report.
7. Full-depth Salmon quantification with Salmon 1.6.0.

The QC branch samples a probabilistic target of 2,000,000 individual reads
(approximately one million read pairs), using seed 12345. The full-depth
branch uses a 64-bit maximum target and therefore copies all reads. The
QC branch's passing matrix/MultiQC gate is an explicit dependency of the
full-depth branch. Salmon tasks depend on the full-depth branch gate.

This port does **not yet reproduce all source BaseJumper publication
artifacts**. In particular, it does not implement the source consensus QC
scoring/report process, every source plot/table/metadata output, or the
separate final-results dashboard and downstream analysis suite. Those are
present in the accepted production results and described in the top-level
RFC. Do not substitute this port's simplified sample metrics for the
BaseJumper consensus score or use them to exclude samples.

## Required setup

1. Install Sprocket with its LSF+Apptainer backend on a host configured to
   submit LSF jobs.
2. Set `image_cache_dir` in `sprocket/sprocket.toml` to a shared writable
   directory. Confirm the `standard` queue name and resource limits against
   the target cluster.
3. Build and publish the custom task image from the workflow root. The Docker
   engine used to build/push this image is not currently available on this
   cluster:

   ```bash
   docker build \
     -f sprocket/containers/rna-workflow.Dockerfile \
     -t ghcr.io/<your-organization>/resolveome-rna:1.0.0 .
   docker push ghcr.io/<your-organization>/resolveome-rna:1.0.0
   ```

   The Dockerfile uses the checked-in `setup/environment.yml`, Salmon 1.6.0,
   and downloads checksum-pinned AMOR, ohchibi, and harrietr source archives at
   build time. Record the final image digest before production use.
4. Create a four-column sample sheet with `sample_id`, `r1`, `r2`, and
   `group` columns. After running `setup/setup.sh` with references enabled,
   use the input builder to fill the sample arrays and all reference paths
   from `runtime_config.json`:

   ```bash
   RUNTIME="/path/to/wdl_runtime"
   "$RUNTIME/env/resolveome-rna/bin/python" wdl/scripts/build_run_inputs.py \
     --workflow sprocket \
     --sample-sheet /path/to/rna_cells.tsv \
     --runtime-config "$RUNTIME/runtime_config.json" \
     --sprocket-container docker://ghcr.io/<your-organization>/resolveome-rna:1.0.0 \
     --sprocket-output /path/to/new/sprocket_inputs.json
   ```

   The builder validates the configured references and preserves sample-sheet
   order. Use a new, empty output directory for the workflow; do not point at
   the completed `RNA_results` tree.

## Example invocation

The following is an example only; it cannot run until Sprocket and the
container prerequisites above are available:

```bash
sprocket run \
  --config sprocket/sprocket.toml \
  sprocket/main.wdl \
  @sprocket/example_inputs.json \
  -o /path/to/new/rna_sprocket_run
```

`example_inputs.json` contains only one example sample and placeholder
paths. It is not a production input manifest.

## Validation before production use

- Run `miniwdl check sprocket/main.wdl` after any WDL edits.
- Build and inspect the container; verify `STAR`, `fastp`, `seqkit`,
  `samtools`, `htseq-count`, `qualimap`, `fastqc`, `picard`, `geneBody_coverage.py`,
  `multiqc`, `salmon`, R, and Python are available at the expected versions.
- Run a tiny, non-production sample through both depth branches on LSF.
- Compare read counts, STAR command/version and mapping metrics, HTSeq
  assigned counts, Qualimap and RSeQC/Picard outputs to the established
  source workflow.
- Implement and validate the missing BaseJumper report/scoring and downstream
  dashboard stages before claiming full-method equivalence.
