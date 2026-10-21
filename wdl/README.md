# RNA workflow: WDL/Cromwell orchestration

`main.wdl` orchestrates the QC-depth and full-depth branches, gates downstream
stages on audited outputs, scatters Salmon quantification by sample, and runs
the final analyses and organizer. The scientific BaseJumper stages remain
Nextflow tasks behind adapters; this is not a pure WDL implementation.

The workflow needs a separate, authorized BaseJumper source checkout. It also
needs the software and reference runtime from `setup/setup.sh`, Cromwell with
Java, and an LSF cluster configured with the custom runtime attributes in
`cromwell_lsf.conf`. The `memory_mb`, `lsf_queue`, and `lsf_walltime` runtime
attributes are Cromwell-specific and are not recognized by generic WDL
engines. `miniwdl check` reports those backend keys as unknown; that is not an
LSF execution test.

## Prepare inputs and launch

Set absolute paths to the installed Cromwell jar, a sample sheet, and the
authorized BaseJumper source. Prepare a **new** results directory; the builder refuses to reuse an existing
one and symlinks FASTQs without copying them.

```bash
RESULTS="/path/to/new/rna_run"
SAMPLES="/path/to/rna_cells.tsv"
BASEJ_SOURCE="/authorized/path/to/basej-rnaqc"
CROMWELL_JAR="/path/to/cromwell.jar"
RUNTIME="/path/to/wdl_runtime"

bash setup/setup.sh --runtime-dir "$RUNTIME"
"$RUNTIME/env/resolveome-rna/bin/python" wdl/scripts/build_run_inputs.py \
  --workflow cromwell \
  --sample-sheet "$SAMPLES" \
  --runtime-config "$RUNTIME/runtime_config.json" \
  --results-dir "$RESULTS" \
  --basej-source-dir "$BASEJ_SOURCE"

java -Dconfig.file=wdl/cromwell_lsf.conf \
  -jar "$CROMWELL_JAR" run main.wdl \
  -i "$RESULTS/wdl_inputs.json" \
  -m "$RESULTS/cromwell_metadata.json"
```

The sample sheet must have exactly `sample_id`, `r1`, `r2`, and `group`
columns, unique safe sample IDs, existing paired FASTQs, and at least five
samples for the source cohort plots. See
`examples/rna_cells.example.tsv`. The runtime config must point to references
that exist and have built STAR/Salmon indexes. The builder creates
`$RESULTS/wdl_inputs.json` and prepares the isolated run directory before
writing it. It refuses to reuse an existing results directory or overwrite an
inputs JSON. `nextflow_threads` and
`nextflow_memory_mb` are adjustable WDL inputs; by default, per-cell
Nextflow/LSF processes request 8 CPUs and 64 GB. Set the Cromwell queue with
the workflow `queue` input.

To generate both Cromwell and Sprocket inputs from one sample sheet, add
`--workflow both` and `--sprocket-container docker://registry/image:tag` to
the builder command. It writes `sprocket_inputs.json` alongside the Cromwell
JSON. For Sprocket alone, choose a new output path and provide an accessible
container image:

```bash
"$RUNTIME/env/resolveome-rna/bin/python" wdl/scripts/build_run_inputs.py \
  --workflow sprocket \
  --sample-sheet "$SAMPLES" \
  --runtime-config "$RUNTIME/runtime_config.json" \
  --sprocket-container docker://registry/image:tag \
  --sprocket-output /path/to/new/sprocket_inputs.json
```

The builder only prepares run directories and JSON; it never submits or starts
a workflow. Sprocket remains a prototype and is not production-validated.

## Gates and results

The WDL checks its run tree and references before starting. It then runs QC
preprocessing, alignment, counts/QC, matrix, and reports; only after those
audits pass does it start full-depth processing. Salmon quantification is
scattered over the cohort. Salmon and full-depth report audits gate final
analysis and results organization. Run products, `final_results/`, logs, and
Cromwell metadata are retained under the new results directory. Never point
this workflow at a completed production tree.

Each stage adapter loads the validated paths from
`step2_reference_audit.json` in the run directory. This keeps later alignment,
counting, and matrix tasks on the same configured shared-reference bundle
checked by preflight; portable setup defaults are used only when no run audit
is configured.

## Validation status

Static WDL parsing is available with `miniwdl check main.wdl`, but its generic
checker warns about the expected Cromwell-specific runtime keys. The pipeline
is not certified on LSF or end-to-end from a clean checkout. A 47-cell
production attempt passed preflight and QC preprocessing but failed at QC
alignment because its stage adapter used portable reference defaults instead
of the configured shared-reference paths. The adapter now reads the validated
reference audit, but this fix still requires a new workflow attempt. The
setup/index build, full WDL path, and remaining stage adapters still require
cluster validation. Do not treat this README as production sign-off.
