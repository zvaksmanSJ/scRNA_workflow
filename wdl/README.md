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
authorized BaseJumper source. Prepare a **new** results directory; the helper
refuses to reuse an existing one and symlinks FASTQs without copying them.

```bash
RESULTS="/path/to/new/rna_run"
SAMPLES="/path/to/rna_cells.tsv"
BASEJ_SOURCE="/authorized/path/to/basej-rnaqc"
CROMWELL_JAR="/path/to/cromwell.jar"

bash setup/setup.sh --runtime-dir /path/to/wdl_runtime
python wdl/scripts/prepare_wdl_run.py \
  --results-dir "$RESULTS" --sample-sheet "$SAMPLES"
python wdl/scripts/build_inputs_json.py \
  --results-dir "$RESULTS" --basej-source-dir "$BASEJ_SOURCE" \
  --runtime-config /path/to/wdl_runtime/runtime_config.json \
  --output "$RESULTS/wdl_inputs.json"

java -Dconfig.file=wdl/cromwell_lsf.conf \
  -jar "$CROMWELL_JAR" run main.wdl \
  -i "$RESULTS/wdl_inputs.json" \
  -m "$RESULTS/cromwell_metadata.json"
```

The sample sheet must have exactly `sample_id`, `r1`, `r2`, and `group`
columns, unique safe sample IDs, existing paired FASTQs, and at least five
samples for the source cohort plots. See
`examples/rna_cells.example.tsv`. The runtime config must point to references
that exist and have built STAR/Salmon indexes. `nextflow_threads` and
`nextflow_memory_mb` are adjustable WDL inputs; by default, per-cell
Nextflow/LSF processes request 8 CPUs and 64 GB. Set the Cromwell queue with
the workflow `queue` input.

## Gates and results

The WDL checks its run tree and references before starting. It then runs QC
preprocessing, alignment, counts/QC, matrix, and reports; only after those
audits pass does it start full-depth processing. Salmon quantification is
scattered over the cohort. Salmon and full-depth report audits gate final
analysis and results organization. Run products, `final_results/`, logs, and
Cromwell metadata are retained under the new results directory. Never point
this workflow at a completed production tree.

## Validation status

Static WDL parsing is available with `miniwdl check main.wdl`, but its generic
checker warns about the expected Cromwell-specific runtime keys. The pipeline
has not been executed end-to-end from a clean checkout or certified on LSF.
The setup/index build, full WDL path, and all stage adapters still require
cluster validation. Do not treat this README as production sign-off.
