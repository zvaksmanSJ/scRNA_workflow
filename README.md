# workflow_RNA/

The ResolveOME single-cell RNA workflow. This folder mirrors the compact
entry-point layout of `workflow_DNA/`; supporting code, setup, task modules,
tests, and historical implementations live in subfolders.

For the complete method, audit history, and known limitations, see
[README_RFC_RNA_Workflow.md](./README_RFC_RNA_Workflow.md).
For a component-by-component comparison of BaseJumper methods and project
additions, plus the QC scoring thresholds, see
[README_WORKFLOW_COMPONENTS_AND_QC.md](./README_WORKFLOW_COMPONENTS_AND_QC.md).

## Folder layout

| Path | What's in it |
|---|---|
| **`main.wdl`** | General-purpose WDL entry point. Uses Cromwell/LSF for orchestration and invokes the established BaseJumper Nextflow stages through audit-gated adapters. It is not a tool-by-tool WDL reimplementation. |
| **`README_RFC_RNA_Workflow.md`** | Methods, execution paths, expected outputs, validation status, and limitations. |
| **`README_WORKFLOW_COMPONENTS_AND_QC.md`** | Component-by-component workflow map, code entry points, BaseJumper/add-on comparison, QC thresholds, and release boundaries. |
| `wdl/` | WDL task modules, helper scripts, Cromwell configurations, and launch instructions. |
| `sprocket/` | Separate Sprocket-native prototype. It is not yet validated or output-equivalent to the BaseJumper-based workflow. |
| `setup/` | Portable conda setup and public Ensembl 104 reference/index preparation. |
| `wdl_runtime/` | Local-only installation target for tools, package downloads, and references; ignored by Git. |
| `python_lsf_pipeline/` | FASTQ preparation, merging, and integrity-audit utilities, including LSF workers. |
| `legacy_nextflow/` | Nextflow/BaseJumper adapters and historical launch helpers used by the WDL stages. |
| `extra_methods/` | Salmon quantification, downstream analyses, dashboard, and results organization. |
| `analysis_resources/` | Publicly distributable pathway resources and local-only branding assets. |
| `examples/` | Generic example sample sheet. |
| `tests/` | Unit and opt-in synthetic smoke tests. |
| `docs/` | Workflow diagrams and design aids. |
| `Backup_code/`, `vendor/`, `references/`, `logs/`, `data/` | Local history, private/unlicensed source, large references, logs, and cohort-specific files. These are excluded from a public source release where appropriate. |

## Which workflow should I use?

- **Cromwell/LSF candidate:** [main.wdl](./main.wdl) provides staged QC and
  full-depth analysis, audit gates, Salmon scatter/gather, and final result
  organization. It still depends on an authorized external BaseJumper source
  checkout, Cromwell, LSF, and a suitable compute environment. A fresh
  end-to-end run has not been certified from a clean checkout.
- **Sprocket prototype:** See [sprocket/README.md](./sprocket/README.md).
  Its tool/container and output-equivalence validation is incomplete.
- **Historical implementation:** See [legacy_nextflow/README.md](./legacy_nextflow/README.md).
  Do not point historical launchers at completed results.

## Getting started (setup only)

On Linux, install conda or mamba, then run:

```bash
bash setup/setup.sh --skip-references
# Or omit --skip-references to download Ensembl 104 and build STAR/Salmon indexes.
```

The reference build requires substantial disk space and at least 48 GB RAM
(64 GB recommended). Setup downloads public Reactome and Ensembl resources.
See
[setup/README.md](./setup/README.md) for dependencies and expected paths.

To prepare a new run, use a four-column TSV with `sample_id`, `r1`, `r2`, and
`group` columns (see
[examples/rna_cells.example.tsv](./examples/rna_cells.example.tsv)). Do not
reuse or overwrite completed production results. Exact WDL preparation and
launch instructions are in [wdl/README.md](./wdl/README.md).

## Quick Start

These examples are for a **new, isolated run**. They do not start or resume any
existing cohort run. Replace every placeholder, review the command plan and
resource requirements, and verify that each results directory is new.

### 1. Set up software and references

Install conda or mamba on Linux. The full setup downloads the pinned software
dependencies, the checksum-verified Reactome resource, public Ensembl 104
references, and builds the STAR and Salmon indexes. Building the STAR index
requires at least 48 GB RAM (64 GB recommended).

```bash
export RNA_RUNTIME_DIR="$PWD/wdl_runtime"
bash setup/setup.sh \
  --runtime-dir "$RNA_RUNTIME_DIR" \
  --threads 8 \
  --memory-mb 64000
```

To install software without downloading/building Ensembl indexes, add
`--skip-references`; that mode is not sufficient to run the workflows.
Prepare a sample sheet with the exact `sample_id`, `r1`, `r2`, and `group`
columns and absolute paths to paired FASTQs.

### 2. Sprocket (direct-tool prototype)

Sprocket must be installed separately with its LSF/Apptainer backend. Build
and publish the task image as described in
[`sprocket/README.md`](./sprocket/README.md). Copy the generic example inputs
and replace its sample arrays, paths, reference indexes, and image with values
for your run:

```bash
cp sprocket/example_inputs.json /path/to/rna_inputs.json
# Edit /path/to/rna_inputs.json; use all samples in identical array order.
sprocket run \
  --config sprocket/sprocket.toml \
  sprocket/main.wdl \
  @/path/to/rna_inputs.json \
  -o /path/to/new_sprocket_run
```

This is a prototype: it has only passed static WDL validation and has not been
run on LSF or certified equivalent to the BaseJumper outputs. Do not use it
for production without completing the validation steps in its README.

### 3. Cromwell/LSF (audit-gated WDL using Nextflow stages)

This route runs the audit-gated WDL entry point. The WDL invokes the
scientific BaseJumper stages as Nextflow tasks, so obtain an authorized
BaseJumper source checkout separately. Cromwell, Java, and a configured LSF
cluster are also required. The helper refuses to reuse an existing results
directory.

```bash
RUNTIME="/path/to/wdl_runtime"
RESULTS="/path/to/new_rna_run"
SAMPLES="/path/to/rna_cells.tsv"
BASEJ_SOURCE="/authorized/path/to/basej-rnaqc"
CROMWELL_JAR="/path/to/cromwell.jar"
export RNA_RUNTIME_DIR="$RUNTIME"

bash setup/setup.sh --runtime-dir "$RUNTIME" --threads 8 --memory-mb 64000
"$RUNTIME/env/resolveome-rna/bin/python" wdl/scripts/prepare_wdl_run.py \
  --results-dir "$RESULTS" --sample-sheet "$SAMPLES"
"$RUNTIME/env/resolveome-rna/bin/python" wdl/scripts/build_inputs_json.py \
  --results-dir "$RESULTS" --basej-source-dir "$BASEJ_SOURCE" \
  --runtime-config "$RUNTIME/runtime_config.json" \
  --output "$RESULTS/wdl_inputs.json"

java -Dconfig.file=wdl/cromwell_lsf.conf \
  -jar "$CROMWELL_JAR" run main.wdl \
  -i "$RESULTS/wdl_inputs.json" \
  -m "$RESULTS/cromwell_metadata.json"
```

See [`wdl/README.md`](./wdl/README.md) for sample-sheet requirements,
configuration, gates, and output locations.

### 4. Direct Nextflow driver (advanced)

This legacy route calls Nextflow directly and bypasses Cromwell orchestration.
It still requires the separately authorized BaseJumper source, a prepared
runtime, references, and LSF. First inspect the command plan; only submit after
checking all paths and selecting a new results directory:

```bash
RUNTIME="/path/to/wdl_runtime"
RESULTS="/path/to/new_nextflow_run"
SAMPLES="/path/to/rna_cells.tsv"
BASEJ_SOURCE="/authorized/path/to/basej-rnaqc"

RNA_RUNTIME_DIR="$RUNTIME" \
  bash legacy_nextflow/run_all_resolveome_rna.sh \
  --plan --cells "$SAMPLES" --outdir "$RESULTS" \
  --source-dir "$BASEJ_SOURCE"

# Explicit LSF submission; the driver performs its own prerequisite checks.
RNA_RUNTIME_DIR="$RUNTIME" \
  bash legacy_nextflow/run_all_resolveome_rna.sh \
  --submit --cells "$SAMPLES" --outdir "$RESULTS" \
  --source-dir "$BASEJ_SOURCE"
```

Use Cromwell/LSF when you want the documented WDL audit gates and stage
orchestration. Use the direct Nextflow route only when you specifically need
the standalone driver and understand its separate execution path.

## Release readiness and distribution

The workflow is published at the public repository
[`zvaksmanSJ/scRNA_workflow`](https://github.com/zvaksmanSJ/scRNA_workflow).
The source-only release excludes cohort data and manifests, local results,
reference indexes, runtime installations, downloaded assets, institutional
branding, backup copies, and the private BaseJumper source.
The private BaseJumper source has no confirmed redistribution license and is
excluded; users must obtain an authorized source separately. The St. Jude
logo is also excluded pending permission, with a text-only dashboard fallback.
No project-level open-source license has been selected or added.
Public visibility does not grant permission to reuse or redistribute the
project code; obtain institutional approval and add an appropriate license
before inviting reuse.
The conda/reference setup, WDL tasks, LSF configuration, and clean-checkout
execution still require validation on a suitable system before claiming
portable production readiness. This source package does not start or modify
any existing production runs.

## RNA-only research extensions

The current validated analysis remains descriptive gene/transcript expression
and QC. A separate, optional RNA-only research-add-on runner is implemented
for expanded QC, exon-bin counts, descriptive PSI, and a fusion-specific
STAR/Arriba branch. It writes to a separate results subfolder and does not
alter BaseJumper BAMs or QC. The dashboard has output-type tabs and displays
the add-on QC statuses/reasons, figures, and evidence when a completed
add-on run is present. These additions have not been run against a cohort or
validated with synthetic fusion/splicing controls, are not yet scatter-based
Sprocket/WDL stages, and must not be described as production-validated. See
`extra_methods/README.md` for the gating, setup, and planned execution steps.
They do not require or consume DNA data.

RNA–DNA identity checks and DNA-assisted RNA variant interpretation are
explicitly deferred until matched DNA/RNA samples are available. Clinical
validation and clinical interpretation are out of scope for this workflow.
