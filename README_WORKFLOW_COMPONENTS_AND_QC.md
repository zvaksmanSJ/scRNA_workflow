# RNA workflow components and QC criteria

This guide compares the BaseJumper-based RNA workflow with additions in this
ResolveOME RNA project. “Our addition” distinguishes changes that preserve the
canonical workflow from optional research analyses. Cohort-specific audit
records and results are intentionally excluded from this source repository.
The optional research runner is implemented but has not been run or validated
against a cohort. Unless stated otherwise, results are descriptive and do not
establish clinical findings.

## How execution flows through the code

The normal Cromwell route and the optional research-add-on route are separate.
The normal route does not launch the optional add-ons.

| Order / area | Entry point or module | What it does |
|---|---|---|
| FASTQ merge plan and execution | [step0_prepare_raw_fastqs.py](./python_lsf_pipeline/step0_prepare_raw_fastqs.py): `plan`, `write_plan`, `merge_cell`, `main` | Groups FASTQs by normalized cell and paired-end mate, creates merge plans, submits/executes merges, and records input/merge provenance. |
| Sample-sheet validation and source manifest | [step1_prepare_samples.py](./python_lsf_pipeline/step1_prepare_samples.py): `read_samples`, `discover`, `write_manifest`, `validate_pairs`, `write_source_csv`, `main` | Validates IDs and paired FASTQs, discovers/records inputs, checks read-file structure, and writes the source manifest used by the BaseJumper adapter. |
| Expanded-cohort validation | [audit_expanded_fastqs.py](./python_lsf_pipeline/audit_expanded_fastqs.py): `audit_one`, `combine`, `main` | Audits each sample and combines the per-sample results into a cohort-level audit; it does not run RNA analysis. |
| Environment/reference setup | [setup.sh](./setup/setup.sh), [step2_prepare_references.py](./setup/step2_prepare_references.py): `check_references`, `create_tx2gene`, `create_genebody_bed`, `create_ref_flat`, `combine_transcript_fastas` | Installs pinned environments and prepares/checks the GRCh38/Ensembl 104 derived reference artifacts, including transcript mapping, gene-body BED, and refFlat files. |
| WDL orchestration | [main.wdl](./main.wdl) | Calls preflight, runs QC-depth tasks, gates full depth on QC audits, scatters Salmon, validates results, then runs analysis and organization. |
| WDL stage adapter tasks | [rna_stages.wdl](./wdl/tasks/rna_stages.wdl): `RnaPreflight`, `RunRnaStage`, `RecordWdlSalmonManifest`, `ValidateSalmon`, `RunFinalAnalysis`, `OrganizeRnaResults` | Defines preflight, stage runner/checker, audit, final-analysis, and result-organization tasks. A stage's audit is the dependency gate for downstream calls. |
| Per-sample WDL Salmon | [rna_salmon.wdl](./wdl/tasks/rna_salmon.wdl) | Quantifies full-depth transcript abundance per sample; WDL scatters this task and gathers outputs for manifest and validation. |
| Stand-alone Salmon submit/validation route | [quantify_transcripts.py](./extra_methods/quantify_transcripts.py): `preflight`, `write_job_manifest`, `submit`, `record_wdl_completions`, `validate_outputs` | Supports plan/submit/validate operations outside the WDL path; gates on full-depth inputs and records LSF submissions/output validation. |
| Legacy QC branch dashboard | [step4_rna_dashboard.py](./extra_methods/step4_rna_dashboard.py): `audit_branch`, `dashboard`, `main` | Audits and presents the two-depth branch outputs and branch-level status; it is distinct from the downstream cohort-analysis dashboard. |
| Final descriptive analysis and dashboard | [run_rna_downstream.py](./extra_methods/run_rna_downstream.py): prerequisite validation, size-factor normalization, Salmon/junction loading, pathway analysis, plotting, `build_dashboard`, `run` | Validates audited prerequisites; generates normalized gene/transcript summaries, Reactome enrichment, junction summaries, figures, QC tables, and the tabbed dashboard. Refresh mode rebuilds the dashboard from completed package outputs. |
| Post-Salmon LSF wrapper | [run_rna_downstream_after_salmon.sh](./extra_methods/run_rna_downstream_after_salmon.sh) | Waits on the submitted per-sample Salmon jobs, validates the cohort outputs, then launches final analysis and organization; it is not the production alignment/counting runner. |
| Results organization | [organize_rna_results.py](./extra_methods/organize_rna_results.py): move/link preflight, duplicate-link checks, move planning, README generation, apply | Plans and checks moves/links, guards against duplicate paths or overwrites, and writes the final folder guide. |
| Optional research add-ons | [run_rna_research_addons.py](./extra_methods/run_rna_research_addons.py): `run_qc_summaries`, `collect_duplication_summary`, `collect_exon_counts`, `collect_psi`, `summarize_fusions`, `build_plan`, `execute_plan` | Emits a command plan, gates on completed full-depth audits, and optionally runs TIN/positional duplication, DEXSeq, SUPPA2, and a separate STAR/Arriba fusion path. This route is not part of `main.wdl`. |
| Sprocket-native prototype | [sprocket/main.wdl](./sprocket/main.wdl), [rna_core.wdl](./sprocket/tasks/rna_core.wdl) | Defines a separate direct-tool prototype. It is not validated as equivalent to the BaseJumper-backed `main.wdl` route. |

The standard WDL call chain is: `RnaPreflight` → QC preprocessing → QC
alignment → QC counts/QC → QC matrix → QC reports → full-depth preprocessing
→ full-depth alignment → full-depth counts/QC → full-depth matrix → full-depth
reports → per-sample Salmon → Salmon manifest/audit → final analysis → result
organization. The user-facing WDL delegates the canonical scientific tasks to
the separately authorized BaseJumper Nextflow source.

## Workflow component comparison

| Workflow component | Function | BaseJumper method / baseline | Our addition |
|---|---|---|---|
| Workflow orchestration | Orders tasks, manages dependencies, and records run results. | BaseJumper's source workflow is a Nextflow process graph. | The general-purpose `main.wdl` provides Cromwell/LSF audit gates around the existing BaseJumper stages, with Salmon scatter and downstream reporting. A separate Sprocket-native WDL prototype exists but is not validated or output-equivalent. |
| References and annotation | Ensures alignment/counting use a compatible genome and gene model. | Uses the genome, STAR index, transcript references, and annotation configured for the BaseJumper run. | Records GRCh38/Ensembl 104 compatibility checks; the private complete STAR index is kept separate from the shared index. References and audit provenance are not modified by the research add-on runner. |
| Sample sheet and FASTQ preparation | Defines each biological sample and creates paired sample-level FASTQs from lane/run inputs. | Uses sample metadata and merges the FASTQ inputs assigned to each sample. | Audits FASTQ identity, pairing, lane/run provenance, gzip integrity, and merged-read totals; retains source files and merge provenance. |
| Read-depth QC branch | Provides a consistent, lower-depth branch for rapid QC. | SeqKit subsampling targets approximately 2,000,000 individual reads (about 1,000,000 paired-end pairs) with seed `12345`; it is not 2,000,000 read pairs. | Keeps the QC-depth branch separate from full depth and records audits so the QC result is not mistaken for the full-depth expression result. |
| Full-depth branch | Processes all available reads for downstream expression and reporting. | Runs the same core preprocessing/alignment/counting method without QC subsampling. | Uses a separate full-depth output branch, with audit gates and independent results; the completed outputs are not overwritten by optional add-ons. |
| Adapter trimming and read QC | Trims adapters and reports read-quality/filter metrics. | Paired-end fastp adapter detection/trimming and its standard read-quality and filtering summaries. | Preserves the source scientific settings; reports the fastp estimates and does not remove duplicates. |
| Genome alignment | Maps reads to the reference genome and provides the canonical alignment used for gene counts. | STAR two-pass alignment with the source RNA settings, followed by the workflow's BAM processing/indexing. | Uses the canonical BaseJumper-compatible alignment unchanged. Fusion discovery, when run, uses a separate chimeric STAR alignment and does not replace the expression BAMs. |
| Gene-level quantification | Produces gene-assigned counts used by the BaseJumper QC and count matrix. | HTSeq uses BAM input, position sorting, unstranded mode, exon features, `gene_id`, union mode, and ignores secondary alignments; the workflow builds the count matrix. | Validates count outputs and adds descriptive normalized-expression summaries. No differential-expression test is reported for the one-condition cohort. |
| Alignment and library QC metrics | Summarizes mapping, read quality, feature composition, RNA metrics, and library characteristics. | BaseJumper reports fastp, STAR, HTSeq, Qualimap, and related report metrics/plots, including mitochondrial and ribosomal fields where available. | Preserves the source metrics in the dashboard, including MT-assigned counts and percentages. MT percent is based on HTSeq-assigned gene counts, not all sequenced or aligned reads. Adds optional RSeQC TIN and positional-duplication summaries. |
| BaseJumper composite QC | Converts five measured criteria into a per-sample 0–5 score. | Applies the five strict thresholds in the QC table below; scores 4–5 are PASS, 3 is Borderline, and 0–2 is FAIL. | Reports the score, each measured criterion, and criterion-specific reason(s) without changing the BaseJumper result. |
| Transcript abundance | Estimates transcript-level abundance. | The BaseJumper source stage set used for this workflow centers on STAR/HTSeq counts and the BaseJumper QC reports; Salmon is not part of the source QC composite. | Adds separately version-matched Salmon transcript TPM/estimated counts and transcript-to-gene summaries. |
| Cohort expression summaries | Visualizes expression and cohort-level patterns. | BaseJumper produces its QC PCA/heatmaps and related reports; these are not a condition-contrast model. | Adds descriptive normalized-count PCA/correlation, variable-gene and transcript summaries, and exploratory Reactome enrichment. It does not infer differential expression or treatment effects. |
| Exon usage and alternative splicing | Quantifies exon features and splice events. | Canonical STAR junction tables provide junction evidence; they are not themselves event-level PSI or differential-splicing calls. | Optional DEXSeq flattened exon-bin counts and SUPPA2 per-event PSI estimates, with figures/tables. No differential test is run; raw event-level inclusion/exclusion counts and PSI uncertainty are not currently provided. |
| Fusion discovery | Identifies candidate chimeric transcripts. | Not part of the standard BaseJumper QC/expression output. | Optional separate STAR chimeric alignment plus Arriba. Reports native confidence/filter labels, breakpoints, split/discordant-read evidence, and a chromosome map. Candidates require research review, not clinical interpretation. |
| Dashboard and deliverables | Presents QC, expression, and run outputs. | BaseJumper provides MultiQC and source QC figures/tables. | Adds an output-tab dashboard with BaseJumper and enhanced QC status/reason tables, MT-read summary, figures, methods, and links to detailed data. Add-on tabs populate only after a successful add-on run and dashboard refresh. |

### Software versions recorded in the workflow configuration

These are the versions pinned by the current setup definitions. The research
add-on environment has not been solved or installed in this session.

| Use | Tools and configured versions |
|---|---|
| Core preprocessing/alignment/counting/reporting | SeqKit 2.13; fastp 0.20.1; STAR 2.7.6a; samtools 1.21; HTSeq 2.0.5; Qualimap 2.2.2d; FastQC 0.12.1; Picard 3.2.0; RSeQC 5.0.4; MultiQC 1.33 |
| Transcript quantification addition | Salmon 1.6.0 |
| Optional RNA research add-ons | Arriba 2.5.1; DEXSeq 1.46.0; SUPPA2 2.3; RSeQC 5.0.4; STAR 2.7.6a |

## BaseJumper QC thresholds

Each of the five criteria contributes one point to the BaseJumper composite
score. The inequalities are strict: a value exactly at a threshold does not
pass that criterion. The proportions are the Qualimap-derived values supplied
to the BaseJumper QC scorer.

| BaseJumper criterion | PASS criterion (one point) | Criterion does not pass | Measurement note |
|---|---:|---:|---|
| Mappability proportion | `> 0.70` (>70%) | `≤ 0.70` | Qualimap-derived mappability proportion. |
| Exonic proportion | `> 0.70` (>70%) | `≤ 0.70` | Calculated from Qualimap raw feature counts; the genomic denominator includes exonic, intronic, intergenic, and overlapping-exon counts. This can differ from a displayed regional percentage. |
| Intergenic proportion | `< 0.10` (<10%) | `≥ 0.10` | Qualimap-derived intergenic proportion. |
| Mitochondrial assigned-count proportion | `< 0.10` (<10%) | `≥ 0.10` | Counts assigned to mitochondrial genes divided by all HTSeq-assigned gene counts. It is not the fraction of all FASTQ reads. |
| Detected protein-coding genes | `> 500` | `≤ 500` | Number of detected protein-coding genes. |

| Composite score | BaseJumper status |
|---:|---|
| 4–5 | **PASS** |
| 3 | **Borderline** |
| 0–2 | **FAIL** |
| Missing/invalid score | **FAIL** in the adapted status reporting |

A failed individual criterion is shown as a reason, including the observed
value or a missing-metric note. The five-point score and the BaseJumper status
are retained as the authoritative BaseJumper label.

## Enhanced QC thresholds added by this project

The enhanced label combines (but does not replace) the BaseJumper status with a
cohort-relative RSeQC transcript-integrity check. It uses the median TIN per
sample and a modified robust z-score:

```text
robust_z = 0.6745 × (sample_median_TIN − cohort_median_TIN) / MAD
MAD = median absolute deviation of sample median TIN values
```

There is **no asserted universal absolute TIN cutoff**. These are research
outlier flags chosen to complement the BaseJumper criteria; they are not
clinical thresholds.

| Enhanced-QC condition | Enhanced status | Reason/reporting behavior |
|---|---|---|
| BaseJumper status is FAIL | **FAIL** | BaseJumper failure is retained regardless of TIN. |
| TIN is missing or invalid | **FAIL** | Required enhanced metric could not be evaluated. |
| Fewer than 8 samples and BaseJumper is not FAIL | **Borderline** | Cohort-relative TIN scoring is not considered sufficiently supported; the sample is retained. |
| TIN file/value missing, malformed, or non-finite | **Add-on validation failure** | The runner stops with an explicit error rather than publishing a success-shaped per-cell status table. |
| At least 8 samples; robust z `≤ -3` | **FAIL** | Severe low-TIN cohort outlier. |
| At least 8 samples; robust z `> -3` and `≤ -2` | **Borderline** | Low-TIN cohort outlier. |
| At least 8 samples; robust z `> -2`, and BaseJumper is PASS | **PASS** | No low-TIN outlier flag. |
| At least 8 samples; robust z `> -2`, and BaseJumper is Borderline | **Borderline** | BaseJumper Borderline status is retained. |

If the cohort MAD is zero, the implementation represents a lower-than-median
TIN as a negative-infinite robust score (severe low outlier) and a
higher-than-median TIN as positive infinity (not a low outlier). No sample is
automatically excluded under either QC system.

### Other add-on evidence is not a cell PASS/FAIL rule

| Added measurement | Reporting rule | Interpretation |
|---|---|---|
| RSeQC positional duplication | Report the position-occurrence estimate; no pass/fail cutoff. | Non-UMI reads at the same mapping position are not necessarily PCR duplicates or distinct molecule counts. No deduplication is performed. |
| DEXSeq exon-bin counts / SUPPA2 PSI | Report counts and descriptive PSI; no pass/fail cutoff or differential test. | These are exploratory feature/event summaries, not proof of differential exon usage or splicing. |
| Arriba fusion candidates | Preserve Arriba confidence/filter fields and read-support counts. The dashboard tiers use Arriba's native labels; an unfiltered candidate must also have at least one reported supporting read to be marked HIGH or REVIEW. | No additional minimum read-count cutoff is invented. All candidates require review; none are clinical calls. |

## Scope and status

- BaseJumper QC status and enhanced QC status are reported separately; neither
  status automatically removes a cell from the analysis.
- RNA–DNA sample identity, RNA variant interpretation, and clinical validation
  remain deferred.
- The optional research-add-on runner is stand-alone, sequential, and not yet
  wired as scatter-based Sprocket/WDL stages. Its environment, tool outputs,
  and synthetic fusion/splicing behavior still require validation.
- For setup, invocation, and gating details, see
  [extra_methods/README.md](./extra_methods/README.md). For full scientific
  rationale and cohort history, see
  [README_RFC_RNA_Workflow.md](./README_RFC_RNA_Workflow.md).

## GitHub source-package boundaries

The intended repository should contain the workflow-owned source, WDL/task
definitions, setup recipes, generic example inputs, tests, and documentation.
The setup scripts are responsible for obtaining or building software and
references at installation time; those generated/downloaded assets should not
be committed.

| Include as source | Exclude from the repository |
|---|---|
| Workflow WDLs, task/helper scripts, source-owned Python/R/shell code, environment YAML definitions, generic examples, tests, and READMEs. | Installed environments (`wdl_runtime/`), downloaded/build references (`references/`), generated run outputs (`RNA_results*/`), tool caches (`.nextflow/`), logs, and generated environment snapshots. |
| Setup scripts that download pinned public inputs or build indexes, plus their checksums and usage documentation. | The downloaded Reactome GMT archive is not committed and must be provisioned separately with the configured checksum; the St. Jude logo is excluded pending redistribution permission. |
| The BaseJumper adapter interfaces and descriptions of required external dependencies. | `vendor/basej-rnaqc/`, its bundled input/test data, and other BaseJumper source until redistribution rights are confirmed; also exclude duplicate `Backup_code/` and private legacy copies. |
| Synthetic sample-sheet templates with fictional IDs and placeholder paths. | `data/rna_cells.local.tsv`, reference manifests with institutional absolute paths, production FASTQs/results, cluster logs, and any credentials or key material. |

The directory `.gitignore` encodes the local-data, generated-asset, and
credential exclusions above. Before the first push, review the exact staged
file list and scan staged text for institutional paths and secrets; a private
repository is not a substitute for that review.
