# ResolveOME Single-Cell RNA Workflow RFC

**Status:** Draft
**RFC number:** 0000 (placeholder; see [RFC process note](#rfc-process-note))
**Workflow version:** 1.0.0

- [Introduction](#introduction)
- [Motivation](#motivation)
- [Discussion](#discussion)
  - [Full-depth alignment and QC](#full-depth-alignment-and-qc)
  - [Cell selection and interpretation of QC](#cell-selection-and-interpretation-of-qc)
  - [Expression and cohort analysis](#expression-and-cohort-analysis)
  - [Sprocket-native WDL port](#sprocket-native-wdl-port)
  - [WDL/Cromwell orchestration](#wdlcromwell-orchestration)
  - [Reference compatibility and limitations](#reference-compatibility-and-limitations)
- [Specification](#specification)
  - [Dependencies](#dependencies)
  - [Reference files](#reference-files)
  - [Workflow (RNA production pipeline)](#workflow-rna-production-pipeline)
  - [Tool-by-tool commands, in execution order](#tool-by-tool-commands-in-execution-order)
  - [Workflow (general-purpose WDL pipelines)](#workflow-general-purpose-wdl-pipelines)
  - [Output organization](#output-organization)
  - [Execution and failure handling](#execution-and-failure-handling)
- [Appendix](#appendix)
  - [Verification and current implementation status](#verification-and-current-implementation-status)
  - [Known issues and fixes baked into this workflow](#known-issues-and-fixes-baked-into-this-workflow)
  - [RFC process note](#rfc-process-note)

# Introduction

This RFC describes a paired-end, full-length single-cell RNA-sequencing
workflow that preserves the relevant BaseJumper RNA processing and reporting
methods while providing separately audited QC-depth and full-depth branches.
It also documents transcript quantification, cohort-level descriptive
analysis, and two WDL orchestration options with different levels of
completeness.

The workflow has three distinct execution surfaces:

1. **The audited RNA production workflow** prepares samples and references,
   then uses adapted BaseJumper Nextflow process definitions with native
   Linux tools and LSF orchestration. It supports an approximately two-million
   individual-read QC branch and a separate full-depth expression branch.
2. **The legacy WDL/Cromwell workflow** orchestrates the existing
   Nextflow-backed stages behind explicit audit gates, then scatters
   full-depth Salmon quantification and runs downstream analysis and result
   organization. This is an HPC-specific orchestration layer, not a
   tool-by-tool reimplementation of BaseJumper.
3. **The Sprocket-native WDL port** invokes the RNA tools directly in task
   containers. It is an incomplete prototype and is not output-equivalent to
   the audited production workflow.

This local draft is modeled on the St. Jude Cloud RFC process. It is not an
accepted RFC or a submission to that repository. The public source repository
does not include the upstream BaseJumper source checkout because
redistribution permission has not been confirmed.

## Current source-package status

The compact source layout and portable setup scripts are published in the
RNA workflow repository, but this source package is **not**
production-certified from a clean public clone. The Cromwell WDL orchestrates
audited stages that still invoke the BaseJumper Nextflow source; that source
is excluded because its redistribution license has not been confirmed. The
St. Jude logo is also excluded pending permission. The Sprocket port remains
incomplete and unvalidated. See the top-level [README](./README.md) and
[WDL guide](./wdl/README.md) for setup, run prerequisites, and blockers.

The cohort audit and results history below documents the already completed
local analyses; it does not constitute certification of the newly reorganized
source tree or authorize relaunching those analyses. The active 47-cell
workflow status is recorded in the Appendix and is separate from the
completed 16-cell cohort.

## Historical cohort validation status

The local 16-sample cohort completed both the approximately two-million-
individual-read QC branch and an independent full-depth branch. The full-depth
matrix contained 40,019 genes across 16 samples. Its BaseJumper-derived QC
labels were 14 PASS and 2 Borderline; every sample was retained. Full-depth
alignment used the documented quick audit (quick checks for all samples and
deep streaming of three), followed by cohort-wide count, matrix, and report
audits. The separate downstream analysis package passed its recorded audits.
Do not resubmit those completed local jobs.

A separate expanded 31-sample cohort had 62 merged FASTQs and 31 provenance
records; its exhaustive merge audit passed for 2,119,874,132 read pairs. No
RNA analysis was launched from that expanded cohort. The data, audit reports,
FASTQ paths, and cluster job records are local results, not part of this
portable source package.

# Motivation

- **Preserve two scientifically distinct read-depth views.** The QC branch
  uses a reproducible, approximately two-million-individual-read subset for
  efficient quality review. Full-depth expression uses all available merged
  read pairs and is not inferred from the QC subset.
- **Make sample preparation and processing auditable.** The workflow records
  source FASTQs, merge provenance, input manifests, tool and reference checks,
  task traces, LSF job IDs, stage state, and output audits. Failed or missing
  tasks must not be reported as successful stages.
- **Retain BaseJumper RNA methods while exposing adaptations.** Source
  settings and reporting are retained where practical; runtime substitutions,
  parser fixes, resource settings, and extra QC reports are identified
  explicitly rather than being presented as bitwise source equivalence.
- **Keep biological QC decisions visible and non-destructive.** QC labels,
  borderline metrics, and limitations are reported for every cell. Cells
  are not silently dropped based on DNA labels, a single metric, or a
  success-shaped fallback.
- **Separate measurement from optional interpretation.** Gene counts,
  transcript quantification, splice-junction support, descriptive
  comparisons, and optional research extensions are distinguished from
  differential expression, differential splicing, variant calling, and
  clinical conclusions that the available design does not support.
- **Support portability without overstating readiness.** Public references,
  pinned environments, sample-sheet input generation, and WDL definitions
  are documented, while the licensed/unavailable upstream source,
  cluster-specific Cromwell/LSF requirements, and unvalidated Sprocket
  implementation remain explicit prerequisites or limitations.

# Discussion

## Full-depth alignment and QC

The workflow validates paired FASTQ inputs and preserves the original files.
For each cell, the QC-depth branch samples both mates with SeqKit seed
12345 to a probabilistic target of approximately two million individual
reads (approximately one million read pairs), then runs raw FastQC, fastp,
trimmed FastQC, STAR two-pass alignment, BAM processing, HTSeq counts, and
RNA QC/reporting. The independent full-depth branch uses all merged read
pairs, without subsampling, and produces separate preprocessing, alignment,
count, matrix, and report outputs.

fastp uses paired-end adapter detection and source-default automatic poly-G
detection; forced poly-G trimming is an explicit optional flag. RNA duplicate
removal is not performed: these libraries have no UMIs, and identical read
sequences alone do not establish PCR duplicates. STAR uses the source
two-pass Basic configuration and source unmapped/primary/chimeric flags.
samtools filtering with `-F 256` excludes secondary alignments but retains
supplementary alignments. HTSeq counts unstranded exon/gene_id features using
the union mode, with secondary alignments ignored.

The audit gates distinguish technical completeness from biological quality.
FASTQ integrity, mate concordance, task return status, count agreement,
matrix schema, BAM sorting/indexing and required reports are validated. A
cell with a low mapping or exonic fraction remains visible for review rather
than being automatically removed.

## Cell selection and interpretation of QC

This is full-length paired-end RNA sequencing, not a droplet/UMI workflow.
The source-derived RNA consensus score uses the established metrics and
thresholds: mapping, exonic and intergenic fractions, mitochondrial
(`MT-` symbol) fraction, and detected protein-coding genes. The original
source score labels are retained (PASS, Borderline, or FAIL) with their
component values; a label is not an automatic inclusion/exclusion rule.
For the completed 16-cell cohort all cells were retained, including the two
Borderline cells.

Qualimap's displayed feature-region fractions and the source scoring
proportions are not identical measures: the source exonic denominator
includes overlapping-exon counts. This distinction is preserved in reports.
FastQC warnings, fastp filtering/duplication estimates, STAR mapping, Picard
RNA metrics, RSeQC gene-body coverage, and mitochondrial counts are
complementary evidence and should not be collapsed into a single pass/fail
claim.

The WDL and Sprocket prototypes do not yet reproduce all source consensus
QC-reporting behavior. Their presence does not replace the production
workflow's audited RNA score or justify automatic cell exclusion.

## Expression and cohort analysis

The full-depth branch preserves the source count matrix and report products.
Salmon is run separately against the matching transcript index/version and
its outputs are audited before the cohort analysis package is produced.
Downstream outputs include descriptive normalized expression, PCA,
correlations, variable-gene summaries, transcript abundance, gene-aggregated
TPM, and observed STAR splice-junction support.

The completed original cohort contains 16 biological replicates from one
SY5Y condition. Its analyses are descriptive; no condition contrast,
differential-expression, differential-transcript-usage, or
differential-splicing test is scientifically defined. RNA-only research
extensions for enhanced library QC, exon-level expression/splicing, and
fusion discovery are separately identified as optional and not yet cohort
validated. RNA variant calling, DNA concordance, ancestry inference, and
clinical interpretation are out of scope.

## Sprocket-native WDL port

A separate WDL 1.2 workflow is available at
[`sprocket/main.wdl`](./sprocket/main.wdl). It invokes the RNA tools directly
inside task containers and does not wrap Nextflow or Cromwell. Its QC-depth
branch gates the full-depth branch, and full-depth Salmon quantification is
scattered per sample.

The port currently covers FASTQ QC/subsampling, fastp, STAR, samtools,
HTSeq, Qualimap, Picard, RSeQC gene-body coverage, a count matrix, MultiQC,
and Salmon. It does not yet reproduce the complete source BaseJumper consensus
QC-reporting process, all source plots/metadata artifacts, or the downstream
analysis dashboard. It must not be described as output-equivalent to the
accepted production run.

Static WDL validation passes with `miniwdl check sprocket/main.wdl`. The
cluster does not have Sprocket or Docker installed, so the image is unbuilt
and the workflow has not been smoke-tested or production-run. See
[`sprocket/README.md`](./sprocket/README.md) for setup prerequisites and the
validation plan.

## WDL/Cromwell orchestration

A retained, legacy WDL entry point is available at [main.wdl](./main.wdl). It
orchestrates the existing BaseJumper Nextflow-backed stages behind audit gates,
scatters full-depth Salmon quantification, and then runs downstream analyses
and result organization. It is not a tool-by-tool reimplementation of
BaseJumper in WDL. Before any WDL execution, create a new isolated results tree;
never point it at this completed `RNA_results` dataset. Preparation and
Cromwell instructions are in [wdl/README.md](./wdl/README.md). No Cromwell
run has been completed end-to-end from a clean checkout. A 47-cell
Cromwell/LSF run is currently in its preflight gate; its exact status and
execution boundary are recorded in the Appendix. Prefer the Sprocket-native
port for new work only after its missing reporting/downstream stages and
container are complete and validated.

### Results organization

`RNA_results/README.md` is the top-level navigation guide. Inputs are under
`RNA_results/inputs/`, the QC-depth and full-depth steps are under
`RNA_results/steps/qc_depth/` and `RNA_results/steps/full_depth/`, canonical
audits are under `RNA_results/audits/`, and logs are under `RNA_results/logs/`.
The user-facing dashboard and analysis package remains
`RNA_results/final_results/`. `organize_rna_results.py --plan` previews the
completed organization and `--apply` performs same-filesystem moves with
compatibility symlinks; it never copies stage trees or removes Nextflow caches.

## Safety and status

**Step-by-step review is now in effect.** Automatic production launch from
[validate_then_run.sh](../RNA_results/validate_then_run.sh) is disabled.
[run_step1_prepare.sh](../RNA_results/run_step1_prepare.sh) submits only input
validation and per-cell merging, with no downstream jobs. Its outputs must
be inspected before alignment is submitted. The early-read inspection in
`RNA_results/step1_input_inspection.tsv` checks 1000 pairs per raw lane input;
it is not a substitute for full gzip/pair validation or full-library FastQC.

Step 1 has now completed for all 16 real RNA cells: 32 merged FASTQs and
16 provenance records verified, no partial files, and exact source/merged
pair-count and compressed-size agreement.
Step 2 compatibility checks passed; the reproducible checker is
[check_step2_references.py](../RNA_results/check_step2_references.py) and its
result is `RNA_results/step2_reference_audit.json`. It checks all GTF intervals,
refFlat exon/transcript membership, countable gene mapping after source-standard
version normalization, housekeeping BED intervals and STAR/FASTA contig lengths.
The preprocessing-only driver **330085142** runs
[run_step3_preprocessing.py](../RNA_results/run_step3_preprocessing.py), reusing
the adapted source process definitions but wiring only raw FastQC, SeqKit
subsampling, fastp and trimmed FastQC. Its products and exact task trace are
under `RNA_results/stage3_qc_preprocessing/`; logs are under `logs_step3/`.
It deliberately does not invoke STAR or submit a dependent downstream driver.
Sampled and trimmed reads are explicitly retained for inspection and the
subsequent alignment stage.

Preprocessing completed and passed the full retained-read audit for all 16
cells: 64 successful tasks, 32 sampled FASTQs, 32 trimmed FASTQs, 16 read-count
files, 16 fastp JSONs and 64 FastQC HTML/ZIP report pairs. Sampling yielded
1,996,604--1,998,424 individual reads per cell (probabilistic target 2M).
Trimming retained 84.6--95.9% of sampled reads; post-trim Q30 ranged
88.7--93.4%. Adapter-content failures disappeared, with one remaining warning.
Duplication, composition and overrepresented-sequence flags remain and must
be interpreted alongside mapping and expression, not silently treated as passes.
Cell 16 has the largest trimming loss, primarily too-short and low-quality reads.
The reproducible audit is
[check_step3_outputs.py](../RNA_results/check_step3_outputs.py), producing
`RNA_results/step3_output_audit.json`.

Alignment-only driver **330096776** runs
[run_step4_alignment.py](../RNA_results/run_step4_alignment.py), reusing source
STAR and samtools processes on the inspected QC-depth trimmed reads. It retains
STAR outputs under `stage4_qc_alignment/star/` and final BAM/BAI under
`stage4_qc_alignment/bam/`. It does not run HTSeq, Qualimap or reports.
Each STAR task requests 8 slots and 64000 MB on one host. Nextflow does not
impose a per-process concurrency cap; LSF schedules tasks according to available
resources. These are QC-depth alignments; the separate full-depth branch is
described in the current progress section below.

Alignment has completed and passed
[check_step4_outputs.py](../RNA_results/check_step4_outputs.py): all 32 tasks,
16 final BAM/BAI pairs, coordinate sorting/sample tags, absence of secondary
records, indexed region retrieval and STAR input/trimmed-count agreement.
Results are in `RNA_results/step4_output_audit.json`. Unique mapping ranges
67.81--83.47%; Cell 13 is lowest and remains in the cohort pending expression
and regional QC, not automatically excluded.

Per-cell counting/QC driver **330097917** runs
[run_step5_counts_qc.py](../RNA_results/run_step5_counts_qc.py) on the audited
QC-depth BAMs: HTSeq, Qualimap, Picard RNA metrics/flagstat and RSeQC gene-body
coverage only. Products are under `stage5_qc_counts/tool_outputs/` and
`stage5_qc_counts/extras/`, with exact tasks in `trace.tsv` and driver logs
under `logs_step5/`. Cohort matrices, PCA, scoring and full-depth analysis
have not been submitted by this staged driver.

Step 5 completed and passed
[check_step5_outputs.py](../RNA_results/check_step5_outputs.py), recorded in
`RNA_results/step5_output_audit.json`: 60664 nonnegative gene counts per cell,
all five HTSeq special categories, Qualimap PDF/text/coverage products,
Picard regional fractions and flagstat, and 100-percentile RSeQC coverage/PDF.
Qualimap was invoked with source PDF output, so an HTML report is not expected.
Detected genes range 4935--8091. Cell 6 has only 46.15% exonic reads in
Qualimap; this is a biological QC flag, not a missing-output failure. No cell
is automatically removed.

Matrix-only driver **330098388** runs
[run_step6_matrix.py](../RNA_results/run_step6_matrix.py): per-cell source
gene/biotype/mitochondrial summaries, their merged tables, and source wide
matrix/housekeeping products/count Parquet. Products are under
`stage6_qc_matrix/summaries/`, `merged/`, and `matrix/`; driver logs under
`logs_step6/`. Sample/gene/count agreement must be audited before PCA,
consensus scoring or cohort reporting is run.

Step 6 completed and passed
[check_step6_outputs.py](../RNA_results/check_step6_outputs.py): 18188 detected
genes across 16 samples, 291008 long-format Parquet rows, with every matrix
and Parquet count matching the original HTSeq count and all assigned counts
preserved. The source matrix retains the union of genes detected in any cell;
the original per-cell HTSeq files retain all 60664 countable genes.
Audit: `RNA_results/step6_output_audit.json`.

Report-only driver **330098435** runs
[run_step7_reports.py](../RNA_results/run_step7_reports.py), using the source
PCA/heatmap, Qualimap summary, consolidated QC/Parquet/consensus and MultiQC
processes on inspected earlier outputs. Products are under
`stage7_qc_reports/pca/`, `qualimap_summary/`, `consensus/`, and `multiqc/`;
logs under `logs_step7/`. Final sample membership and QC labels must be
verified before the QC-depth branch is considered complete.
The first report driver stopped during input discovery because Nextflow
`channel.fromPath` defaults to files, not directories. Qualimap outputs were
present and verified; the staged driver now explicitly selects `type: 'dir'`.
After confirming no active report tasks remained, it was resumed as
**330098437**. The staged report driver accepts `--resume` (default off) only
after inspecting the failed run; it does not rerun input preparation or counting.

The report run finished, but output inspection found a source HTSeq report
parser bug: with `--additional-attr=gene_name`, HTSeq emits three columns,
whereas the source parser treated the second (gene symbol) as a count and
silently reported zeros. The runtime adapter now reads the last column as
an integer and raises on parsing errors. A regression test executes the
adapted Python block on three-column counts and malformed values.
Corrected report-only job **330098618** writes to `stage7_qc_reports_v2/`,
preserving the earlier reports. This stage must still pass the final metrics
audit; upstream HTSeq files and matrices were correct and have not been changed.
The staged driver accepts `--output-name stage7_qc_reports_v2` for this explicit
new report iteration; default remains `stage7_qc_reports`.

**QC-depth staged validation is now complete.** Corrected report job 330098618
finished and passed
[check_step7_outputs.py](../RNA_results/check_step7_outputs.py). All 16 assigned
counts and reported HTSeq categories match original counts; all five source
score components, statuses, and QC Parquet metrics agree. Mandatory plots,
per-cell metric files and MultiQC are present. The results are 14 PASS and
2 Borderline (Cells 6 and 13); no cell is removed. Audit:
`RNA_results/step7_output_audit.json`. The accepted report is
[multiqc_report.html](../RNA_results/stage7_qc_reports_v2/multiqc/multiqc_report.html),
not the superseded report in `stage7_qc_reports/`.
The source exonic scoring denominator includes overlapping-exon counts as
implemented in `parse_qualimap.R`, so scoring proportions differ from the
rounded percentages printed directly by Qualimap; that behavior is preserved.

Full-depth preprocessing-only driver **330098666** runs
[run_step3_preprocessing.py](../RNA_results/run_step3_preprocessing.py)
with `--depth full` (default `qc`). It uses all merged pairs, source fastp
settings and raw/trimmed FastQC, writing separately to
`stage3_full_preprocessing/` with logs in `logs_full_step3/`. The full-depth
audit uses
[check_step3_outputs.py](../RNA_results/check_step3_outputs.py) `--depth full`,
which additionally requires retained preparation pair counts equal the
original merged counts and writes `step3_full_output_audit.json`.
No full-depth alignment, counting or reports are automatically submitted.
Full-depth preprocessing job 330098666 has finished. All 16 fastp JSONs are
present; preliminary retention ranges from 84.6% to 95.9%. Full-file audit job **330107009** completed deep checks for 11 cells before
being stopped at the user's request to shorten audits. A quick all-cohort audit
passed and is recorded separately in `RNA_results/step3_full_output_audit_quick.json`.
It reuses those 11 completed deep checks and checks the first 1,000 pairs in both
prepared and trimmed files for every cell, copy-all compressed sizes, reported
read totals against merge provenance, filtering accounting, all 64 successful
tasks and all FastQC archives. Retention is 84.6-95.9%; post-trim Q30 is 88.7-93.5%.
Full gzip integrity and independently counted totals remain unverified for the
other five cells. Audit logs are under `RNA_results/logs_full_step3_audit/`.
Full-depth alignment driver **330116427** ran
`run_step4_alignment.py --depth full --accept-quick-audit`; its outputs are in
`stage4_full_alignment/`. The final alignment audit passed in quick mode:
all 16 cells passed STAR/BAM, header, sorting, indexing, and input checks, with
three BAMs deeply streamed. The earlier generated configuration used a
four-task cap; this affected only that completed alignment run.

The continuation [run_full_staged_continuation.sh](../RNA_results/run_full_staged_continuation.sh)
ran serially behind that alignment job. Its Step 5 counts/QC audit, Step 6
matrix/count-concordance audit, and Step 7 reports audit all passed for the
full 16-cell cohort. Accepted products are in
`stage5_full_counts/`, `stage6_full_matrix/`, and `stage7_full_reports/`;
audit JSONs are `step5_full_output_audit.json`,
`step6_full_output_audit.json`, and `step7_full_output_audit.json`.
The job logs and `full_staged_status.txt` preserve the execution trace.

The annotated standalone Obsidian canvas is
[basejumper_scRNA_tools_annotated.canvas](./docs/basejumper_scRNA_tools_annotated.canvas).
It covers tools, references, scientific settings, QC/full branches, output
paths, audit gates and limitations. The Mac reference canvas and destination
were inaccessible from this Linux session; transfer the file to your Mac's
`RNA_Workflow` folder. Canvas status is a snapshot, not a live dashboard.
The successful staged QC and full-depth runs do not by themselves certify the
single-launch workflow's revised BAM-map publication; the synthetic integration
failure described below remains a separate, outstanding check.

This workflow is separate from DNA. It uses the new **resolveome-rna** conda
environment; it does not modify `resolveome-align` or DNA results.
**Real RNA QC-depth and full-depth staged processing are complete for all 16
cells.** Every staged output audit passed. The alignment audit remains explicitly
quick mode as described above; this is not a claim that every full-depth BAM
received a complete end-to-end read scan.

Historical validation: the complete human reference build finished
successfully, and a native Nextflow/LSF probe verified the RNA environment,
one-host placement and a 16-GB total allocation for a 2-slot task. The first
end-to-end synthetic run reached source plotting and exposed a missing `gt` R
package; it was added to the environment and preflight. Both synthetic
branches completed, but final BAM publication audit failed. LSF job
**330078337** exited with code 1 after running
[validate_then_run.sh](../RNA_results/validate_then_run.sh): complete preflight,
both synthetic branches and the mandatory output audit.
Logs are under `RNA_results/logs_validation/`. The BAM publication adapter
has been changed to explicitly publish BAM/BAI fields from the source map;
that correction has not been certified by a rerun of the integrated synthetic
launcher. The successful staged production runs do not bypass that separate
publication check, and the integrated launcher should not be assumed
certified from the staged results.

The completed cohort was prepared from paired raw RNA FASTQs using a
cohort-specific LSF merge and audit outside this source package. Input
inventories record physical run/lane identities; original FASTQs remain
unchanged. Each RNA sample was retained independently of any DNA QC labels.
The selected analysis included both the approximately two-million-individual-
read QC branch and a separate full-depth expression branch, with independent
outputs and passing stage audits. Run-specific input paths and submission
scripts are intentionally not distributed; use the portable sample-sheet and
WDL preparation instructions for a new cohort.

# Specification

## Dependencies

The audited production workflow uses the dedicated `resolveome-rna` Conda
environment and native Linux executables, with Nextflow and LSF available on
the cluster. The separate Salmon environment pins Salmon 1.6.0 to match the
existing shared index. Pinned production tools include STAR 2.7.6a, fastp
0.20.1, SeqKit 2.13, samtools 1.21, HTSeq 2.0.5, Qualimap 2.2.2d, FastQC
0.12.1, Picard 3.2.0, RSeQC 5.0.4, MultiQC 1.33, Python 3.11, and R 4.3.
R plotting dependencies and their versions are listed in the environment
definition and runtime manifests.

The Cromwell route additionally requires Java, Cromwell 82 or a compatible
version, Womtool for validation, and an LSF backend configured for the
cluster. Its custom runtime attributes are not portable to an arbitrary WDL
engine. The Sprocket route requires a working Sprocket installation and
container runtime; those are not available in the documented cluster
environment, and its task images have not been built.

The portable package definitions are
[setup/environment.yml](./setup/environment.yml),
[setup/analysis_salmon_environment.yml](./setup/analysis_salmon_environment.yml),
and [setup/setup.sh](./setup/setup.sh). They describe software setup, but do
not bundle the upstream BaseJumper source or local production data.

## Reference files

The reference assets are a coordinated GRCh38/Ensembl release 104 set. The
portable installer places generated files under
`wdl_runtime/references/Ensembl104/` by default:

| Resource | Path or role |
|---|---|
| Primary-assembly FASTA | `GRCh38.primary.fa` |
| Ensembl 104 GTF | `Homo_sapiens.GRCh38.104.chr.gtf` |
| Transcript-to-gene mapping | `Ensembl104.tx2gene.tsv` |
| Transcript BED12 | `Ensembl104.transcripts.bed` |
| STAR genome index | `STAR/` |
| Salmon transcript index | `salmon_index/` |
| Picard refFlat | `Ensembl104.refFlat.txt` |

The production installation uses shared project GRCh38 references under
`pipeline_resources/genomes/Homo_sapiens/NCBI/GRCh38/Sequence/`, with a
private complete STAR index under `workflow_RNA/references/STAR_Ensembl104`.
The STAR index was built with STAR 2.7.6a and `sjdbOverhang=100`. The matching
Ensembl 104 GTF, tx2gene mapping, housekeeping BED, and refFlat are audited
together. The source runtime normalizes versioned transcript identifiers
where required; shared annotation files are not rewritten.

The reference bundle passed the documented compatibility audit (194 FASTA/
STAR contigs, 47 GTF contigs, 60,664 genes, 1,498,409 exons, 236,920 refFlat
transcripts, 258,710 tx2gene rows, and 3,802 housekeeping BED records).
Reference files and indices are not copied into the public source repository.
A different assembly or annotation release requires a matching rebuilt
bundle and a new compatibility audit.

## Workflow (RNA production pipeline)

The production route is a staged, audit-gated RNA workflow, not a single
unconditional launcher:

1. **Prepare inputs.** Validate the sample manifest and paired FASTQs;
   merge lane/flowcell inputs without modifying originals and retain merge
   provenance. Stage immutable input links and record sample membership.
2. **Preflight references and runtime.** Verify tool versions, reference
   compatibility, Nextflow/LSF availability, and configured paths before
   analysis submission.
3. **QC-depth branch.** For each cell, sample approximately two million
   individual reads, run raw and trimmed FastQC around fastp trimming, align
   with STAR, create and audit BAM/BAI, count with HTSeq, collect Qualimap/
   Picard/RSeQC metrics, and produce the source matrix, QC scores, figures,
   and MultiQC report.
4. **Full-depth branch.** Independently use every merged read pair, then
   repeat preprocessing, alignment, counting, matrix generation, and
   reporting under separate output paths. Each stage is submitted only after
   its predecessor's required audit passes.
5. **Transcript abundance and cohort reporting.** After full-depth
   alignment/count audits pass, scatter Salmon quantification per sample,
   validate each `quant.sf`, and generate the descriptive downstream
   expression package.
6. **Organize results.** Build a navigable final-results/dashboard view
   without deleting Nextflow work directories, audit records, or source
   outputs.

Sample-sheet and JSON generation, Sprocket examples, and Cromwell commands
are documented in [README.md](./README.md) and [wdl/README.md](./wdl/README.md).
Production execution must use a new, isolated results root and the documented
stage-by-stage review process.

## Tool-by-tool commands, in execution order

The ordered commands and representative invocation settings are detailed in
[Ordered tool commands](#ordered-tool-commands). The commands are templates;
the exact commands, task traces, and runtime configuration for a run are
retained with its execution outputs. The sequence is input validation and
sampling/copy-all, FastQC/fastp, STAR, samtools, HTSeq, Qualimap and
supplemental metrics, source expression matrix/plots, and MultiQC/final
audit. Salmon and optional downstream analyses are described separately
under [Analyses after this workflow](#analyses-after-this-workflow).

## Workflow (general-purpose WDL pipelines)

Two WDL execution models are documented:

- **Cromwell/LSF:** `main.wdl` imports task modules under `wdl/tasks/` and
  invokes the adapted BaseJumper Nextflow stages through audit-gated tasks.
  Reference paths are passed as strings so Cromwell does not localize large
  shared references without their adjacent indices. Salmon is scattered per
  sample; downstream cohort analysis and result organization follow its
  validation. The workflow is cluster-specific and is not certified on a
  generic WDL engine.
- **Sprocket-native:** `sprocket/main.wdl` invokes tools directly inside task
  containers. Its current task set covers QC/subsampling, fastp, STAR,
  samtools, HTSeq, Qualimap, Picard, RSeQC, matrices, MultiQC, and Salmon,
  but omits complete source consensus reporting and the downstream analysis
  dashboard. Static validation is not evidence of runtime or output parity.

Neither WDL definition should be described as a complete, independently
validated replacement for the audited production workflow. See the
individual [WDL guide](./wdl/README.md) and
[Sprocket guide](./sprocket/README.md) for JSON schemas, launch requirements,
and validation limitations.

## Output organization

Run-specific output trees are separate from source code and preserve raw
inputs by reference. The historical 16-cell cohort uses `RNA_results/`;
the distinct combined cohort is staged under `RNA_results_Expended_full47/`.
The portable WDL/Sprocket package does not include those run data.

The audited result layout keeps input manifests and merged FASTQs, QC-depth
and full-depth stage outputs, canonical audit JSON files, task/driver logs,
Nextflow work caches, Salmon results, and the final analysis/dashboard
separate. The `RNA_results/README.md` navigation guide and
`organize_rna_results.py` document the exact local layout and non-destructive
organization behavior. The top-level
[Outputs and validation](#outputs-and-validation) section below records
required products and resume safeguards.

## Execution and failure handling

Preflight and planning are separate from production submission. The
step-by-step launch procedure disables automatic production launch; the
preparation step does not submit downstream analysis. Before advancing a
stage, inspect its exit status, required files, sample membership, and
audit report. A successful process exit alone is not sufficient evidence of
a valid stage.

Failures in source report parsing or required cohort products must be
surfaced, not converted into zero counts, missing-success defaults, or
partial success. Each retry uses a distinct output iteration or is explicitly
resumed only after checking that no related driver/task is still active.
Preserve run manifests, task logs, audit reports, inputs, and Nextflow
caches; do not delete caches or reuse a completed result directory to force
a rerun. Any change to input membership, references, code, or scientific
settings requires a fresh isolated results root.

## Detailed implementation sections

## Steps and source parity

| Stage | Implementation |
|---|---|
| 0. Raw preparation | `python_lsf_pipeline/`: optional FASTQ discovery, merging, and integrity audit utilities; input data are not modified |
| 1. Samples | Validate the four-column sample sheet, then stage an isolated run manifest with immutable FASTQ symlinks |
| 2. References | `setup/`: public Ensembl 104 FASTA/GTF/transcripts, derived annotation helpers, STAR and Salmon indexes |
| 3. QC-depth branch | Approximately 2 million individual reads per sample, followed by trimming, alignment, counting, QC, matrices and reports |
| 4. Full-depth branch | Same principal methods without read subsampling; each stage depends on the previous passing audit |
| 5. Transcript abundance | One full-depth Salmon task per sample, then cohort-wide validation |
| 6. Final analyses | Cohort summaries, figures, dashboard, and organized results |

The WDL implementation preserves BaseJumper scientific process definitions
through Nextflow adapters; it is not a pure-WDL reimplementation or an
untouched stock BaseJumper run. The upstream source checkout is not included
because redistribution permission has not been established. Obtain it
separately from an authorized source. The completed local runs and their audit
summaries below are historical evidence, not a guarantee that a fresh public
checkout is production-certified or bitwise identical.

### Preserved source settings

- Paired-read subsampling proportion = target **individual reads** / (R1+R2).
  Thus 2 million reads means approximately 1 million pairs, not 2 million pairs.
  SeqKit uses seed 12345 for both mates; sampling is probabilistic, not an exact
  fixed-count draw.
- fastp 0.20.1: `--detect_adapter_for_pe`; automatic poly-G detection is retained.
  Optional `--trim-poly-g` forces poly-G trimming. No RNA duplicate removal.
- STAR 2.7.6a: two-pass Basic, sorted BAM, source unmapped/primary/chimeric flags.
- samtools `-F 256`: excludes secondary alignments; **does not exclude**
  supplementary alignments. SM/PL read-group tags are set as in the source.
- HTSeq union/exon/gene_id counting, unstranded by default, secondary ignored.
- Source Ensembl104 biotypes and `MT-` gene-symbol mitochondrial definition.
- Source normalization, variable-gene plots, housekeeping genes, consensus
  scores, per-cell metrics, expression/QC Parquet schemas and MultiQC layouts.

### Explicit adaptations

1. Native conda tools on Linux x86/LSF replace source cloud/custom Docker images.
   Source R scripts and assets are copied into the runtime, not rewritten.
2. Gzip CRC, four-line record structure, mate names and counts are checked;
   mismatches/corruption are fatal. Input files are never repaired in place.
3. Source subsampled `FINAL_READS` counted only R1, whereas `TOTAL_READS` counted
   both mates. The runtime reports **R1+R2 consistently** and revalidates sampled pairs.
4. Failures in source matrix/PCA/report tasks terminate the run instead of being
   silently ignored. Source scripts may still have internal zero-input behaviors;
   the final audit verifies every expected cell and mandatory product.
5. STAR BAM sorting gets an explicit maximum 8 GB (or one-quarter of task RAM,
   if smaller). Qualimap Java heap leaves 8 GB for native/non-heap overhead.
6. Raw/trimmed FastQC, Picard RNA metrics and RSeQC are additions. The upstream
   README mentions gene-body coverage, but its current main workflow does not
   wire that process; this adapter does.
7. HTSeq **2.0.5** replaces source container 0.13.5 to work with Python 3.11;
   its three-column output is parsed using the final numeric count column.
   R/package versions also differ from the source custom images. Numerical or
   pixel-perfect identity is not claimed.

## Environment

Definition: [setup/environment.yml](./setup/environment.yml) and
[setup/analysis_salmon_environment.yml](./setup/analysis_salmon_environment.yml).
Installer: [setup/setup.sh](./setup/setup.sh). The default runtime location is
`wdl_runtime/`; set `RNA_RUNTIME_DIR` or pass `--runtime-dir` to relocate it.

Pinned tools include STAR 2.7.6a, fastp 0.20.1, SeqKit 2.13, samtools 1.21,
HTSeq 2.0.5, Qualimap 2.2.2d, FastQC 0.12.1, Picard 3.2.0, RSeQC 5.0.4,
MultiQC 1.33, Python 3.11, R 4.3, ggplot2 3.5 and ggtree 3.10.
The source plotting dependencies AMOR and ohchibi are installed from pinned
archives; harrietr is pinned from CRAN. Runtime environment manifests and
package archives stay under `wdl_runtime/`, outside the source package.
Cromwell, Java, and an LSF installation are separate requirements. No
cluster-specific conda prefix or project reference directory is required.

## References

The portable setup downloads Ensembl release 104 under
`wdl_runtime/references/Ensembl104/` (or the configured runtime directory).

| Resource | Generated path |
|---|---|
| Primary-assembly FASTA | `GRCh38.primary.fa` |
| Ensembl release-104 GTF | `Homo_sapiens.GRCh38.104.chr.gtf` |
| Five-column tx2gene | `Ensembl104.tx2gene.tsv` |
| Transcript BED12 | `Ensembl104.transcripts.bed` |
| STAR index | `STAR/` |
| Salmon transcript index | `salmon_index/` |
| Picard refFlat | `Ensembl104.refFlat.txt` |

Annotation helpers are derived from the downloaded GTF, preventing mismatched
assemblies. Downloads are validated against Ensembl `CHECKSUMS`. STAR index
construction requires at least 48 GB RAM (64 GB recommended); plan for about
100 GB free disk for the full software/reference build.

```bash
# Install software and public references/indexes.
bash setup/setup.sh --threads 8 --memory-mb 64000
```

See [wdl/README.md](./wdl/README.md) for sample-sheet validation, isolated run
preparation, Cromwell inputs JSON, and launch commands. A fresh results
directory is mandatory; completed results are never a default target.

## Portable setup and WDL launch

The default install builds tools and public references under `wdl_runtime/`.
For software only, use `bash setup/setup.sh --skip-references`; otherwise
allow at least 48 GB RAM and approximately 100 GB free disk for references
and indexes.

Prepare a sample sheet with exactly `sample_id`, `r1`, `r2`, and `group`
columns and at least five samples for cohort plotting. Set up a **new** results
directory, then follow the complete, current Cromwell/LSF procedure in
[wdl/README.md](./wdl/README.md). The workflow accepts queue and resource
overrides including `nextflow_threads` and `nextflow_memory_mb`; default
per-cell tasks request 8 CPUs and 64 GB. These are cluster resource requests,
not total workflow limits.

The Cromwell configuration uses custom `memory_mb`, `lsf_queue`, and
`lsf_walltime` runtime keys. The WDL is not certified on a generic WDL engine,
and a clean-checkout LSF run has not been completed. `miniwdl check` reports
these Cromwell-specific keys as unknown while still parsing the WDL.

The alternative in `sprocket/` is still a prototype. It does not yet reproduce
the complete BaseJumper consensus QC report and downstream dashboard, and its
container has not been built or executed.

For source tests, use the RNA environment and run
`python -m unittest discover -s tests -q`; the opt-in
`tests/test_real_rna_smoke.py` uses synthetic reads and is separate from a
production cohort run.

## Ordered tool commands

Templates below show default settings; actual source commands and output paths
are retained in each branch's Nextflow work-directory `.command.sh` files.
Let `SID`, `R1`, `R2`, `BAM`, `GTF`, `INDEX`, `TX2GENE`, `BED`, `REFFLAT`
represent the current sample and configured resource paths.

### 1. Validate and subsample/copy

```bash
python python_lsf_pipeline/step1_prepare_samples.py --validate-pair "$R1" "$R2"
seqkit sample -p "$PROPORTION" -s 12345 -j 2 -o "${SID}_subsampled_R1.fastq.gz" "$R1"
seqkit sample -p "$PROPORTION" -s 12345 -j 2 -o "${SID}_subsampled_R2.fastq.gz" "$R2"
```

QC `PROPORTION=2000000/(2*validated_pair_count)`. If total is below the target,
or in full-depth mode, copy both mates instead of sampling. The latter uses a
target of `9223372036854775807` to select the source's no-sampling branch.
Raw FastQC runs on the original input; raw and trimmed files are never confused.

### 2. fastp and FastQC

```bash
fastqc -t 2 -o . "$R1" "$R2"
fastp --thread 2 --in1 "$SAMPLED_R1" --in2 "$SAMPLED_R2" \
  --out1 "${SID}_R1_trim.fastq.gz" --out2 "${SID}_R2_trim.fastq.gz" \
  --json "${SID}_fastp.json" --html "${SID}_fastp.html" --detect_adapter_for_pe
fastqc -t 2 -o . "${SID}_R1_trim.fastq.gz" "${SID}_R2_trim.fastq.gz"
```

Optional forced poly-G appends `--trim_poly_g`. Version defaults: Q15 qualified
base, 40% unqualified limit, at most 5 Ns, minimum read length 15, minimum poly-G
length 10. No UMI extraction, 3'-tag counting or duplicate removal is added.

### 3. STAR

```bash
STAR --genomeDir "$INDEX" --runThreadN 8 --limitBAMsortRAM 8000000000 \
  --outSAMunmapped Within KeepPairs --twopassMode Basic --outReadsUnmapped None \
  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "${SID}_" \
  --readFilesCommand zcat --readFilesIn "$TRIM_R1" "$TRIM_R2" \
  --outSAMstrandField intronMotif --outSAMattributes All \
  --outFilterMultimapNmax 10 --outSAMprimaryFlag OneBestScore \
  --chimSegmentMin 12 --chimJunctionOverhangMin 8 --chimOutJunctionFormat 1 \
  --outSAMattrRGline "ID:$SID"
```

### 4. samtools primary filter, tag update, sort/index

```bash
samtools view -F 256 -h "$STAR_BAM" \
  | sed -e 's/\tSM:[^\t]*//' -e 's/\tPL:[^\t]*//' \
  | sed "s/^\\(@RG\\tID:[^\\t]*\\)/\\1\\tSM:${SID}\\tPL:Illumina/" \
  | samtools view -h -b - | samtools sort -O BAM -o "${SID}.bam"
samtools index "${SID}.bam"
samtools quickcheck -v "${SID}.bam"
samtools flagstat "${SID}.bam" > "${SID}.flagstat.txt"
```

### 5. HTSeq counts and source annotation

```bash
htseq-count -f bam -r pos -s no -t exon -i gene_id \
  --additional-attr=gene_name -m union --secondary-alignments=ignore \
  "$BAM" "$GTF" > "${SID}.htseq_counts.tsv"
Rscript "$RUNTIME/scripts/htseq_to_gene_df.R" "${SID}.htseq_counts.tsv" "$TX2GENE"
```

Per-cell gene, biotype and mitochondrial tables are merged into the source
`df_gene_counts_starhtseq.tsv`, `df_gene_types_detected_summary_starhtseq.tsv`,
`df_mt_gene_counts_starhtseq.tsv`.

### 6. Qualimap and supplemental RNA metrics

```bash
qualimap --java-mem-size=54G rnaseq -gtf "$GTF" -bam "$BAM" --sorted \
  -outfile "${SID}.report.pdf" -oc "${SID}.oc" -outdir "qualimap_outdir_${SID}"
picard -Xmx8g CollectRnaSeqMetrics I="$BAM" O="${SID}.rna_metrics" \
  REF_FLAT="$REFFLAT" STRAND_SPECIFICITY=NONE
geneBody_coverage.py -r "$BED" -i "$BAM" -o gene_body
```

Qualimap heap is `task.memory.toGiga()-8` (rounded Nextflow GiB, minimum 1).
The 54G template reflects a 64000 MB request in binary GiB accounting.
Picard lacks a matching rRNA interval list here: its ribosomal metric is not
claimed available. Source HTSeq rRNA/Mt_rRNA biotype fractions are retained.
RSeQC uses the provided source housekeeping BED12, not a different whole-transcript
annotation or an invented 5'/3' ratio.

### 7. Source expression matrix and plots

```bash
Rscript "$RUNTIME/scripts/htseq_summarydf_to_matrix.R" GRCh38
Rscript "$RUNTIME/scripts/pca_heatmap_qc_salmon_htseq_nometadata.R"
Rscript "$RUNTIME/scripts/parse_qualimap.R"
Rscript "$RUNTIME/scripts/rna_qc_plot.R" \
  --matrix_file "$MATRIX" --metrics_file combined_selected_metrics.csv \
  --metadata_file metadata_enriched.csv
```

R consensus strict comparisons: mappable >0.70, exonic >0.70,
intergenic <0.10, mitochondrial <0.10, detected protein-coding genes >500.
Scores 4/5 PASS; 3 Borderline; 0–2 FAIL. Full-depth cells are not filtered out
based on these QC labels.

### 8. Parquet, MultiQC and final audit

Source embedded Python creates count and per-biosample QC Parquet files.
Source MultiQC uses the source-generated config/logo:

```bash
multiqc . -n multiqc_report.html -c multiqc_config_runtime.yaml --force
multiqc "$BRANCH/extras" "$BRANCH/tool_outputs" --force \
  --outdir "$BRANCH/reports" --filename multiqc_supplemental.html
python extra_methods/step4_rna_dashboard.py --outdir "$OUTDIR" --cells rna_cells.tsv
```

## Outputs and validation

Default new root: `P/RNA_results`.

Real-data preparation adds `raw_inputs.json`, `rna_cells.tsv`, `fastqs/`
(merged mates and per-cell merge provenance), and `logs_prepare/`.
The standalone raw preparer accepts `--raw-dir DIR --outdir DIR` to plan
(no FASTQ merge), or `--manifest FILE --merge-index N` to merge the
1-based Nth planned cell. These flags have no default paths and must be supplied.

- `qc_2M/`: source-comparison QC, BAMs, raw/trimmed FastQC, RNA metrics,
  gene-body coverage, expression and QC Parquet, source figures/MultiQC.
- `full_depth/`: corresponding full-depth counts and outputs.
- `workflow_dashboard.html`: links both depths without overwriting source reports.
- `output_audit.json`: expected sample membership and verified outputs.
- `runtime/basej-rnaqc/`: runtime source and provenance manifest.
- `scripts/rna_native.config`: actual per-process resource/environment settings.
- `rna_input.csv`: exact BaseJumper input snapshot.
- `tool_preflight.json`, `workflow_state.json`, `workflow.lock`.
- `logs_driver/`: branch logs, trace, execution report and timeline.
- `work_qc_2M/`, `work_full_depth/`: Nextflow cache and exact task commands.

Never delete work caches to attempt a resume. A changed sample manifest,
scientific parameters, references or code requires a new output root.
An interrupted stage refuses automatic resubmission unless `--resume` is explicitly
provided after inspecting lingering drivers/tasks.

Regression tests:

```bash
conda run -n resolveome-rna python -m unittest \
  test_rna_workflow test_rna_downstream test_quantify_transcripts
```

Opt-in real-tool synthetic integration:

```bash
conda run -n resolveome-rna python test_real_rna_smoke.py \
  --directory /path/to/new_synthetic_fixture --mode both
```

## Analyses after this workflow

The separate full-depth cohort analysis uses the already audited count matrix,
full-depth STAR junction tables and full-depth trimmed FASTQs. It does not
rerun the production workflow or alter its stage outputs. The cohort contains
16 biological replicates of one SY5Y condition, so these results are
descriptive; no condition-based differential test is scientifically defined.

### RNA-only research extensions

These are optional, separate analyses of RNA outputs. They must not change the
BaseJumper-derived counts, BAMs, or QC decisions. The standard analysis remains
usable without them. No DNA data, matched-DNA identity checks, or variant
calling are part of this RNA-only branch.

#### Expanded library QC

The current workflow already reports FastQC/fastp quality, STAR mapping,
Qualimap feature-region composition, Picard RNA metrics, RSeQC gene-body
coverage, mitochondrial counts/fraction, and a fastp duplication estimate.
Keep the no-deduplication policy: these libraries have no UMIs, so identical
read sequences are not sufficient evidence of PCR duplicates.

The optional `extra_methods/run_rna_research_addons.py` runner implements:

- per-sample median RSeQC TIN with cohort-relative robust-z outlier flags
  (borderline at z ≤−2, fail at z ≤−3, minimum cohort size 8);
- an RSeQC positional duplication estimate and figure, explicitly labeled as
  a mapping-position estimate rather than UMI molecule counting;
- separate BaseJumper and enhanced-QC status/reason columns, with no automatic
  sample exclusion.

The BaseJumper thresholds remain the authoritative first label; enhanced
status adds TIN outlier review. These thresholds are research flags, not
clinical cutoffs. Ribosomal fraction, insert-size/read-length summaries, and
documented depth-saturation curves remain future improvements. RNA-to-DNA
genotype concordance/sample identity is deferred to a later matched DNA/RNA
add-on. The optional module is implemented but has not been run or validated on
a cohort.

#### Exon-level expression and splicing

The optional runner adds DEXSeq-style flattened exon-bin counts and SUPPA2
event PSI estimates from Salmon TPM, with heatmaps and per-sample tables. It:

- derives exon-bin counts from audited BAMs with the Ensembl 104 annotation;
- creates event definitions and descriptive event PSI estimates, including
  skipped-exon events where SUPPA2 can define them;
- does not currently preserve event-level raw inclusion/exclusion counts or
  uncertainty estimates alongside PSI, so PSI is exploratory and must not be
  interpreted as a differential-splicing result;
- runs no differential exon usage/splicing test, including for the one-
  condition SY5Y cohort.

The add-on environment pins DEXSeq 1.46.0, SUPPA2 2.3, and the related tools.
The runner is not yet a Sprocket/WDL scatter stage and has not been run or
validated with synthetic or real data. For the existing one-condition SY5Y
cohort, report descriptive exon and junction evidence only; do not manufacture
groups or report a differential-splicing result.

#### RNA fusion discovery

The existing STAR junction tables are not fusion calls, and the current
alignment configuration does not provide the complete fusion-caller evidence
required by Arriba. Add a distinct, full-depth fusion branch that:

- realigns the already trimmed full-depth FASTQs to the same GRCh38/Ensembl
  104 reference with STAR chimeric alignments enabled in Arriba-compatible
  mode, without replacing or modifying the standard expression BAMs;
- calls candidate fusions with Arriba and retains passed and discarded
  candidates, supporting-read evidence, STAR parameters, references, and logs;
- summarizes recurrent/high-confidence candidates for research review.

The optional runner pins Arriba 2.5.1 and uses the existing GRCh38/Ensembl 104
inputs plus the caller-supplied blacklist, known-fusion, and protein-domain
resources. The dashboard summarizes breakpoints and Arriba confidence/filter
labels alongside split/discordant support; no unsupported local read-count
threshold is asserted. The analysis uses an isolated fusion-specific STAR
alignment and retains discarded calls. The runner has not been validated with
synthetic fusion-positive/negative controls and is not yet a Sprocket/WDL
scatter stage. Fusion predictions are research candidates, not clinically
reportable findings.

#### Explicitly deferred

- RNA-to-DNA sample identity and genotype concordance are deferred until the
  matched DNA/RNA sample map and both validated workflows are available.
- RNA variant calling and DNA-assisted variant interpretation are deferred
  until matched DNA/RNA samples are ready. The RNA-only workflow will not
  report SNV/indel calls.
- Clinical assay validation, interpretation, and clinical reporting are out
  of scope pending discussion with the clinical group. No clinical workflow
  or diagnostic claim is being added here.

#### Selected method references

- Luecken and Theis, [Current best practices in single-cell RNA-seq
  analysis](https://doi.org/10.15252/msb.20188746), for QC and analysis
  principles.
- Piskol et al., [Reliable identification of genomic variants from RNA-seq
  data](https://doi.org/10.1016/j.ajhg.2013.08.008), for RNA-specific variant
  interpretation caveats; variant analysis remains deferred here.
- Uhrig et al., [Accurate and efficient detection of gene fusions from RNA
  sequencing data](https://doi.org/10.1101/gr.257246.119), for Arriba and
  STAR-based fusion evidence.
- Anders and Reyes, [DEXSeq: Inference of differential exon usage in RNA-Seq](https://doi.org/10.18129/B9.bioc.DEXSeq),
  for exon-bin analysis.
- Shen et al., [rMATS: robust and flexible detection of differential
  alternative splicing](https://doi.org/10.1073/pnas.1419161111), for
  replicate-aware event-level testing.
- Huang and Sanguinetti, [BRIE2: computational identification of splicing
  phenotypes from single-cell transcriptomic experiments](https://doi.org/10.1186/s13059-021-02461-5),
  as a method to consider only when its single-cell assay assumptions fit.

### Matched-version Salmon transcript abundance

The shared Salmon index records version 1.6.0. A separate environment keeps
that exact quantifier version isolated from the ResolveOME production
environment. From `workflow_RNA/`:

```bash
conda run -n resolveome-rna python extra_methods/quantify_transcripts.py --plan
conda run -n resolveome-rna python extra_methods/quantify_transcripts.py --submit
# After all 16 LSF jobs finish:
conda run -n resolveome-rna python extra_methods/quantify_transcripts.py --validate
conda run -n resolveome-rna python extra_methods/run_rna_downstream.py
```

The submitter gates on the full-depth alignment audit, validates the 16
full-depth trimmed pairs and shared index version, refuses existing outputs,
and persists each LSF submission before proceeding. Validation checks every
`quant.sf` schema, estimates and TPM total. The completed execution used
`run_rna_downstream_after_salmon.sh`, an LSF job with `ended()` dependencies on
all 16 Salmon jobs; it validated Salmon outputs and built the analysis package
automatically. Driver and Salmon logs are retained under
`RNA_results/logs/downstream/`; the former package log paths remain symlinks.
The analysis status is PASS. Do not resubmit the completed jobs.

### Cohort characterization and organized results

After the Salmon output audit passes:

```bash
conda run -n resolveome-rna python extra_methods/run_rna_downstream.py
```

This creates `RNA_results/final_results/` with:

- descriptive median-ratio normalized expression, PCA, correlations and
  abundance/variability summaries;
- exploratory Reactome over-representation among variable genes and
  within-cohort gene-wise expression scores;
- Salmon transcript TPM/NumReads and gene-aggregated TPM;
- STAR unique splice-junction support, without PSI or differential-splicing
  claims;
- per-replicate expression outliers, explicit QC/provenance and a dashboard.

See `RNA_results/final_results/README.md` for the folder map and limitations.
The dashboard includes the requested method/design credit and an explicit
non-endorsement statement for the St. Jude logo. The BAM folder is a symlink
to the existing audited alignments, not a duplicate copy.

The dashboard also includes BaseJumper's original QC figures and all 68
per-sample fields from `rnaqc_all_metrics.tsv`, including mitochondrial read
counts and percentages. MT percentage is calculated from HTSeq counts assigned
to `MT-` genes divided by all assigned gene counts. The fastp duplication
estimate is reported but not removed. To refresh only the dashboard from
completed analysis outputs without rerunning the analyses:

```bash
conda run -n resolveome-rna python extra_methods/run_rna_downstream.py --refresh-dashboard
```

### Conditional analyses

Differential expression, differential transcript usage, differential
expression, and contrast-specific pathway enrichment require biological
group labels and a valid design. They must not be inferred by splitting these
single-condition biological replicates. If a future experiment adds groups,
use a pre-specified design with replicate and batch metadata and a
count-modeling method appropriate to that design. The delivered pathway
scores are not GSVA/ssGSEA, and the per-replicate enriched-gene table is not
a validated cell-line marker analysis.

## Reference compatibility and limitations

- This is full-length paired-end RNA, not a droplet/UMI pipeline.
- RNA-only research extensions have an optional stand-alone runner but are not
  yet scatter-based Sprocket/WDL stages and have not been run or validated.
- Source plotting assumes a multi-cell, sufficiently expressed cohort.
  Full launcher requires at least 5 samples; low expression/insufficient variable
  genes can still fail source R plots and will be surfaced rather than hidden.
- No inferred strand setting: default unstranded matches BaseJumper.
- No invented RNA somatic/lineage or differential-expression conclusions.
  Counts/PCA/heatmaps do not establish ancestry or treatment effects.
- No RNA variant calls, DNA concordance, or clinical interpretation are
  included.
- Shared reference files are not duplicated, and public/source-code version
  equivalence is documented rather than asserted.
- Synthetic readiness is not a claim of completed production RNA results.

# Appendix

## Verification and current implementation status

**Completed 16-cell cohort:** the QC-depth and full-depth staged workflows
completed with passing per-stage audits. The full-depth expression matrix
contained 40,019 genes across 16 samples. Source-derived QC labels were 14
PASS and 2 Borderline; all 16 cells were retained. The full-depth alignment
audit used quick checks for all samples and deeper streaming checks for three,
so it does not claim a complete independent read scan of every BAM. Salmon
quantification and the descriptive downstream analysis package passed their
recorded audits. These results are historical and must not be resubmitted.

**Combined 47-cell Cromwell/LSF workflow:** a separate run tree at
`RNA_results_Expended_full47/` was prepared from 31 expanded-cohort cells and
the 16 original RNA cells. The 94 FASTQ paths are symlinks to the preserved
inputs. The combined sample manifest, generated Cromwell JSON, tool
preflight, and reference audit are retained there. Womtool 82 validates the
WDL and input JSON; the local tool preflight and reference audit passed.

Cromwell workflow `f69178cf-652b-4856-8d8a-8fd3b1cceaeb` passed preflight
(LSF job `334440049`) and the 47-cell QC preprocessing audit (job
`334440397`). QC alignment job `334443961` then exited with code 1 before
submitting STAR tasks: the stage driver used portable default references
instead of the configured shared-reference paths, and could not find
`wdl_runtime/references/Ensembl104/GRCh38.primary.fa`. The reference handoff
is corrected in the source, but not yet validated in a new workflow attempt.
No other QC or full-depth analysis stages passed in this attempt. Its logs,
execution directory, and audited preprocessing products are preserved. Do
not reuse this results tree: preflight intentionally refuses existing stage
outputs. Consult the run manifest and current Cromwell/LSF state before
deciding how to continue.

**Sprocket port:** `miniwdl check sprocket/main.wdl` passes static validation,
but Sprocket and Docker are unavailable in the documented cluster setup.
Images are unbuilt, and no runtime smoke test or production run has been
performed. The Sprocket workflow is not output-equivalent to the production
workflow.

A workflow is considered complete only when the expected stage manifests,
job outcomes, per-sample membership, required output checks, and downstream
reports have been audited. Static WDL validation or a successful preflight
does not establish successful cohort processing.

## Known issues and fixes baked into this workflow

The following issues were discovered during staged RNA workflow development.
The corresponding fixes are present in the documented runtime/adapters;
preserve the distinction between completed stage validation and unresolved
integration checks.

- **Three-column HTSeq report parsing.** With `--additional-attr=gene_name`,
  HTSeq writes gene ID, symbol, and count. The source report parser treated
  the symbol column as the count and could silently report zeros. The runtime
  adapter now parses the final numeric column and raises on malformed or
  negative counts; a regression test covers valid and invalid rows.
- **Cromwell reference localization.** Passing shared FASTA assets as WDL
  `File` inputs localized the FASTA without its adjacent index. The WDL now
  passes shared reference locations as strings and the configured reference
  paths are checked by preflight.
- **Cromwell JSON-file parsing.** `write_json()` returns a file path, not
  inline JSON. The adapter now reads the generated file before extracting
  sample identifiers.
- **Preflight shell variables and configured reference paths.** Early
  attempts exposed unset `basej_source_dir`/`runtime_dir` values and a helper
  that assumed setup-local references. The task initializes its paths,
  activates the RNA environment before checks, and the reference audit reads
  the actual preflight report.
- **Groovy integer division in task resource calculations.** A `/` expression
  produced a `BigDecimal` where an integer was required. The task uses
  integer division for the corresponding resource calculation.
- **Nextflow discovery of directory outputs.** `channel.fromPath` defaults
  to files; the staged report driver now explicitly selects directories for
  Qualimap results. A corrected report iteration was written to a separate
  output directory, preserving the earlier run.
- **Integrated synthetic BAM publication audit remains unverified.** A
  synthetic integrated launcher completed its branches but failed the final
  BAM publication audit. The adapter was changed to publish BAM/BAI fields
  explicitly from the source map, but the integrated synthetic launcher has
  not been rerun to certify that fix. Successful staged production audits do
  not close this separate integration check.
- **Audit scope is explicit.** The accepted full-depth alignment audit for
  the 16-cell cohort was quick mode with three deeply streamed BAMs. For the
  47-cell run, preflight completion alone is not evidence that any analysis
  task passed.

## RFC process note

The St. Jude Cloud RFC repository describes a draft, stakeholder-discussion,
and broader-community-discussion process. This local document is not an
accepted RFC. Before formal submission, review it with the repository
maintainers; if accepted for discussion, the maintainers assign the RFC
number and the placeholder `0000` should be updated in the title, filename,
and cross-references as appropriate.
