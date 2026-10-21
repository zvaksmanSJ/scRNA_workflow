# Python/LSF FASTQ utilities

This directory contains sample-sheet validation, raw FASTQ planning and
per-cell merging, and cohort/per-cell integrity audits. Inputs are read-only;
merged FASTQs and audit reports are written to paths supplied by the operator.
Worker scripts use `RNA_PYTHON` or the runtime installed under
`RNA_RUNTIME_DIR`.

## Sample sheet

The canonical tab-separated columns are `sample_id`, `r1`, `r2`, and `group`.
Sample identifiers must be unique and shell-safe, and both FASTQ paths must
exist. A generic template is in `../examples/rna_cells.example.tsv`.

## Explicit LSF audit submission

The audit helper requires the sample manifest, raw-input plan, and a new audit
directory. It is an LSF-specific convenience, not part of the Cromwell WDL
stage graph.

```bash
bash python_lsf_pipeline/submit_fastq_audit.sh --submit \
  --manifest /path/to/rna_cells.tsv \
  --raw-plan /path/to/raw_inputs.json \
  --audit-dir /path/to/new/audit/fastq \
  --queue standard
```

After all array elements complete, combine and validate their reports with
`audit_expanded_fastqs.py --combine-report`. A successful submission is not an
audit pass; inspect every task and require a passing aggregate report before
downstream processing.
