# RNA analyses and result utilities

This directory contains the Salmon quantification planner/validator, the
downstream cohort analysis and dashboard builder, the result-tree organizer,
and the BaseJumper-output audit/dashboard utility. These are post-alignment
analyses and reporting helpers, not the low-level Sprocket task definitions.

## RNA-only research extensions

`run_rna_research_addons.py` implements a separately gated runner for expanded
QC, exon-bin counts, descriptive event PSI, and fusion discovery. The module has
not been run against a cohort or validated with a fusion-positive/negative
fixture; the optional environment has not been installed or solved. Do not infer
that its outputs exist or are production-validated.

The runner is separate from the standard BaseJumper outputs:

- retain current BaseJumper-derived QC and expression outputs unchanged;
- write to `final_results/10_rna_research_addons/` and refuse to overwrite an
  existing add-on directory;
- use RSeQC TIN and positional duplication (an estimate, not a UMI molecule
  count), DEXSeq flattened exon-bin counts, and SUPPA2 event PSI; no
  differential test is run for the one-condition cohort;
- use an independent STAR fusion alignment plus Arriba, preserving passed and
  discarded calls and logs without changing canonical expression BAMs;
- refresh the dashboard only after a successful add-on run to populate its
  Enhanced QC, Exon & splicing, and Fusions tabs.

First install the optional environment with `setup/setup.sh --research-addons`.
Run the script without `--execute` to print a planned command manifest, then
review the references, inputs, resource needs, and resolved tools. The
execution command requires the completed full-depth results directory, matching
GRCh38/Ensembl 104 references, the STAR index, gene-body BED, and Arriba
blacklist/known-fusion/protein-domain resources:

```bash
conda run -p "${RNA_RUNTIME_DIR:-wdl_runtime}/env/resolveome-rna-research" \
  python extra_methods/run_rna_research_addons.py \
  --results-dir /path/to/completed/RNA_results \
  --gtf /path/to/Ensembl104.gtf \
  --fasta /path/to/GRCh38.fa \
  --star-index /path/to/STAR_Ensembl104 \
  --genebody-bed /path/to/GRCh38_gene_body.bed \
  --arriba-blacklist /path/to/arriba_blacklist.tsv \
  --arriba-known-fusions /path/to/known_fusions.tsv \
  --arriba-protein-domains /path/to/protein_domains.gff3 \
  --threads 8 --execute
```

This stand-alone runner currently executes its stages sequentially; it is not
yet a scatter-based Sprocket/WDL stage. Do not submit it as a production cohort
run until the optional environment, resource profile, and synthetic fusion and
splicing outputs have been validated. To display completed add-on results,
refresh the dashboard without rerunning the standard analysis:

```bash
conda run -p "${RNA_RUNTIME_DIR:-wdl_runtime}/env/resolveome-rna" \
  python extra_methods/run_rna_downstream.py --refresh-dashboard
```

No RNA-to-DNA identity comparison, RNA variant calling, clinical validation,
or clinical interpretation is included. Those require separately agreed
matched DNA/RNA inputs or clinical-group review.

The Cromwell WDL invokes `quantify_transcripts.py`,
`run_rna_downstream.py`, and `organize_rna_results.py`. Their paths are
resolved from `RESOLVEOME_RNA_WORKFLOW`, `RNA_RUNTIME_DIR`,
`RESOLVEOME_RNA_RESULTS`, and the configured reference tree. The optional
`run_rna_downstream_after_salmon.sh` requires an explicit completed results
directory through `RESOLVEOME_RNA_RESULTS`; it does not default to a
historical cohort.

For new runs, use an isolated results tree and require passing upstream
audits. Do not rerun analysis against completed results except for a
specifically read-only audit or dashboard refresh.
