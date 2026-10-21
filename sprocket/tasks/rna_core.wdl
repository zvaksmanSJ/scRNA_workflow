version 1.2

# Prototype tasks for the Sprocket-native DAG. These are not yet output-
# equivalent to the BaseJumper-backed main.wdl workflow; validate and compare
# each output before treating this port as a replacement.
task ValidateInputs {
  input {
    Array[String] sample_ids
    Array[File] r1_fastqs
    Array[File] r2_fastqs
    String image
  }

  command <<<
    set -euo pipefail
    python3 - <<'PY'
import json
import re

sample_ids = json.load(open(r'''~{write_json(sample_ids)}'''))
r1 = json.load(open(r'''~{write_json(r1_fastqs)}'''))
r2 = json.load(open(r'''~{write_json(r2_fastqs)}'''))
if not sample_ids or len(sample_ids) != len(r1) or len(r1) != len(r2):
    raise SystemExit("Provide equal, nonempty sample, R1, and R2 arrays")
if len(sample_ids) != len(set(sample_ids)):
    raise SystemExit("Sample IDs must be unique")
if any(not re.fullmatch(r"[A-Za-z0-9_.-]+", name) for name in sample_ids):
    raise SystemExit("Sample IDs may contain only letters, digits, dot, underscore, and hyphen")
for path in r1 + r2:
    if not path:
        raise SystemExit("FASTQ paths must not be empty")
with open("input_audit.txt", "w") as report:
    report.write(f"PASS\\t{len(sample_ids)} samples\\n")
PY
  >>>

  output {
    File audit = "input_audit.txt"
  }

  requirements {
    container: image
    cpu: 1
    memory: "2 GB"
  }
}

# QC/full-depth paired-end preprocessing, including the configured depth
# branch behavior and fastp outputs.
task Preprocess {
  input {
    String sample_id
    File r1
    File r2
    Int target_individual_reads
    Int seed = 12345
    Int threads = 2
    String image
    File? upstream_audit
  }

  command <<<
    set -euo pipefail
    mkdir -p sampled trimmed qc/fastqc_raw qc/fastqc_trimmed
    gzip -t "~{r1}" "~{r2}"

    audit_path="~{if defined(upstream_audit) then select_first([upstream_audit]) else ""}"
    if [ -n "$audit_path" ]; then
      grep -q '^PASS' "$audit_path" || {
        echo "Upstream audit did not pass" >&2
        exit 2
      }
    fi

    fastqc -t 2 -o qc/fastqc_raw "~{r1}" "~{r2}"
    n1=$(zcat "~{r1}" | awk 'END {if (NR % 4) exit 2; print NR / 4}')
    n2=$(zcat "~{r2}" | awk 'END {if (NR % 4) exit 2; print NR / 4}')
    if [ "$n1" -ne "$n2" ]; then
      echo "Input mate read counts differ: R1=$n1 R2=$n2" >&2
      exit 3
    fi
    total=$((n1 + n2))
    if [ "$total" -le "~{target_individual_reads}" ]; then
      cp "~{r1}" "sampled/~{sample_id}_R1.fastq.gz"
      cp "~{r2}" "sampled/~{sample_id}_R2.fastq.gz"
      final=$total
    else
      proportion=$(awk -v target="~{target_individual_reads}" -v total="$total" \
        'BEGIN {printf "%.18f", target/total}')
      seqkit sample -p "$proportion" -s "~{seed}" -j 2 \
        -o "sampled/~{sample_id}_R1.fastq.gz" "~{r1}"
      seqkit sample -p "$proportion" -s "~{seed}" -j 2 \
        -o "sampled/~{sample_id}_R2.fastq.gz" "~{r2}"
      sampled_pairs=$(zcat "sampled/~{sample_id}_R1.fastq.gz" | awk 'END {if (NR % 4) exit 2; print NR / 4}')
      sampled_pairs_r2=$(zcat "sampled/~{sample_id}_R2.fastq.gz" | awk 'END {if (NR % 4) exit 2; print NR / 4}')
      if [ "$sampled_pairs" -ne "$sampled_pairs_r2" ]; then
        echo "Sampled mate read counts differ" >&2
        exit 4
      fi
      final=$((2 * sampled_pairs))
    fi

    fastp --thread "~{threads}" \
      --in1 "sampled/~{sample_id}_R1.fastq.gz" \
      --in2 "sampled/~{sample_id}_R2.fastq.gz" \
      --out1 "trimmed/~{sample_id}_R1_trim.fastq.gz" \
      --out2 "trimmed/~{sample_id}_R2_trim.fastq.gz" \
      --json "qc/~{sample_id}_fastp.json" \
      --html "qc/~{sample_id}_fastp.html" \
      --detect_adapter_for_pe \
      2> "qc/~{sample_id}_fastp.log"
    gzip -t "trimmed/~{sample_id}_R1_trim.fastq.gz" \
      "trimmed/~{sample_id}_R2_trim.fastq.gz"
    fastqc -t 2 -o qc/fastqc_trimmed \
      "trimmed/~{sample_id}_R1_trim.fastq.gz" \
      "trimmed/~{sample_id}_R2_trim.fastq.gz"
    printf '%s\n%s\n' "$total" "$final" > "qc/~{sample_id}_read_counts.txt"
  >>>

  output {
    File sampled_r1 = "sampled/~{sample_id}_R1.fastq.gz"
    File sampled_r2 = "sampled/~{sample_id}_R2.fastq.gz"
    File trimmed_r1 = "trimmed/~{sample_id}_R1_trim.fastq.gz"
    File trimmed_r2 = "trimmed/~{sample_id}_R2_trim.fastq.gz"
    File fastp_json = "qc/~{sample_id}_fastp.json"
    File fastp_html = "qc/~{sample_id}_fastp.html"
    File read_counts = "qc/~{sample_id}_read_counts.txt"
    Array[File] raw_fastqc_html = glob("qc/fastqc_raw/*_fastqc.html")
    Array[File] raw_fastqc_zip = glob("qc/fastqc_raw/*_fastqc.zip")
    Array[File] trimmed_fastqc_html = glob("qc/fastqc_trimmed/*_fastqc.html")
    Array[File] trimmed_fastqc_zip = glob("qc/fastqc_trimmed/*_fastqc.zip")
  }

  requirements {
    container: image
    cpu: threads
    memory: "8 GB"
  }
}

# Creates the canonical STAR RNA alignment outputs consumed by counting and QC.
task Align {
  input {
    String sample_id
    File r1
    File r2
    Directory star_index
    Int threads = 8
    String image
  }

  command <<<
    set -euo pipefail
    mkdir -p star bam
    STAR \
      --genomeDir "~{star_index}" \
      --runThreadN "~{threads}" \
      --limitBAMsortRAM 8000000000 \
      --outSAMunmapped Within KeepPairs \
      --twopassMode Basic \
      --outReadsUnmapped None \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "star/~{sample_id}_" \
      --readFilesCommand zcat \
      --readFilesIn "~{r1}" "~{r2}" \
      --outSAMstrandField intronMotif \
      --outSAMattributes All \
      --outFilterMultimapNmax 10 \
      --outSAMprimaryFlag OneBestScore \
      --chimSegmentMin 12 \
      --chimJunctionOverhangMin 8 \
      --chimOutJunctionFormat 1 \
      --outSAMattrRGline "ID:~{sample_id}"
    bam_source="star/~{sample_id}_Aligned.sortedByCoord.out.bam"
    samtools view -F 256 -h "$bam_source" \
      | sed -e 's/\tSM:[^\t]*//' -e 's/\tPL:[^\t]*//' \
            -e 's/^\(@RG\tID:[^\t]*\)/\1\tSM:~{sample_id}\tPL:ILLUMINA/' \
      | samtools view -h -b - \
      | samtools sort -@ 2 -O BAM -o "bam/~{sample_id}.bam" -
    samtools index -@ 2 "bam/~{sample_id}.bam"
    samtools quickcheck -v "bam/~{sample_id}.bam"
    samtools flagstat "bam/~{sample_id}.bam" > "bam/~{sample_id}.flagstat.txt"
  >>>

  output {
    File bam = "bam/~{sample_id}.bam"
    File bai = "bam/~{sample_id}.bam.bai"
    File flagstat = "bam/~{sample_id}.flagstat.txt"
    File star_log_final = "star/~{sample_id}_Log.final.out"
    File chimeric_junctions = "star/~{sample_id}_Chimeric.out.junction"
  }

  requirements {
    container: image
    cpu: threads
    memory: "64 GB"
  }
}

# Produces per-sample gene counts and supporting alignment/QC metrics.
task CountAndQc {
  input {
    String sample_id
    File bam
    File bai
    File gtf
    File ref_flat
    File genebody_bed
    Int threads = 2
    String image
  }

  command <<<
    set -euo pipefail
    mkdir -p counts qualimap picard gene_body
    samtools quickcheck -v "~{bam}"
    htseq-count \
      -f bam -r pos -s no -t exon -i gene_id --additional-attr=gene_name -m union \
      --secondary-alignments=ignore \
      "~{bam}" "~{gtf}" > "counts/~{sample_id}.htseq_counts.tsv"
    unset DISPLAY
    mkdir -p tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
    qualimap --java-mem-size=16G rnaseq \
      -gtf "~{gtf}" -bam "~{bam}" --sorted \
      -outfile "~{sample_id}.report.pdf" -oc "~{sample_id}.oc" \
      -outdir "qualimap/~{sample_id}"
    samtools flagstat "~{bam}" > "picard/~{sample_id}.flagstat.txt"
    picard -Xmx8g CollectRnaSeqMetrics \
      I="~{bam}" O="picard/~{sample_id}.rna_metrics" \
      REF_FLAT="~{ref_flat}" STRAND_SPECIFICITY=NONE
    geneBody_coverage.py -r "~{genebody_bed}" -i "~{bam}" -o "gene_body/~{sample_id}"
  >>>

  output {
    File htseq_counts = "counts/~{sample_id}.htseq_counts.tsv"
    Directory qualimap_output = "qualimap/~{sample_id}"
    File picard_metrics = "picard/~{sample_id}.rna_metrics"
    File flagstat = "picard/~{sample_id}.flagstat.txt"
    File gene_body_coverage = "gene_body/~{sample_id}.geneBodyCoverage.txt"
    File gene_body_plot = "gene_body/~{sample_id}.geneBodyCoverage.curves.pdf"
  }

  requirements {
    container: image
    cpu: threads
    memory: "24 GB"
  }
}

# Gathers per-sample counts into the cohort matrix expected by downstream QC.
task BuildMatrix {
  input {
    Array[String] sample_ids
    Array[File] htseq_counts
    String image
  }

  command <<<
    set -euo pipefail
    python3 - <<'PY'
import csv
import json
from pathlib import Path
import pandas as pd

samples = json.load(open(r'''~{write_json(sample_ids)}'''))
files = json.load(open(r'''~{write_json(htseq_counts)}'''))
if len(samples) != len(files) or not samples:
    raise SystemExit("Sample IDs and HTSeq files must be nonempty and aligned")
columns = {}
symbols = {}
for sample, path in zip(samples, files):
    values = {}
    with open(path, newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 3:
                raise SystemExit(f"Malformed HTSeq row in {path}: {row!r}")
            gene_id, symbol, raw_count = row[0], row[1], row[-1]
            if gene_id.startswith("__"):
                continue
            try:
                count = int(raw_count)
            except ValueError as error:
                raise SystemExit(f"Invalid HTSeq count in {path}: {row!r}") from error
            if count < 0:
                raise SystemExit(f"Negative HTSeq count in {path}: {row!r}")
            values[gene_id] = count
            symbols[gene_id] = symbol
    columns[sample] = values
matrix = pd.DataFrame(columns).fillna(0).astype("int64")
matrix = matrix.loc[matrix.sum(axis=1) > 0]
matrix.insert(0, "gene_name", [symbols.get(gene, "") for gene in matrix.index])
matrix.insert(0, "gene_id", matrix.index)
matrix.to_csv("matrix_gene_counts_starhtseq.txt", sep="\t", index=False)
long = matrix.melt(id_vars=["gene_id", "gene_name"], var_name="sample_id", value_name="count")
long.to_parquet("matrix_gene_counts_starhtseq.parquet", index=False)
metrics = []
for sample in samples:
    total = int(pd.to_numeric(matrix[sample]).sum())
    mt = int(sum(int(matrix.loc[i, sample]) for i in matrix.index
                 if str(matrix.loc[i, "gene_name"]).upper().startswith("MT-")))
    metrics.append({"sample_id": sample, "assigned_counts": total,
                    "mitochondrial_counts": mt,
                    "mitochondrial_fraction": mt / total if total else 0.0})
pd.DataFrame(metrics).to_csv("sample_expression_metrics.tsv", sep="\t", index=False)
Path("audit.txt").write_text(f"PASS\\t{len(samples)} samples\\t{len(matrix)} genes\\n")
PY
  >>>

  output {
    File matrix = "matrix_gene_counts_starhtseq.txt"
    File matrix_parquet = "matrix_gene_counts_starhtseq.parquet"
    File sample_metrics = "sample_expression_metrics.tsv"
    File audit = "audit.txt"
  }

  requirements {
    container: image
    cpu: 2
    memory: "16 GB"
  }
}

# Collects tool reports into a run-level MultiQC report.
task MultiqcReport {
  input {
    Array[File] fastp_jsons
    Array[File] fastp_htmls
    Array[File] fastqc_html
    Array[File] fastqc_zip
    Array[File] star_logs
    Array[Directory] qualimap_dirs
    String image
  }

  command <<<
    set -euo pipefail
    mkdir -p reports
    multiqc --force --outdir reports --filename multiqc_report.html \
      ~{sep(" ", fastp_jsons)} ~{sep(" ", fastp_htmls)} \
      ~{sep(" ", fastqc_html)} ~{sep(" ", fastqc_zip)} \
      ~{sep(" ", star_logs)} ~{sep(" ", qualimap_dirs)}
  >>>

  output {
    File report = "reports/multiqc_report.html"
  }

  requirements {
    container: image
    cpu: 2
    memory: "8 GB"
  }
}

# Quantifies transcripts per sample from full-depth trimmed reads.
task SalmonQuant {
  input {
    String sample_id
    File r1
    File r2
    Directory salmon_index
    Int threads = 8
    String image
    File upstream_gate
  }

  command <<<
    set -euo pipefail
    grep -q '^PASS' "~{upstream_gate}"
    salmon quant -i "~{salmon_index}" -l A \
      --validateMappings --gcBias \
      -1 "~{r1}" -2 "~{r2}" -p "~{threads}" \
      -o "salmon/~{sample_id}"
    test -s "salmon/~{sample_id}/quant.sf"
  >>>

  output {
    File quant = "salmon/~{sample_id}/quant.sf"
    File metadata = "salmon/~{sample_id}/aux_info/meta_info.json"
  }

  requirements {
    container: image
    cpu: threads
    memory: "16 GB"
  }
}

# Stops dependent tasks when an upstream branch audit is absent or not PASS.
task BranchGate {
  input {
    File audit
    File multiqc_report
    String image
  }

  command <<<
    set -euo pipefail
    grep -q '^PASS' "~{audit}"
    test -s "~{multiqc_report}"
    printf 'PASS\\tbranch audit and MultiQC report\\n' > branch_gate.txt
  >>>

  output {
    File gate = "branch_gate.txt"
  }

  requirements {
    container: image
    cpu: 1
    memory: "2 GB"
  }
}
