#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

################################################################################
# 12_align_longrna.sh
# Step 12: Align piRNA-rescue–unmapped reads to human longRNA (Bowtie1)
#
# Usage:
#   bash scripts/12_align_longrna.sh [--qc]
#
# Description:
#   - Aligns reads unmapped after piRNA Bowtie2 rescue to long RNA reference
#   - Reports all best alignments (--best --strata -a)
#   - Adds NH tags for multimappers (fractional assignment)
#   - Runs featureCounts (-M --fraction -O)
#   - Optional FastQC + MultiQC on longRNA-unmapped reads
#   - Uses a single per-sample log combining alignment + featureCounts output
#
# Inputs:
#   unmapped/pirna_bowtie2/*_after_pirna_bowtie2.fastq.gz
#
# Outputs:
#   aligned/longrna_bowtie1/*.bam
#   unmapped/longrna_bowtie1/*_after_longrna.fastq.gz
#   counts/longrna_bowtie1/*_longrna_counts.txt
#   counts/longrna_bowtie1/longrna_counts_bowtie1.csv
#   log/longrna_bowtie1/*.log
#   qc/longrna_bowtie1/{fastqc,multiqc}
################################################################################

THREADS=10
DO_QC=false
[[ "${1:-}" == "--qc" ]] && DO_QC=true

# --- References (HUMAN) ---
REF_INDEX="reference/longrna/human_longRNA"
SAF_FILE="reference/longrna/human_longRNA.saf"
ANNOT_FILE="reference/longrna/human_longRNA_annotation.tsv"

# --- Directories ---
INPUT_DIR="unmapped/pirna_bowtie2"
ALIGN_DIR="aligned/longrna_bowtie1"
UNMAP_DIR="unmapped/longrna_bowtie1"
COUNTS_DIR="counts/longrna_bowtie1"
LOG_DIR="log/longrna_bowtie1"
QC_DIR="qc/longrna_bowtie1"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$LOG_DIR"
$DO_QC && mkdir -p "$QC_FASTQC" "$QC_MULTIQC"

inputs=( "$INPUT_DIR"/*_after_pirna_bowtie2.fastq.gz )
[[ ${#inputs[@]} -gt 0 ]] || { echo "❌ No FASTQ files found in $INPUT_DIR"; exit 1; }

echo "============================================================"
echo "🎯 Step 12: Aligning human longRNA (Bowtie1)"
echo "QC enabled: $DO_QC"
echo "============================================================"

################################################################################
# Main loop
################################################################################
for fq in "${inputs[@]}"; do
  sample=$(basename "$fq" "_after_pirna_bowtie2.fastq.gz")
  sample="${sample#trimmed_}"

  bam_out="${ALIGN_DIR}/${sample}_longrna.bam"
  unmapped_fq="${UNMAP_DIR}/${sample}_after_longrna.fastq.gz"
  log_file="${LOG_DIR}/${sample}.log"

  echo "--------------------------------------------"
  echo "🧩 Processing sample: $sample"
  echo "--------------------------------------------"

  if [[ ! -s "$bam_out" ]]; then
    echo "🔧 Bowtie1 alignment…" | tee -a "$log_file"

    bowtie -S -p "$THREADS" -a --best --strata -v 1 \
      -x "$REF_INDEX" -q "$fq" \
      --un >(gzip > "$unmapped_fq") 2>> "$log_file" \
    | samtools view -h - \
    | awk 'BEGIN{OFS="\t"} /^@/ {print; next} {if(NF>=11) print}' \
    | samtools view -bS - \
    | samtools sort -@ "$THREADS" -o "$bam_out" - 2>> "$log_file"

    if [[ ! -s "$bam_out" ]]; then
      echo "⚠️ No alignments for $sample — skipping." | tee -a "$log_file"
      continue
    fi

    echo "🔖 Adding NH tags…" | tee -a "$log_file"

    # ---- NH tagging (QNAME-grouped) ----
    samtools view -h "$bam_out" | \
      awk 'BEGIN{OFS="\t"}
           /^@/ {print; next}
           {q=$1; c[q]++; buf[q]=buf[q]"\n"$0}
           END{
             for(q in buf){
               n=c[q];
               split(buf[q], a, "\n");
               for(i=2;i<=length(a);i++){
                 line=a[i];
                 if(line=="") continue;
                 gsub(/\tNH:i:[0-9]+/,"",line);
                 print line "\tNH:i:" n;
               }
             }
           }' \
    | samtools view -bS - \
    | samtools sort -@ "$THREADS" -o "${bam_out%.bam}_NH.bam" -

    mv "${bam_out%.bam}_NH.bam" "$bam_out"
    samtools index "$bam_out"

    echo "✅ BAM ready: $bam_out" | tee -a "$log_file"
  else
    echo "⏩ BAM exists → skipping alignment." | tee -a "$log_file"
  fi

  # ------------------------------------------------------------------
  # featureCounts
  # ------------------------------------------------------------------
  fc_out="${COUNTS_DIR}/${sample}_longrna_counts.txt"
  if [[ ! -s "$fc_out" && -s "$bam_out" ]]; then
    echo "🧮 featureCounts…" | tee -a "$log_file"

    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
      -M --fraction -O -s 0 \
      -o "$fc_out" "$bam_out" &>> "$log_file" \
      || echo "⚠️ featureCounts failed for $sample" | tee -a "$log_file"
  fi

  # ------------------------------------------------------------------
  # QC
  # ------------------------------------------------------------------
  if $DO_QC; then
    qc_zip="${QC_FASTQC}/${sample}_after_longrna_fastqc.zip"
    if [[ -s "$unmapped_fq" && ! -s "$qc_zip" ]]; then
      echo "🔍 FastQC…" | tee -a "$log_file"
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq" || \
        echo "⚠️ FastQC failed" | tee -a "$log_file"
    fi
  fi
done

################################################################################
# Combine counts and annotate
################################################################################
combined_csv="${COUNTS_DIR}/longrna_counts_bowtie1.csv"
echo "📊 Combining counts → $combined_csv"

Rscript - <<EOF
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

files <- list.files("${COUNTS_DIR}", pattern="_longrna_counts.txt$", full.names=TRUE)
if (length(files)==0) quit(status=0)

ann <- read_tsv("${ANNOT_FILE}", show_col_types=FALSE)

df_list <- lapply(files, function(f){
  sample <- sub("_longrna_counts.txt","",basename(f))
  dat <- read_tsv(f, comment="#", show_col_types=FALSE)
  count_col <- names(dat)[ncol(dat)]
  dat %>%
    rename(genename = Geneid) %>%
    select(genename, Length, !!sym(count_col)) %>%
    rename(fc_length = Length,
           !!sample := !!sym(count_col))
})

fc <- Reduce(function(x,y) full_join(x,y, by=c("genename","fc_length")), df_list)

merged <- fc %>%
  left_join(ann, by="genename") %>%
  mutate(length = coalesce(length, fc_length)) %>%
  select(genename, geneid, biotype_1, biotype_2,
         chromosome, start, end, strand, length, source,
         everything(), -fc_length)

names(merged) <- sub("_trimmed$", "", names(merged))
write_csv(merged, "${combined_csv}")
EOF

################################################################################
# MultiQC
################################################################################
if $DO_QC; then
  echo "📈 Running MultiQC…"
  multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force || \
    echo "⚠️ MultiQC warnings"
fi

echo "🎉 Human longRNA Bowtie1 alignment complete"