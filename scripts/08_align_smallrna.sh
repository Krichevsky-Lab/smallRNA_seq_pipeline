#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 8_align_smallrna.sh
#
# Align tRNA-rescue–unmapped reads to human small RNA references using Bowtie1.
#
# Usage:
#   bash scripts/8_align_smallrna.sh [--qc]
#
# Inputs:
#   - unmapped/trna_bowtie2/*_after_trna_bowtie2.fastq.gz
#
# Outputs:
#   - aligned/smallrna_bowtie1/*.bam
#   - unmapped/smallrna_bowtie1/*_after_smallrna.fastq.gz
#   - counts/smallrna_bowtie1/*_counts.txt
#   - counts/smallrna_bowtie1/smallrna_counts_bowtie1.csv
#   - log/smallrna_bowtie1/*_bowtie1.log
#   - qc/smallrna_bowtie1/** (optional)
#
# Notes:
#   - Bowtie1 best/strata alignment
#   - NH tags added for multimappers (fractional counting)
#   - featureCounts: -M --fraction -O
#   - Optional FastQC + MultiQC
# ------------------------------------------------------------------------------

THREADS=10
DO_QC=false
[[ "${1:-}" == "--qc" ]] && DO_QC=true

# --- References ---
REF_INDEX="reference/smallrna/human_smallRNA"
SAF_FILE="reference/smallrna/human_smallRNA.saf"
ANNOT_FILE="reference/smallrna/human_smallRNA_annotation.tsv"

# --- Directories ---
INPUT_DIR="unmapped/trna_bowtie2"
ALIGN_DIR="aligned/smallrna_bowtie1"
UNMAP_DIR="unmapped/smallrna_bowtie1"
COUNTS_DIR="counts/smallrna_bowtie1"
LOG_DIR="log/smallrna_bowtie1"
QC_DIR="qc/smallrna_bowtie1"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$LOG_DIR" "$QC_FASTQC" "$QC_MULTIQC"

echo "============================================================"
echo "🎯 Step 8: Align to human smallRNA (Bowtie1)"
echo "============================================================"

for fq in "$INPUT_DIR"/*_after_trna_bowtie2.fastq.gz; do
  [[ -e "$fq" ]] || { echo "❌ No FASTQ files found"; exit 1; }

  sample=$(basename "$fq" _after_trna_bowtie2.fastq.gz)
  sample=${sample#trimmed_}
  bam_out="${ALIGN_DIR}/${sample}_smallrna.bam"
  unmapped_fq="${UNMAP_DIR}/${sample}_after_smallrna.fastq.gz"
  log_file="${LOG_DIR}/${sample}_bowtie1.log"

  # --- Bowtie alignment ---
  if [[ -s "$bam_out" ]]; then
    echo "⏩ $sample BAM exists"
  else
    echo "🔧 Aligning $sample..."
    bowtie -S -p "$THREADS" -a --best --strata -v 1 \
      -x "$REF_INDEX" -q "$fq" \
      --un >(gzip > "$unmapped_fq") 2> "$log_file" \
      | samtools view -h - \
      | awk 'BEGIN{OFS="\t"} /^@/ {print; next} {if($0==""||$1=="*"||NF<11) next; print}' \
      | samtools view -bS - \
      | samtools sort -@ "$THREADS" -o "$bam_out" - 2>> "$log_file"

    if [[ ! -s "$bam_out" || $(stat -c%s "$bam_out") -lt 500 ]]; then
      echo "⚠️ Empty BAM for $sample" | tee -a "$log_file"
      rm -f "$bam_out"
      continue
    fi

    echo "🔖 Adding NH tags..."
    tmp_sam="${LOG_DIR}/${sample}_tmp.sam"
    tmp_best="${LOG_DIR}/${sample}_tmp_best.sam"

    samtools sort -@ "$THREADS" -n -O SAM "$bam_out" > "$tmp_sam"

    awk -v OFS="\t" '
      /^@/ {print; next}
      {
        if ($1 != q && q != "") {
          for (i=1;i<=c;i++) print buf[i] "\tNH:i:" c
          delete buf; c=0
        }
        q=$1; buf[++c]=$0
      }
      END {
        for (i=1;i<=c;i++) print buf[i] "\tNH:i:" c
      }' "$tmp_sam" > "$tmp_best"

    samtools sort -@ "$THREADS" -O BAM -o "$bam_out" "$tmp_best"
    samtools index "$bam_out"
    rm -f "$tmp_sam" "$tmp_best"
  fi

  # --- featureCounts ---
  fc_out="${COUNTS_DIR}/${sample}_smallrna_counts.txt"
  if [[ ! -s "$fc_out" && -s "$bam_out" ]]; then
    echo "🧮 Counting $sample..."
    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
      -M --fraction -O -s 0 -o "$fc_out" "$bam_out" &>> "$log_file"
  fi

  # --- QC ---
  if $DO_QC; then
    qc_zip="${QC_FASTQC}/${sample}_after_smallrna_fastqc.zip"
    if [[ -s "$unmapped_fq" && ! -s "$qc_zip" ]]; then
      fastqc -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq"
    fi
  fi
done

# --- MultiQC ---
if $DO_QC; then
  if compgen -G "${QC_FASTQC}/*fastqc.zip" > /dev/null; then
    multiqc -o "$QC_MULTIQC" "$QC_FASTQC"
  fi
fi

# --- Combine counts ---
combined_csv="${COUNTS_DIR}/smallrna_counts_bowtie1.csv"
if [[ ! -s "$combined_csv" ]]; then
  echo "🧩 Combining counts..."

  Rscript - <<EOF
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

files <- list.files("${COUNTS_DIR}", pattern="_smallrna_counts.txt$", full.names=TRUE)
if (length(files)==0) quit("no", status=0)

ann <- read_tsv("${ANNOT_FILE}", show_col_types=FALSE)

df_list <- lapply(files, function(f){
  sample <- sub("_smallrna_counts.txt","",basename(f))
  # remove leading "trimmed_" if present
  sample <- str_remove(sample,"^trimmed_")
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
  select(genename, geneid, biotype_1, biotype_2, chromosome,
         start, end, strand, length, source,
         everything(), -fc_length)

# ---- REMOVE _trimmed FROM SAMPLE COLUMNS ----
clean_names <- names(merged)
clean_names <- sub("_trimmed$", "", clean_names)
names(merged) <- clean_names

write_csv(merged, "${combined_csv}")
EOF
fi

echo "✅ Human smallRNA Bowtie1 alignment complete"
