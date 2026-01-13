#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 6_align_trna.sh
# Align rRNA-Bowtie2-unmapped reads to human tRNA using Bowtie1.
#
# Usage:
#   bash scripts/6_align_trna.sh [--qc]
#
# Inputs:
#   - unmapped/rrna_bowtie2/*_after_rrna_bowtie2.fastq.gz
#
# Outputs:
#   - aligned/trna_bowtie1/*.bam
#   - unmapped/trna_bowtie1/*_after_trna.fastq.gz
#   - counts/trna_bowtie1/*_counts.txt
#   - counts/trna_bowtie1/trna_counts_bowtie1.csv
#   - counts/trna_fragment/**
#   - log/trna_bowtie1/*_bowtie1.log
#   - qc/trna_bowtie1/** (optional)
#
# Notes:
#   - Bowtie1: -a --best --strata
#   - NH tags added for multimappers
#   - featureCounts: -M --fraction -O
#   - tRNA fragment classification performed post-alignment
# ------------------------------------------------------------------------------

THREADS=10
DO_QC=false
[[ "${1:-}" == "--qc" ]] && DO_QC=true

# --- References (HUMAN) ---
REF_INDEX="reference/trna/human_tRNA"
SAF_FILE="reference/trna/human_tRNA.saf"
ANNOT_FILE="reference/trna/human_tRNA_annotation.tsv"

# --- Directories ---
INPUT_DIR="unmapped/rrna_bowtie2"
ALIGN_DIR="aligned/trna_bowtie1"
UNMAP_DIR="unmapped/trna_bowtie1"
COUNTS_DIR="counts/trna_bowtie1"
TRF_DIR="counts/trna_fragment"
LOG_DIR="log/trna_bowtie1"
QC_DIR="qc/trna_bowtie1"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$TRF_DIR" \
         "$LOG_DIR" "$QC_FASTQC" "$QC_MULTIQC"

echo "============================================================"
echo "🎯 Step 6: HUMAN tRNA alignment (Bowtie1)"
echo "============================================================"

################################################################################
# Main loop
################################################################################
for fq in "$INPUT_DIR"/*_after_rrna_bowtie2.fastq.gz; do
  [[ -e "$fq" ]] || { echo "❌ No FASTQs found in $INPUT_DIR"; exit 1; }

  base="$(basename "$fq" _after_rrna_bowtie2.fastq.gz)"
  sample="${base#trimmed_}"
  sample_base="${sample%%_after_*}"

  bam_out="${ALIGN_DIR}/${sample_base}_trna_bowtie1.bam"
  unmapped_fq="${UNMAP_DIR}/${sample_base}_after_trna.fastq.gz"
  log_file="${LOG_DIR}/${sample_base}_bowtie1.log"
  fc_out="${COUNTS_DIR}/${sample_base}_trna_bowtie1_counts.txt"

  echo "🧩 Processing sample: $sample_base" | tee -a "$log_file"

  # ------------------------------------------------------------------
  # Bowtie1 alignment
  # ------------------------------------------------------------------
  if [[ -s "$bam_out" ]]; then
    echo "⏩ BAM exists — skipping alignment." | tee -a "$log_file"
  else
    echo "🔧 Bowtie1 alignment..." | tee -a "$log_file"

    ( [[ "$fq" =~ \.gz$ ]] && zcat "$fq" || cat "$fq" ) \
      | bowtie -S -p "$THREADS" -a --best --strata -v 1 \
          -x "$REF_INDEX" -q - \
          --un >(gzip > "$unmapped_fq") 2>> "$log_file" \
      | samtools view -h -@ "$THREADS" - \
      | samtools sort -n -@ "$THREADS" -O SAM - \
      | awk -v OFS="\t" '
          /^@/ { print; next }
          {
            if ($1 != q && q != "") {
              for (i=1; i<=c; i++) print buf[i] "\tNH:i:" c;
              delete buf; c=0;
            }
            q=$1; c++; buf[c]=$0;
          }
          END {
            for (i=1; i<=c; i++) print buf[i] "\tNH:i:" c
          }
        ' \
      | samtools sort -@ "$THREADS" -O BAM -o "$bam_out" -

    samtools index "$bam_out"
    echo "✅ BAM created: $bam_out" | tee -a "$log_file"
  fi

  # ------------------------------------------------------------------
  # featureCounts
  # ------------------------------------------------------------------
  if [[ ! -s "$fc_out" && -s "$bam_out" ]]; then
    echo "🧮 featureCounts (fractional)..." | tee -a "$log_file"

    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
      -M --fraction -O -s 0 -o "$fc_out" "$bam_out" &>> "$log_file" \
      || echo "⚠️ featureCounts failed" | tee -a "$log_file"
  fi

  # ------------------------------------------------------------------
  # Optional QC
  # ------------------------------------------------------------------
  if $DO_QC; then
    qc_zip="${QC_FASTQC}/${sample_base}_after_trna_fastqc.zip"
    if [[ -s "$unmapped_fq" && ! -s "$qc_zip" ]]; then
      echo "🔍 FastQC..." | tee -a "$log_file"
      fastqc -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq" \
        || echo "⚠️ FastQC failed" | tee -a "$log_file"
    fi
  fi
done

################################################################################
# MultiQC
################################################################################
if $DO_QC; then
  if compgen -G "${QC_FASTQC}/*fastqc.zip" > /dev/null; then
    echo "📊 Running MultiQC..."
    multiqc -o "$QC_MULTIQC" "$QC_FASTQC" \
      || echo "⚠️ MultiQC warnings, continuing..."
  else
    echo "⚠️ No FastQC files — skipping MultiQC."
  fi
fi

################################################################################
# Combine counts
################################################################################
combined_csv="${COUNTS_DIR}/trna_counts_bowtie1.csv"
if [[ ! -s "$combined_csv" ]]; then
  echo "🧩 Combining per-sample counts → ${combined_csv}"

  Rscript - <<EOF
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

files <- list.files("${COUNTS_DIR}", pattern="_trna_bowtie1_counts.txt$", full.names=TRUE)
if (length(files)==0) quit(status=0)

ann <- read_tsv("${ANNOT_FILE}", show_col_types=FALSE)

df_list <- lapply(files, function(f){
  sample <- sub("_trna_bowtie1_counts.txt","",basename(f))
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

names(merged) <- sub("_trimmed$", "", names(merged))
write_csv(merged, "${combined_csv}")
EOF
fi

################################################################################
# tRF classification
################################################################################
echo "🔎 Running tRF classification..."
Rscript scripts/classify_tRF.R \
  "$ALIGN_DIR" "$SAF_FILE" "$ANNOT_FILE" \
  "${TRF_DIR}/tRF_counts_bowtie1.csv" "${TRF_DIR}/per_sample"

echo "🎉 HUMAN tRNA Bowtie1 alignment + tRF classification complete!"
