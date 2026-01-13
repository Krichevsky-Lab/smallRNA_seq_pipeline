#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 7_align_trna_rescue.sh
#
# Rescue tRNA-unmapped reads using Bowtie2 (local alignment).
#
# Usage:
#   bash scripts/7_align_trna_rescue.sh [--qc]
#
# Inputs:
#   - unmapped/trna_bowtie1/*_after_trna.fastq.gz
#
# Outputs:
#   - aligned/trna_bowtie2/*.bam
#   - unmapped/trna_bowtie2/*_after_trna_bowtie2.fastq.gz
#   - counts/trna_bowtie2/*_counts.txt
#   - counts/trna_bowtie2/trna_counts_bowtie2.csv
#   - counts/trna_fragment/**
#   - log/trna_bowtie2/*_bowtie2.log
#   - qc/trna_bowtie2/** (optional)
#
# Notes:
#   - Bowtie2 local, very-sensitive-local
#   - Retains best AS:i alignments (ties kept)
#   - NH tags added for fractional counting
#   - featureCounts: -M --fraction -O
#   - tRF classification performed post-alignment
# ------------------------------------------------------------------------------

THREADS=10
DO_QC=false
[[ "${1:-}" == "--qc" ]] && DO_QC=true

# ------------------------------------------------------------------
# References (HUMAN)
# ------------------------------------------------------------------
REF_INDEX="reference/trna/human_tRNA"
SAF_FILE="reference/trna/human_tRNA.saf"
ANNOT_FILE="reference/trna/human_tRNA_annotation.tsv"

# ------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------
INPUT_DIR="unmapped/trna_bowtie1"
ALIGN_DIR="aligned/trna_bowtie2"
UNMAP_DIR="unmapped/trna_bowtie2"
COUNTS_DIR="counts/trna_bowtie2"
LOG_DIR="log/trna_bowtie2"
QC_DIR="qc/trna_bowtie2"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"
TRF_DIR="counts/trna_fragment"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$LOG_DIR" \
         "$QC_FASTQC" "$QC_MULTIQC" "$TRF_DIR"

echo "============================================================"
echo "🎯 Step 7: HUMAN tRNA rescue alignment (Bowtie2 local)"
echo "============================================================"

################################################################################
# Main loop
################################################################################
for fq in "$INPUT_DIR"/*_after_trna.fastq.gz; do
  [[ -e "$fq" ]] || { echo "❌ No FASTQ files found in $INPUT_DIR"; exit 1; }

  base="$(basename "$fq" _after_trna.fastq.gz)"
  sample="${base#trimmed_}"
  sample_base="${sample%%_after_*}"

  bam_out="${ALIGN_DIR}/${sample_base}_trna_bowtie2.bam"
  unmapped_fq="${UNMAP_DIR}/${sample_base}_after_trna_bowtie2.fastq.gz"
  log_file="${LOG_DIR}/${sample_base}_bowtie2.log"
  tmp_sam="${ALIGN_DIR}/${sample_base}_tmp.sam"
  tmp_best="${ALIGN_DIR}/${sample_base}_tmp_best.sam"

  echo "🧩 Processing sample: $sample_base" | tee -a "$log_file"

  # ------------------------------------------------------------------
  # Bowtie2 rescue alignment
  # ------------------------------------------------------------------
  if [[ -s "$bam_out" ]]; then
    echo "⏩ BAM exists — skipping Bowtie2 rescue." | tee -a "$log_file"
  else
    echo "🔧 Bowtie2 rescue alignment (local)..." | tee -a "$log_file"

    bowtie2 -p "$THREADS" \
      --local --very-sensitive-local -L 15 -N 1 \
      -D 20 -R 3 -i S,1,0.60 \
      --ma 2 --mp 2,2 --rdg 3,2 --rfg 3,2 \
      --score-min L,0,0.8 \
      --no-unal --all \
      -x "$REF_INDEX" -U "$fq" --un-gz "$unmapped_fq" \
      -S "$tmp_sam" 2>> "$log_file"

    if [[ ! -s "$tmp_sam" ]]; then
      echo "⚠️ No alignments — skipping sample." | tee -a "$log_file"
      continue
    fi

    # --------------------------------------------------------------
    # Retain best AS:i alignments + add NH tags
    # --------------------------------------------------------------
    awk -v OFS="\t" '
      BEGIN { q=""; c=0; maxAS=-999999 }
      /^@/ { print; next }
      {
        if (match($0,/AS:i:([0-9-]+)/,m)) ascore=m[1]; else ascore=-999999;
        if (q != "" && $1 != q) {
          for (i=1;i<=c;i++) if (as[i]==maxAS) best[i]=1;
          nBest=0; for (i=1;i<=c;i++) if (best[i]==1) nBest++;
          for (i=1;i<=c;i++) if (best[i]==1) print buf[i] "\tNH:i:" nBest;
          delete buf; delete as; delete best; c=0; maxAS=-999999;
        }
        q=$1; c++; buf[c]=$0; as[c]=ascore;
        if (ascore>maxAS) maxAS=ascore;
      }
      END {
        for (i=1;i<=c;i++) if (as[i]==maxAS) best[i]=1;
        nBest=0; for (i=1;i<=c;i++) if (best[i]==1) nBest++;
        for (i=1;i<=c;i++) if (best[i]==1) print buf[i] "\tNH:i:" nBest;
      }
    ' "$tmp_sam" > "$tmp_best"

    samtools view -@ "$THREADS" -bS "$tmp_best" \
      | samtools sort -@ "$THREADS" -o "$bam_out"

    samtools index "$bam_out"
    rm -f "$tmp_sam" "$tmp_best"

    echo "✅ BAM created: $bam_out" | tee -a "$log_file"
  fi

  # ------------------------------------------------------------------
  # featureCounts
  # ------------------------------------------------------------------
  fc_out="${COUNTS_DIR}/${sample_base}_trna_bowtie2_counts.txt"
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
    qc_zip="${QC_FASTQC}/${sample_base}_after_trna_bowtie2_fastqc.zip"
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
      || echo "⚠️ MultiQC encountered issues."
  else
    echo "⚠️ No FastQC files — skipping MultiQC."
  fi
fi

################################################################################
# Combine counts and annotate
################################################################################
combined_csv="${COUNTS_DIR}/trna_counts_bowtie2.csv"
if [[ ! -s "$combined_csv" ]]; then
  echo "🧩 Combining per-sample rescue counts → ${combined_csv}"

  Rscript - <<EOF
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

files <- list.files("${COUNTS_DIR}", pattern="_trna_bowtie2_counts.txt$", full.names=TRUE)
if (length(files)==0) quit(status=0)

ann <- read_tsv("${ANNOT_FILE}", show_col_types=FALSE)

df_list <- lapply(files, function(f){
  sample <- sub("_trna_bowtie2_counts.txt","",basename(f))
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
  "${TRF_DIR}/tRF_counts_bowtie2.csv" "${TRF_DIR}/per_sample_bowtie2"

echo "🎉 HUMAN tRNA Bowtie2 rescue alignment + tRF classification complete!"
