#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 4_align_rrna.sh
# Align UniVec-unmapped reads to human rRNA using Bowtie1.
#
# Usage:
#   bash scripts/4_align_rrna.sh [--qc]
#
# Inputs:
#   - unmapped/univec/*_after_univec.fq.gz
#
# Outputs:
#   - aligned/rrna_bowtie1/*.bam
#   - unmapped/rrna_bowtie1/*.fastq.gz
#   - counts/rrna_bowtie1/*_counts.txt
#   - counts/rrna_bowtie1/rrna_counts_bowtie1.csv
#   - log/rrna_bowtie1/*_bowtie1.log
#   - log/rrna_bowtie1/read_counts.csv
#   - qc/rrna_bowtie1/** (optional)
#
# Notes:
#   - Bowtie1: -a --best --strata
#   - NH tags added for multimappers
#   - featureCounts: -M --fraction -O
# ------------------------------------------------------------------------------

THREADS=10
DO_QC=false
[[ "${1:-}" == "--qc" ]] && DO_QC=true

# --- References (HUMAN) ---
REF_INDEX="reference/rrna/human_rRNA"
SAF_FILE="reference/rrna/human_rRNA.saf"
ANNOT_FILE="reference/rrna/human_rRNA_annotation.tsv"

# --- Directories ---
INPUT_DIR="unmapped/univec"
ALIGN_DIR="aligned/rrna_bowtie1"
UNMAP_DIR="unmapped/rrna_bowtie1"
COUNTS_DIR="counts/rrna_bowtie1"
LOG_DIR="log/rrna_bowtie1"
QC_DIR="qc/rrna_bowtie1"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" \
         "$LOG_DIR" "$QC_FASTQC" "$QC_MULTIQC"

inputs=( "$INPUT_DIR"/*_after_univec.fq.gz )
[[ ${#inputs[@]} -gt 0 ]] || {
  echo "❌ No FASTQ found in $INPUT_DIR (expected *_after_univec.fq.gz)"
  exit 1
}

echo "============================================================"
echo "🎯 Step 4: Aligning to human rRNA (Bowtie1)"
echo "============================================================"

################################################################################
# Function: process_sample
################################################################################
process_sample() {
  local in_file="$1"
  local sample="${in_file##*/}"
  sample="${sample%_after_univec.fq.gz}"

  local bam_out="${ALIGN_DIR}/${sample}_rrna_bowtie1.bam"
  local unmapped_fq="${UNMAP_DIR}/${sample}_after_rrna.fastq.gz"
  local log_file="${LOG_DIR}/${sample}_bowtie1.log"
  local fc_out="${COUNTS_DIR}/${sample}_rrna_bowtie1_counts.txt"

  echo "🧩 Processing sample: $sample" | tee -a "$log_file"

  # --- Bowtie1 alignment ---
  if [[ ! -s "$bam_out" ]]; then
    echo "🔧 Bowtie1 alignment..." | tee -a "$log_file"

    if ! ( ( [[ "$in_file" =~ \.gz$ ]] && zcat "$in_file" || cat "$in_file" ) \
      | bowtie -S -q -v 2 -a --best --strata -p "$THREADS" \
          --un >(gzip > "$unmapped_fq") "$REF_INDEX" - 2>> "$log_file" \
      | samtools view -@ "$THREADS" -h - \
      | samtools sort -@ "$THREADS" -n -O SAM - \
      | awk -v OFS="\t" '
          /^@/ {print; next}
          {
            if($1 != q && q != "") {
              for(i=1; i<=c; i++) print buf[i] "\tNH:i:" c;
              delete buf; c=0;
            }
            q=$1; c++; buf[c]=$0;
          }
          END {for(i=1; i<=c; i++) print buf[i] "\tNH:i:" c}
        ' \
      | samtools sort -@ "$THREADS" -O BAM -o "$bam_out" - ); then
      echo "❌ Bowtie1 failed" | tee -a "$log_file"
      return 1
    fi

    samtools index "$bam_out"
    echo "✅ BAM created: $bam_out" | tee -a "$log_file"
  else
    echo "⏩ BAM exists — skipping alignment." | tee -a "$log_file"
  fi

  # --- featureCounts ---
  if [[ ! -s "$fc_out" && -s "$bam_out" ]]; then
    echo "🧮 featureCounts (fractional)..." | tee -a "$log_file"

    if ! featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
          -M --fraction -O -s 0 -o "$fc_out" "$bam_out" &>> "$log_file"; then
      echo "⚠️ featureCounts failed" | tee -a "$log_file"
    fi
  fi

  # --- QC ---
  if $DO_QC; then
    local qc_zip="${QC_FASTQC}/${sample}_after_rrna_fastqc.zip"
    if [[ -s "$unmapped_fq" && ! -s "$qc_zip" ]]; then
      echo "🔍 FastQC..." | tee -a "$log_file"
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq"
    fi
  fi
}

################################################################################
# Run all samples
################################################################################
rc_all=0
for f in "${inputs[@]}"; do
  process_sample "$f" || rc_all=1
done

################################################################################
# Combine counts + annotation
################################################################################
combined_csv="${COUNTS_DIR}/rrna_counts_bowtie1.csv"
echo "📊 Combining counts → ${combined_csv}"

Rscript - <<EOF
suppressPackageStartupMessages({
  library(dplyr); library(readr)
})

counts_dir <- "${COUNTS_DIR}"
ann <- read_tsv("${ANNOT_FILE}", show_col_types=FALSE)
files <- list.files(counts_dir, pattern="_rrna_bowtie1_counts.txt$", full.names=TRUE)
if (length(files)==0) quit(status=0)

df_list <- lapply(files, function(f){
  sample <- sub("_rrna_bowtie1_counts.txt","",basename(f))
  dat <- read_tsv(f, comment="#", show_col_types=FALSE)
  count_col <- names(dat)[ncol(dat)]
  dat %>%
    rename(genename=Geneid) %>%
    select(genename, Length, !!sym(count_col)) %>%
    rename(fc_length = Length,
           !!sample := !!sym(count_col))
})

fc <- Reduce(function(x,y) full_join(x,y,by=c("genename","fc_length")), df_list)

merged <- fc %>%
  left_join(ann, by="genename") %>%
  mutate(length = coalesce(length, fc_length)) %>%
  select(genename, geneid, biotype_1, biotype_2, chromosome,
         start, end, strand, length, source,
         everything(), -fc_length)

names(merged) <- sub("_trimmed$", "", names(merged))
write_csv(merged, "${combined_csv}")
EOF

################################################################################
# MultiQC
################################################################################
if $DO_QC; then
  echo "📈 Running MultiQC..."
  multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force || true
fi

################################################################################
# Read-count summary (FIXED)
################################################################################
echo "📄 Summarizing reads processed..."
summary="${LOG_DIR}/read_counts.csv"
echo "Sample,Reads_Processed" > "$summary"

for f in "${LOG_DIR}"/*_bowtie1.log; do
  raw_sample=$(basename "$f" _bowtie1.log)
  sample="${raw_sample%_trimmed}"

  reads=$(
    awk '
      /# reads processed:/ {print $4; exit}
      /reads; of these:/  {print $1; exit}
    ' "$f"
  )

  [[ -n "$reads" ]] && echo "${sample},${reads}" >> "$summary"
done

echo "✅ Saved read summary: $summary"
echo "🎉 Human rRNA Bowtie1 alignment + counting complete"

exit $rc_all
