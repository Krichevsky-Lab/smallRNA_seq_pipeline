#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 5_align_rrna_rescue.sh
# Rescue rRNA-unmapped reads using Bowtie2 (local alignment).
#
# Usage:
#   bash scripts/5_align_rrna_rescue.sh [--qc]
#
# Inputs:
#   - unmapped/rrna_bowtie1/*.fastq.gz
#
# Outputs:
#   - aligned/rrna_bowtie2/*.bam
#   - unmapped/rrna_bowtie2/*.fastq.gz
#   - counts/rrna_bowtie2/*_counts.txt
#   - counts/rrna_bowtie2/rrna_counts_bowtie2.csv
#   - log/rrna_bowtie2/*_bowtie2.log
#   - qc/rrna_bowtie2/** (optional)
#
# Notes:
#   - Bowtie2 local, very-sensitive-local
#   - Retains only best AS:i score per read
#   - NH tags added for multimappers
#   - featureCounts: -M --fraction -O
# ------------------------------------------------------------------------------

THREADS=10
DO_QC=false
[[ "${1:-}" == "--qc" ]] && DO_QC=true

# -------------------------
# References (HUMAN)
# -------------------------
REF_INDEX="reference/rrna/human_rRNA"
SAF_FILE="reference/rrna/human_rRNA.saf"
ANNOT_FILE="reference/rrna/human_rRNA_annotation.tsv"

# -------------------------
# Directories (GLOBAL)
# -------------------------
INPUT_DIR="unmapped/rrna_bowtie1"
ALIGN_DIR="aligned/rrna_bowtie2"
UNMAP_DIR="unmapped/rrna_bowtie2"
COUNTS_DIR="counts/rrna_bowtie2"
LOG_DIR="log/rrna_bowtie2"
QC_DIR="qc/rrna_bowtie2"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" \
         "$LOG_DIR" "$QC_FASTQC" "$QC_MULTIQC"

inputs=( "$INPUT_DIR"/*.fastq.gz )
[[ ${#inputs[@]} -gt 0 ]] || {
  echo "❌ No FASTQ found in $INPUT_DIR"
  exit 1
}

echo "============================================================"
echo "🎯 Step 5: HUMAN rRNA rescue alignment (Bowtie2 local)"
echo "============================================================"

################################################################################
# Function: process_sample
################################################################################
process_sample() {
  local in_file="$1"
  local sample="${in_file##*/}"
  sample="${sample%.fastq.gz}"
  local sample_base="${sample%%_after_*}"

  local bam_out="${ALIGN_DIR}/${sample_base}_rrna_bowtie2.bam"
  local unmapped_fq="${UNMAP_DIR}/${sample_base}_after_rrna_bowtie2.fastq.gz"
  local log_file="${LOG_DIR}/${sample_base}_bowtie2.log"
  local fc_out="${COUNTS_DIR}/${sample_base}_rrna_bowtie2_counts.txt"

  echo "🧩 Processing: $sample_base" | tee -a "$log_file"

  # -------------------------
  # Bowtie2 alignment
  # -------------------------
  if [[ ! -s "$bam_out" ]]; then
    echo "🔧 Bowtie2 local alignment (best AS:i)..." | tee -a "$log_file"

    if ! bowtie2 -p "$THREADS" --local --very-sensitive-local --no-unal --all \
        -x "$REF_INDEX" -U "$in_file" --un-gz "$unmapped_fq" 2>> "$log_file" \
      | samtools view -h -@ "$THREADS" - \
      | samtools sort -n -@ "$THREADS" -O SAM - \
      | awk -v OFS="\t" '
          /^@/ { print; next }
          {
            if (match($0,/AS:i:([0-9-]+)/,m)) ascore=m[1]; else ascore=-999999;
            if ($1 != q && q != "") {
              for (i=1;i<=c;i++) if (as[i]==maxAS) best[i]=1;
              nBest=0; for(i=1;i<=c;i++) if(best[i]==1) nBest++;
              for (i=1;i<=c;i++) if (best[i]==1) print buf[i] "\tNH:i:" nBest;
              delete buf; delete as; delete best; c=0; maxAS=-999999;
            }
            q=$1; c++; buf[c]=$0; as[c]=ascore;
            if (ascore>maxAS) maxAS=ascore;
          }
          END {
            for (i=1;i<=c;i++) if (as[i]==maxAS) best[i]=1;
            nBest=0; for(i=1;i<=c;i++) if(best[i]==1) nBest++;
            for (i=1;i<=c;i++) if (best[i]==1) print buf[i] "\tNH:i:" nBest;
          }
        ' \
      | samtools sort -@ "$THREADS" -O BAM -o "$bam_out" -; then
      echo "❌ Bowtie2 failed for $sample_base" | tee -a "$log_file"
      return 1
    fi

    samtools index "$bam_out"
    echo "✅ BAM created: $bam_out" | tee -a "$log_file"
  else
    echo "⏩ BAM exists — skipping alignment." | tee -a "$log_file"
  fi

  # -------------------------
  # featureCounts
  # -------------------------
  if [[ ! -s "$fc_out" && -s "$bam_out" ]]; then
    echo "🧮 featureCounts (fractional)..." | tee -a "$log_file"

    if ! featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
        -M --fraction -O -s 0 -o "$fc_out" "$bam_out" &>> "$log_file"; then
      echo "⚠️ featureCounts failed" | tee -a "$log_file"
    fi
  fi

  # -------------------------
  # Optional QC
  # -------------------------
  if $DO_QC; then
    local qc_zip="${QC_FASTQC}/${sample_base}_after_rrna_bowtie2_fastqc.zip"
    if [[ -s "$unmapped_fq" && ! -s "$qc_zip" ]]; then
      echo "🔍 FastQC..." | tee -a "$log_file"
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq" || \
        echo "⚠️ FastQC failed" | tee -a "$log_file"
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
# Combine counts
################################################################################
combined_csv="${COUNTS_DIR}/rrna_counts_bowtie2.csv"
echo "📊 Combining counts → $combined_csv"

Rscript - <<EOF
suppressPackageStartupMessages({
  library(dplyr); library(readr)
})

counts_dir <- "${COUNTS_DIR}"
ann <- read_tsv("${ANNOT_FILE}", show_col_types=FALSE)
files <- list.files(counts_dir, pattern="_rrna_bowtie2_counts.txt$", full.names=TRUE)
if (length(files)==0) quit(status=0)

df_list <- lapply(files, function(f){
  sample <- sub("_rrna_bowtie2_counts.txt","",basename(f))
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

################################################################################
# MultiQC
################################################################################
if $DO_QC; then
  echo "📈 Running MultiQC..."
  multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force || true
fi

echo "🎉 HUMAN rRNA Bowtie2 rescue alignment complete"
exit $rc_all
