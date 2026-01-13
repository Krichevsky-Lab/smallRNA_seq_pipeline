#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

################################################################################
# 13_align_longrna_rescue.sh
# Step 13: Rescue unmapped longRNA reads using Bowtie2 (local mode)
#
# Usage:
#   bash scripts/13_align_longrna_rescue.sh [--qc]
#
# Description:
#   - Re-aligns reads unmapped after longRNA Bowtie1 using Bowtie2 local mode
#   - Retains only best AS:i score alignments (best stratum + ties)
#   - Adds NH tags for fractional counting
#   - Runs featureCounts (-M --fraction -O)
#   - Optional FastQC + MultiQC on longRNA-unmapped reads
#   - Uses a single per-sample log combining alignment + featureCounts
#
# Inputs:
#   unmapped/longrna_bowtie1/*_after_longrna.fastq.gz
#
# Outputs:
#   aligned/longrna_bowtie2/*.bam
#   unmapped/longrna_bowtie2/*_after_longrna_bowtie2.fastq.gz
#   counts/longrna_bowtie2/*_longrna_bowtie2_counts.txt
#   counts/longrna_bowtie2/longrna_counts_bowtie2.csv
#   log/longrna_bowtie2/*.log
#   qc/longrna_bowtie2/{fastqc,multiqc}
################################################################################

THREADS=10
DO_QC=false
rc_all=0
[[ "${1:-}" == "--qc" ]] && DO_QC=true

# --- References (HUMAN) ---
REF_INDEX="reference/longrna/human_longRNA"
SAF_FILE="reference/longrna/human_longRNA.saf"
ANNOT_FILE="reference/longrna/human_longRNA_annotation.tsv"

# --- Directories ---
INPUT_DIR="unmapped/longrna_bowtie1"
ALIGN_DIR="aligned/longrna_bowtie2"
UNMAP_DIR="unmapped/longrna_bowtie2"
COUNTS_DIR="counts/longrna_bowtie2"
LOG_DIR="log/longrna_bowtie2"
QC_DIR="qc/longrna_bowtie2"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$LOG_DIR"
$DO_QC && mkdir -p "$QC_FASTQC" "$QC_MULTIQC"

inputs=( "$INPUT_DIR"/*_after_longrna.fastq.gz )
[[ ${#inputs[@]} -gt 0 ]] || { echo "❌ No FASTQ files found in $INPUT_DIR"; exit 1; }

echo "============================================================"
echo "🎯 Step 13: Human longRNA rescue alignment (Bowtie2 local)"
echo "QC enabled: $DO_QC"
echo "============================================================"

################################################################################
# Main loop
################################################################################
for fq in "${inputs[@]}"; do
  sample=$(basename "$fq" "_after_longrna.fastq.gz")
  sample="${sample#trimmed_}"

  bam_out="${ALIGN_DIR}/${sample}_longrna_bowtie2.bam"
  unmapped_fq="${UNMAP_DIR}/${sample}_after_longrna_bowtie2.fastq.gz"
  log_file="${LOG_DIR}/${sample}.log"

  echo "--------------------------------------------"
  echo "🧩 Processing sample: $sample"
  echo "--------------------------------------------"

  # ------------------------------------------------------------------
  # Bowtie2 local rescue alignment
  # ------------------------------------------------------------------
  if [[ ! -s "$bam_out" ]]; then
    echo "🔧 Bowtie2 local alignment…" | tee -a "$log_file"

    tmp_sam="${LOG_DIR}/${sample}_tmp.sam"
    tmp_best="${LOG_DIR}/${sample}_tmp_best.sam"

    if ! bowtie2 -p "$THREADS" \
        --local --very-sensitive-local -L 20 -N 1 --score-min L,0,0.9 \
        --no-unal --all \
        -x "$REF_INDEX" -U "$fq" \
        --un-gz "$unmapped_fq" \
        -S "$tmp_sam" 2>> "$log_file"; then
      echo "❌ Bowtie2 failed for $sample" | tee -a "$log_file"
      continue
    fi

    # Select best AS:i alignments + NH tags
    awk -v OFS="\t" '
      /^@/ { print; next }
      {
        if (match($0,/AS:i:([0-9-]+)/,m)) ascore=m[1]; else ascore=-999999;
        if ($1 != q && q != "") {
          for (i=1;i<=c;i++) if (as[i]==maxAS) best[i]=1;
          nBest=0; for(i=1;i<=c;i++) if(best[i]==1) nBest++;
          for (i=1;i<=c;i++) if(best[i]==1) print buf[i] "\tNH:i:" nBest;
          delete buf; delete as; delete best; c=0; maxAS=-999999;
        }
        q=$1; c++; buf[c]=$0; as[c]=ascore;
        if (ascore>maxAS) maxAS=ascore;
      }
      END {
        for (i=1;i<=c;i++) if (as[i]==maxAS) best[i]=1;
        nBest=0; for(i=1;i<=c;i++) if(best[i]==1) nBest++;
        for (i=1;i<=c;i++) if(best[i]==1) print buf[i] "\tNH:i:" nBest;
      }
    ' "$tmp_sam" > "$tmp_best"

    samtools sort -@ "$THREADS" -O BAM -o "$bam_out" "$tmp_best"
    samtools index "$bam_out"
    rm -f "$tmp_sam" "$tmp_best"

    echo "✅ BAM created: $bam_out" | tee -a "$log_file"
  else
    echo "⏩ BAM exists → skipping alignment" | tee -a "$log_file"
  fi

  # ------------------------------------------------------------------
  # featureCounts (fractional)
  # ------------------------------------------------------------------
  fc_out="${COUNTS_DIR}/${sample}_longrna_bowtie2_counts.txt"
  if [[ ! -s "$fc_out" && -s "$bam_out" ]]; then
    echo "🧮 featureCounts…" | tee -a "$log_file"
    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
      -M --fraction -O -s 0 \
      -o "$fc_out" "$bam_out" &>> "$log_file" \
      || echo "⚠️ featureCounts failed" | tee -a "$log_file"
  fi

  # ------------------------------------------------------------------
  # QC
  # ------------------------------------------------------------------
  if $DO_QC; then
    qc_zip="${QC_FASTQC}/${sample}_after_longrna_bowtie2_fastqc.zip"
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
combined_csv="${COUNTS_DIR}/longrna_counts_bowtie2.csv"
echo "📊 Combining per-sample counts → $combined_csv"

Rscript - <<EOF
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

files <- list.files("${COUNTS_DIR}", pattern="_longrna_bowtie2_counts.txt$", full.names=TRUE)
if (length(files)==0) quit(status=0)

ann <- read_tsv("${ANNOT_FILE}", show_col_types=FALSE)

df_list <- lapply(files, function(f){
  sample <- sub("_longrna_bowtie2_counts.txt","",basename(f))
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

echo "🎉 Human longRNA Bowtie2 rescue alignment complete"
exit $rc_all
