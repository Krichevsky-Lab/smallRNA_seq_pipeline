#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

################################################################################
# Step 9: Rescue smallRNA-unmapped reads using Bowtie2 (local mode)
#
# Usage:
#   bash scripts/9_align_smallrna_rescue.sh [--qc]
#
# Description:
#   - Re-aligns reads unmapped from smallRNA Bowtie1 step using Bowtie2 (local)
#   - Optimized for short RNAs (miRNA-sized, ~18–25 nt)
#   - Retains best AS:i alignments only (best stratum + ties)
#   - Adds NH tags for fractional counting
#   - Runs featureCounts (-M --fraction -O)
#   - Optionally runs FastQC + MultiQC
#
# Inputs:
#   - unmapped/smallrna_bowtie1/*_after_smallrna.fastq.gz
#
# Outputs:
#   - aligned/smallrna_bowtie2/*.bam
#   - unmapped/smallrna_bowtie2/*_after_smallrna_bowtie2.fastq.gz
#   - counts/smallrna_bowtie2/*_counts.txt
#   - counts/smallrna_bowtie2/smallrna_counts_bowtie2.csv
#   - log/smallrna_bowtie2/**
#   - qc/smallrna_bowtie2/** (optional)
#
################################################################################

THREADS=10
DO_QC=false

if [[ "${1:-}" == "--qc" ]]; then
  DO_QC=true
fi

# --- HUMAN References ---
REF_INDEX="reference/smallrna/human_smallRNA"
SAF_FILE="reference/smallrna/human_smallRNA.saf"
ANNOT_FILE="reference/smallrna/human_smallRNA_annotation.tsv"

# --- Directories (NO TISSUE) ---
INPUT_DIR="unmapped/smallrna_bowtie1"
ALIGN_DIR="aligned/smallrna_bowtie2"
UNMAP_DIR="unmapped/smallrna_bowtie2"
COUNTS_DIR="counts/smallrna_bowtie2"
LOG_DIR="log/smallrna_bowtie2"
QC_DIR="qc/smallrna_bowtie2"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$LOG_DIR" \
         "$QC_FASTQC" "$QC_MULTIQC"

inputs=( "$INPUT_DIR"/*_after_smallrna.fastq.gz )
[[ ${#inputs[@]} -gt 0 ]] || { echo "❌ No FASTQ found in $INPUT_DIR"; exit 1; }

echo "============================================================"
echo "🎯 Step 9: Human smallRNA Bowtie2 rescue alignment + counting"
echo "============================================================"

################################################################################
# Function: process_sample
################################################################################
process_sample() {
  local in_file="$1"
  local sample="${in_file##*/}"
  sample="${sample%_after_smallrna.fastq.gz}"
  local sample_base="${sample#trimmed_}"

  local bam_out="${ALIGN_DIR}/${sample_base}_smallrna_bowtie2.bam"
  local unmapped_fq="${UNMAP_DIR}/${sample_base}_after_smallrna_bowtie2.fastq.gz"
  local log_file="${LOG_DIR}/${sample_base}.log"
  local fc_out="${COUNTS_DIR}/${sample_base}_smallrna_bowtie2_counts.txt"

  echo "-------------------------------------------"
  echo "🧩 Processing sample: $sample_base"

  if [[ -s "$bam_out" && -s "$fc_out" ]]; then
    echo "⏩ Skipping ${sample_base} — outputs exist"
    return 0
  fi

  local tmp_sam="${LOG_DIR}/${sample_base}_tmp.sam"
  local tmp_best="${LOG_DIR}/${sample_base}_tmp_best.sam"

  {
    echo "🔧 Bowtie2 local alignment"
    bowtie2 -p "$THREADS" \
      --local --very-sensitive-local -L 15 -N 1 --score-min L,0,0.8 \
      --no-unal --all \
      -x "$REF_INDEX" -U "$in_file" --un-gz "$unmapped_fq" \
      -S "$tmp_sam"

    echo "🧮 Selecting best AS:i alignments + adding NH tags"
    awk -v OFS="\t" '
      /^@/ { print; next }
      {
        if (match($0,/AS:i:([0-9-]+)/,m)) ascore = m[1]; else ascore = -999999;
        if ($1 != q && q != "") {
          for (i=1;i<=c;i++) if (as[i] == maxAS) best[i]=1;
          nBest=0; for(i=1;i<=c;i++) if(best[i]) nBest++;
          for(i=1;i<=c;i++) if(best[i]) print buf[i] "\tNH:i:" nBest;
          delete buf; delete as; delete best; c=0; maxAS=-999999;
        }
        q=$1; c++; buf[c]=$0; as[c]=ascore;
        if (ascore > maxAS) maxAS = ascore;
      }
      END {
        for (i=1;i<=c;i++) if (as[i]==maxAS) best[i]=1;
        nBest=0; for(i=1;i<=c;i++) if(best[i]) nBest++;
        for (i=1;i<=c;i++) if(best[i]) print buf[i] "\tNH:i:" nBest;
      }
    ' "$tmp_sam" > "$tmp_best"

    samtools sort -@ "$THREADS" -O BAM -o "$bam_out" "$tmp_best"
    samtools index "$bam_out"
    rm -f "$tmp_sam" "$tmp_best"

    echo "🧮 featureCounts (fractional)"
    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
      -M --fraction -O -s 0 \
      -o "$fc_out" "$bam_out" 2>&1

    echo "✅ Done: ${sample_base}"
  } > "$log_file" 2>&1

  if $DO_QC; then
    local qc_zip="${QC_FASTQC}/${sample_base}_after_smallrna_bowtie2_fastqc.zip"
    if [[ -s "$unmapped_fq" && ! -s "$qc_zip" ]]; then
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq" || true
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
# Combine counts and annotate
################################################################################

combined_csv="${COUNTS_DIR}/smallrna_counts_bowtie2.csv"
echo "📊 Combining per-sample counts → ${combined_csv}"

Rscript - <<EOF
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

counts_dir <- "${COUNTS_DIR}"
ann_file   <- "${ANNOT_FILE}"

files <- list.files(counts_dir, pattern="_smallrna_bowtie2_counts.txt$", full.names=TRUE)
if (length(files)==0) quit(status=0)

ann <- read_tsv(ann_file, show_col_types=FALSE)

df_list <- lapply(files, function(f){
  sample <- sub("_smallrna_bowtie2_counts.txt","",basename(f))
  dat <- read_tsv(f, comment="#", show_col_types=FALSE)
  count_col <- names(dat)[ncol(dat)]
  dat %>%
    rename(genename = Geneid) %>%
    select(genename, Length, !!sym(count_col)) %>%
    rename(fc_length = Length,
           !!sample := !!sym(count_col))
})

fc <- Reduce(function(a,b) full_join(a,b, by=c("genename","fc_length")), df_list)

merged <- fc %>%
  left_join(ann, by="genename") %>%
  mutate(length = coalesce(length, fc_length)) %>%
  select(genename, geneid, biotype_1, biotype_2, chromosome,
         start, end, strand, length, source,
         everything(), -fc_length)

clean_names <- names(merged)
clean_names <- sub("_trimmed$", "", clean_names)
names(merged) <- clean_names

write_csv(merged, "${combined_csv}")
message("✅ Saved: ${combined_csv}")
EOF

################################################################################
# MultiQC
################################################################################
if $DO_QC; then
  echo "📈 Running MultiQC…"
  multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force || true
fi

################################################################################
# Done
################################################################################

echo "✅ Human smallRNA Bowtie2 rescue alignment complete"
exit $rc_all
