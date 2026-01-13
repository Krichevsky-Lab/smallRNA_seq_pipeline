#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

################################################################################
# 11_align_pirna_rescue.sh
# Step 11: Rescue unmapped piRNA reads using Bowtie2 (local mode)
#
# Usage:
#   bash scripts/11_align_pirna_rescue.sh [--qc]
#
# Description:
#   - Aligns reads unmapped from the piRNA Bowtie1 step using Bowtie2 (local)
#   - Optimized for piRNA lengths (24–33 nt)
#   - Retains only best AS-score alignments (best stratum)
#   - Adds NH tags for fractional assignment
#   - Runs featureCounts (-M --fraction -O)
#   - Optionally runs FastQC + MultiQC
#
# Inputs:
#   - unmapped/pirna_bowtie1/*_after_pirna.fastq.gz
#
# Outputs:
#   - aligned/pirna_bowtie2/*.bam
#   - unmapped/pirna_bowtie2/*_after_pirna_bowtie2.fastq.gz
#   - counts/pirna_bowtie2/*_counts.txt
#   - counts/pirna_bowtie2/pirna_counts_bowtie2.csv
#   - log/pirna_bowtie2/**
#   - qc/pirna_bowtie2/** (optional)
################################################################################

THREADS=10
DO_QC=false
[[ "${1:-}" == "--qc" ]] && DO_QC=true

# --- References (HUMAN) ---
REF_INDEX="reference/pirna/human_piRNA"
SAF_FILE="reference/pirna/human_piRNA.saf"
ANNOT_FILE="reference/pirna/human_piRNA_annotation.tsv"

# --- Directories (no tissue) ---
INPUT_DIR="unmapped/pirna_bowtie1"
ALIGN_DIR="aligned/pirna_bowtie2"
UNMAP_DIR="unmapped/pirna_bowtie2"
COUNTS_DIR="counts/pirna_bowtie2"
LOG_DIR="log/pirna_bowtie2"
ALIGN_LOG_DIR="${LOG_DIR}/align"
FC_LOG_DIR="${LOG_DIR}/featurecounts"
QC_DIR="qc/pirna_bowtie2"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" \
         "$ALIGN_LOG_DIR" "$FC_LOG_DIR" "$QC_FASTQC" "$QC_MULTIQC"

inputs=( "$INPUT_DIR"/*_after_pirna.fastq.gz )
[[ ${#inputs[@]} -gt 0 ]] || { echo "❌ No FASTQ files found in $INPUT_DIR"; exit 1; }

echo "============================================================"
echo "🎯 Step 11: Human piRNA rescue alignment (Bowtie2 local)"
echo "============================================================"

################################################################################
# Function: process_sample
################################################################################
process_sample() {
  local in_file="$1"
  local sample="${in_file##*/}"
  sample="${sample%_after_pirna.fastq.gz}"
  local sample_base="${sample#trimmed_}"

  echo "🧩 Processing sample: $sample_base"

  local bam_out="${ALIGN_DIR}/${sample_base}_pirna_bowtie2.bam"
  local unmapped_fq="${UNMAP_DIR}/${sample_base}_after_pirna_bowtie2.fastq.gz"
  local bowtie_log="${ALIGN_LOG_DIR}/${sample_base}.log"
  local fc_out="${COUNTS_DIR}/${sample_base}_pirna_bowtie2_counts.txt"
  local fc_log="${FC_LOG_DIR}/${sample_base}.log"

  ##############################################################################
  # Step 1. Bowtie2 alignment (local mode, relaxed for piRNA length)
  ##############################################################################
  if [[ ! -s "$bam_out" ]]; then
    echo "🔧 Aligning ${sample_base} (Bowtie2 local, relaxed)…"

    local tmp_sam="${ALIGN_LOG_DIR}/${sample_base}_tmp.sam"
    local tmp_best="${ALIGN_LOG_DIR}/${sample_base}_tmp_best.sam"

    if ! bowtie2 -p "$THREADS" \
      --local --very-sensitive-local -L 15 -N 1 --score-min L,0,0.8 \
      --no-unal --all \
      -x "$REF_INDEX" -U "$in_file" --un-gz "$unmapped_fq" \
      -S "$tmp_sam" 2> "$bowtie_log"; then
      echo "❌ Bowtie2 failed for ${sample_base}. See ${bowtie_log}"
      return 1
    fi

    # Best-score selection + NH tag assignment
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
        if (ascore > maxAS) maxAS = ascore;
      }
      END {
        for (i=1;i<=c;i++) if (as[i]==maxAS) best[i]=1;
        nBest=0; for(i=1;i<=c;i++) if (best[i]==1) nBest++;
        for (i=1;i<=c;i++) if (best[i]==1) print buf[i] "\tNH:i:" nBest;
      }
    ' "$tmp_sam" > "$tmp_best"

    samtools sort -@ "$THREADS" -O BAM -o "$bam_out" "$tmp_best"
    samtools index "$bam_out"
    rm -f "$tmp_sam" "$tmp_best"

    echo "✅ Created BAM: $bam_out"
  else
    echo "⏩ Skipping: BAM exists."
  fi

  ##############################################################################
  # Step 2. featureCounts
  ##############################################################################
  if [[ ! -s "$fc_out" && -s "$bam_out" ]]; then
    echo "🧮 Running featureCounts for ${sample_base}..."
    if ! featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
         -M --fraction -O -s 0 -o "$fc_out" "$bam_out" > "$fc_log" 2>&1; then
      echo "⚠️ featureCounts failed for ${sample_base}. See $fc_log"
    fi
  else
    echo "⏩ Skipping featureCounts."
  fi

  ##############################################################################
  # Step 3. QC
  ##############################################################################
  if $DO_QC; then
    local qc_zip="${QC_FASTQC}/${sample_base}_after_pirna_bowtie2_fastqc.zip"
    if [[ -s "$unmapped_fq" && ! -s "$qc_zip" ]]; then
      echo "🔍 Running FastQC for ${sample_base}..."
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq" || \
        echo "⚠️ FastQC failed for ${sample_base}"
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

combined_csv="${COUNTS_DIR}/pirna_counts_bowtie2.csv"
echo "📊 Combining per-sample counts → ${combined_csv}"

Rscript - <<EOF
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

counts_dir <- "${COUNTS_DIR}"
ann_file   <- "${ANNOT_FILE}"

files <- list.files(counts_dir, pattern="_pirna_bowtie2_counts.txt$", full.names=TRUE)
if (length(files)==0) quit(status=0)

ann <- read_tsv(ann_file, show_col_types=FALSE)

df_list <- lapply(files, function(f){
  sample <- sub("_pirna_bowtie2_counts.txt","",basename(f))
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

# ---- REMOVE _trimmed FROM SAMPLE NAMES ----
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
  echo "📈 Running MultiQC..."
  multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force || \
    echo "⚠️ MultiQC returned warnings."
fi

echo "✅ Human piRNA Bowtie2 rescue alignment complete"
exit $rc_all
