#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

################################################################################
# 10_align_pirna.sh
# Step 10: Align smallRNA-rescue–unmapped reads to human piRNA (Bowtie1)
#
# Usage:
#   bash scripts/10_align_pirna.sh [--qc]
#
# Description:
#   - Aligns reads unmapped after smallRNA Bowtie2 rescue to piRNA reference
#   - Reports all best alignments (--best --strata -a)
#   - Adds NH tags for multimappers (fractional assignment)
#   - Runs featureCounts (-M --fraction -O)
#   - Optionally runs FastQC + MultiQC
#
# Inputs:
#   - unmapped/smallrna_bowtie2/*_after_smallrna_bowtie2.fastq.gz
#
# Outputs:
#   - aligned/pirna_bowtie1/*.bam
#   - unmapped/pirna_bowtie1/*_after_pirna.fastq.gz
#   - counts/pirna_bowtie1/*_counts.txt
#   - counts/pirna_bowtie1/pirna_counts_bowtie1.csv
#   - log/pirna_bowtie1/**
#   - qc/pirna_bowtie1/** (optional)
################################################################################

THREADS=10
DO_QC=false
[[ "${1:-}" == "--qc" ]] && DO_QC=true

# --- References (HUMAN) ---
REF_INDEX="reference/pirna/human_piRNA"
SAF_FILE="reference/pirna/human_piRNA.saf"
ANNOT_FILE="reference/pirna/human_piRNA_annotation.tsv"

# --- Directories (NO TISSUE) ---
INPUT_DIR="unmapped/smallrna_bowtie2"
ALIGN_DIR="aligned/pirna_bowtie1"
UNMAP_DIR="unmapped/pirna_bowtie1"
COUNTS_DIR="counts/pirna_bowtie1"
LOG_DIR="log/pirna_bowtie1"
QC_DIR="qc/pirna_bowtie1"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$LOG_DIR" "$QC_FASTQC" "$QC_MULTIQC"

echo "============================================================"
echo "🎯 Step 10: Aligning to **human piRNA** (Bowtie1)"
echo "============================================================"

for fq in "$INPUT_DIR"/*_after_smallrna_bowtie2.fastq.gz; do
  [[ -e "$fq" ]] || { echo "❌ No FASTQ files in $INPUT_DIR"; exit 1; }

  sample=$(basename "$fq" _after_smallrna_bowtie2.fastq.gz)
  sample=${sample#trimmed_}
  bam_out="${ALIGN_DIR}/${sample}_pirna.bam"
  unmapped_fq="${UNMAP_DIR}/${sample}_after_pirna.fastq.gz"
  log_file="${LOG_DIR}/${sample}_bowtie1.log"

  # --- Skip if BAM exists ---
  if [[ -s "$bam_out" ]]; then
    echo "⏩ $sample BAM exists → skipping."
  else
    echo "🔧 Aligning $sample..."

    bowtie -S -p "$THREADS" -a --best --strata -v 1 \
      -x "$REF_INDEX" -q "$fq" \
      --un >(gzip > "$unmapped_fq") 2> "$log_file" \
      | samtools view -h - \
      | awk 'BEGIN{OFS="\t"} /^@/ {print; next} {
               if($0=="" || $1=="*" || NF<11) next;
               print
             }' \
      | samtools view -bS - \
      | samtools sort -@ "$THREADS" -o "$bam_out" - 2>> "$log_file"

    # --- Check output ---
    if [[ ! -s "$bam_out" || $(stat -c%s "$bam_out") -lt 1000 ]]; then
      echo "⚠️ No alignments for $sample (empty BAM)" | tee -a "$log_file"
      rm -f "$bam_out"
      continue
    fi

    echo "✅ Created BAM: $bam_out"

    # --- Add NH tags for multimappers ---
    echo "🔖 Adding NH tags..."
    samtools view -h "$bam_out" | \
      awk 'BEGIN{OFS="\t"}
           /^@/ {print; next}
           {
             q=$1; count[q]++; store[q]=store[q]"\n"$0;
           }
           END{
             for(q in store){
               n=count[q];
               split(store[q], lines, "\n");
               for(i=2;i<=length(lines);i++){
                 line=lines[i];
                 if(line=="") continue;
                 gsub(/\tNH:i:[0-9]+/, "", line);
                 print line "\tNH:i:" n;
               }
             }
           }' \
      | samtools view -bS -o "${bam_out%.bam}_NH.bam" -

    mv "${bam_out%.bam}_NH.bam" "$bam_out"
    echo "✅ NH tags added to $sample"
  fi


  # --- featureCounts ---
fc_out="${COUNTS_DIR}/${sample}_pirna_counts.txt"
if [[ ! -s "$fc_out" && -s "$bam_out" ]]; then
  echo "🧮 Counting reads for $sample..." >> "$log_file"
  featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
    -M --fraction -O -s 0 -o "$fc_out" "$bam_out" &>> "$log_file" || \
    echo "⚠️ featureCounts failed for $sample" >> "$log_file"
fi


  # --- QC ---
  if $DO_QC; then
    qc_zip="${QC_FASTQC}/${sample}_after_pirna_fastqc.zip"
    if [[ -s "$unmapped_fq" ]]; then
      [[ -s "$qc_zip" ]] && echo "⏩ FastQC exists → skipping." || \
        fastqc -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq" || \
        echo "⚠️ FastQC failed for $sample"
    fi
  fi

done

# --- MultiQC ---
if $DO_QC; then
  if compgen -G "${QC_FASTQC}/*fastqc.zip" > /dev/null; then
    echo "📊 Running MultiQC..."
    multiqc -o "$QC_MULTIQC" "$QC_FASTQC" || \
      echo "⚠️ MultiQC encountered issues."
  else
    echo "⚠️ No FastQC files → skipping MultiQC."
  fi
fi

# --- Combine, annotate ---
combined_csv="${COUNTS_DIR}/pirna_counts_bowtie1.csv"
if [[ ! -s "$combined_csv" ]]; then
  echo "🧩 Combining per-sample piRNA counts…"
  Rscript - <<EOF
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

files <- list.files("${COUNTS_DIR}", pattern="_pirna_counts.txt$", full.names=TRUE)
if (length(files)==0) quit("no", status=0)

ann <- read_tsv("${ANNOT_FILE}", show_col_types=FALSE)

df_list <- lapply(files, function(f){
  sample <- sub("_pirna_counts.txt","",basename(f))
  sample <- str_remove(sample, "^trimmed_")
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
message("✅ Saved: ${combined_csv}")
EOF
fi

echo "✅ Human piRNA Bowtie1 alignment complete"
