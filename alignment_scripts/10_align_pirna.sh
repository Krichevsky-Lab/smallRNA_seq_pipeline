#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 10_align_pirna.sh
# Align smallRNA-rescue–unmapped reads to piRNA reference using Bowtie1.
#
# Usage:
#   bash alignment_scripts/10_align_pirna.sh --species <mouse|human> [--qc]
#   SPECIES=mouse bash alignment_scripts/10_align_pirna.sh [--qc]
#
# Inputs:
#   - alignment_output/unmapped/smallrna_bowtie2/*_after_smallrna_bowtie2.fastq.gz
#
# Outputs:
#   - alignment_output/aligned/pirna_bowtie1/*.bam
#   - alignment_output/unmapped/pirna_bowtie1/*_after_pirna.fastq.gz
#   - alignment_output/counts/pirna_bowtie1/*_counts.txt
#   - alignment_output/counts/pirna_bowtie1/pirna_counts_bowtie1.csv
#   - alignment_output/log/pirna_bowtie1/**
#   - alignment_output/qc/pirna_bowtie1/** (optional)
#
# Notes:
#   - Bowtie1 best/strata; NH tags for multimappers
#   - featureCounts: -M --fraction -O
# ------------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common.sh"

SPECIES="${SPECIES:-}"
DO_QC=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --species) SPECIES="$2"; shift ;;
        --qc)      DO_QC=true ;;
        *)         log_fail "Unknown argument: $1"; exit 1 ;;
    esac
    shift
done
[[ -n "$SPECIES" ]] || { log_fail "--species required (mouse or human). Set via arg or \$SPECIES env."; exit 1; }
[[ "$SPECIES" == "mouse" || "$SPECIES" == "human" ]] || { log_fail "Unknown species '$SPECIES'"; exit 1; }

REF_INDEX="reference/$SPECIES/pirna/${SPECIES}_piRNA"
SAF_FILE="reference/$SPECIES/pirna/${SPECIES}_piRNA.saf"
ANNOT_FILE="reference/$SPECIES/pirna/${SPECIES}_piRNA_annotation.tsv"

require_index "$REF_INDEX"
require_file  "$SAF_FILE"
require_file  "$ANNOT_FILE"

INPUT_DIR="alignment_output/unmapped/smallrna_bowtie2"
ALIGN_DIR="alignment_output/aligned/pirna_bowtie1"
UNMAP_DIR="alignment_output/unmapped/pirna_bowtie1"
COUNTS_DIR="alignment_output/counts/pirna_bowtie1"
LOG_DIR="alignment_output/log/pirna_bowtie1"
QC_DIR="alignment_output/qc/pirna_bowtie1"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$LOG_DIR" "$QC_FASTQC" "$QC_MULTIQC"

step_banner 10 "piRNA alignment (Bowtie1)"
log_info "Species: $SPECIES"

inputs=( "$INPUT_DIR"/*_after_smallrna_bowtie2.fastq.gz )
[[ ${#inputs[@]} -gt 0 ]] || { log_fail "No FASTQ files in $INPUT_DIR"; exit 1; }
n=${#inputs[@]}
i=0
rc_all=0
for fq in "${inputs[@]}"; do
  i=$((i+1))
  sample=$(basename "$fq" _after_smallrna_bowtie2.fastq.gz)
  bam_out="${ALIGN_DIR}/${sample}_pirna.bam"
  unmapped_fq="${UNMAP_DIR}/${sample}_after_pirna.fastq.gz"
  log_file="${LOG_DIR}/${sample}_bowtie1.log"

  echo "  Processing: $sample ($i/$n)"

  if check_bam "$bam_out"; then
    log_skip "$sample — BAM exists"
  else
    log_info "$sample — Bowtie1 alignment..."

    bowtie -S -p "$THREADS" -a --best --strata -v 1 \
      -x "$REF_INDEX" -q "$fq" \
      --un >(gzip > "$unmapped_fq") 2> "$log_file" \
    | samtools view -@ "$THREADS" -h - \
    | samtools sort -@ "$THREADS" -n -O SAM - \
    | awk -v OFS="\t" '
        /^@/ { print; next }
        {
          if ($1 != q && q != "") {
            for (i=1;i<=c;i++) print buf[i] "\tNH:i:" c
            delete buf; c=0
          }
          q=$1; buf[++c]=$0
        }
        END { for (i=1;i<=c;i++) print buf[i] "\tNH:i:" c }
      ' \
    | samtools sort -@ "$THREADS" -O BAM -o "$bam_out" -

    if ! check_bam "$bam_out"; then
      log_warn "$sample — no alignments (empty BAM)"
      continue
    fi

    samtools index "$bam_out"
    log_done "$sample — BAM created"
  fi

  fc_out="${COUNTS_DIR}/${sample}_pirna_counts.txt"
  if ! check_tabular "$fc_out" && check_bam "$bam_out"; then
    log_info "$sample — featureCounts..."
    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
      -M --fraction -O -s 0 -o "$fc_out" "$bam_out" &>> "$log_file" \
      || log_warn "featureCounts failed for $sample"
  fi

  if $DO_QC; then
    qc_zip="${QC_FASTQC}/${sample}_after_pirna_fastqc.zip"
    if check_gz "$unmapped_fq" && [[ ! -f "$qc_zip" ]]; then
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq" \
        || log_warn "FastQC failed for $sample"
    fi
  fi
done

if $DO_QC; then
  if compgen -G "${QC_FASTQC}/*fastqc.zip" > /dev/null; then
    log_info "Running MultiQC..."
    multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force || log_warn "MultiQC encountered issues."
  fi
fi

combined_csv="${COUNTS_DIR}/pirna_counts_bowtie1.csv"
if ! check_tabular "$combined_csv"; then
  log_info "Combining counts -> ${combined_csv}"
  Rscript "$SCRIPT_DIR/combine_counts.R" \
    "$COUNTS_DIR" "_pirna_counts\\.txt$" "_pirna_counts.txt" \
    "$ANNOT_FILE" "$combined_csv"
fi

log_done "Step 10 — piRNA Bowtie1 alignment complete"
exit $rc_all
