#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 06_align_trna.sh
# Align rRNA-rescue–unmapped reads to tRNA using Bowtie1.
#
# Usage:
#   bash alignment_scripts/06_align_trna.sh --species <mouse|human> [--qc]
#   SPECIES=mouse bash alignment_scripts/06_align_trna.sh [--qc]
#
# Inputs:
#   - alignment_output/unmapped/rrna_bowtie2/*_after_rrna_bowtie2.fastq.gz
#
# Outputs:
#   - alignment_output/aligned/trna_bowtie1/*.bam
#   - alignment_output/unmapped/trna_bowtie1/*_after_trna.fastq.gz
#   - alignment_output/counts/trna_bowtie1/*_counts.txt
#   - alignment_output/counts/trna_bowtie1/trna_counts_bowtie1.csv
#   - alignment_output/counts/trna_fragment/tRF_counts_bowtie1.csv
#   - alignment_output/counts/trna_fragment/per_sample_bowtie1/
#   - alignment_output/log/trna_bowtie1/*_bowtie1.log
#   - alignment_output/qc/trna_bowtie1/** (optional)
#
# Notes:
#   - Bowtie1: -a --best --strata
#   - NH tags added for multimappers
#   - featureCounts: -M --fraction -O
#   - tRNA fragment classification via classify_tRF.R
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

REF_INDEX="reference/$SPECIES/trna/${SPECIES}_tRNA"
SAF_FILE="reference/$SPECIES/trna/${SPECIES}_tRNA.saf"
ANNOT_FILE="reference/$SPECIES/trna/${SPECIES}_tRNA_annotation.tsv"

require_index "$REF_INDEX"
require_file  "$SAF_FILE"
require_file  "$ANNOT_FILE"

INPUT_DIR="alignment_output/unmapped/rrna_bowtie2"
ALIGN_DIR="alignment_output/aligned/trna_bowtie1"
UNMAP_DIR="alignment_output/unmapped/trna_bowtie1"
COUNTS_DIR="alignment_output/counts/trna_bowtie1"
LOG_DIR="alignment_output/log/trna_bowtie1"
QC_DIR="alignment_output/qc/trna_bowtie1"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" \
         "$LOG_DIR" "$QC_FASTQC" "$QC_MULTIQC"

step_banner 6 "tRNA alignment (Bowtie1)"
log_info "Species: $SPECIES"

inputs=( "$INPUT_DIR"/*_after_rrna_bowtie2.fastq.gz )
[[ ${#inputs[@]} -gt 0 ]] || { log_fail "No FASTQs found in $INPUT_DIR"; exit 1; }
n=${#inputs[@]}
i=0
rc_all=0
for fq in "${inputs[@]}"; do
  i=$((i+1))
  base="$(basename "$fq" _after_rrna_bowtie2.fastq.gz)"
  sample_base="${base%%_after_*}"

  bam_out="${ALIGN_DIR}/${sample_base}_trna_bowtie1.bam"
  unmapped_fq="${UNMAP_DIR}/${sample_base}_after_trna.fastq.gz"
  log_file="${LOG_DIR}/${sample_base}_bowtie1.log"
  fc_out="${COUNTS_DIR}/${sample_base}_trna_bowtie1_counts.txt"

  echo "  Processing: $sample_base ($i/$n)"

  if check_bam "$bam_out"; then
    log_skip "$sample_base — BAM exists"
  else
    log_info "$sample_base — Bowtie1 alignment..."

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
    log_done "$sample_base — BAM created"
  fi

  if ! check_tabular "$fc_out" && check_bam "$bam_out"; then
    log_info "$sample_base — featureCounts..."
    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
      -M --fraction -O -s 0 -o "$fc_out" "$bam_out" &>> "$log_file" \
      || log_warn "featureCounts failed for $sample_base"
  fi

  if $DO_QC; then
    qc_zip="${QC_FASTQC}/${sample_base}_after_trna_fastqc.zip"
    if check_gz "$unmapped_fq" && [[ ! -f "$qc_zip" ]]; then
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq" \
        || log_warn "FastQC failed for $sample_base"
    fi
  fi
done

if $DO_QC; then
  if compgen -G "${QC_FASTQC}/*fastqc.zip" > /dev/null; then
    log_info "Running MultiQC..."
    multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force || log_warn "MultiQC warnings, continuing..."
  fi
fi

combined_csv="${COUNTS_DIR}/trna_counts_bowtie1.csv"
if ! check_tabular "$combined_csv"; then
  log_info "Combining counts -> ${combined_csv}"
  Rscript "$SCRIPT_DIR/combine_counts.R" \
    "$COUNTS_DIR" "_trna_bowtie1_counts\\.txt$" "_trna_bowtie1_counts.txt" \
    "$ANNOT_FILE" "$combined_csv"
fi

TRF_DIR="alignment_output/counts/trna_fragment"
mkdir -p "$TRF_DIR"
trf_csv="${TRF_DIR}/tRF_counts_bowtie1.csv"
if ! check_tabular "$trf_csv"; then
  log_info "Running tRF classification (Bowtie1)..."
  Rscript "$SCRIPT_DIR/classify_tRF.R" \
    "$ALIGN_DIR" "$SAF_FILE" "$ANNOT_FILE" \
    "$trf_csv" "${TRF_DIR}/per_sample_bowtie1"
fi

log_done "Step 6 — tRNA Bowtie1 alignment + tRF classification complete"
exit $rc_all
