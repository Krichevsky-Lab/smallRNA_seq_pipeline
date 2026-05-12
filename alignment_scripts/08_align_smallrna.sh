#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 08_align_smallrna.sh
# Align tRNA-rescue–unmapped reads to small RNA references using Bowtie1.
#
# Usage:
#   bash alignment_scripts/08_align_smallrna.sh --species <mouse|human> [--qc]
#   SPECIES=mouse bash alignment_scripts/08_align_smallrna.sh [--qc]
#
# Inputs:
#   - alignment_output/unmapped/trna_bowtie2/*_after_trna_bowtie2.fastq.gz
#
# Outputs:
#   - alignment_output/aligned/smallrna_bowtie1/*.bam
#   - alignment_output/unmapped/smallrna_bowtie1/*_after_smallrna.fastq.gz
#   - alignment_output/counts/smallrna_bowtie1/*_counts.txt
#   - alignment_output/counts/smallrna_bowtie1/smallrna_counts_bowtie1.csv
#   - alignment_output/log/smallrna_bowtie1/*_bowtie1.log
#   - alignment_output/qc/smallrna_bowtie1/** (optional)
#
# Notes:
#   - Bowtie1 best/strata alignment
#   - NH tags added for multimappers (fractional counting)
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

REF_INDEX="reference/$SPECIES/smallrna/${SPECIES}_smallRNA"
SAF_FILE="reference/$SPECIES/smallrna/${SPECIES}_smallRNA.saf"
ANNOT_FILE="reference/$SPECIES/smallrna/${SPECIES}_smallRNA_annotation.tsv"

require_index "$REF_INDEX"
require_file  "$SAF_FILE"
require_file  "$ANNOT_FILE"

INPUT_DIR="alignment_output/unmapped/trna_bowtie2"
ALIGN_DIR="alignment_output/aligned/smallrna_bowtie1"
UNMAP_DIR="alignment_output/unmapped/smallrna_bowtie1"
COUNTS_DIR="alignment_output/counts/smallrna_bowtie1"
LOG_DIR="alignment_output/log/smallrna_bowtie1"
QC_DIR="alignment_output/qc/smallrna_bowtie1"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$LOG_DIR" "$QC_FASTQC" "$QC_MULTIQC"

step_banner 8 "smallRNA alignment (Bowtie1)"
log_info "Species: $SPECIES"

inputs=( "$INPUT_DIR"/*_after_trna_bowtie2.fastq.gz )
[[ ${#inputs[@]} -gt 0 ]] || { log_fail "No FASTQ files found in $INPUT_DIR"; exit 1; }
n=${#inputs[@]}
i=0
rc_all=0
for fq in "${inputs[@]}"; do
  i=$((i+1))
  sample=$(basename "$fq" _after_trna_bowtie2.fastq.gz)
  bam_out="${ALIGN_DIR}/${sample}_smallrna.bam"
  unmapped_fq="${UNMAP_DIR}/${sample}_after_smallrna.fastq.gz"
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
      log_warn "$sample — empty BAM, skipping"
      continue
    fi

    samtools index "$bam_out"
    log_done "$sample — BAM created"
  fi

  fc_out="${COUNTS_DIR}/${sample}_smallrna_counts.txt"
  if ! check_tabular "$fc_out" && check_bam "$bam_out"; then
    log_info "$sample — featureCounts..."
    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
      -M --fraction -O -s 0 -o "$fc_out" "$bam_out" &>> "$log_file"
  fi

  if $DO_QC; then
    qc_zip="${QC_FASTQC}/${sample}_after_smallrna_fastqc.zip"
    if check_gz "$unmapped_fq" && [[ ! -f "$qc_zip" ]]; then
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq"
    fi
  fi
done

if $DO_QC; then
  if compgen -G "${QC_FASTQC}/*fastqc.zip" > /dev/null; then
    log_info "Running MultiQC..."
    multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force || log_warn "MultiQC warnings"
  fi
fi

combined_csv="${COUNTS_DIR}/smallrna_counts_bowtie1.csv"
if ! check_tabular "$combined_csv"; then
  log_info "Combining counts -> ${combined_csv}"
  Rscript "$SCRIPT_DIR/combine_counts.R" \
    "$COUNTS_DIR" "_smallrna_counts\\.txt$" "_smallrna_counts.txt" \
    "$ANNOT_FILE" "$combined_csv"
fi

log_done "Step 8 — smallRNA Bowtie1 alignment complete"
exit $rc_all
