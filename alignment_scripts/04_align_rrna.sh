#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 04_align_rrna.sh
# Align UniVec-unmapped reads to rRNA using Bowtie1.
#
# Usage:
#   bash alignment_scripts/04_align_rrna.sh --species <mouse|human> [--qc]
#   SPECIES=mouse bash alignment_scripts/04_align_rrna.sh [--qc]
#
# Inputs:
#   - alignment_output/unmapped/univec/*_after_univec.fq.gz
#
# Outputs:
#   - alignment_output/aligned/rrna_bowtie1/*.bam
#   - alignment_output/unmapped/rrna_bowtie1/*.fastq.gz
#   - alignment_output/counts/rrna_bowtie1/*_counts.txt
#   - alignment_output/counts/rrna_bowtie1/rrna_counts_bowtie1.csv
#   - alignment_output/log/rrna_bowtie1/*_bowtie1.log
#   - alignment_output/qc/rrna_bowtie1/** (optional)
#
# Notes:
#   - Bowtie1: -v 1 -a --best --strata
#   - NH tags added for multimappers
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

REF_INDEX="reference/$SPECIES/rrna/${SPECIES}_rRNA"
SAF_FILE="reference/$SPECIES/rrna/${SPECIES}_rRNA.saf"
ANNOT_FILE="reference/$SPECIES/rrna/${SPECIES}_rRNA_annotation.tsv"

require_index "$REF_INDEX"
require_file  "$SAF_FILE"
require_file  "$ANNOT_FILE"

INPUT_DIR="alignment_output/unmapped/univec"
ALIGN_DIR="alignment_output/aligned/rrna_bowtie1"
UNMAP_DIR="alignment_output/unmapped/rrna_bowtie1"
COUNTS_DIR="alignment_output/counts/rrna_bowtie1"
LOG_DIR="alignment_output/log/rrna_bowtie1"
QC_DIR="alignment_output/qc/rrna_bowtie1"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" \
         "$LOG_DIR" "$QC_FASTQC" "$QC_MULTIQC"

inputs=( "$INPUT_DIR"/*_after_univec.fq.gz )
[[ ${#inputs[@]} -gt 0 ]] || {
  log_fail "No FASTQ found in $INPUT_DIR (expected *_after_univec.fq.gz)"
  exit 1
}

step_banner 4 "rRNA alignment (Bowtie1)"
log_info "Species: $SPECIES | Samples: ${#inputs[@]}"

process_sample() {
  local in_file="$1"
  local idx="$2"
  local total="$3"
  local sample="${in_file##*/}"
  sample="${sample%_after_univec.fq.gz}"

  local bam_out="${ALIGN_DIR}/${sample}_rrna_bowtie1.bam"
  local unmapped_fq="${UNMAP_DIR}/${sample}_after_rrna.fastq.gz"
  local log_file="${LOG_DIR}/${sample}_bowtie1.log"
  local fc_out="${COUNTS_DIR}/${sample}_rrna_bowtie1_counts.txt"

  echo "  Processing: $sample ($idx/$total)"

  if ! check_bam "$bam_out"; then
    log_info "$sample — Bowtie1 alignment..."

    if ! ( ( [[ "$in_file" =~ \.gz$ ]] && zcat "$in_file" || cat "$in_file" ) \
      | bowtie -S -q -v 1 -a --best --strata -p "$THREADS" \
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
      log_fail "Bowtie1 failed for $sample"
      return 1
    fi

    samtools index "$bam_out"
    log_done "$sample — BAM created"
  else
    log_skip "$sample — BAM exists"
  fi

  if ! check_tabular "$fc_out" && check_bam "$bam_out"; then
    log_info "$sample — featureCounts..."
    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
          -M --fraction -O -s 0 -o "$fc_out" "$bam_out" &>> "$log_file" \
      || log_warn "featureCounts failed for $sample"
  fi

  if $DO_QC; then
    local qc_zip="${QC_FASTQC}/${sample}_after_rrna_fastqc.zip"
    if check_gz "$unmapped_fq" && [[ ! -f "$qc_zip" ]]; then
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq"
    fi
  fi
}

n=${#inputs[@]}
i=0
rc_all=0
for f in "${inputs[@]}"; do
  i=$((i+1))
  process_sample "$f" "$i" "$n" || rc_all=1
done

combined_csv="${COUNTS_DIR}/rrna_counts_bowtie1.csv"
log_info "Combining counts -> ${combined_csv}"
Rscript "$SCRIPT_DIR/combine_counts.R" \
  "$COUNTS_DIR" "_rrna_bowtie1_counts\\.txt$" "_rrna_bowtie1_counts.txt" \
  "$ANNOT_FILE" "$combined_csv"

if $DO_QC; then
  log_info "Running MultiQC..."
  multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force || true
fi

summary="${LOG_DIR}/read_counts.csv"
echo "Sample,Reads_Processed" > "$summary"
for f in "${LOG_DIR}"/*_bowtie1.log; do
  raw_sample=$(basename "$f" _bowtie1.log)
  sample="${raw_sample%_trimmed}"
  reads=$(awk '/# reads processed:/ {print $4; exit} /reads; of these:/ {print $1; exit}' "$f")
  [[ -n "$reads" ]] && echo "${sample},${reads}" >> "$summary"
done

log_done "Step 4 — rRNA Bowtie1 alignment complete"
exit $rc_all
