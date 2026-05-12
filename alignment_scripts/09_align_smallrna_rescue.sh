#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 09_align_smallrna_rescue.sh
# Rescue smallRNA-unmapped reads using Bowtie2 (local mode).
#
# Usage:
#   bash alignment_scripts/09_align_smallrna_rescue.sh --species <mouse|human> [--qc]
#   SPECIES=mouse bash alignment_scripts/09_align_smallrna_rescue.sh [--qc]
#
# Inputs:
#   - alignment_output/unmapped/smallrna_bowtie1/*_after_smallrna.fastq.gz
#
# Outputs:
#   - alignment_output/aligned/smallrna_bowtie2/*.bam
#   - alignment_output/unmapped/smallrna_bowtie2/*_after_smallrna_bowtie2.fastq.gz
#   - alignment_output/counts/smallrna_bowtie2/*_counts.txt
#   - alignment_output/counts/smallrna_bowtie2/smallrna_counts_bowtie2.csv
#   - alignment_output/log/smallrna_bowtie2/**
#   - alignment_output/qc/smallrna_bowtie2/** (optional)
#
# Notes:
#   - Bowtie2 local, very-sensitive-local; optimized for miRNA length (~18-25 nt)
#   - Retains best AS:i alignments only
#   - NH tags added for fractional counting
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

INPUT_DIR="alignment_output/unmapped/smallrna_bowtie1"
ALIGN_DIR="alignment_output/aligned/smallrna_bowtie2"
UNMAP_DIR="alignment_output/unmapped/smallrna_bowtie2"
COUNTS_DIR="alignment_output/counts/smallrna_bowtie2"
LOG_DIR="alignment_output/log/smallrna_bowtie2"
QC_DIR="alignment_output/qc/smallrna_bowtie2"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$LOG_DIR" \
         "$QC_FASTQC" "$QC_MULTIQC"

inputs=( "$INPUT_DIR"/*_after_smallrna.fastq.gz )
[[ ${#inputs[@]} -gt 0 ]] || { log_fail "No FASTQ found in $INPUT_DIR"; exit 1; }

step_banner 9 "smallRNA rescue alignment (Bowtie2)"
log_info "Species: $SPECIES | Samples: ${#inputs[@]}"

process_sample() {
  local in_file="$1"
  local idx="$2"
  local total="$3"
  local sample="${in_file##*/}"
  sample="${sample%_after_smallrna.fastq.gz}"

  local bam_out="${ALIGN_DIR}/${sample}_smallrna_bowtie2.bam"
  local unmapped_fq="${UNMAP_DIR}/${sample}_after_smallrna_bowtie2.fastq.gz"
  local log_file="${LOG_DIR}/${sample}.log"
  local fc_out="${COUNTS_DIR}/${sample}_smallrna_bowtie2_counts.txt"

  echo "  Processing: $sample ($idx/$total)"

  if check_bam "$bam_out" && check_tabular "$fc_out"; then
    log_skip "$sample — outputs exist"
    return 0
  fi

  local tmp_sam="${LOG_DIR}/${sample}_tmp.sam"
  local tmp_best="${LOG_DIR}/${sample}_tmp_best.sam"

  if ! check_bam "$bam_out"; then
    log_info "$sample — Bowtie2 local alignment..."

    bowtie2 -p "$THREADS" \
      --local --very-sensitive-local -L 15 -N 1 --score-min L,0,0.8 \
      --no-unal --all \
      -x "$REF_INDEX" -U "$in_file" --un-gz "$unmapped_fq" \
      -S "$tmp_sam" 2>> "$log_file"

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
    log_done "$sample — BAM created"
  fi

  if ! check_tabular "$fc_out" && check_bam "$bam_out"; then
    log_info "$sample — featureCounts..."
    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
      -M --fraction -O -s 0 \
      -o "$fc_out" "$bam_out" &>> "$log_file"
  fi

  if $DO_QC; then
    local qc_zip="${QC_FASTQC}/${sample}_after_smallrna_bowtie2_fastqc.zip"
    if check_gz "$unmapped_fq" && [[ ! -f "$qc_zip" ]]; then
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq" || true
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

combined_csv="${COUNTS_DIR}/smallrna_counts_bowtie2.csv"
log_info "Combining counts -> ${combined_csv}"
Rscript "$SCRIPT_DIR/combine_counts.R" \
  "$COUNTS_DIR" "_smallrna_bowtie2_counts\\.txt$" "_smallrna_bowtie2_counts.txt" \
  "$ANNOT_FILE" "$combined_csv"

if $DO_QC; then
  log_info "Running MultiQC..."
  multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force || true
fi

log_done "Step 9 — smallRNA Bowtie2 rescue complete"
exit $rc_all
