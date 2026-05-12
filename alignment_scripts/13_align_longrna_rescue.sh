#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 13_align_longrna_rescue.sh
# Rescue unmapped longRNA reads using Bowtie2 (local mode).
#
# Usage:
#   bash alignment_scripts/13_align_longrna_rescue.sh --species <mouse|human> [--qc]
#   SPECIES=mouse bash alignment_scripts/13_align_longrna_rescue.sh [--qc]
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
#
# Notes:
#   - Bowtie2 local; best AS:i retained + NH tags
#   - featureCounts: -M --fraction -O
# ------------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common.sh"

SPECIES="${SPECIES:-}"
DO_QC=false
rc_all=0
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

REF_INDEX="reference/$SPECIES/longrna/${SPECIES}_longRNA"
SAF_FILE="reference/$SPECIES/longrna/${SPECIES}_longRNA.saf"
ANNOT_FILE="reference/$SPECIES/longrna/${SPECIES}_longRNA_annotation.tsv"

require_index "$REF_INDEX"
require_file  "$SAF_FILE"
require_file  "$ANNOT_FILE"

INPUT_DIR="alignment_output/unmapped/longrna_bowtie1"
ALIGN_DIR="alignment_output/aligned/longrna_bowtie2"
UNMAP_DIR="alignment_output/unmapped/longrna_bowtie2"
COUNTS_DIR="alignment_output/counts/longrna_bowtie2"
LOG_DIR="alignment_output/log/longrna_bowtie2"
QC_DIR="alignment_output/qc/longrna_bowtie2"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$LOG_DIR"
$DO_QC && mkdir -p "$QC_FASTQC" "$QC_MULTIQC"

inputs=( "$INPUT_DIR"/*_after_longrna.fastq.gz )
[[ ${#inputs[@]} -gt 0 ]] || { log_fail "No FASTQ files found in $INPUT_DIR"; exit 1; }

step_banner 13 "longRNA rescue alignment (Bowtie2)"
log_info "Species: $SPECIES | QC: $DO_QC"

n=${#inputs[@]}
i=0
for fq in "${inputs[@]}"; do
  i=$((i+1))
  sample=$(basename "$fq" "_after_longrna.fastq.gz")

  bam_out="${ALIGN_DIR}/${sample}_longrna_bowtie2.bam"
  unmapped_fq="${UNMAP_DIR}/${sample}_after_longrna_bowtie2.fastq.gz"
  log_file="${LOG_DIR}/${sample}.log"

  echo "  Processing: $sample ($i/$n)"

  if ! check_bam "$bam_out"; then
    log_info "$sample — Bowtie2 local alignment..."

    tmp_sam="${LOG_DIR}/${sample}_tmp.sam"
    tmp_best="${LOG_DIR}/${sample}_tmp_best.sam"

    if ! bowtie2 -p "$THREADS" \
        --local --very-sensitive-local -L 20 -N 1 --score-min L,0,0.9 \
        --no-unal --all \
        -x "$REF_INDEX" -U "$fq" \
        --un-gz "$unmapped_fq" \
        -S "$tmp_sam" 2>> "$log_file"; then
      log_fail "Bowtie2 failed for $sample"
      rc_all=1
      continue
    fi

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
    log_done "$sample — BAM created"
  else
    log_skip "$sample — BAM exists"
  fi

  fc_out="${COUNTS_DIR}/${sample}_longrna_bowtie2_counts.txt"
  if ! check_tabular "$fc_out" && check_bam "$bam_out"; then
    log_info "$sample — featureCounts..."
    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
      -M --fraction -O -s 0 \
      -o "$fc_out" "$bam_out" &>> "$log_file" \
      || log_warn "featureCounts failed for $sample"
  fi

  if $DO_QC; then
    qc_zip="${QC_FASTQC}/${sample}_after_longrna_bowtie2_fastqc.zip"
    if check_gz "$unmapped_fq" && [[ ! -f "$qc_zip" ]]; then
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq" \
        || log_warn "FastQC failed for $sample"
    fi
  fi
done

combined_csv="${COUNTS_DIR}/longrna_counts_bowtie2.csv"
log_info "Combining counts -> $combined_csv"
Rscript "$SCRIPT_DIR/combine_counts.R" \
  "$COUNTS_DIR" "_longrna_bowtie2_counts\\.txt$" "_longrna_bowtie2_counts.txt" \
  "$ANNOT_FILE" "$combined_csv"

if $DO_QC; then
  log_info "Running MultiQC..."
  multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force \
    || log_warn "MultiQC warnings"
fi

log_done "Step 13 — longRNA Bowtie2 rescue complete"
exit $rc_all
