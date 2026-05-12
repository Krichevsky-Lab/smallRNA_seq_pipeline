#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 07_align_trna_rescue.sh
# Rescue tRNA-unmapped reads using Bowtie2 (local alignment, relaxed mismatch/gap
# penalties to recover reads with modification-induced RT mismatches).
#
# Usage:
#   bash alignment_scripts/07_align_trna_rescue.sh --species <mouse|human> [--qc]
#   SPECIES=mouse bash alignment_scripts/07_align_trna_rescue.sh [--qc]
#
# Inputs:
#   - alignment_output/unmapped/trna_bowtie1/*_after_trna.fastq.gz
#
# Outputs:
#   - alignment_output/aligned/trna_bowtie2/*.bam
#   - alignment_output/unmapped/trna_bowtie2/*_after_trna_bowtie2.fastq.gz
#   - alignment_output/counts/trna_bowtie2/*_counts.txt
#   - alignment_output/counts/trna_bowtie2/trna_counts_bowtie2.csv
#   - alignment_output/counts/trna_fragment/tRF_counts_bowtie2.csv
#   - alignment_output/counts/trna_fragment/per_sample_bowtie2/
#   - alignment_output/log/trna_bowtie2/*_bowtie2.log
#   - alignment_output/qc/trna_bowtie2/** (optional)
#
# Notes:
#   - Bowtie2 local, very-sensitive-local
#   - Reduced --mp/--rdg/--rfg penalties to recover reads with RT mismatches from
#     tRNA modifications (m1A, m3C, pseudouridine, inosine) that Bowtie1 (-v 1) missed
#   - Retains best AS:i alignments (ties kept)
#   - NH tags added for fractional counting
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

INPUT_DIR="alignment_output/unmapped/trna_bowtie1"
ALIGN_DIR="alignment_output/aligned/trna_bowtie2"
UNMAP_DIR="alignment_output/unmapped/trna_bowtie2"
COUNTS_DIR="alignment_output/counts/trna_bowtie2"
LOG_DIR="alignment_output/log/trna_bowtie2"
QC_DIR="alignment_output/qc/trna_bowtie2"
QC_FASTQC="${QC_DIR}/fastqc"
QC_MULTIQC="${QC_DIR}/multiqc"

mkdir -p "$ALIGN_DIR" "$UNMAP_DIR" "$COUNTS_DIR" "$LOG_DIR" \
         "$QC_FASTQC" "$QC_MULTIQC"

step_banner 7 "tRNA rescue alignment (Bowtie2)"
log_info "Species: $SPECIES"

inputs=( "$INPUT_DIR"/*_after_trna.fastq.gz )
[[ ${#inputs[@]} -gt 0 ]] || { log_fail "No FASTQ files found in $INPUT_DIR"; exit 1; }
n=${#inputs[@]}
i=0
for fq in "${inputs[@]}"; do
  i=$((i+1))
  base="$(basename "$fq" _after_trna.fastq.gz)"
  sample_base="${base%%_after_*}"

  bam_out="${ALIGN_DIR}/${sample_base}_trna_bowtie2.bam"
  unmapped_fq="${UNMAP_DIR}/${sample_base}_after_trna_bowtie2.fastq.gz"
  log_file="${LOG_DIR}/${sample_base}_bowtie2.log"
  tmp_sam="${ALIGN_DIR}/${sample_base}_tmp.sam"
  tmp_best="${ALIGN_DIR}/${sample_base}_tmp_best.sam"

  echo "  Processing: $sample_base ($i/$n)"

  if check_bam "$bam_out"; then
    log_skip "$sample_base — BAM exists"
  else
    log_info "$sample_base — Bowtie2 local alignment..."

    # Reduced mismatch/gap penalties to recover reads carrying tRNA modifications
    # (m1A, m3C, pseudouridine, inosine) that cause RT misincorporations rejected by Bowtie1.
    bowtie2 -p "$THREADS" \
      --local --very-sensitive-local -L 15 -N 1 \
      -D 20 -R 3 -i S,1,0.60 \
      --ma 2 --mp 2,2 --rdg 3,2 --rfg 3,2 \
      --score-min L,0,0.8 \
      --no-unal --all \
      -x "$REF_INDEX" -U "$fq" --un-gz "$unmapped_fq" \
      -S "$tmp_sam" 2>> "$log_file"

    if [[ ! -s "$tmp_sam" ]]; then
      log_warn "$sample_base — no alignments, skipping"
      continue
    fi

    awk -v OFS="\t" '
      BEGIN { q=""; c=0; maxAS=-999999 }
      /^@/ { print; next }
      {
        if (match($0,/AS:i:([0-9-]+)/,m)) ascore=m[1]; else ascore=-999999;
        if (q != "" && $1 != q) {
          for (i=1;i<=c;i++) if (as[i]==maxAS) best[i]=1;
          nBest=0; for (i=1;i<=c;i++) if (best[i]==1) nBest++;
          for (i=1;i<=c;i++) if (best[i]==1) print buf[i] "\tNH:i:" nBest;
          delete buf; delete as; delete best; c=0; maxAS=-999999;
        }
        q=$1; c++; buf[c]=$0; as[c]=ascore;
        if (ascore>maxAS) maxAS=ascore;
      }
      END {
        for (i=1;i<=c;i++) if (as[i]==maxAS) best[i]=1;
        nBest=0; for (i=1;i<=c;i++) if (best[i]==1) nBest++;
        for (i=1;i<=c;i++) if (best[i]==1) print buf[i] "\tNH:i:" nBest;
      }
    ' "$tmp_sam" > "$tmp_best"

    samtools view -@ "$THREADS" -bS "$tmp_best" \
      | samtools sort -@ "$THREADS" -o "$bam_out"

    samtools index "$bam_out"
    rm -f "$tmp_sam" "$tmp_best"
    log_done "$sample_base — BAM created"
  fi

  fc_out="${COUNTS_DIR}/${sample_base}_trna_bowtie2_counts.txt"
  if ! check_tabular "$fc_out" && check_bam "$bam_out"; then
    log_info "$sample_base — featureCounts..."
    featureCounts -T "$THREADS" -F SAF -a "$SAF_FILE" \
      -M --fraction -O -s 0 -o "$fc_out" "$bam_out" &>> "$log_file" \
      || log_warn "featureCounts failed for $sample_base"
  fi

  if $DO_QC; then
    qc_zip="${QC_FASTQC}/${sample_base}_after_trna_bowtie2_fastqc.zip"
    if check_gz "$unmapped_fq" && [[ ! -f "$qc_zip" ]]; then
      fastqc -q -t "$THREADS" -o "$QC_FASTQC" "$unmapped_fq" \
        || log_warn "FastQC failed for $sample_base"
    fi
  fi
done

if $DO_QC; then
  if compgen -G "${QC_FASTQC}/*fastqc.zip" > /dev/null; then
    log_info "Running MultiQC..."
    multiqc "$QC_FASTQC" -o "$QC_MULTIQC" --force || log_warn "MultiQC encountered issues."
  fi
fi

combined_csv="${COUNTS_DIR}/trna_counts_bowtie2.csv"
if ! check_tabular "$combined_csv"; then
  log_info "Combining counts -> ${combined_csv}"
  Rscript "$SCRIPT_DIR/combine_counts.R" \
    "$COUNTS_DIR" "_trna_bowtie2_counts\\.txt$" "_trna_bowtie2_counts.txt" \
    "$ANNOT_FILE" "$combined_csv"
fi

TRF_DIR="alignment_output/counts/trna_fragment"
mkdir -p "$TRF_DIR"
trf_csv="${TRF_DIR}/tRF_counts_bowtie2.csv"
if ! check_tabular "$trf_csv"; then
  log_info "Running tRF classification (Bowtie2)..."
  Rscript "$SCRIPT_DIR/classify_tRF.R" \
    "$ALIGN_DIR" "$SAF_FILE" "$ANNOT_FILE" \
    "$trf_csv" "${TRF_DIR}/per_sample_bowtie2"
fi

log_done "Step 7 — tRNA Bowtie2 rescue + tRF classification complete"
