#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ------------------------------------------------------------------------------
# 3_align_univec.sh
# Filter UniVec contaminants using Bowtie1 (single-end).
#
# Usage:
#   bash alignment_scripts/03_align_univec.sh --species <mouse|human> [--qc]
#   SPECIES=mouse bash alignment_scripts/03_align_univec.sh [--qc]
#
# Inputs:
#   - alignment_output/trimmed/*.fastq.gz
#
# Outputs:
#   - alignment_output/unmapped/univec/*_after_univec.fq.gz
#   - alignment_output/log/univec/*.log
#   - alignment_output/qc/univec/fastqc/ (optional)
#   - alignment_output/qc/univec/multiqc/ (optional)
#
# Notes:
#   - Uses Bowtie1 with UniVec reference
#   - QC mode runs FastQC + MultiQC
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

REF_INDEX="reference/$SPECIES/univec/${SPECIES}_univec"

require_index "$REF_INDEX"

INPUT_DIR="alignment_output/trimmed"
OUT_DIR="alignment_output/unmapped/univec"
LOG_DIR="alignment_output/log/univec"

QC_BASE="alignment_output/qc/univec"
FASTQC_DIR="${QC_BASE}/fastqc"
MULTIQC_DIR="${QC_BASE}/multiqc"

mkdir -p "$OUT_DIR" "$LOG_DIR"
[[ "$DO_QC" == true ]] && mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR"

step_banner 3 "UniVec contaminant filtering"
log_info "Species: $SPECIES | Input: $INPUT_DIR | Output: $OUT_DIR | QC: $DO_QC"

# -----------------------------
# Collect FASTQs safely
# -----------------------------
fqs=("$INPUT_DIR"/*.fastq.gz)

if [[ ${#fqs[@]} -eq 0 ]]; then
  log_fail "No FASTQ files found in $INPUT_DIR"
  exit 1
fi

# -----------------------------
# Bowtie filtering
# -----------------------------
n=${#fqs[@]}
i=0
for fq in "${fqs[@]}"; do
  i=$((i+1))
  sample="$(basename "$fq" .fastq.gz)"
  out_fq="${OUT_DIR}/${sample}_after_univec.fq.gz"
  log_file="${LOG_DIR}/${sample}_univec.log"

  echo ""
  echo "  Processing: $sample ($i/$n)"

  if check_gz "$out_fq"; then
    log_skip "$sample — already filtered"
    continue
  fi

  log_info "Running Bowtie1..."
  bowtie -p "$THREADS" \
    -x "$REF_INDEX" \
    -q "$fq" \
    --un >(gzip > "$out_fq") \
    /dev/null \
    > "$log_file" 2>&1 \
    || log_warn "Bowtie failed for $sample (see $log_file)"

  log_done "$sample — UniVec filtering done"

  if [[ "$DO_QC" == true ]]; then
    qc_zip="${FASTQC_DIR}/${sample}_after_univec_fastqc.zip"
    if [[ -f "$qc_zip" ]]; then
      log_skip "FastQC exists for $sample"
    else
      fastqc -q -t "$THREADS" -o "$FASTQC_DIR" "$out_fq"
    fi
  fi
done

if [[ "$DO_QC" == true ]]; then
  log_info "Running MultiQC..."
  multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" --force \
    || log_warn "MultiQC returned warnings — continuing."
  log_done "Report: $MULTIQC_DIR/multiqc_report.html"
fi

echo ""
log_done "Step 3 — UniVec filtering complete"

exit 0
