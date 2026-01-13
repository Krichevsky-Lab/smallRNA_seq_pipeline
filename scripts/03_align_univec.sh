#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ------------------------------------------------------------------------------
# 3_align_univec.sh
# Filter UniVec contaminants using Bowtie1 (single-end).
#
# Usage:
#   bash scripts/3_align_univec.sh [--qc]
#
# Inputs:
#   - trimmed/*.fastq.gz
#
# Outputs:
#   - unmapped/univec/*_after_univec.fq.gz
#   - log/univec/*.log
#   - qc/univec/fastqc/ (optional)
#   - qc/univec/multiqc/ (optional)
#
# Notes:
#   - Uses Bowtie1 with UniVec reference
#   - QC mode runs FastQC + MultiQC
# ------------------------------------------------------------------------------

THREADS=10
DO_QC=false
[[ "${1:-}" == "--qc" ]] && DO_QC=true

REF_INDEX="reference/univec/UniVec"
INPUT_DIR="trimmed"
OUT_DIR="unmapped/univec"
LOG_DIR="log/univec"

QC_BASE="qc/univec"
FASTQC_DIR="${QC_BASE}/fastqc"
MULTIQC_DIR="${QC_BASE}/multiqc"

mkdir -p "$OUT_DIR" "$LOG_DIR"
[[ "$DO_QC" == true ]] && mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR"

echo "============================================================"
echo " 🎯 Step 3: Filter UniVec contaminants"
echo " Input:   $INPUT_DIR"
echo " Output:  $OUT_DIR"
echo " QC:      $DO_QC"
echo "============================================================"

# -----------------------------
# Collect FASTQs safely
# -----------------------------
fqs=("$INPUT_DIR"/*.fastq.gz)

if [[ ${#fqs[@]} -eq 0 ]]; then
  echo "❌ No FASTQ files found in $INPUT_DIR"
  exit 1
fi

# -----------------------------
# Bowtie filtering
# -----------------------------
for fq in "${fqs[@]}"; do

  sample="$(basename "$fq" .fastq.gz)"
  out_fq="${OUT_DIR}/${sample}_after_univec.fq.gz"
  log_file="${LOG_DIR}/${sample}_univec.log"

  echo ""
  echo "------------------------------------------------"
  echo " 🔧 Processing: $sample"
  echo "------------------------------------------------"

  if [[ -s "$out_fq" ]]; then
    echo " ⏩ Output exists — skipping."
    continue
  fi

  echo " 🔨 Running Bowtie1..."
  bowtie -p "$THREADS" \
    -x "$REF_INDEX" \
    -q "$fq" \
    --un >(gzip > "$out_fq") \
    /dev/null \
    > "$log_file" 2>&1 \
    || echo "⚠️ Bowtie failed for $sample (see $log_file)"

  echo " ✅ UniVec filtering complete"

  # -----------------------------
  # Optional per-sample QC
  # -----------------------------
  if [[ "$DO_QC" == true ]]; then
    qc_zip="${FASTQC_DIR}/${sample}_after_univec_fastqc.zip"
    if [[ -s "$qc_zip" ]]; then
      echo " ⏩ FastQC exists — skipping."
    else
      fastqc -t "$THREADS" -o "$FASTQC_DIR" "$out_fq"
    fi
  fi
done

# -----------------------------
# MultiQC (optional)
# -----------------------------
if [[ "$DO_QC" == true ]]; then
  echo ""
  echo "------------------------------------------------"
  echo " 📊 Running MultiQC..."
  echo "------------------------------------------------"

  multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" --force \
    || echo "⚠️ MultiQC returned warnings — continuing."

  echo " 📄 Report: $MULTIQC_DIR/multiqc_report.html"
fi

echo ""
echo "============================================================"
echo " ✅ UniVec filtering complete"
echo "============================================================"

exit 0
