#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------------
# 2_trim.sh
# Trim TruSeq Small RNA adapters using Cutadapt.
#
# Adapter:
#   TGGAATTCTCGGGTGCCAAGG
#
# Usage:
#   bash scripts/2_trim.sh [--qc]
#
# Inputs:
#   - raw_data/*.fastq.gz
#
# Outputs:
#   - trimmed/*_trimmed.fastq.gz
#   - log/trim/*.log
#   - qc/trimmed/fastqc/ (optional)
#   - qc/trimmed/multiqc/ (optional)
#
# Notes:
#   - Enforces length >=15 nt
#   - Filters low-quality and ambiguous reads
#   - QC mode runs FastQC + MultiQC
# ------------------------------------------------------------------------------

RAW_DIR="raw_data"
TRIM_DIR="trimmed"
LOG_DIR="log/trim"
QC_BASE="qc/trimmed"
FASTQC_DIR="${QC_BASE}/fastqc"
MULTIQC_DIR="${QC_BASE}/multiqc"

ADAPTER_3P="TGGAATTCTCGGGTGCCAAGG"
MIN_LEN=15

# -----------------------------
# Optional QC mode
# -----------------------------
DO_QC="false"
if [[ "${1:-}" == "--qc" ]]; then
  DO_QC="true"
fi

# -----------------------------
# Create folders
# -----------------------------
mkdir -p "$TRIM_DIR" "$LOG_DIR"
if [[ "$DO_QC" == "true" ]]; then
  mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR"
fi

echo "================================================"
echo " Step 2: TruSeq Small RNA adapter trimming"
echo " Adapter:     ${ADAPTER_3P}"
echo " Min length:  >${MIN_LEN} nt"
echo " QC enabled:  ${DO_QC}"
echo "================================================"

# -----------------------------
# Loop over samples
# -----------------------------
for fq in ${RAW_DIR}/*.fastq.gz; do
  [[ -e "$fq" ]] || { echo "❌ No FASTQ files found in ${RAW_DIR}"; exit 1; }

  fname="$(basename "$fq")"
  sample="${fname%.fastq.gz}"

  out="${TRIM_DIR}/${sample}_trimmed.fastq.gz"
  log="${LOG_DIR}/${sample}.log"

  echo ""
  echo "------------------------------------------------"
  echo "[`date '+%Y-%m-%d %H:%M:%S'`] Processing: ${sample}"
  echo "------------------------------------------------"

  if [[ -f "$out" ]]; then
    echo " ➤ Already trimmed — skipping."
    continue
  fi

  echo " ➤ Running Cutadapt..."
  cutadapt \
    -a "$ADAPTER_3P" \
    -O 5 \
    -e 0.1 \
    -q 20 \
    --max-n 0 \
    -m "$MIN_LEN" \
    -o "$out" \
    "$fq" \
    > "$log" 2>&1

  echo " ➤ Trim complete"
  echo " ➤ Log: $log"

  if [[ "$DO_QC" == "true" ]]; then
    echo " ➤ Running FastQC..."
    fastqc -o "$FASTQC_DIR" "$out" > "${FASTQC_DIR}/${sample}_fastqc.log" 2>&1
  fi
done

# -----------------------------
# MultiQC
# -----------------------------
if [[ "$DO_QC" == "true" ]]; then
  echo ""
  echo "------------------------------------------------"
  echo " Running MultiQC..."
  echo "------------------------------------------------"

  multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" > "${MULTIQC_DIR}/multiqc.log" 2>&1
  echo " ➤ MultiQC report generated"
fi

echo ""
echo "===================================="
echo " TruSeq trimming — ALL DONE"
echo "===================================="
