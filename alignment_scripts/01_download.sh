#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------------
# 01_download.sh
# Download raw FASTQ files from NCBI SRA using input/runlist.txt.
#
# Usage:
#   bash alignment_scripts/01_download.sh [--qc]
#
# Inputs:
#   - input/runlist.txt   (one SRR accession per line)
#
# Outputs:
#   - alignment_output/raw_data/*.fastq.gz
#   - alignment_output/log/download/download.log
#   - alignment_output/qc/raw/fastqc/   (optional)
#   - alignment_output/qc/raw/multiqc/  (optional)
#
# Notes:
#   - Skips files already present in raw_data/
#   - --qc runs FastQC + MultiQC on raw_data/ after download
# ------------------------------------------------------------------------------

SPECIES="${SPECIES:-}"
DO_QC=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --species) SPECIES="$2"; shift ;;
        --qc)      DO_QC=true ;;
        *)         echo "[FAIL]  Unknown argument: $1"; exit 1 ;;
    esac
    shift
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common.sh"

RUNLIST="input/runlist.txt"

if [[ ! -f "$RUNLIST" ]]; then
  echo "[FAIL]  Not found: $RUNLIST"
  exit 1
fi

# ---- Directories ----
RAW_DIR="alignment_output/raw_data"
LOG_DIR="alignment_output/log/download"
LOG_FILE="${LOG_DIR}/download.log"
QC_DIR="alignment_output/qc/raw"
FASTQC_DIR="${QC_DIR}/fastqc"
MULTIQC_DIR="${QC_DIR}/multiqc"

mkdir -p "$RAW_DIR" "$LOG_DIR"

{
  echo "============================================================"
  echo "[INFO]  Step 1: Download raw FASTQs from SRA"
  echo "Runlist: $RUNLIST"
  echo "Output:  $RAW_DIR"
  echo "QC:      $DO_QC"
  echo "============================================================"
} | tee "$LOG_FILE"

# ==============================================================================
# Download loop
# ==============================================================================
while read -r SRR || [[ -n "$SRR" ]]; do
  SRR="$(echo "$SRR" | tr -d '\r\n ')"
  [[ -z "$SRR" ]] && continue

  OUTFILE="${RAW_DIR}/${SRR}.fastq.gz"

  if [[ -s "$OUTFILE" ]]; then
    echo "  [SKIP]  $SRR already exists — skipping." | tee -a "$LOG_FILE"
    continue
  fi

  echo "  [INFO]  Downloading $SRR ..." | tee -a "$LOG_FILE"
  wget -q --show-progress -O "$OUTFILE" \
    "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=${SRR}" || true

  if [[ -s "$OUTFILE" ]]; then
    echo "  [DONE]  $SRR downloaded." | tee -a "$LOG_FILE"
  else
    echo "  [FAIL]  Failed to download $SRR" | tee -a "$LOG_FILE"
    rm -f "$OUTFILE"
  fi
done < "$RUNLIST"

# ==============================================================================
# Optional QC
# ==============================================================================
if [[ "$DO_QC" == true ]]; then
  shopt -s nullglob
  raw_files=( "$RAW_DIR"/*.fastq.gz )

  if [[ ${#raw_files[@]} -eq 0 ]]; then
    echo "[WARN]  No FASTQ files in $RAW_DIR — skipping QC." | tee -a "$LOG_FILE"
  else
    mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR"

    echo "[INFO]  Running FastQC on ${#raw_files[@]} file(s)..." | tee -a "$LOG_FILE"
    fastqc -q -t "$THREADS" -o "$FASTQC_DIR" "${raw_files[@]}" 2>&1 | tee -a "$LOG_FILE"

    echo "[INFO]  Running MultiQC..." | tee -a "$LOG_FILE"
    multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" --force 2>&1 | tee -a "$LOG_FILE"

    echo "[INFO]  Report: $MULTIQC_DIR/multiqc_report.html" | tee -a "$LOG_FILE"
  fi
fi

{
  echo "============================================================"
  echo "[DONE]  Step 1 complete — logs: $LOG_FILE"
  echo "============================================================"
} | tee -a "$LOG_FILE"
