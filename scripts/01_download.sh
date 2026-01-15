#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------------
# 1_download.sh
# Download raw FASTQ files from NCBI SRA using an SRR runlist.
#
# Usage:
#   bash 1_download.sh
#
# Inputs:
#   - input/runlist.txt
#
# Outputs:
#   - raw_data/*.fastq.gz
#   - log/download.log
#
# Notes:
#   - Skips existing FASTQs
#   - Intended to be run as Step 1 of the pipeline
# ------------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RUNLIST="${SCRIPT_DIR}/../input/runlist.txt"

if [[ ! -f "$RUNLIST" ]]; then
  echo "❌ runlist.txt not found at $RUNLIST"
  exit 1
fi

# ---- Directories ----
RAW_DIR="raw_data"
LOG_DIR="log/download"
LOG_FILE="${LOG_DIR}/download.log"

mkdir -p "$RAW_DIR" "$LOG_DIR"


{
  echo "============================================================"
  echo "📥 Step 1: Downloading raw FASTQs"
  echo "Runlist: $RUNLIST"
  echo "Output:  $RAW_DIR"
  echo "============================================================"
} | tee "$LOG_FILE"

# ---- Download loop ----
while read -r SRR || [[ -n "$SRR" ]]; do
  SRR="$(echo "$SRR" | tr -d '\r\n ')"
  [[ -z "$SRR" ]] && continue

  OUTFILE="${RAW_DIR}/${SRR}.fastq.gz"

  if [[ -f "$OUTFILE" ]]; then
    echo "✅ $SRR already exists — skipping." | tee -a "$LOG_FILE"
    continue
  fi

  echo "🔽 Downloading $SRR ..." | tee -a "$LOG_FILE"
  wget -q --show-progress -O "$OUTFILE" \
    "https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=${SRR}" || true

  if [[ -s "$OUTFILE" ]]; then
    echo "✅ Successfully downloaded $SRR" | tee -a "$LOG_FILE"
  else
    echo "❌ Failed to download $SRR" | tee -a "$LOG_FILE"
    rm -f "$OUTFILE"
  fi
done < "$RUNLIST"

{
  echo "============================================================"
  echo "🎉 Download step complete"
  echo "Logs saved to: $LOG_FILE"
  echo "============================================================"
} | tee -a "$LOG_FILE"
