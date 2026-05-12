#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------------
# 02_trim.sh
# Trim Illumina TruSeq Small RNA adapters using Cutadapt.
#
# Adapter (canonical):
#   TGGAATTCTCGGGTGCCAAGG
#
# Usage:
#   bash alignment_scripts/02_trim.sh [--species <mouse|human>] [--qc]
#
# Inputs:
#   - alignment_output/raw_data/*.fastq.gz OR *.fq.gz
#
# Outputs:
#   - alignment_output/trimmed/*_trimmed.fastq.gz
#   - alignment_output/log/trim/*.log
#   - alignment_output/qc/trimmed/fastqc/ (optional)
#   - alignment_output/qc/trimmed/multiqc/ (optional)
#
# Notes:
#   - Handles partial/truncated adapters
#   - Enforces minimum length and quality from config.sh
#   - Does NOT discard untrimmed reads
# ------------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common.sh"

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

RAW_DIR="alignment_output/raw_data"
TRIM_DIR="alignment_output/trimmed"
LOG_DIR="alignment_output/log/trim"
QC_BASE="alignment_output/qc/trimmed"
FASTQC_DIR="${QC_BASE}/fastqc"
MULTIQC_DIR="${QC_BASE}/multiqc"

ADAPTERS=(
  "TGGAATTCTCGGGTGCCAAGG"
  "GGAATTCTCGGGTGCCAAGG"
  "AATTCTCGGGTGCCAAGG"
  "ATTCTCGGGTGCCAAGG"
)

mkdir -p "$TRIM_DIR" "$LOG_DIR"
[[ "$DO_QC" == true ]] && mkdir -p "$FASTQC_DIR" "$MULTIQC_DIR"

shopt -s nullglob
FASTQS=( "$RAW_DIR"/*.fastq.gz "$RAW_DIR"/*.fq.gz )

if [[ ${#FASTQS[@]} -eq 0 ]]; then
  log_fail "No FASTQ files found in ${RAW_DIR}"
  exit 1
fi

step_banner 2 "Adapter trimming (Cutadapt)"
log_info "Adapters: ${ADAPTERS[*]}"
log_info "Min length: ${MIN_READ_LEN} nt  |  Quality: Phred>=${PHRED_QUAL}  |  QC: ${DO_QC}"

for fq in "${FASTQS[@]}"; do
  fname="$(basename "$fq")"
  sample="${fname%.fastq.gz}"
  sample="${sample%.fq.gz}"

  out="${TRIM_DIR}/${sample}_trimmed.fastq.gz"
  log="${LOG_DIR}/${sample}.log"

  echo ""
  echo "  Processing: ${sample}"

  if check_gz "$out"; then
    log_skip "${sample} — already trimmed"
    continue
  fi

  ADAPT_ARGS=()
  for a in "${ADAPTERS[@]}"; do
    ADAPT_ARGS+=( -a "$a" )
  done

  cutadapt \
    "${ADAPT_ARGS[@]}" \
    -a "A{10}" \
    -O 5 \
    -e 0.1 \
    --trim-n \
    --max-n 0 \
    -q "$PHRED_QUAL" \
    -m "$MIN_READ_LEN" \
    -o "$out" \
    "$fq" \
    > "$log" 2>&1

  log_done "trimmed: $out"

  if [[ "$DO_QC" == true ]]; then
    log_info "FastQC on ${sample}..."
    fastqc -t "$THREADS" -o "$FASTQC_DIR" "$out" > "${FASTQC_DIR}/${sample}_fastqc.log" 2>&1
  fi
done

if [[ "$DO_QC" == true ]]; then
  echo ""
  log_info "Running MultiQC..."
  multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR" --force > "${MULTIQC_DIR}/multiqc.log" 2>&1
  log_done "MultiQC report: ${MULTIQC_DIR}/multiqc_report.html"
fi

echo ""
log_done "Step 2 — adapter trimming complete"
