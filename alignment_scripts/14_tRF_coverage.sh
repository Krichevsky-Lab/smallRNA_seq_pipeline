#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# ------------------------------------------------------------------------------
# 14_tRF_coverage.sh
# Compute per-interval tRNA fragment coverage from Bowtie1 + Bowtie2 BAMs.
#
# Usage:
#   bash alignment_scripts/14_tRF_coverage.sh --species <mouse|human>
#   SPECIES=mouse bash alignment_scripts/14_tRF_coverage.sh
# ------------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common.sh"

SPECIES="${SPECIES:-}"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --species) SPECIES="$2"; shift ;;
        --qc)      ;;  # accepted but unused
        *)         log_fail "Unknown argument: $1"; exit 1 ;;
    esac
    shift
done
[[ -n "$SPECIES" ]] || { log_fail "--species required (mouse or human). Set via arg or \$SPECIES env."; exit 1; }
[[ "$SPECIES" == "mouse" || "$SPECIES" == "human" ]] || { log_fail "Unknown species '$SPECIES'"; exit 1; }

DIR1="alignment_output/aligned/trna_bowtie1"
DIR2="alignment_output/aligned/trna_bowtie2"
ANNOT="reference/$SPECIES/trna/${SPECIES}_tRNA.saf"

require_file "$ANNOT"

mkdir -p "analysis/data"

OUTFILE="analysis/data/tRNA_intervals.tsv"
TMP=$(mktemp)
LOOKUP=$(mktemp)
trap 'rm -f "$TMP" "$LOOKUP"' EXIT

step_banner 14 "tRNA fragment (tRF) coverage"
log_info "Species: $SPECIES"

# =========================
# BUILD LOOKUP TABLE
# =========================
log_info "Building tRNA coordinate lookup..."

awk 'BEGIN{OFS="\t"}
NR>1 {
  trna_start[$1] = $3
  trna_end[$1]   = $4
}
END {
  for (k in trna_start) {
    print k, trna_start[k], trna_end[k]
  }
}' "$ANNOT" > "$LOOKUP"

# =========================
# FIND SAMPLES
# =========================
bam_files=( "$DIR1"/*.bam )
[[ ${#bam_files[@]} -gt 0 ]] || { log_fail "No BAM files found in $DIR1"; exit 1; }

n=${#bam_files[@]}
log_info "Found $n samples"

# =========================
# HEADER
# =========================
echo -e "sample\ttRNA\tstart\tend\tstrand\tlength\tfragment\tweight" > "$TMP"

# =========================
# PROCESS EACH SAMPLE
# =========================
i=0
for bam1 in "${bam_files[@]}"; do
  i=$((i+1))
  fname=$(basename "$bam1")
  sample="${fname%_trna_bowtie1.bam}"
  sample_clean="${sample%_trimmed}"

  bam2="${DIR2}/${sample}_trna_bowtie2.bam"

  if [[ ! -f "$bam2" ]]; then
    log_warn "Skipping $sample_clean (missing Bowtie2 BAM: $bam2)"
    continue
  fi

  log_info "Processing: $sample_clean ($i/$n)"

  (
    samtools view -F 0x900 "$bam1"
    samtools view -F 0x900 "$bam2"
  ) | awk -v SAMPLE="$sample_clean" -v LOOKUP="$LOOKUP" '

  BEGIN {
    OFS="\t"

    while ((getline < LOOKUP) > 0) {
      trna_start[$1] = $2
      trna_end[$1]   = $3
    }
  }

  {
    flag  = $2
    tRNA  = $3
    start = $4
    cigar = $6

    if (tRNA == "*" || cigar == "*") next
    if (!(tRNA in trna_start)) next

    # Parse CIGAR for reference-consuming length (POSIX awk compatible)
    aligned_len = 0
    c = cigar

    while (match(c, /[0-9]+[MIDNSHP=X]/)) {
      tok = substr(c, RSTART, RLENGTH)
      op  = substr(tok, RLENGTH, 1)
      len = substr(tok, 1, RLENGTH - 1) + 0
      if (op == "M" || op == "D" || op == "=" || op == "X") {
        aligned_len += len
      }
      c = substr(c, RSTART + RLENGTH)
    }

    if (aligned_len == 0) next

    end = start + aligned_len - 1

    strand = "+"
    if (and(flag,16)) strand = "-"

    nh = 1
    if (match($0, /NH:i:([0-9]+)/, a)) nh = a[1]
    weight = 1 / nh

    tstart = trna_start[tRNA]
    tend   = trna_end[tRNA]
    tlen   = tend - tstart + 1

    rel_start = start - tstart
    rel_end   = end   - tstart

    near5 = (rel_start <= 10)
    near3 = (rel_end >= (tlen - 10))

    if (near5 && near3) fragment = "FL-tRNA"
    else if (near5)     fragment = "5prime-tRF"
    else if (near3)     fragment = "3prime-tRF"
    else                fragment = "internal-tRF"

    print SAMPLE, tRNA, start, end, strand, aligned_len, fragment, weight
  }' >> "$TMP"

done

# =========================
# COLLAPSE IDENTICAL FRAGMENTS
# =========================
log_info "Collapsing identical intervals..."

{
  echo -e "sample\ttRNA\tstart\tend\tstrand\tlength\tfragment\tcount"
  tail -n +2 "$TMP" | \
  awk 'BEGIN{OFS="\t"}
  {
    key = $1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7
    sum[key] += $8
  }
  END {
    for (k in sum) {
      print k, sum[k]
    }
  }'
} > "$OUTFILE"

log_done "Step 14 — tRF coverage: ${OUTFILE}"
exit 0
