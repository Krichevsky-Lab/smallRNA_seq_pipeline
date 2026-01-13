#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

################################################################################
# 00_make_reference.sh
#
# Build reference indices and annotation files from *_all.fa FASTA files.
#
# UniVec:
#   - Bowtie1 index ONLY
#   - No FASTA / SAF / annotation
#
# All other references:
#   - Simplified FASTA (>genename)
#   - SAF
#   - Annotation TSV (source = <folder>_alignment)
#   - Bowtie1 + Bowtie2 indices
#
# Manifest (CSV):
#   reference_name, folder, source, total_entries, biotype_1, biotype_1_count
#
# Usage:
#   bash scripts/00_make_reference.sh
################################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REFERENCE_DIR="$PROJECT_ROOT/reference"
MANIFEST="$REFERENCE_DIR/reference_manifest.csv"

# ------------------ init manifest ------------------
echo "reference_name,folder,source,total_entries,biotype_1,biotype_1_count" \
  > "$MANIFEST"

LAST_CLASS=""

# ------------------------------------------------------------------------------
# Loop over references
# ------------------------------------------------------------------------------
for INPUT in "$REFERENCE_DIR"/**/*_all.fa "$REFERENCE_DIR"/**/*_all.fasta; do
  [[ -e "$INPUT" ]] || continue

  INPUT="$(realpath "$INPUT")"
  OUTDIR="$(dirname "$INPUT")"
  BASE="$(basename "$INPUT")"

  BASENAME="${BASE%_all.fa}"
  BASENAME="${BASENAME%_all.fasta}"
  REF_FOLDER="$(basename "$OUTDIR")"
  SOURCE="${REF_FOLDER}_alignment"
  INDEX_PREFIX="${OUTDIR}/${BASENAME}"

  # Print RNA class header once
  if [[ "$REF_FOLDER" != "$LAST_CLASS" ]]; then
    echo "🧬 Making reference for ${REF_FOLDER}"
    LAST_CLASS="$REF_FOLDER"
  fi

  # Detect UniVec
  IS_UNIVEC=false
  [[ "$OUTDIR" == *"/univec"* || "$BASE" =~ [Uu]ni[Vv]ec ]] && IS_UNIVEC=true

  # ============================================================================
  # UniVec (Bowtie1 only, no annotation → still counted)
  # ============================================================================
  if [[ "$IS_UNIVEC" == true ]]; then
    [[ -f "${INDEX_PREFIX}.1.ebwt" ]] || \
      bowtie-build "$INPUT" "$INDEX_PREFIX" > /dev/null 2>&1

    total_n=$(grep -c "^>" "$INPUT")

    # UniVec has no biotypes → record as NA
    echo "$BASENAME,$REF_FOLDER,univec_alignment,$total_n,NA,NA" \
      >> "$MANIFEST"

    echo "✅ Completed: $BASENAME"
    continue
  fi

  FA_SIMPLE="${OUTDIR}/${BASENAME}.fa"
  SAF_OUT="${OUTDIR}/${BASENAME}.saf"
  ANN_OUT="${OUTDIR}/${BASENAME}_annotation.tsv"

  # ------------------ Simplified FASTA ------------------
  awk -F"|" '
    /^>/ { gsub(/^>/,"",$1); print ">"$1; next }
    { print }
  ' "$INPUT" > "$FA_SIMPLE"

  dups=$(grep "^>" "$FA_SIMPLE" | sort | uniq -d)
  if [[ -n "$dups" ]]; then
    echo "❌ Duplicate genenames detected in $FA_SIMPLE"
    exit 1
  fi

  # ------------------ SAF ------------------
  echo -e "GeneID\tChr\tStart\tEnd\tStrand" > "$SAF_OUT"
  awk '
    /^>/ {
      if (seq != "") print name"\t"name"\t1\t"length(seq)"\t*"
      name=substr($1,2); seq=""; next
    }
    { seq=seq$0 }
    END {
      if (seq != "") print name"\t"name"\t1\t"length(seq)"\t*"
    }
  ' "$FA_SIMPLE" >> "$SAF_OUT"

  # ------------------ Annotation ------------------
  echo -e "genename\tgeneid\tbiotype_1\tbiotype_2\tchromosome\tstart\tend\tstrand\tlength\tsource" \
    > "$ANN_OUT"

  awk -F"|" -v SRC="$SOURCE" '
    /^>/ {
      gsub(/^>/,"",$1)
      genename=$1; rawid=$2; biotype1=$3; biotype2=$3; header=$4
      geneid=(rawid ~ /^ENS|^NM_|^NR_|^XM_|^XR_/)?rawid:genename

      if (match(header,/^([^:]+):([0-9]+)-([0-9]+)\(([+-])\)$/,m)) {
        chr=m[1]; start=m[2]; end=m[3]; strand=m[4]
      } else if (match(header,/^([^:]+):([0-9]+)-([0-9]+)$/,m)) {
        chr=m[1]; start=m[2]; end=m[3]; strand="*"
      } else {
        chr=header; start="NA"; end="NA"; strand="*"
      }

      flen=(start~/^[0-9]+$/ && end~/^[0-9]+$/)?end-start+1:"NA"
      print genename"\t"geneid"\t"biotype1"\t"biotype2"\t"chr"\t"start"\t"end"\t"strand"\t"flen"\t"SRC
    }
  ' "$INPUT" >> "$ANN_OUT"

  # ------------------ Indices (silent) ------------------
  [[ -f "${INDEX_PREFIX}.1.ebwt" ]] || \
    bowtie-build "$FA_SIMPLE" "$INDEX_PREFIX" > /dev/null 2>&1

  [[ -f "${INDEX_PREFIX}.1.bt2" ]] || \
    bowtie2-build --threads 8 "$FA_SIMPLE" "$INDEX_PREFIX" > /dev/null 2>&1

  # ------------------ Counts ------------------
  full_n=$(grep -c "^>" "$INPUT")
  simple_n=$(grep -c "^>" "$FA_SIMPLE")
  saf_n=$(( $(wc -l < "$SAF_OUT") - 1 ))
  ann_n=$(( $(wc -l < "$ANN_OUT") - 1 ))

  echo "📊 Entry counts:"
  echo "   Input FASTA:      $full_n"
  echo "   Simplified FASTA: $simple_n"
  echo "   SAF entries:      $saf_n"
  echo "   Annotation rows:  $ann_n"
  echo "✅ Completed: $BASENAME"

  # ------------------ Manifest: biotype counts ------------------
  total_entries="$ann_n"

  awk -F'\t' -v NAME="$BASENAME" -v FOLDER="$REF_FOLDER" -v SRC="$SOURCE" -v TOTAL="$total_entries" '
    NR > 1 {
      counts[$3]++
    }
    END {
      for (b in counts) {
        print NAME","FOLDER","SRC","TOTAL","b","counts[b]
      }
    }
  ' "$ANN_OUT" >> "$MANIFEST"

done
