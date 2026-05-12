#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

################################################################################
# 00_make_reference.sh
#
# Split <species>_reference.fa (single combined file) by RNA class and build
# per-class subfolders inside reference/<species>/ with indices and annotation.
#
# Header format: >genename|accession|biotype|coordinates
# Entries with unknown biotype (UniVec) go to univec/ → Bowtie1 only.
# All other classes → simplified FASTA, SAF, annotation TSV, Bowtie1+2 indices.
#
# Output structure:
#   reference/<species>/
#     <species>_reference.fa       ← input (must exist before running)
#     pirna/      <species>_piRNA.fa  .saf  _annotation.tsv  + bowtie indices
#     longrna/    <species>_longRNA.*
#     smallrna/   <species>_smallRNA.*
#     trna/       <species>_tRNA.*
#     rrna/       <species>_rRNA.*
#     univec/     <species>_univec.*  (Bowtie1 only)
#     reference_manifest.csv
#
# Manifest (CSV):
#   species, reference_name, folder, source, total_entries, biotype_1, biotype_1_count
#
# Usage:
#   bash alignment_scripts/00_make_reference.sh mouse
#   bash alignment_scripts/00_make_reference.sh human
################################################################################

# ------------------ species argument ------------------
if [[ $# -ne 1 ]]; then
    echo "Usage: bash alignment_scripts/00_make_reference.sh <species>"
    echo "  e.g. bash alignment_scripts/00_make_reference.sh mouse"
    echo "       bash alignment_scripts/00_make_reference.sh human"
    exit 1
fi

SPECIES="$1"
if [[ "$SPECIES" != "mouse" && "$SPECIES" != "human" ]]; then
    echo "[FAIL]  Unknown species '$SPECIES'. Expected 'mouse' or 'human'."
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REFERENCE_DIR="$(dirname "$SCRIPT_DIR")/reference/$SPECIES"
INPUT="$REFERENCE_DIR/${SPECIES}_reference.fa"
MANIFEST="$REFERENCE_DIR/reference_manifest.csv"

[[ -f "$INPUT" ]] || { echo "[FAIL]  Not found: $INPUT"; exit 1; }

echo "[INFO]  Species: $SPECIES"
echo "[INFO]  Reference dir: $REFERENCE_DIR"

SPLIT_DIR="$(mktemp -d)"
trap 'rm -rf "$SPLIT_DIR"' EXIT

# Folder name → output file prefix (species-prefixed)
declare -A FOLDER_PREFIX=(
    [pirna]="${SPECIES}_piRNA"
    [longrna]="${SPECIES}_longRNA"
    [smallrna]="${SPECIES}_smallRNA"
    [trna]="${SPECIES}_tRNA"
    [rrna]="${SPECIES}_rRNA"
    [univec]="${SPECIES}_univec"
)

# Processing order (univec last)
FOLDERS=(pirna longrna smallrna trna rrna univec)

# ------------------ init manifest ------------------
echo "species,reference_name,folder,source,total_entries,biotype_1,biotype_1_count" > "$MANIFEST"

# ==============================================================================
# Step 1: Split input into per-folder temp FASTA files
# ==============================================================================
echo "[INFO]  Splitting $INPUT by RNA class..."

for folder in "${FOLDERS[@]}"; do
    mkdir -p "$REFERENCE_DIR/$folder"
    > "$SPLIT_DIR/${folder}.fa"
done

awk -F"|" -v SDIR="$SPLIT_DIR" '
BEGIN {
    bf["piRNA"]          = "pirna"
    bf["lncRNA"]         = "longrna"
    bf["protein_coding"] = "longrna"
    bf["mRNA_fragment"]  = "longrna"
    bf["miRNA"]          = "smallrna"
    bf["miRNA_hairpin"]  = "longrna"
    bf["miRNA_isomiR"]   = "smallrna"
    bf["snoRNA"]         = "smallrna"
    bf["snRNA"]          = "smallrna"
    bf["scaRNA"]         = "smallrna"
    bf["sRNA"]           = "smallrna"
    bf["misc_RNA"]       = "smallrna"
    bf["RNase_P_RNA"]    = "smallrna"
    bf["tRNA"]           = "trna"
    bf["Mt_tRNA"]        = "trna"
    bf["rRNA"]           = "rrna"
    bf["Mt_rRNA"]        = "rrna"
    cur = ""
}
/^>/ {
    folder = (NF == 4) ? (($3 in bf) ? bf[$3] : "univec") : "univec"
    cur    = SDIR "/" folder ".fa"
}
cur != "" { print >> cur }
' "$INPUT"

echo "[DONE]  Split complete"

# ==============================================================================
# Step 2: Process each folder
# ==============================================================================
for folder in "${FOLDERS[@]}"; do
    prefix="${FOLDER_PREFIX[$folder]}"
    OUTDIR="$REFERENCE_DIR/$folder"
    FULL_FA="$SPLIT_DIR/${folder}.fa"
    SOURCE="${folder}_alignment"
    INDEX_PREFIX="$OUTDIR/$prefix"

    [[ -s "$FULL_FA" ]] || { echo "[WARN]  Skipping $folder (empty)"; continue; }

    echo "[INFO]  Making reference for $folder"

    # --------------------------------------------------------------------------
    # UniVec: Bowtie1 only
    # --------------------------------------------------------------------------
    if [[ "$folder" == "univec" ]]; then
        [[ -f "${INDEX_PREFIX}.1.ebwt" ]] || \
            bowtie-build "$FULL_FA" "$INDEX_PREFIX" > /dev/null 2>&1

        total_n=$(grep -c "^>" "$FULL_FA")
        echo "$SPECIES,$prefix,$folder,univec_alignment,$total_n,NA,NA" >> "$MANIFEST"
        echo "[DONE]  Completed: $prefix"
        continue
    fi

    FA_SIMPLE="${OUTDIR}/${prefix}.fa"
    SAF_OUT="${OUTDIR}/${prefix}.saf"
    ANN_OUT="${OUTDIR}/${prefix}_annotation.tsv"

    # ------------------ Simplified FASTA ------------------
    awk -F"|" '
      /^>/ { gsub(/^>/,"",$1); print ">"$1; next }
      { print }
    ' "$FULL_FA" > "$FA_SIMPLE"

    dups=$(grep "^>" "$FA_SIMPLE" | sort | uniq -d)
    if [[ -n "$dups" ]]; then
        echo "[FAIL]  Duplicate genenames in $FA_SIMPLE:"
        echo "$dups"
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
      END { if (seq != "") print name"\t"name"\t1\t"length(seq)"\t*" }
    ' "$FA_SIMPLE" >> "$SAF_OUT"

    # ------------------ Annotation ------------------
    echo -e "genename\tgeneid\tbiotype_1\tbiotype_2\tchromosome\tstart\tend\tstrand\tlength\tsource" \
        > "$ANN_OUT"

    awk -F"|" -v SRC="$SOURCE" '
      /^>/ {
        gsub(/^>/,"",$1)
        genename=$1; rawid=$2; biotype1=$3; biotype2=$3; header=$4
        geneid=(rawid ~ /^ENS|^NM_|^NR_|^XM_|^XR_/) ? rawid : genename

        if (match(header,/^([^:]+):([0-9]+)-([0-9]+)\(([+-])\)$/,m)) {
            chr=m[1]; start=m[2]; end=m[3]; strand=m[4]
        } else if (match(header,/^([^:]+):([0-9]+)-([0-9]+)$/,m)) {
            chr=m[1]; start=m[2]; end=m[3]; strand="*"
        } else {
            chr=header; start="NA"; end="NA"; strand="*"
        }

        flen=(start~/^[0-9]+$/ && end~/^[0-9]+$/) ? end-start+1 : "NA"
        print genename"\t"geneid"\t"biotype1"\t"biotype2"\t"chr"\t"start"\t"end"\t"strand"\t"flen"\t"SRC
      }
    ' "$FULL_FA" >> "$ANN_OUT"

    # ------------------ Indices ------------------
    [[ -f "${INDEX_PREFIX}.1.ebwt" ]] || \
        bowtie-build "$FA_SIMPLE" "$INDEX_PREFIX" > /dev/null 2>&1

    [[ -f "${INDEX_PREFIX}.1.bt2" ]] || \
        bowtie2-build --threads 8 "$FA_SIMPLE" "$INDEX_PREFIX" > /dev/null 2>&1

    # ------------------ Counts ------------------
    full_n=$(grep -c "^>" "$FULL_FA")
    simple_n=$(grep -c "^>" "$FA_SIMPLE")
    saf_n=$(( $(wc -l < "$SAF_OUT") - 1 ))
    ann_n=$(( $(wc -l < "$ANN_OUT") - 1 ))

    echo "[INFO]  Entry counts: input=$full_n  simplified=$simple_n  SAF=$saf_n  annotation=$ann_n"
    echo "[DONE]  Completed: $prefix"

    # ------------------ Manifest ------------------
    awk -F'\t' -v SP="$SPECIES" -v NAME="$prefix" -v FOLDER="$folder" -v SRC="$SOURCE" -v TOTAL="$ann_n" '
      NR > 1 { counts[$3]++ }
      END { for (b in counts) print SP","NAME","FOLDER","SRC","TOTAL","b","counts[b] }
    ' "$ANN_OUT" >> "$MANIFEST"

done
