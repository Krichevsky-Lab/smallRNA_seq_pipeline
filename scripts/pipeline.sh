#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

################################################################################
# pipeline.sh — Sequential controller for human small RNA-seq pipeline
#
# DESCRIPTION:
#   Orchestrates execution of the AD_pipeline from reference building through
#   final count matrix merging. Steps are executed sequentially unless a subset
#   is specified. This script dispatches step-specific scripts located in scripts/.
#
#   Each processing step is implemented as a standalone script and can be run
#   independently, provided its required inputs and references are available.
#
# INPUTS:
#   - scripts/runlist.txt                  : List of SRA accessions
#   - reference/                           : Reference FASTA files provided
#                                           (e.g., <class>_full.fa); indices and
#                                           auxiliary files are built by the pipeline
#   - raw FASTQ files (downloaded or provided locally)
#
# OUTPUTS:
#   - trimmed/                             : Adapter-trimmed FASTQ files
#   - aligned/<class>_<aligner>/           : BAM files per RNA class and aligner
#   - unmapped/<class>_<aligner>/          : Unmapped FASTQ files passed to next step
#   - counts/<class>_<aligner>/            : Per-sample and merged count matrices
#   - log/<class>_<aligner>/               : Per-sample logs and read summaries
#   - qc/<class>_<aligner>/                : FastQC and MultiQC reports (when --qc)
#
# QC:
#   - When invoked with `--qc`, supported steps run FastQC on intermediate FASTQ
#     files and generate per-step MultiQC reports.
#   - QC outputs are written to `qc/<step_name>/`.
#
# tRF CLASSIFICATION:
#   - Steps 5 (rRNA rescue) and 6 (tRNA alignment) invoke `classify_tRF.R`
#     to classify tRNA-derived fragments (tRFs) after alignment and counting.
#
# SUPPORTED STEPS:
#   0   Build references
#   1   Download raw data
#   2   Adapter trimming (Cutadapt)
#   3   UniVec filtering
#   4   rRNA alignment
#   5   rRNA rescue (Bowtie2) + tRF classification (classify_tRF.R)
#   6   tRNA alignment            + tRF classification (classify_tRF.R)
#   7   tRNA rescue (Bowtie2)
#   8   smallRNA alignment
#   9   smallRNA rescue (Bowtie2)
#   10  piRNA alignment
#   11  piRNA rescue (Bowtie2)
#   12  longRNA alignment
#   13  longRNA rescue (Bowtie2)
#   14  Merge and annotate count matrices (R)
#
# USAGE:
#   bash scripts/pipeline.sh [--qc] [--steps N-M | n1,n2,n3]
################################################################################

### --------------------
### Defaults
### --------------------
DO_QC=false
STEP_SELECTION=""

### --------------------
### Parse arguments
### --------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --qc)
            DO_QC=true
            ;;
        --steps)
            STEP_SELECTION="$2"
            shift
            ;;
        *)
            echo "❌ Unknown argument: $1"
            exit 1
            ;;
    esac
    shift
done

### --------------------
### Define step scripts
### --------------------
declare -A STEP_SCRIPTS=(
    [0]="00_make_reference.sh"
    [1]="01_download.sh"
    [2]="02_trim.sh"
    [3]="03_align_univec.sh"
    [4]="04_align_rrna.sh"
    [5]="05_align_rrna_rescue.sh"
    [6]="06_align_trna.sh"
    [7]="07_align_trna_rescue.sh"
    [8]="08_align_smallrna.sh"
    [9]="09_align_smallrna_rescue.sh"
    [10]="10_align_pirna.sh"
    [11]="11_align_pirna_rescue.sh"
    [12]="12_align_longrna.sh"
    [13]="13_align_longrna_rescue.sh"
    [14]="14_merge.R"
)

ALL_STEPS=({0..14})

### --------------------
### Step selection parser
### --------------------
parse_steps() {
    local sel="$1"
    local parsed=()

    if [[ -z "$sel" ]]; then
        echo "${ALL_STEPS[@]}"
        return
    fi

    if [[ "$sel" =~ ^([0-9]+)-([0-9]+)$ ]]; then
        for ((i=${BASH_REMATCH[1]}; i<=${BASH_REMATCH[2]}; i++)); do
            parsed+=("$i")
        done
        echo "${parsed[@]}"
        return
    fi

    if [[ "$sel" =~ , ]]; then
        IFS=',' read -ra parsed <<< "$sel"
        echo "${parsed[@]}"
        return
    fi

    echo "$sel"
}

STEPS_TO_RUN=($(parse_steps "$STEP_SELECTION"))

### --------------------
### Resolve paths
### --------------------
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${PIPELINE_DIR}/.." && pwd)"
cd "$ROOT_DIR"

### --------------------
### Runtime logging
### --------------------
RUNTIME_DIR="log"
RUNTIME_CSV="${RUNTIME_DIR}/pipeline_runtime.csv"
mkdir -p "$RUNTIME_DIR"

if [[ ! -f "$RUNTIME_CSV" ]]; then
    echo "timestamp,step,script,elapsed_seconds,status" > "$RUNTIME_CSV"
fi

### --------------------
### Run pipeline
### --------------------
echo "============================================================"
echo " RUNNING AD PIPELINE"
echo " STEPS: ${STEPS_TO_RUN[*]}"
echo " QC MODE: $DO_QC"
echo "============================================================"

for step in "${STEPS_TO_RUN[@]}"; do
    script="${STEP_SCRIPTS[$step]}"

    [[ -z "$script" ]] && { echo "❌ Unknown step: $step"; exit 1; }

    echo
    echo "------------------------------------------------------------"
    echo "▶ STEP $step: $script"
    echo "------------------------------------------------------------"

    start_ts=$(date +%s)
    start_human=$(date +"%Y-%m-%d %H:%M:%S")
    status="OK"

    if [[ "$script" == *.R ]]; then
        Rscript "$PIPELINE_DIR/$script" || status="FAIL"
    else
        if $DO_QC; then
            bash "$PIPELINE_DIR/$script" --qc || status="FAIL"
        else
            bash "$PIPELINE_DIR/$script" || status="FAIL"
        fi
    fi

    end_ts=$(date +%s)
    elapsed=$((end_ts - start_ts))

    echo "${start_human},${step},${script},${elapsed},${status}" >> "$RUNTIME_CSV"

    [[ "$status" == "FAIL" ]] && {
        echo "❌ Step $step failed — stopping pipeline"
        exit 1
    }

    echo "✔ Completed Step $step (${elapsed}s)"
done

echo
echo "============================================================"
echo " 🎉 PIPELINE COMPLETE"
echo " Runtime log: ${RUNTIME_CSV}"
echo "============================================================"
