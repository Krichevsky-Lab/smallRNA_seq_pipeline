#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

################################################################################
# pipeline.sh — Sequential controller for small RNA-seq pipeline
#
# DESCRIPTION:
#   Orchestrates execution of the pipeline from raw data download through
#   final count matrix merging. Steps are executed sequentially unless a subset
#   is specified. This script does not perform data processing itself; it
#   dispatches step-specific scripts located in the alignment_scripts/ directory.
#
# SUPPORTED STEPS:
#   1   Download raw data
#   2   Adapter trimming (Cutadapt)
#   3   UniVec filtering
#   4   rRNA alignment (Bowtie1)
#   5   rRNA rescue (Bowtie2)
#   6   tRNA alignment (Bowtie1)
#   7   tRNA rescue (Bowtie2)
#   8   smallRNA alignment
#   9   smallRNA rescue (Bowtie2)
#   10  piRNA alignment (Bowtie1)
#   11  piRNA rescue (Bowtie2)
#   12  longRNA alignment (Bowtie1)
#   13  longRNA rescue (Bowtie2)
#   14  tRF coverage (per-interval tRNA fragment extraction)
#   15  Merge and annotate count matrices (R)
#
# USAGE:
#   bash alignment_scripts/pipeline.sh <mouse|human> [OPTIONS]
#
# EXAMPLES:
#   bash alignment_scripts/pipeline.sh mouse
#   bash alignment_scripts/pipeline.sh human --qc
#   bash alignment_scripts/pipeline.sh mouse --steps 4-6
#   bash alignment_scripts/pipeline.sh mouse --steps 2,5,9
#   bash alignment_scripts/pipeline.sh mouse --status
#   bash alignment_scripts/pipeline.sh mouse --force --steps 6-7
#
# OPTIONS:
#   mouse|human  Required. Species to process (positional first argument)
#   --qc         Enable per-step QC (FastQC + MultiQC) where supported
#   --steps      Subset of steps to run (range N-M or comma-separated list)
#   --status     Print done/pending status for each step, then exit
#   --force      Ignore completion stamps; re-run steps even if already done
#
# OUTPUTS:
#   - alignment_output/log/pipeline_runtime.csv        (updated after each completed step)
#   - alignment_output/log/checkpoints/<species>_step_<N>.done  (stamp files)
#
# NOTE:
#   SPECIES is exported as an environment variable so all dispatched step
#   scripts can read it directly without needing to parse it themselves.
#
################################################################################

SPECIES="${1:-}"
shift || true

DO_QC=false
STEP_SELECTION=""
FORCE=false
SHOW_STATUS=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --qc)      DO_QC=true ;;
        --steps)   STEP_SELECTION="$2"; shift ;;
        --force)   FORCE=true ;;
        --status)  SHOW_STATUS=true ;;
        *)         echo "[FAIL]  Unknown argument: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$SPECIES" ]]; then
    echo "[FAIL]  Species is required as first argument (mouse or human)"
    echo "Usage: bash alignment_scripts/pipeline.sh <mouse|human> [--qc] [--steps N-M] [--status] [--force]"
    exit 1
fi

if [[ "$SPECIES" != "mouse" && "$SPECIES" != "human" ]]; then
    echo "[FAIL]  Unknown species '$SPECIES'. Expected 'mouse' or 'human'."
    exit 1
fi

export SPECIES

### --------------------
### Step definitions
### --------------------
declare -A STEP_SCRIPTS=(
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
    [14]="14_tRF_coverage.sh"
    [15]="15_merge.R"
)

declare -A STEP_NAMES=(
    [1]="Download raw data"
    [2]="Adapter trimming"
    [3]="UniVec filtering"
    [4]="rRNA alignment (Bowtie1)"
    [5]="rRNA rescue (Bowtie2)"
    [6]="tRNA alignment (Bowtie1)"
    [7]="tRNA rescue (Bowtie2)"
    [8]="smallRNA alignment (Bowtie1)"
    [9]="smallRNA rescue (Bowtie2)"
    [10]="piRNA alignment (Bowtie1)"
    [11]="piRNA rescue (Bowtie2)"
    [12]="longRNA alignment (Bowtie1)"
    [13]="longRNA rescue (Bowtie2)"
    [14]="tRF coverage"
    [15]="Merge and annotate count matrices"
)

ALL_STEPS=({1..15})

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
        IFS=',' read -ra parts <<< "$sel"
        echo "${parts[@]}"
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
### Checkpoint directory
### --------------------
CHECKPOINT_DIR="alignment_output/log/checkpoints"
mkdir -p "$CHECKPOINT_DIR"

stamp_path() { echo "${CHECKPOINT_DIR}/${SPECIES}_step_${1}.done"; }

step_is_done() {
    local stamp
    stamp=$(stamp_path "$1")
    [[ -f "$stamp" ]]
}

mark_done() {
    local stamp
    stamp=$(stamp_path "$1")
    date '+%Y-%m-%d %H:%M:%S' > "$stamp"
}

### --------------------
### --status mode
### --------------------
if [[ "$SHOW_STATUS" == true ]]; then
    echo ""
    echo "  Pipeline status — species: $SPECIES"
    echo "  -----------------------------------------------"
    printf "  %-5s  %-8s  %s\n" "STEP" "STATUS" "NAME"
    echo "  -----------------------------------------------"
    for step in "${ALL_STEPS[@]}"; do
        if step_is_done "$step"; then
            done_ts=$(cat "$(stamp_path "$step")")
            printf "  %-5s  %-8s  %s  [%s]\n" "$step" "DONE" "${STEP_NAMES[$step]}" "$done_ts"
        else
            printf "  %-5s  %-8s  %s\n" "$step" "pending" "${STEP_NAMES[$step]}"
        fi
    done
    echo "  -----------------------------------------------"
    echo ""
    exit 0
fi

### --------------------
### Runtime logging setup
### --------------------
RUNTIME_DIR="alignment_output/log"
RUNTIME_CSV="${RUNTIME_DIR}/pipeline_runtime.csv"

mkdir -p "$RUNTIME_DIR"

if [[ ! -f "$RUNTIME_CSV" ]]; then
    echo "timestamp,species,step,script,elapsed_seconds,status" > "$RUNTIME_CSV"
fi

### --------------------
### Run steps
### --------------------
echo ""
echo "============================================================"
echo "  PIPELINE START"
echo "  Species: $SPECIES"
echo "  Steps:   ${STEPS_TO_RUN[*]}"
echo "  QC:      $DO_QC"
echo "  Force:   $FORCE"
echo "============================================================"

for step in "${STEPS_TO_RUN[@]}"; do
    script="${STEP_SCRIPTS[$step]:-}"

    if [[ -z "$script" ]]; then
        echo "[FAIL]  Unknown step: $step"
        exit 1
    fi

    if [[ "$FORCE" == false ]] && step_is_done "$step"; then
        done_ts=$(cat "$(stamp_path "$step")")
        echo ""
        echo "  Step $step — ${STEP_NAMES[$step]}"
        echo "[SKIP]  Already completed at ${done_ts} — use --force to re-run"
        continue
    fi

    start_ts=$(date +%s)
    start_time=$(date +"%Y-%m-%d %H:%M:%S")
    status="OK"

    if [[ "$script" == *.R ]]; then
        if ! Rscript "$PIPELINE_DIR/$script" "$SPECIES"; then
            status="FAIL"
        fi
    else
        args=(--species "$SPECIES")
        [[ "$DO_QC" == true ]] && args+=(--qc)
        if ! bash "$PIPELINE_DIR/$script" "${args[@]}"; then
            status="FAIL"
        fi
    fi

    end_ts=$(date +%s)
    elapsed=$((end_ts - start_ts))

    echo "${start_time},${SPECIES},${step},${script},${elapsed},${status}" >> "$RUNTIME_CSV"

    if [[ "$status" == "FAIL" ]]; then
        echo "[FAIL]  Step $step failed after ${elapsed}s — stopping pipeline"
        exit 1
    fi

    mark_done "$step"
    echo "[DONE]  Step $step completed in ${elapsed}s"
done

echo ""
echo "============================================================"
echo "  PIPELINE COMPLETE — runtime log: ${RUNTIME_CSV}"
echo "============================================================"
