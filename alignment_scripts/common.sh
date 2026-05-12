#!/usr/bin/env bash
# common.sh — shared settings and functions, sourced by all step scripts

THREADS=10
MIN_READ_LEN=15
PHRED_QUAL=20

log_info() { echo "[INFO]  $*"; }
log_skip() { echo "[SKIP]  $*"; }
log_done() { echo "[DONE]  $*"; }
log_warn() { echo "[WARN]  $*"; }
log_fail() { echo "[FAIL]  $*" >&2; }

step_banner() {
    local step="$1" name="$2"
    local ts; ts=$(date '+%Y-%m-%d %H:%M')
    echo ""
    echo "------------------------------------------------------------"
    printf "  Step %-2s | %-6s | %s\n" "$step" "${SPECIES:-unknown}" "$ts"
    echo "  $name"
    echo "------------------------------------------------------------"
}

# Exit with error if a required file is missing
require_file() {
    local f="$1"
    if [[ ! -f "$f" ]]; then
        log_fail "Required file not found: $f"
        exit 1
    fi
}

# Exit with error if a Bowtie1 or Bowtie2 index is missing
require_index() {
    local idx="$1"
    if [[ ! -f "${idx}.1.ebwt" && ! -f "${idx}.1.bt2" && ! -f "${idx}.1.bt2l" ]]; then
        log_fail "Alignment index not found: ${idx}.1.ebwt (or .bt2)"
        log_fail "Run 00_make_reference.sh first to build the index."
        exit 1
    fi
}

# BAM: samtools quickcheck verifies header + EOF block (catches truncation)
check_bam() {
    local f="$1"
    [[ -f "$f" ]] || return 1
    if samtools quickcheck "$f" 2>/dev/null; then
        return 0
    fi
    log_warn "Truncated BAM removed: $f"
    rm -f "$f" "${f}.bai"
    return 1
}

# FASTQ.gz: checks file is non-empty and first record is readable
check_gz() {
    local f="$1"
    [[ -s "$f" ]] || return 1
    if zcat "$f" 2>/dev/null | head -1 | grep -q '^[@>]'; then
        return 0
    fi
    log_warn "Truncated/unreadable gz removed: $f"
    rm -f "$f"
    return 1
}

# Tabular (TSV/CSV/TXT): must have header + at least one data row
check_tabular() {
    local f="$1"
    [[ -f "$f" ]] || return 1
    if [[ $(wc -l < "$f") -ge 2 ]]; then
        return 0
    fi
    log_warn "Empty or header-only file removed: $f"
    rm -f "$f"
    return 1
}
