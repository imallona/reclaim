#!/usr/bin/env bash
# Pull derived single-cell outputs from vm into barbara for integrative analysis.
# Run on barbara. Pass --dry-run to preview.
#
# Excludes raw FASTQs, BAMs, and kallisto bus intermediates.

set -euo pipefail

SRC_HOST=vm
SRC_ROOT=/mnt/src/repeats/results
DEST_ROOT="$HOME/src/repeats/results"

DRY_RUN=()
if [[ "${1-}" == "--dry-run" ]]; then
    DRY_RUN=(--dry-run)
fi

mkdir -p "$DEST_ROOT"

RSYNC_OPTS=(-avh --partial --prune-empty-dirs "${DRY_RUN[@]}")

FILTERS=(
    --filter='+ */'
    --filter='+ Solo.out/***'
    --filter='+ Log.final.out'
    --filter='+ counts.mtx'
    --filter='+ counts.barcodes.txt'
    --filter='+ counts.genes.txt'
    --filter='+ run_info.json'
    --filter='+ counts/***'
    --filter='+ report_rds/***'
    --filter='+ ruv_rds/***'
    --filter='+ evaluation/***'
    --filter='+ metadata/***'
    --filter='+ logs/***'
    --filter='+ benchmarks/***'
    --filter='+ *.html'
    --filter='- data/'
    --filter='- tmp/'
    --filter='- *.bam'
    --filter='- *.bam.bai'
    --filter='- *.bus'
    --filter='- matrix.ec'
    --filter='- transcripts.txt'
    --filter='- _STARgenome/'
    --filter='- *'
)

sync_tree() {
    local sub="$1"
    echo "=== syncing ${sub} ==="
    rsync "${RSYNC_OPTS[@]}" "${FILTERS[@]}" \
        "${SRC_HOST}:${SRC_ROOT}/${sub}/" \
        "${DEST_ROOT}/${sub}/"
}

sync_tree gse230647_sc
sync_tree paper/tdp43

echo "Done. Destination: ${DEST_ROOT}"
