#!/usr/bin/env bash
# Push lightweight downstream-ready artefacts from barbara to vm.
# Run on barbara. Pass --dry-run to preview.
#
# Ships count matrices, report RDS, metadata, logs, benchmarks, HTML reports,
# and simulation ground truth + evaluation.
# Skips raw FASTQs, BAMs, kallisto bus/h5, salmon aux, aligner indices,
# and the entire shared/ tree (vm has its own).

set -euo pipefail

DST_HOST=vm
SRC_ROOT="$HOME/src/repeats/results"
DST_ROOT=/mnt/src/repeats/results

DRY_RUN=()
if [[ "${1-}" == "--dry-run" ]]; then
    DRY_RUN=(--dry-run)
fi

ssh "$DST_HOST" "mkdir -p ${DST_ROOT}"

RSYNC_OPTS=(-avh --partial --prune-empty-dirs "${DRY_RUN[@]}")

# Ordered include/exclude filters.
# Trailing "/***" means "this directory and everything inside it".
BULK_FILTERS=(
    --filter='+ */'
    --filter='+ counts/***'
    --filter='+ report_rds/***'
    --filter='+ metadata/***'
    --filter='+ logs/***'
    --filter='+ benchmarks/***'
    --filter='+ *.html'
    --filter='+ star/bulk/*/Log.final.out'
    --filter='+ star/bulk/*/SJ.out.tab'
    --filter='+ star/bulk/*/*_counts.tsv'
    --filter='+ star/bulk/*/*_counts.tsv.summary'
    --filter='+ kallisto/*/bulk/*/abundance.tsv'
    --filter='+ kallisto/*/bulk/*/run_info.json'
    --filter='+ alevin/*/bulk/*/quant.sf'
    --filter='+ alevin/*/bulk/*/cmd_info.json'
    --filter='+ alevin/*/bulk/*/lib_format_counts.json'
    --filter='- data/'
    --filter='- tmp/'
    --filter='- *.bam'
    --filter='- *.bam.bai'
    --filter='- *.h5'
    --filter='- output.bus'
    --filter='- output.sorted.bus'
    --filter='- matrix.ec'
    --filter='- transcripts.txt'
    --filter='- aux_info/'
    --filter='- libParams/'
    --filter='- *'
)

SIM_FILTERS=(
    --filter='+ */'
    --filter='+ counts/***'
    --filter='+ evaluation/***'
    --filter='+ simulations/***'
    --filter='+ metadata/***'
    --filter='+ logs/***'
    --filter='+ benchmarks/***'
    --filter='+ *.html'
    --filter='- tmp/'
    --filter='- *.bam'
    --filter='- *.bam.bai'
    --filter='- *'
)

push_tree() {
    local sub="$1"
    shift
    local filters=("$@")
    echo "=== syncing ${sub} ==="
    rsync "${RSYNC_OPTS[@]}" "${filters[@]}" \
        "${SRC_ROOT}/${sub}/" \
        "${DST_HOST}:${DST_ROOT}/${sub}/"
}

push_tree gse126543_bulk "${BULK_FILTERS[@]}"
push_tree gse230647_bulk "${BULK_FILTERS[@]}"

for sim in simulation_smartseq2 simulation_chromium \
           simulation_smartseq2_noise_0pct simulation_smartseq2_noise_1pct \
           simulation_smartseq2_noise_5pct simulation_smartseq2_noise_10pct \
           simulation_chromium_noise_0pct simulation_chromium_noise_1pct \
           simulation_chromium_noise_5pct simulation_chromium_noise_10pct; do
    if ssh "$DST_HOST" "true" && [[ -d "${SRC_ROOT}/${sim}" ]]; then
        push_tree "${sim}" "${SIM_FILTERS[@]}"
    fi
done

echo "Done. Destination: ${DST_HOST}:${DST_ROOT}"
