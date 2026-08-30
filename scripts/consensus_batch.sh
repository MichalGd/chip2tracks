#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

output_root="${OUTPUT_DIR}/05_peaks/consensus"
status_dir="$output_root/.cohort_status"
mkdir -p "$output_root" "$status_dir"
rm -f -- "$status_dir"/*.tsv 2>/dev/null || true

consensus_worker() {
    local cohort="$1"
    local -a args=(
        --sample-manifest "$SAMPLE_MANIFEST"
        --cohort-manifest "$COHORT_MANIFEST"
        --output-root "$output_root"
        --minimum-support "$CONSENSUS_MIN_BIOLOGICAL_SAMPLES"
        --cohort-id "$cohort"
        --status-file "$status_dir/${cohort}.tsv"
    )
    is_true "$ALLOW_SINGLE_SAMPLE_CONSENSUS" && args+=(--allow-single)
    is_true "$REQUIRE_ALL_ENABLED_TRACKS" && args+=(--require-all)
    run_logged python3 "${SCRIPT_DIR}/build_consensus.py" "${args[@]}"
}

parallel_pool_init "$MERGE_PARALLEL_JOBS"
while IFS=$'\t' read -r cohort rest; do
    [[ "$cohort" == "cohort_id" ]] && continue
    parallel_pool_submit "consensus:$cohort" consensus_worker "$cohort"
done < "$COHORT_MANIFEST"
pool_failed=false
if ! parallel_pool_wait_all; then pool_failed=true; fi

status="$output_root/consensus_status.tsv"
printf 'cohort_id\tstatus\ttotal_samples\tsuccessful_peak_samples\texcluded_samples\tregions\treason\n' > "$status"
for status_file in "$status_dir"/*.tsv; do
    [[ -s "$status_file" ]] || continue
    tail -n +2 "$status_file" >> "$status"
done
if is_true "$pool_failed"; then
    exit 1
fi
exit 0
