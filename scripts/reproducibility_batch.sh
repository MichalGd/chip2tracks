#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
require_config

# Build the ATACseq2tracks-style biological-support consensus for every peak class.
bash "${SCRIPT_DIR}/consensus_batch.sh"

root="${OUTPUT_DIR}/05_peaks/reproducibility"
mkdir -p "$root" "${OUTPUT_DIR}/logs/reproducibility"
status="$root/status.tsv"
printf 'cohort_id\tcomparison\tmethod\tstatus\tpassing_regions\toutput\treason\n' > "$status"

if ! is_true "$RUN_IDR"; then
    printf '.\t.\tidr\tSKIPPED\t0\t.\tRUN_IDR=false\n' >> "$status"
    exit 0
fi

while IFS=$'\t' read -r cohort_id cohort_key genome assay factor antibody layout target_class duplicate_policy primary_caller primary_class biological_samples sample_keys conditions; do
    [[ "$cohort_id" == "cohort_id" ]] && continue
    if [[ "$primary_caller" != "macs3" || "$primary_class" != "narrow" ]]; then
        printf '%s\t.\tidr\tNOT_APPLICABLE\t0\t.\tbroad_or_non_macs3_primary\n' "$cohort_id" >> "$status"
        continue
    fi
    IFS=',' read -r -a samples <<< "$sample_keys"
    if (( ${#samples[@]} < 2 )); then
        printf '%s\t.\tidr\tSKIPPED\t0\t.\tfewer_than_two_biological_replicates\n' "$cohort_id" >> "$status"
        continue
    fi
    for ((i=0; i<${#samples[@]}-1; i++)); do
        for ((j=i+1; j<${#samples[@]}; j++)); do
            left="${samples[$i]}"; right="${samples[$j]}"
            comparison="${left}__vs__${right}"
            pair_dir="${root}/${cohort_id}"
            raw="${pair_dir}/${comparison}.idr.narrowPeak"
            passing="${pair_dir}/${comparison}.idr${IDR_SOFT_THRESHOLD}.bed"
            log="${OUTPUT_DIR}/logs/reproducibility/${comparison}.idr.log"
            left_peaks="${OUTPUT_DIR}/05_peaks/per_sample/${left}/macs3/${left}.macs3.idr_peaks.narrowPeak"
            right_peaks="${OUTPUT_DIR}/05_peaks/per_sample/${right}/macs3/${right}.macs3.idr_peaks.narrowPeak"
            mkdir -p "$pair_dir"
            if [[ ! -s "$left_peaks" || ! -s "$right_peaks" ]]; then
                printf '%s\t%s\tidr\tERROR\t0\t.\tmissing_ranked_narrowPeak_input\n' "$cohort_id" "$comparison" >> "$status"
                continue
            fi
            if run_logged "$IDR_COMMAND" --samples "$left_peaks" "$right_peaks" \
                --input-file-type narrowPeak --rank "$IDR_RANK" \
                --soft-idr-threshold "$IDR_SOFT_THRESHOLD" --output-file "$raw" >"$log" 2>&1; then
                cutoff="$(awk -v p="$IDR_SOFT_THRESHOLD" 'BEGIN{print -log(p)/log(10)}')"
                awk -v cutoff="$cutoff" 'BEGIN{OFS="\t"} NF>=12 && $12>=cutoff {print $1,$2,$3}' \
                    "$raw" > "$passing"
                count="$(awk 'NF>=3 {n++} END{print n+0}' "$passing")"
                printf '%s\t%s\tidr\tSUCCESS\t%s\t%s\t.\n' "$cohort_id" "$comparison" "$count" "$passing" >> "$status"
            else
                printf '%s\t%s\tidr\tERROR\t0\t.\tidr_command_failed\n' "$cohort_id" "$comparison" >> "$status"
            fi
        done
    done
done < "$COHORT_MANIFEST"

if awk -F '\t' 'NR>1 && $4=="ERROR" {found=1} END{exit !found}' "$status"; then
    warn "one or more IDR comparisons failed; see $status"
fi
