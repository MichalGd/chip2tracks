#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

root="${OUTPUT_DIR}/04_tracks/spikein"
control_root="${OUTPUT_DIR}/04_tracks/spikein_control"
qc_root="${OUTPUT_DIR}/06_qc/spikein"
count_parts="${qc_root}/count_parts"
rm -rf -- "$root" "$control_root" "$qc_root"
mkdir -p "$root" "$control_root" "$count_parts" \
    "${OUTPUT_DIR}/03_alignment/spikein/filtered" "${OUTPUT_DIR}/logs/spikein"
if [[ "$SPIKEIN_MODE" == "none" ]]; then
    printf '{"status":"SKIPPED","reason":"SPIKEIN_MODE=none"}\n' > "${root}/SKIPPED.json"
    exit 0
fi

process_spike_bam() {
    local sample_key="$1" layout="$2" cohort_id="$3" spike_stage="$4" spike_lot="$5" ratio="$6"
    local raw marked prefiltered filtered metrics mask host_count spike_count
    raw="${OUTPUT_DIR}/03_alignment/spikein/spike/${sample_key}.${SPIKEIN_MODE}.bam"
    marked="${OUTPUT_DIR}/03_alignment/spikein/filtered/${sample_key}.${SPIKEIN_MODE}.marked.bam"
    prefiltered="${OUTPUT_DIR}/03_alignment/spikein/filtered/${sample_key}.${SPIKEIN_MODE}.preblacklist.bam"
    filtered="${OUTPUT_DIR}/03_alignment/spikein/filtered/${sample_key}.${SPIKEIN_MODE}.analysis.bam"
    metrics="${OUTPUT_DIR}/03_alignment/metrics/${sample_key}.${SPIKEIN_MODE}.duplicate_metrics.txt"
    run_logged "$PICARD_COMMAND" MarkDuplicates I="$raw" O="$marked" M="$metrics" REMOVE_DUPLICATES=false \
        ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
        >"${OUTPUT_DIR}/logs/spikein/${sample_key}.picard.log" 2>&1
    if [[ "$layout" == "PE" ]]; then
        [[ "$SPIKEIN_DUPLICATE_POLICY" == "remove" ]] && mask=3852 || mask=2828
        samtools view -@ "$THREADS_SAMTOOLS" -b -q "$SPIKEIN_MIN_MAPQ" -f 2 -F "$mask" -o "$prefiltered" "$marked"
    else
        [[ "$SPIKEIN_DUPLICATE_POLICY" == "remove" ]] && mask=3844 || mask=2820
        samtools view -@ "$THREADS_SAMTOOLS" -b -q "$SPIKEIN_MIN_MAPQ" -F "$mask" -o "$prefiltered" "$marked"
    fi
    if [[ -n "$SPIKEIN_BLACKLIST" && "$SPIKEIN_BLACKLIST" != "." ]]; then
        bedtools intersect -v -abam "$prefiltered" -b "$SPIKEIN_BLACKLIST" > "$filtered"
        rm -f -- "$prefiltered"
    else
        mv "$prefiltered" "$filtered"
    fi
    samtools quickcheck "$filtered"
    samtools index -@ "$THREADS_SAMTOOLS" "$filtered"
    host_count="$(signal_count "$(analysis_bam_path "$sample_key")" "$layout")"
    spike_count="$(signal_count "$filtered" "$layout")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample_key" "$layout" "$cohort_id" \
        "$spike_stage" "$spike_lot" "$ratio" "$host_count" "$spike_count" \
        > "${count_parts}/${sample_key}.tsv"
}

parallel_pool_init "$SPIKEIN_PARALLEL_JOBS"
while IFS=$'\t' read -r \
    sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
    is_control control_type control_id control_key duplicate_policy blacklist ratio spike_stage spike_lot batch donor output_prefix \
    technical_units fastq_1_list fastq_2_list cohort_id cohort_key primary_caller primary_class; do
    [[ "$sample_key" == "sample_key" ]] && continue
    parallel_pool_submit "$sample_key" process_spike_bam "$sample_key" "$layout" "$cohort_id" \
        "$spike_stage" "$spike_lot" "$ratio"
done < "$SAMPLE_MANIFEST"
parallel_pool_wait_all

counts="${qc_root}/spikein_counts.tsv"
printf 'sample_key\tlayout\tcohort_id\tspikein_stage\tspikein_lot\tspikein_to_host_ratio\thost_observations\tspike_observations\n' > "$counts"
find "$count_parts" -maxdepth 1 -type f -name '*.tsv' -print0 | sort -z | xargs -0 -r cat >> "$counts"
rm -rf -- "$count_parts"

qc_args=(--scale-target "$SPIKEIN_SCALE_TARGET" --fail-below "$SPIKEIN_MIN_OBSERVATIONS_FAIL"
    --warn-below "$SPIKEIN_MIN_OBSERVATIONS_WARN" --warn-low-fraction "$SPIKEIN_WARN_LOW_FRACTION"
    --warn-high-fraction "$SPIKEIN_WARN_HIGH_FRACTION")
is_true "$ALLOW_FAILED_SPIKEIN" && qc_args+=(--allow-failed)
run_logged python3 "${SCRIPT_DIR}/spikein_qc.py" "$counts" "${qc_root}/spikein_scaling.tsv" "${qc_args[@]}"

write_spike_tracks() {
    local sample_key="$1" layout="$2" cohort_id="$3" host_scale="$4" spike_scale="$5"
    local genome chrom_sizes host_bam spike_bam host_dir host_tmp spike_tmp
    genome="$(awk -F '\t' -v key="$sample_key" 'NR>1 && $1==key {print $5; exit}' "$SAMPLE_MANIFEST")"
    chrom_sizes="$(reference_value CHROM_SIZES "$genome")"
    host_bam="$(analysis_bam_path "$sample_key")"
    spike_bam="${OUTPUT_DIR}/03_alignment/spikein/filtered/${sample_key}.${SPIKEIN_MODE}.analysis.bam"
    host_dir="${root}/${cohort_id}"
    mkdir -p "$host_dir"
    host_tmp="${host_dir}/${sample_key}.SpikeInScaled.${SPIKEIN_MODE}.bedGraph"
    set_bamcoverage_signal_args "$sample_key" "$layout"
    host_args=(--bam "$host_bam" --outFileName "$host_tmp" --outFileFormat bedgraph --binSize "$TRACK_BIN_SIZE"
        --scaleFactor "$host_scale" --numberOfProcessors "$THREADS_BAMCOVERAGE")
    host_args+=("${BAMCOVERAGE_SIGNAL_ARGS[@]}")
    run_logged bamCoverage "${host_args[@]}"
    bedtools sort -faidx "$chrom_sizes" -i "$host_tmp" > "${host_tmp}.sorted"
    is_true "$GENERATE_COVERAGE_BIGWIGS" && run_logged bedGraphToBigWig "${host_tmp}.sorted" "$chrom_sizes" \
        "${host_dir}/${sample_key}.SpikeInScaled.${SPIKEIN_MODE}.bw"
    if is_true "$GENERATE_COVERAGE_BEDGRAPHS"; then mv "${host_tmp}.sorted" "$host_tmp"; \
        else rm -f -- "$host_tmp" "${host_tmp}.sorted"; fi

    if is_true "$GENERATE_SPIKEIN_CONTROL_TRACKS"; then
        spike_tmp="${control_root}/${sample_key}.${SPIKEIN_MODE}.CPM.bedGraph"
        set_bamcoverage_signal_args "$sample_key" "$layout"
        spike_args=(--bam "$spike_bam" --outFileName "$spike_tmp" --outFileFormat bedgraph --binSize "$TRACK_BIN_SIZE"
            --scaleFactor "$spike_scale" --numberOfProcessors "$THREADS_BAMCOVERAGE")
        spike_args+=("${BAMCOVERAGE_SIGNAL_ARGS[@]}")
        run_logged bamCoverage "${spike_args[@]}"
        bedtools sort -faidx "$SPIKEIN_CHROM_SIZES" -i "$spike_tmp" > "${spike_tmp}.sorted"
        is_true "$GENERATE_COVERAGE_BIGWIGS" && run_logged bedGraphToBigWig "${spike_tmp}.sorted" \
            "$SPIKEIN_CHROM_SIZES" "${control_root}/${sample_key}.${SPIKEIN_MODE}.CPM.bw"
        if is_true "$GENERATE_COVERAGE_BEDGRAPHS"; then mv "${spike_tmp}.sorted" "$spike_tmp"; \
            else rm -f -- "$spike_tmp" "${spike_tmp}.sorted"; fi
    fi
}

parallel_pool_init "$SPIKEIN_PARALLEL_JOBS"
while IFS=$'\t' read -r sample_key layout cohort_id spike_stage spike_lot ratio host_count spike_count \
        fraction host_scale spike_scale median status failures warnings; do
    [[ "$sample_key" == "sample_key" || "$status" == "FAILED" ]] && continue
    parallel_pool_submit "$sample_key" write_spike_tracks "$sample_key" "$layout" "$cohort_id" "$host_scale" "$spike_scale"
done < "${qc_root}/spikein_scaling.tsv"
parallel_pool_wait_all
printf 'status\nSUCCESS\n' > "${root}/stage_status.tsv"
