#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

mkdir -p "${OUTPUT_DIR}/05_peaks/per_sample" "${OUTPUT_DIR}/logs/peakcalling"

genome_size() {
    local genome="$1" size
    size="$(reference_value EFFECTIVE_GENOME_SIZE "$genome")"
    [[ "$size" =~ ^[0-9]+$ && "$size" -gt 0 ]] || \
        die "EFFECTIVE_GENOME_SIZE_${genome^^} must be a positive integer"
    printf '%s' "$size"
}

call_macs3_class() {
    local key="$1" layout="$2" genome="$3" peak_class="$4" target_bam="$5" control_bam="$6" out="$7"
    mkdir -p "$out"
    local format=() control=() common=()
    [[ "$control_bam" != "." ]] && control=(-c "$control_bam")
    if [[ "$layout" == "PE" ]]; then
        format=(-f BAMPE)
    else
        format=(-f BAM)
        [[ "$MACS3_SE_MODE" == "fixed" ]] && \
            format+=(--nomodel --shift "$MACS3_SE_SHIFT" --extsize "$MACS3_SE_EXTSIZE")
    fi
    common=(callpeak -t "$target_bam" "${control[@]}" "${format[@]}" -g "$(genome_size "$genome")"
        --keep-dup "$MACS3_KEEP_DUP" --outdir "$out")
    is_true "$MACS3_GENERATE_SIGNAL_TRACKS" && common+=(-B --SPMR)
    if [[ "$peak_class" == "narrow" ]]; then
        local summit_args=() narrow="${out}/${key}.macs3.narrow_peaks.narrowPeak"
        rm -f "$narrow"
        is_true "$MACS3_CALL_SUMMITS" && summit_args+=(--call-summits)
        run_logged "$MACS3_COMMAND" "${common[@]}" -q "$MACS3_QVALUE" -n "${key}.macs3.narrow" "${summit_args[@]}" || return $?
        [[ ! -s "$narrow" ]] || cut -f1-3 "$narrow" > "${out}/${key}.macs3.narrow.bed"
        if is_true "$RUN_IDR"; then
            rm -f "${out}/${key}.macs3.idr_peaks.narrowPeak"
            run_logged "$MACS3_COMMAND" "${common[@]}" -p "$IDR_MACS3_PVALUE" \
                -n "${key}.macs3.idr" --call-summits || return $?
        fi
    else
        local broad="${out}/${key}.macs3.broad_peaks.broadPeak"
        rm -f "$broad"
        run_logged "$MACS3_COMMAND" "${common[@]}" -q "$MACS3_QVALUE" -n "${key}.macs3.broad" --broad --broad-cutoff "$MACS3_BROAD_CUTOFF" || return $?
        [[ ! -s "$broad" ]] || cut -f1-3 "$broad" > "${out}/${key}.macs3.broad.bed"
    fi
}

call_epic2() {
    local key="$1" layout="$2" genome="$3" target_bam="$4" control_bam="$5" out="$6"
    mkdir -p "$out"
    local output="${out}/${key}.epic2.broad.tsv" control=() layout_args=() chrom_sizes effective total fraction
    chrom_sizes="$(reference_value CHROM_SIZES "$genome")"
    effective="$(genome_size "$genome")"
    total="$(awk '{n += $2} END {printf "%.0f", n}' "$chrom_sizes")"
    fraction="$(awk -v e="$effective" -v t="$total" 'BEGIN {if(t<=0) exit 1; printf "%.8f", e/t}')"
    [[ "$control_bam" != "." ]] && control=(--control "$control_bam")
    [[ "$layout" == "PE" ]] && layout_args+=(--guess-bampe)
    rm -f "$output" "${out}/${key}.epic2.broad.bed"
    run_logged "$EPIC2_COMMAND" --treatment "$target_bam" "${control[@]}" \
        --chromsizes "$chrom_sizes" --effective-genome-fraction "$fraction" \
        --bin-size "$EPIC2_BIN_SIZE" --gaps-allowed "$EPIC2_GAP_SIZE" \
        --fragment-size "$EPIC2_FRAGMENT_SIZE" \
        --false-discovery-rate-cutoff "$EPIC2_FDR" --keep-duplicates --mapq 0 \
        "${layout_args[@]}" --output "$output" || return $?
    [[ ! -s "$output" ]] || awk 'BEGIN{OFS="\t"} $0 !~ /^(#|track|browser)/ && NF>=3 {print $1,$2,$3}' \
        "$output" > "${out}/${key}.epic2.broad.bed"
}

worker() {
    local key="$1" layout="$2" genome="$3" target_class="$4" control_key="$5" primary_caller="$6" primary_class="$7"
    local target_bam control_bam="." root caller_status_file metadata_file
    local primary primary_status primary_count primary_reason failed_callers="." any_problem=false
    declare -A caller_status=()
    target_bam="$(analysis_bam_path "$key")"
    [[ "$control_key" != "." ]] && control_bam="$(analysis_bam_path "$control_key")"
    root="${OUTPUT_DIR}/05_peaks/per_sample/${key}"
    mkdir -p "$root"
    caller_status_file="${root}/caller_status.tsv"
    metadata_file="${root}/peakcall_metadata.tsv"
    printf 'sample_key\tcaller\tpeak_class\tstatus\tpeak_count\tpeak_file\tlog\treason\n' > "$caller_status_file"

    run_caller() {
        local caller="$1" peak_class="$2" peak_file="$3" log="$4"
        shift 4
        local status count=0 reason command_status=0
        mkdir -p "$(dirname "$peak_file")" "$(dirname "$log")"
        rm -f "$peak_file"
        if ( "$@" ) >"$log" 2>&1; then command_status=0; else command_status=$?; fi
        if (( command_status != 0 )); then
            status="ERROR"; reason="caller_exit_${command_status}"; : > "$peak_file"
        elif [[ -s "$peak_file" ]]; then
            status="SUCCESS"
            count="$(awk 'NF>=3 && $0 !~ /^(#|track|browser)/ {n++} END{print n+0}' "$peak_file")"
            reason="."
        else
            status="EMPTY"; reason="no_peaks"; : > "$peak_file"
        fi
        caller_status["${caller}:${peak_class}"]="$status"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$key" "$caller" "$peak_class" "$status" "$count" "$peak_file" "$log" "$reason" >> "$caller_status_file"
        if [[ "$status" == "ERROR" ]] || { [[ "$status" == "EMPTY" ]] && ! is_true "$ALLOW_EMPTY_PEAKS"; }; then
            any_problem=true
        fi
    }

    if [[ ",$PEAK_CALLERS," == *,macs3,* ]]; then
        if [[ "$target_class" == "narrow" || "$target_class" == "mixed" ]]; then
            run_caller macs3 narrow "${root}/macs3/${key}.macs3.narrow.bed" \
                "${OUTPUT_DIR}/logs/peakcalling/${key}.macs3.narrow.log" \
                call_macs3_class "$key" "$layout" "$genome" narrow "$target_bam" "$control_bam" "${root}/macs3"
        fi
        if [[ "$target_class" == "broad" || "$target_class" == "mixed" ]]; then
            run_caller macs3 broad "${root}/macs3/${key}.macs3.broad.bed" \
                "${OUTPUT_DIR}/logs/peakcalling/${key}.macs3.broad.log" \
                call_macs3_class "$key" "$layout" "$genome" broad "$target_bam" "$control_bam" "${root}/macs3"
        fi
    fi
    if [[ ",$PEAK_CALLERS," == *,epic2,* && ( "$target_class" == "broad" || "$target_class" == "mixed" ) ]]; then
        run_caller epic2 broad "${root}/epic2/${key}.epic2.broad.bed" \
            "${OUTPUT_DIR}/logs/peakcalling/${key}.epic2.broad.log" \
            call_epic2 "$key" "$layout" "$genome" "$target_bam" "$control_bam" "${root}/epic2"
    fi

    primary="${root}/${primary_caller}/${key}.${primary_caller}.${primary_class}.bed"
    primary_status="${caller_status[${primary_caller}:${primary_class}]:-ERROR}"
    primary_count=0
    [[ -s "$primary" ]] && primary_count="$(awk 'NF>=3 && $0 !~ /^(#|track|browser)/ {n++} END{print n+0}' "$primary")"
    case "$primary_status" in
        SUCCESS) primary_reason="." ;;
        EMPTY) primary_reason="primary_caller_produced_no_peaks" ;;
        *) primary_reason="primary_caller_error" ;;
    esac
    failed_callers="$(awk -F '\t' 'NR>1 && $4!="SUCCESS" {print $2 ":" $3 "=" $4}' "$caller_status_file" | paste -sd, -)"
    [[ -n "$failed_callers" ]] || failed_callers="."
    printf 'sample_key\tcontrol_key\tprimary_caller\tprimary_class\tstatus\tprimary_peak_count\tcaller_warnings\treason\n' > "$metadata_file"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$key" "$control_key" "$primary_caller" "$primary_class" "$primary_status" "$primary_count" \
        "$failed_callers" "$primary_reason" >> "$metadata_file"

    if [[ "$PEAKCALL_FAILURE_POLICY" == "fail" ]] && is_true "$any_problem"; then
        echo "ERROR: peak caller failure for $key: $failed_callers" >&2
        return 1
    fi
    if [[ "$primary_status" != "SUCCESS" ]]; then
        warn "excluding $key from consensus peak contribution: $primary_reason"
    elif [[ "$failed_callers" != "." ]]; then
        warn "non-primary peak caller warning for $key: $failed_callers"
    fi
}

parallel_pool_init "$PEAKCALL_PARALLEL_JOBS"
while IFS=$'\t' read -r \
    sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
    is_control control_type control_id control_key duplicate_policy blacklist ratio spike_stage spike_lot batch donor output_prefix \
    technical_units fastq_1_list fastq_2_list cohort_id cohort_key primary_caller primary_class; do
    [[ "$sample_key" == "sample_key" || "$is_control" == "TRUE" ]] && continue
    parallel_pool_submit "$sample_key" worker "$sample_key" "$layout" "$genome" "$target_class" \
        "$control_key" "$primary_caller" "$primary_class"
done < "$SAMPLE_MANIFEST"
parallel_pool_wait_all

if is_true "$MACS3_GENERATE_SIGNAL_TRACKS"; then
    signal_root="${OUTPUT_DIR}/04_tracks/control_normalized"
    mkdir -p "$signal_root"
    signal_status="$signal_root/status.tsv"
    printf 'sample_key\tbackground\tmetric\tbedgraph\tbigwig\n' > "$signal_status"
    while IFS=$'\t' read -r \
        sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
        is_control control_type control_id control_key duplicate_policy blacklist ratio spike_stage spike_lot batch donor output_prefix \
        technical_units fastq_1_list fastq_2_list cohort_id cohort_key primary_caller primary_class; do
        [[ "$sample_key" == "sample_key" || "$is_control" == "TRUE" ]] && continue
        macs_dir="${OUTPUT_DIR}/05_peaks/per_sample/${sample_key}/macs3"
        pileup="$(find "$macs_dir" -maxdepth 1 -type f -name '*_treat_pileup.bdg' -print -quit 2>/dev/null || true)"
        lambda="$(find "$macs_dir" -maxdepth 1 -type f -name '*_control_lambda.bdg' -print -quit 2>/dev/null || true)"
        [[ -s "$pileup" && -s "$lambda" ]] || continue
        raw="$signal_root/${sample_key}.FE.unsorted.bedGraph"
        bedgraph="$signal_root/${sample_key}.FE.bedGraph"
        bigwig="$signal_root/${sample_key}.FE.bw"
        run_logged "$MACS3_COMMAND" bdgcmp -t "$pileup" -c "$lambda" -m FE -o "$raw"
        LC_ALL=C sort -k1,1 -k2,2n "$raw" > "$bedgraph"
        rm -f "$raw"
        if is_true "$GENERATE_COVERAGE_BIGWIGS"; then
            run_logged bedGraphToBigWig "$bedgraph" "$(reference_value CHROM_SIZES "$genome")" "$bigwig"
        else
            bigwig="."
        fi
        if ! is_true "$GENERATE_COVERAGE_BEDGRAPHS"; then
            rm -f -- "$bedgraph"
            bedgraph="."
        fi
        background="$([[ "$control_key" == "." ]] && echo local_lambda_no_matched_control || echo matched_control_and_local_lambda)"
        printf '%s\t%s\tfold_enrichment\t%s\t%s\n' "$sample_key" "$background" "$bedgraph" "$bigwig" >> "$signal_status"
    done < "$SAMPLE_MANIFEST"
fi

summary="${OUTPUT_DIR}/05_peaks/per_sample/peakcall_status.tsv"
printf 'sample_key\tcontrol_key\tprimary_caller\tprimary_class\tstatus\tprimary_peak_count\tcaller_warnings\treason\n' > "$summary"
warning_count=0
target_count=0
primary_success_count=0
while IFS=$'\t' read -r \
    sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
    is_control rest; do
    [[ "$sample_key" == "sample_key" || "$is_control" == "TRUE" ]] && continue
    target_count=$((target_count + 1))
    metadata="${OUTPUT_DIR}/05_peaks/per_sample/${sample_key}/peakcall_metadata.tsv"
    [[ -s "$metadata" ]] || die "peak-call metadata missing for $sample_key"
    tail -n 1 "$metadata" >> "$summary"
    status="$(awk -F '\t' 'NR==2 {print $5}' "$metadata")"
    warnings="$(awk -F '\t' 'NR==2 {print $7}' "$metadata")"
    [[ "$status" == "SUCCESS" && "$warnings" == "." ]] || warning_count=$((warning_count + 1))
    [[ "$status" == "SUCCESS" ]] && primary_success_count=$((primary_success_count + 1))
done < "$SAMPLE_MANIFEST"
if (( target_count > 0 && primary_success_count == 0 )); then
    printf 'status\ttarget_samples\tprimary_successes\twarning_samples\tfailure_policy\nFAILED\t%s\t0\t%s\t%s\n' \
        "$target_count" "$warning_count" "$PEAKCALL_FAILURE_POLICY" \
        > "${OUTPUT_DIR}/05_peaks/per_sample/stage_status.tsv"
    die "all $target_count primary peak calls failed or produced no accepted peaks"
fi
stage_status="$([[ "$warning_count" -eq 0 ]] && echo SUCCESS || echo COMPLETED_WITH_WARNINGS)"
printf 'status\ttarget_samples\tprimary_successes\twarning_samples\tfailure_policy\n%s\t%s\t%s\t%s\t%s\n' \
    "$stage_status" "$target_count" "$primary_success_count" "$warning_count" "$PEAKCALL_FAILURE_POLICY" \
    > "${OUTPUT_DIR}/05_peaks/per_sample/stage_status.tsv"
