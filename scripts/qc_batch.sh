#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

qc_root="${OUTPUT_DIR}/06_qc"
per_sample_root="$qc_root/alignment_and_complexity/.per_sample"
mkdir -p "$qc_root/alignment_and_complexity" "$qc_root/fragment_length_and_periodicity" \
    "$qc_root/frip_and_peak_reproducibility" "$qc_root/correlation_pca_fingerprint" \
    "$qc_root/tss_signal_profile/reference" "$qc_root/controls" "$per_sample_root" \
    "${OUTPUT_DIR}/logs/qc"

library_complexity() {
    local bam="$1" layout="$2" key="$3" output="$4"
    local tmp total distinct once twice nrf pbc1 pbc2
    tmp="$(mktemp -d "$qc_root/alignment_and_complexity/.complexity.XXXXXX")"
    if [[ "$layout" == "PE" ]]; then
        samtools sort -n -@ "$THREADS_SAMTOOLS" -o "$tmp/name.bam" "$bam"
        bedtools bamtobed -bedpe -i "$tmp/name.bam" | \
            awk 'BEGIN{OFS="\t"} $1==$4 {s=($2<$5?$2:$5); e=($3>$6?$3:$6); if(e>s) print $1,s,e}' > "$tmp/units.bed"
    else
        bedtools bamtobed -i "$bam" | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' > "$tmp/units.bed"
    fi
    LC_ALL=C sort "$tmp/units.bed" | uniq -c > "$tmp/multiplicity.txt"
    read -r total distinct once twice < <(
        awk '{total+=$1; distinct++; if($1==1)once++; if($1==2)twice++}
             END{print total+0,distinct+0,once+0,twice+0}' "$tmp/multiplicity.txt"
    )
    nrf="$(awk -v d="$distinct" -v n="$total" 'BEGIN{if(n)printf "%.6f",d/n;else print "NA"}')"
    pbc1="$(awk -v o="$once" -v d="$distinct" 'BEGIN{if(d)printf "%.6f",o/d;else print "NA"}')"
    pbc2="$(awk -v o="$once" -v t="$twice" 'BEGIN{if(t)printf "%.6f",o/t;else print "Inf"}')"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$key" "$layout" "$total" "$distinct" "$once" "$twice" "$nrf" "$pbc1" "$pbc2" > "$output"
    rm -rf -- "$tmp"
}

qc_worker() {
    local sample_key="$1" layout="$2" genome="$3" is_control="$4" control_key="$5"
    local duplicate_policy="$6" cohort_id="$7" primary_caller="$8" primary_class="$9"
    local bam count unit retained_bam fragment_exclude phantompeak_path tagalign
    local consensus tmp total in_peaks frip control_bam
    bam="$(analysis_bam_path "$sample_key")"
    count="$(signal_count "$bam" "$layout" "$duplicate_policy")"
    unit="$([[ "$layout" == "PE" ]] && echo fragment || echo read)"
    printf '%s\t%s\t%s\t%s\n' "$sample_key" "$layout" "$unit" "$count" \
        > "$per_sample_root/${sample_key}.observations.tsv"
    samtools flagstat -@ "$THREADS_SAMTOOLS" "$bam" \
        > "$qc_root/alignment_and_complexity/${sample_key}.flagstat.txt"
    samtools stats -@ "$THREADS_SAMTOOLS" "$bam" \
        > "$qc_root/alignment_and_complexity/${sample_key}.samtools_stats.txt"

    retained_bam="${OUTPUT_DIR}/03_alignment/filtered/q30_dup-retained/${sample_key}.host.q30.dup-retained.bam"
    if is_true "$RUN_LIBRARY_COMPLEXITY"; then
        [[ -s "$retained_bam" ]] || die "duplicate-retained BAM required for complexity QC: $retained_bam"
        library_complexity "$retained_bam" "$layout" "$sample_key" \
            "$per_sample_root/${sample_key}.complexity.tsv"
    fi
    if [[ "$layout" == "PE" ]] && is_true "$RUN_FRAGMENT_QC"; then
        fragment_exclude="$(signal_exclude_mask "$layout" "$duplicate_policy")"
        samtools view -f 66 -F "$fragment_exclude" "$bam" |
            awk -v maximum="$FRAGMENT_PLOT_MAX_BP" 'BEGIN{OFS="\t"} {t=$9; if(t<0)t=-t; if(t>0 && t<=maximum)n[t]++} END{print "fragment_length","count"; for(i=1;i<=maximum;i++)print i,n[i]+0}' \
            > "$qc_root/fragment_length_and_periodicity/${sample_key}.fragment_lengths.tsv"
    fi
    if is_true "$RUN_PRESEQ"; then
        preseq lc_extrap -B -o "$qc_root/alignment_and_complexity/${sample_key}.preseq.txt" "$retained_bam" \
            >"${OUTPUT_DIR}/logs/qc/${sample_key}.preseq.log" 2>&1 || warn "preseq failed for $sample_key"
    fi
    if is_true "$RUN_CROSS_CORRELATION"; then
        phantompeak_path="$(command -v "$PHANTOMPEAK_COMMAND" || true)"
        if [[ -n "$phantompeak_path" ]]; then
            tagalign="$qc_root/fragment_length_and_periodicity/${sample_key}.q30_dup-retained.tagAlign.gz"
            bedtools bamtobed -i "$retained_bam" | \
                awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"N",1000,$6}' | gzip -c > "$tagalign"
            "$phantompeak_path" -c="$tagalign" \
                -savp="$qc_root/fragment_length_and_periodicity/${sample_key}.cross_correlation.pdf" \
                -out="$qc_root/fragment_length_and_periodicity/${sample_key}.phantompeak.tsv" \
                >"${OUTPUT_DIR}/logs/qc/${sample_key}.phantompeak.log" 2>&1 || warn "cross-correlation failed for $sample_key"
        fi
    fi

    if [[ "$is_control" == "FALSE" ]]; then
        consensus="$(find "${OUTPUT_DIR}/05_peaks/consensus/${cohort_id}/${primary_caller}/${primary_class}" \
            -maxdepth 1 -type f -name '*.consensus.bed' -print -quit 2>/dev/null || true)"
        if [[ -n "$consensus" ]]; then
            if [[ "$layout" == "PE" ]]; then
                tmp="$(mktemp -d "$qc_root/frip_and_peak_reproducibility/.frip.XXXXXX")"
                samtools sort -n -@ "$THREADS_SAMTOOLS" -o "$tmp/name.bam" "$bam"
                bedtools bamtobed -bedpe -i "$tmp/name.bam" |
                    awk 'BEGIN{OFS="\t"} $1==$4 {start=($2<$5?$2:$5); end=($3>$6?$3:$6); if(end>start)print $1,start,end}' > "$tmp/fragments.bed"
                total="$(wc -l < "$tmp/fragments.bed")"
                in_peaks="$(bedtools intersect -u -a "$tmp/fragments.bed" -b "$consensus" | wc -l)"
                rm -rf -- "$tmp"
            else
                total="$count"
                in_peaks="$(bedtools intersect -u -abam "$bam" -b "$consensus" | samtools view -c -)"
            fi
            frip="$(awk -v a="$in_peaks" -v n="$total" 'BEGIN{if(n>0)printf "%.8f",a/n; else print "NA"}')"
            printf 'sample_key\tsignal_unit\ttotal\tin_consensus\tfrip\n%s\t%s\t%s\t%s\t%s\n' \
                "$sample_key" "$unit" "$total" "$in_peaks" "$frip" \
                > "$qc_root/frip_and_peak_reproducibility/${sample_key}.frip.tsv"
        fi
        if [[ "$control_key" != "." ]]; then
            control_bam="$(analysis_bam_path "$control_key")"
            fingerprint_args=(-b "$bam" "$control_bam" --labels "$sample_key" "$control_key"
                --plotFile "$qc_root/controls/${sample_key}.target_control_fingerprint.png"
                --outRawCounts "$qc_root/controls/${sample_key}.target_control_fingerprint.tsv"
                --numberOfProcessors "$THREADS_BAMCOVERAGE")
            [[ "$layout" == "PE" ]] && fingerprint_args+=(--samFlagInclude 66 --extendReads)
            plotFingerprint "${fingerprint_args[@]}" >"${OUTPUT_DIR}/logs/qc/${sample_key}.fingerprint.log" 2>&1 || \
                warn "fingerprint failed for $sample_key"
        fi
    fi
}

parallel_pool_init "$QC_SAMPLE_PARALLEL_JOBS"
while IFS=$'\t' read -r \
    sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
    is_control control_type control_id control_key duplicate_policy blacklist ratio spike_stage spike_lot batch donor output_prefix \
    technical_units fastq_1_list fastq_2_list cohort_id cohort_key primary_caller primary_class; do
    [[ "$sample_key" == "sample_key" ]] && continue
    parallel_pool_submit "$sample_key" qc_worker "$sample_key" "$layout" "$genome" "$is_control" \
        "$control_key" "$duplicate_policy" "$cohort_id" "$primary_caller" "$primary_class"
done < "$SAMPLE_MANIFEST"
parallel_pool_wait_all

summary="$qc_root/alignment_and_complexity/observation_counts.tsv"
printf 'sample_key\tlayout\tsignal_unit\tanalysis_observations\n' > "$summary"
complexity_summary="$qc_root/alignment_and_complexity/library_complexity.tsv"
printf 'sample_key\tlayout\ttotal_observations\tdistinct_observations\tonce\ttwice\tNRF\tPBC1\tPBC2\n' > "$complexity_summary"
target_bams=()
target_labels=()
declare -A prepared_tss=()
while IFS=$'\t' read -r \
    sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class condition treatment cell_type \
    is_control rest; do
    [[ "$sample_key" == "sample_key" ]] && continue
    cat "$per_sample_root/${sample_key}.observations.tsv" >> "$summary"
    if is_true "$RUN_LIBRARY_COMPLEXITY"; then
        cat "$per_sample_root/${sample_key}.complexity.tsv" >> "$complexity_summary"
    fi
    if [[ "$is_control" == "FALSE" ]]; then
        target_bams+=("$(analysis_bam_path "$sample_key")")
        target_labels+=("$sample_key")
    fi
    if is_true "$RUN_TSS_SIGNAL_PROFILE" && [[ -z "${prepared_tss[$genome]:-}" ]]; then
        tss="$(optional_reference_value TSS_BED "$genome")"
        if [[ -z "$tss" || "$tss" == "." ]]; then
            tss="$qc_root/tss_signal_profile/reference/${genome}.tss.bed"
            [[ -s "$tss" ]] || python3 "${SCRIPT_DIR}/prepare_tss_bed.py" "$(reference_value GTF "$genome")" "$tss"
        fi
        prepared_tss[$genome]="$tss"
    fi
done < "$SAMPLE_MANIFEST"

if is_true "$RUN_REPLICATE_CORRELATION" && (( ${#target_bams[@]} >= 2 )); then
    matrix="$qc_root/correlation_pca_fingerprint/target_bins.npz"
    raw="$qc_root/correlation_pca_fingerprint/target_bins.tsv"
    multiBamSummary bins --bamfiles "${target_bams[@]}" --labels "${target_labels[@]}" \
        --numberOfProcessors "$THREADS_BAMCOVERAGE" --outFileName "$matrix" --outRawCounts "$raw"
    plotCorrelation --corData "$matrix" --corMethod spearman --whatToPlot heatmap \
        --skipZeros --plotFile "$qc_root/correlation_pca_fingerprint/spearman_heatmap.png" \
        --outFileCorMatrix "$qc_root/correlation_pca_fingerprint/spearman_matrix.tsv"
    plotPCA --corData "$matrix" --plotFile "$qc_root/correlation_pca_fingerprint/pca.png" \
        --outFileNameData "$qc_root/correlation_pca_fingerprint/pca.tsv"
fi

tss_worker() {
    local sample_key="$1" genome="$2" bw tss
    bw="${OUTPUT_DIR}/04_tracks/cpm/${sample_key}.CPM.bw"
    [[ -s "$bw" ]] || return 0
    tss="${prepared_tss[$genome]}"
    computeMatrix reference-point --referencePoint TSS -b "$TSS_PROFILE_UPSTREAM" -a "$TSS_PROFILE_DOWNSTREAM" \
        -R "$tss" -S "$bw" -o "$qc_root/tss_signal_profile/${sample_key}.matrix.gz" \
        --numberOfProcessors "$THREADS_BAMCOVERAGE"
    plotProfile -m "$qc_root/tss_signal_profile/${sample_key}.matrix.gz" \
        -out "$qc_root/tss_signal_profile/${sample_key}.descriptive_TSS_profile.png" \
        --plotTitle "$sample_key descriptive TSS signal"
}

if is_true "$RUN_TSS_SIGNAL_PROFILE"; then
    parallel_pool_init "$QC_SAMPLE_PARALLEL_JOBS"
    while IFS=$'\t' read -r sample_key sample_id replicate layout genome rest; do
        [[ "$sample_key" == "sample_key" ]] && continue
        parallel_pool_submit "tss:$sample_key" tss_worker "$sample_key" "$genome"
    done < "$SAMPLE_MANIFEST"
    parallel_pool_wait_all
fi

printf 'status\tthreshold_mode\tparallel_jobs\nSUCCESS\tdescriptive\t%s\n' \
    "$QC_SAMPLE_PARALLEL_JOBS" > "$qc_root/stage_status.tsv"
