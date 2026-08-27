#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"
# shellcheck source=lib/parallel_jobs.sh
source "${SCRIPT_DIR}/lib/parallel_jobs.sh"
require_config

CPM_ROOT="${OUTPUT_DIR}/04_tracks/cpm"
COVERAGE_LOG_ROOT="${OUTPUT_DIR}/logs/coverage"
rm -rf -- "$CPM_ROOT"
mkdir -p "$CPM_ROOT" "$CPM_ROOT/permissive" "$CPM_ROOT/intermediate" \
    "$CPM_ROOT/stringent" "$CPM_ROOT/.mapping_parts" "$CPM_ROOT/.mapping_cache" "$COVERAGE_LOG_ROOT"

estimate_se_fragment_length() {
    local key="$1" layout="$2" bam="$3" genome="$4" metadata="$CPM_ROOT/${key}.fragment_length.tsv"
    local length="." source="paired_fragment" log="$COVERAGE_LOG_ROOT/${key}.macs3_predictd.log"
    if [[ "$layout" == "PE" ]]; then
        :
    elif [[ "$SE_SIGNAL_MODE" == "read" ]]; then
        source="unextended_read"
        :
    elif [[ "$SE_FRAGMENT_LENGTH" =~ ^[1-9][0-9]*$ ]]; then
        length="$SE_FRAGMENT_LENGTH"
        source="configured"
    else
        local effective_genome_size predictd_dir
        effective_genome_size="$(reference_value EFFECTIVE_GENOME_SIZE "$genome")"
        predictd_dir="$COVERAGE_LOG_ROOT/${key}.macs3_predictd"
        mkdir -p "$predictd_dir"
        if (cd "$predictd_dir" && run_logged "$MACS3_COMMAND" predictd -i "$bam" -f BAM \
                -g "$effective_genome_size") >"$log" 2>&1; then
            length="$(grep -Eio 'predicted fragment length is[[:space:]]+[0-9]+[[:space:]]+bps|#[[:space:]]*d[[:space:]]*=[[:space:]]*[0-9]+' "$log" \
                | grep -Eo '[0-9]+' | tail -n 1 || true)"
        fi
        if [[ "$length" =~ ^[1-9][0-9]*$ ]]; then
            source="macs3_predictd"
        else
            length="$SE_FRAGMENT_LENGTH_FALLBACK"
            source="configured_fallback"
            warn "MACS3 could not estimate an SE fragment length for $key; using $length bp"
        fi
    fi
    printf 'sample_key\tfragment_length_bp\tsource\tlog\n%s\t%s\t%s\t%s\n' \
        "$key" "$length" "$source" "$([[ -s "$log" ]] && echo "$log" || echo .)" > "$metadata"
    printf '%s' "$length"
}

policy_description() {
    local policy="$1" analysis_duplicate_policy="${2:-remove}"
    case "$policy" in
        analysis)
            if [[ "$analysis_duplicate_policy" == "retain" ]]; then
                printf 'q30_dup-retained\t30\tretained'
            else
                printf 'q30_dup-removed\t30\tremoved'
            fi
            ;;
        permissive) printf 'q0_dup-retained\t0\tretained' ;;
        intermediate) printf 'q0_dup-removed\t0\tremoved' ;;
        stringent) printf 'q30_dup-removed\t30\tremoved' ;;
        *) die "unknown CPM policy: $1" ;;
    esac
}

mapping_composition() {
    local bam="$1" layout="$2"
    local view_args=()
    if [[ "$layout" == "PE" ]]; then view_args=(-f 66 -F 3840); else view_args=(-F 3844); fi
    samtools view "${view_args[@]}" "$bam" | awk '
        BEGIN {total=0; mapq0=0; low=0; xs=0}
        {
            total++
            if ($5 == 0) mapq0++
            if ($5 < 30) low++
            for (i=12; i<=NF; i++) {
                if ($i ~ /^XS:[AifZHB]:/) {xs++; break}
            }
        }
        END {
            p0=(total ? 100*mapq0/total : 0)
            plow=(total ? 100*low/total : 0)
            pxs=(total ? 100*xs/total : 0)
            printf "%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f", total, mapq0, p0, low, plow, xs, pxs
        }'
}

write_track() {
    local key="$1" layout="$2" genome="$3" policy="$4" bam="$5" output_dir="$6" fragment_length="$7" analysis_duplicate_policy="$8"
    local count scale chrom_sizes tmp bedgraph bigwig log description branch mapq duplicates composition normalized_source
    [[ -s "$bam" ]] || die "CPM input BAM missing for $key ($policy): $bam"
    count="$(signal_count "$bam" "$layout")"
    (( count > 0 )) || die "zero observations for $key ($policy)"
    scale="$(awk -v n="$count" 'BEGIN {printf "%.15g", 1000000/n}')"
    chrom_sizes="$(reference_value CHROM_SIZES "$genome")"
    mkdir -p "$output_dir"
    tmp="${output_dir}/${key}.CPM.unsorted.bedGraph"
    bedgraph="${output_dir}/${key}.CPM.bedGraph"
    bigwig="${output_dir}/${key}.CPM.bw"
    log="${COVERAGE_LOG_ROOT}/${key}.${policy}.log"
    rm -f -- "$tmp" "$bedgraph" "$bigwig" "${bedgraph}.gz"

    set_bamcoverage_signal_args "$key" "$layout" "$fragment_length"
    local args=(--bam "$bam" --outFileName "$tmp" --outFileFormat bedgraph --binSize "$TRACK_BIN_SIZE"
        --scaleFactor "$scale" --numberOfProcessors "$THREADS_BAMCOVERAGE")
    args+=("${BAMCOVERAGE_SIGNAL_ARGS[@]}")
    local extra=()
    [[ -n "$BAMCOVERAGE_COMMON_ARGS" ]] && read -r -a extra <<< "$BAMCOVERAGE_COMMON_ARGS"
    run_logged bamCoverage "${args[@]}" "${extra[@]}" >"$log" 2>&1
    bedtools sort -faidx "$chrom_sizes" -i "$tmp" > "$bedgraph"
    awk 'BEGIN{ok=1} NF!=4 || $2<0 || $3<=$2 || $4<0 || $4!=$4 {ok=0} END{exit !ok}' "$bedgraph" || \
        die "invalid bedGraph: $bedgraph"
    if is_true "$GENERATE_COVERAGE_BIGWIGS"; then
        run_logged bedGraphToBigWig "$bedgraph" "$chrom_sizes" "$bigwig"
    fi
    is_true "$GENERATE_COVERAGE_BEDGRAPHS" || rm -f -- "$bedgraph"
    rm -f -- "$tmp"

    description="$(policy_description "$policy" "$analysis_duplicate_policy")"
    IFS=$'\t' read -r branch mapq duplicates <<< "$description"
    printf 'sample_key\tpolicy\tbam_branch\tmapq\tduplicates\tsignal_unit\tsignal_count\tscale\tformula\tse_fragment_length_bp\n' \
        > "${output_dir}/${key}.normalization_metadata.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tC*1e6/L\t%s\n' \
        "$key" "$policy" "$branch" "$mapq" "$duplicates" \
        "$([[ "$layout" == "PE" ]] && echo fragment || echo read)" "$count" "$scale" \
        "$fragment_length" >> "${output_dir}/${key}.normalization_metadata.tsv"

    composition="$(awk -F '\t' -v branch="$branch" '$1==branch {sub(/^[^\t]*\t/, ""); print; exit}' "$MAPPING_COMPOSITION_CACHE")"
    if [[ -z "$composition" ]]; then
        composition="$(mapping_composition "$bam" "$layout")"
        printf '%s\t%s\n' "$branch" "$composition" >> "$MAPPING_COMPOSITION_CACHE"
    fi
    normalized_source="deseq2_consensus"
    [[ "$policy" != "analysis" ]] && normalized_source="deseq2_robust_cpm/${policy}"
    printf '%s\t%s\tcpm/%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$key" "$policy" "$policy" "$normalized_source" "$branch" "$mapq" "$duplicates" \
        "$([[ "$layout" == "PE" ]] && echo fragment || echo read)" "$composition" \
        >> "$CPM_ROOT/.mapping_parts/${key}.tsv"
}

worker() {
    local key="$1" layout="$2" genome="$3" analysis_duplicate_policy="$4" analysis_bam fragment_length
    : > "$CPM_ROOT/.mapping_parts/${key}.tsv"
    MAPPING_COMPOSITION_CACHE="$CPM_ROOT/.mapping_cache/${key}.tsv"
    : > "$MAPPING_COMPOSITION_CACHE"
    analysis_bam="$(analysis_bam_path "$key")"
    fragment_length="$(estimate_se_fragment_length "$key" "$layout" "$analysis_bam" "$genome")"
    write_track "$key" "$layout" "$genome" analysis "$analysis_bam" "$CPM_ROOT" "$fragment_length" "$analysis_duplicate_policy"
    if is_true "$GENERATE_CPM_PERMISSIVE_TRACKS"; then
        write_track "$key" "$layout" "$genome" permissive "$(policy_bam_path "$key" permissive)" \
            "$CPM_ROOT/permissive" "$fragment_length" "$analysis_duplicate_policy"
    fi
    if is_true "$GENERATE_CPM_INTERMEDIATE_TRACKS"; then
        write_track "$key" "$layout" "$genome" intermediate "$(policy_bam_path "$key" intermediate)" \
            "$CPM_ROOT/intermediate" "$fragment_length" "$analysis_duplicate_policy"
    fi
    if is_true "$GENERATE_CPM_STRINGENT_TRACKS"; then
        write_track "$key" "$layout" "$genome" stringent "$(policy_bam_path "$key" stringent)" \
            "$CPM_ROOT/stringent" "$fragment_length" "$analysis_duplicate_policy"
    fi
}

if is_true "$GENERATE_CPM_TRACKS"; then
    parallel_pool_init "$TRACK_PARALLEL_JOBS"
    while IFS=$'\t' read -r sample_key sample_id replicate layout genome assay_profile factor antibody_id target_class \
            condition treatment cell_type is_control control_type control_id control_key analysis_duplicate_policy rest; do
        [[ "$sample_key" == "sample_key" ]] && continue
        parallel_pool_submit "$sample_key" worker "$sample_key" "$layout" "$genome" "$analysis_duplicate_policy"
    done < "$SAMPLE_MANIFEST"
    parallel_pool_wait_all
else
    printf '{"status":"SKIPPED","reason":"GENERATE_CPM_TRACKS=false"}\n' > "$CPM_ROOT/SKIPPED.json"
fi

mapping_summary="$CPM_ROOT/mapping_composition.tsv"
printf 'sample_key\tpolicy\tcpm_family\tnormalized_family_source\tbam_branch\tmapq_policy\tduplicates\tsignal_unit\ttotal_observations\tmapq0_observations\tmapq0_percent\tmapq_lt30_observations\tmapq_lt30_percent\txs_tagged_candidate_multimappers\txs_tagged_percent\n' \
    > "$mapping_summary"
find "$CPM_ROOT/.mapping_parts" -maxdepth 1 -type f -name '*.tsv' -print0 | sort -z | xargs -0 -r cat >> "$mapping_summary"
rm -rf -- "$CPM_ROOT/.mapping_parts" "$CPM_ROOT/.mapping_cache"
cat > "$CPM_ROOT/mapping_composition_definitions.tsv" <<'EOF'
metric	definition
mapq0_observations	Coverage signal units whose representative alignment has MAPQ 0; ambiguous but not a universal multimapper definition.
mapq_lt30_observations	Coverage signal units excluded by the stringent MAPQ 30 policy; includes ambiguous and other low-confidence alignments.
xs_tagged_candidate_multimappers	Coverage signal units whose representative Bowtie2 alignment contains an XS tag indicating a reported alternative alignment score.
EOF

printf 'family\tenabled\tmapq\tduplicates\tpath\n' > "$CPM_ROOT/track_family_status.tsv"
printf 'analysis\t%s\tconfigured\tconfigured\t04_tracks/cpm\n' "$GENERATE_CPM_TRACKS" >> "$CPM_ROOT/track_family_status.tsv"
for family in permissive intermediate stringent; do
    case "$family" in
        permissive) configured="$GENERATE_CPM_PERMISSIVE_TRACKS"; mapq=0; duplicates=retained ;;
        intermediate) configured="$GENERATE_CPM_INTERMEDIATE_TRACKS"; mapq=0; duplicates=removed ;;
        stringent) configured="$GENERATE_CPM_STRINGENT_TRACKS"; mapq=30; duplicates=removed ;;
    esac
    if is_true "$GENERATE_CPM_TRACKS" && is_true "$configured"; then enabled=true; else enabled=false; fi
    printf '%s\t%s\t%s\t%s\t04_tracks/cpm/%s\n' "$family" "$enabled" "$mapq" "$duplicates" "$family" \
        >> "$CPM_ROOT/track_family_status.tsv"
done
printf 'status\nSUCCESS\n' > "$CPM_ROOT/stage_status.tsv"
