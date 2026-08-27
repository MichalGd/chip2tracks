#!/usr/bin/env bash
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }
note() { echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"; }

require_config() {
    [[ -n "${C2T_CONFIG:-}" && -f "$C2T_CONFIG" ]] || die "C2T_CONFIG is not a readable file"
    # shellcheck disable=SC1090
    source "$C2T_CONFIG"
    : "${OUTPUT_DIR:?OUTPUT_DIR missing}"
    SAMPLE_MANIFEST="${OUTPUT_DIR}/00_metadata/sample_manifest.tsv"
    COHORT_MANIFEST="${OUTPUT_DIR}/00_metadata/cohort_manifest.tsv"
    [[ -s "$SAMPLE_MANIFEST" ]] || die "sample manifest missing: $SAMPLE_MANIFEST"
}

is_true() { [[ "${1,,}" == "true" ]]; }

reference_value() {
    local prefix="${1:?prefix}" genome="${2:?genome}" key value
    key="${prefix}_$(printf '%s' "$genome" | tr '[:lower:].-' '[:upper:]__')"
    value="${!key:-}"
    [[ -n "$value" ]] || die "missing reference setting $key"
    printf '%s' "$value"
}

optional_reference_value() {
    local prefix="${1:?prefix}" genome="${2:?genome}" key
    key="${prefix}_$(printf '%s' "$genome" | tr '[:lower:].-' '[:upper:]__')"
    printf '%s' "${!key:-}"
}

signal_exclude_mask() {
    local layout="${1:?layout}" duplicate_policy="${2:?duplicate policy}"
    case "$duplicate_policy" in
        retain|retained) duplicate_policy=retain ;;
        remove|removed) duplicate_policy=remove ;;
        *) die "unsupported duplicate policy for signal counting: $duplicate_policy" ;;
    esac
    case "$layout:$duplicate_policy" in
        PE:retain) printf '2828' ;;
        PE:remove) printf '3852' ;;
        SE:retain) printf '2820' ;;
        SE:remove) printf '3844' ;;
        *) die "unsupported layout for signal counting: $layout" ;;
    esac
}

signal_count() {
    local bam="${1:?bam}" layout="${2:?layout}" duplicate_policy="${3:?duplicate policy}"
    local exclude
    exclude="$(signal_exclude_mask "$layout" "$duplicate_policy")"
    if [[ "$layout" == "PE" ]]; then
        # Count one properly paired, primary representative read per fragment.
        samtools view -c -f 66 -F "$exclude" "$bam"
    else
        samtools view -c -F "$exclude" "$bam"
    fi
}

analysis_bam_path() {
    printf '%s/03_alignment/analysis/%s.host.analysis.bam' "$OUTPUT_DIR" "${1:?sample_key}"
}

policy_bam_path() {
    local key="${1:?sample_key}" policy="${2:?policy}"
    case "$policy" in
        analysis) analysis_bam_path "$key" ;;
        permissive) printf '%s/03_alignment/filtered/q0_dup-retained/%s.host.q0.dup-retained.bam' "$OUTPUT_DIR" "$key" ;;
        intermediate) printf '%s/03_alignment/filtered/q0_dup-removed/%s.host.q0.dup-removed.bam' "$OUTPUT_DIR" "$key" ;;
        stringent) printf '%s/03_alignment/filtered/q30_dup-removed/%s.host.q30.dup-removed.bam' "$OUTPUT_DIR" "$key" ;;
        *) die "unknown BAM policy: $policy" ;;
    esac
}

set_bamcoverage_signal_args() {
    local key="${1:?sample_key}" layout="${2:?layout}" resolved_length="${3:-}"
    BAMCOVERAGE_SIGNAL_ARGS=()
    if [[ "$layout" == "PE" ]]; then
        # One properly paired read represents each analysis fragment.
        BAMCOVERAGE_SIGNAL_ARGS+=(--samFlagInclude 66 --extendReads)
        return 0
    fi
    [[ "$layout" == "SE" ]] || die "unsupported layout for coverage: $layout"
    [[ "$SE_SIGNAL_MODE" == "extend" ]] || return 0
    if [[ -z "$resolved_length" ]]; then
        if [[ "$SE_FRAGMENT_LENGTH" =~ ^[1-9][0-9]*$ ]]; then
            resolved_length="$SE_FRAGMENT_LENGTH"
        else
            local metadata="${OUTPUT_DIR}/04_tracks/cpm/${key}.fragment_length.tsv"
            if [[ -s "$metadata" ]]; then
                resolved_length="$(awk -F '\t' 'NR==2 {print $2; exit}' "$metadata")"
            fi
            if [[ ! "$resolved_length" =~ ^[1-9][0-9]*$ ]]; then
                resolved_length="$SE_FRAGMENT_LENGTH_FALLBACK"
                warn "using SE_FRAGMENT_LENGTH_FALLBACK=$resolved_length for $key"
            fi
        fi
    fi
    [[ "$resolved_length" =~ ^[1-9][0-9]*$ ]] || die "invalid SE fragment length for $key: $resolved_length"
    BAMCOVERAGE_SIGNAL_ARGS+=(--extendReads "$resolved_length")
}

record_command() {
    is_true "${WRITE_COMMAND_LOG:-true}" || return 0
    local quoted=() argument
    for argument in "$@"; do quoted+=("$(printf '%q' "$argument")"); done
    printf '%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${quoted[*]}" >> "${OUTPUT_DIR}/00_metadata/commands.log"
}

run_logged() {
    record_command "$@"
    "$@"
}
