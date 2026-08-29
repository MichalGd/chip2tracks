#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT_DIR=""
URL_BASE_OVERRIDE=""
URL_BASE_WAS_SET=false
TRACK_PREFIX_OVERRIDE=""
PYTHON_COMMAND="${PYTHON_COMMAND:-python3}"

usage() {
    echo "Usage: regenerate_ucsc_tracks.sh --output-dir OUTPUT_DIR [--url-base URL] [--track-prefix PREFIX]" >&2
}

while (( $# )); do
    case "$1" in
        --output-dir) OUTPUT_DIR="${2:?missing output directory}"; shift 2 ;;
        --url-base) URL_BASE_OVERRIDE="${2:?missing URL base}"; URL_BASE_WAS_SET=true; shift 2 ;;
        --track-prefix) TRACK_PREFIX_OVERRIDE="${2:?missing track prefix}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done
[[ -n "$OUTPUT_DIR" ]] || { usage; exit 2; }
config="${OUTPUT_DIR}/00_metadata/resolved_config.tsv"
[[ -s "$config" && -d "${OUTPUT_DIR}/04_tracks" ]] || {
    echo "ERROR: not a completed chip2tracks output directory: $OUTPUT_DIR" >&2
    exit 2
}
command -v "$PYTHON_COMMAND" >/dev/null 2>&1 || {
    echo "ERROR: Python command is unavailable: $PYTHON_COMMAND" >&2
    exit 2
}
url_base="$(awk -F '\t' '$1=="UCSC_BIGDATA_URL_BASE" {print $2; exit}' "$config")"
track_prefix="$(awk -F '\t' '$1=="UCSC_TRACK_PREFIX" {print $2; exit}' "$config")"
[[ "$URL_BASE_WAS_SET" == "false" ]] || url_base="$URL_BASE_OVERRIDE"
[[ -z "$TRACK_PREFIX_OVERRIDE" ]] || track_prefix="$TRACK_PREFIX_OVERRIDE"
[[ -n "$track_prefix" ]] || track_prefix=CHIP
"$PYTHON_COMMAND" "${ROOT}/scripts/generate_ucsc_tracks.py" "$OUTPUT_DIR" \
    --url-base "$url_base" --track-prefix "$track_prefix"
echo "UCSC descriptors regenerated: ${OUTPUT_DIR}/09_browser/ucsc"
