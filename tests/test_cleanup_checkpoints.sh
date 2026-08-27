#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf -- "$TMP"' EXIT
OUTPUT="$TMP/output"
mkdir -p "$OUTPUT/00_metadata" "$OUTPUT/02_trimmed_fastq" "$OUTPUT/03_alignment/sorted" \
    "$OUTPUT/.checkpoints" "$OUTPUT/10_reports"
printf 'sample_key\nS1\n' > "$OUTPUT/00_metadata/sample_manifest.tsv"
printf 'report\n' > "$OUTPUT/10_reports/pipeline_report.html"
printf 'fastq\n' > "$OUTPUT/02_trimmed_fastq/S1.fastq.gz"
printf 'bam\n' > "$OUTPUT/03_alignment/sorted/S1.bam"
for stage in preprocess alignment filtering; do printf '{}\n' > "$OUTPUT/.checkpoints/$stage.json"; done
cat > "$TMP/config.sh" <<EOF
OUTPUT_DIR='$OUTPUT'
ENABLE_AUTOMATIC_CLEANUP=true
KEEP_TRIMMED_FASTQ=false
KEEP_RAW_ALIGNMENT_BAMS=false
KEEP_MARKED_BAMS=true
KEEP_FILTERED_BAMS=true
KEEP_SPIKEIN_BAMS=true
EOF
C2T_CONFIG="$TMP/config.sh" bash "$ROOT/scripts/cleanup.sh"
test ! -e "$OUTPUT/02_trimmed_fastq"
test ! -e "$OUTPUT/03_alignment/sorted"
test ! -e "$OUTPUT/.checkpoints/preprocess.json"
test ! -e "$OUTPUT/.checkpoints/alignment.json"
test -e "$OUTPUT/.checkpoints/filtering.json"
grep -q 'invalidated_after_intermediate_cleanup' "$OUTPUT/00_metadata/cleanup_manifest.tsv"
echo "Cleanup checkpoint-invalidation regression test passed"
