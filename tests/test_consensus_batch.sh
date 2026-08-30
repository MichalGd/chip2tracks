#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf -- "$TMP"' EXIT

mkdir -p "$TMP/bin" "$TMP/output/00_metadata"
cat > "$TMP/bin/python3" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
status_file=""
while (( $# > 0 )); do
    case "$1" in
        --status-file) status_file="$2"; shift 2 ;;
        *) shift ;;
    esac
done
if [[ "${MOCK_CONSENSUS_FAIL:-false}" == "true" ]]; then
    exit 7
fi
[[ -n "$status_file" ]]
mkdir -p "$(dirname "$status_file")"
printf 'cohort_id\tstatus\ttotal_samples\tsuccessful_peak_samples\texcluded_samples\tregions\treason\n' \
    > "$status_file"
printf 'C1\tSUCCESS\t2\t2\t0\t1\t.\n' >> "$status_file"
MOCK
chmod +x "$TMP/bin/python3"

printf 'sample_key\tis_control\tcohort_id\nA\tFALSE\tC1\nB\tFALSE\tC1\n' \
    > "$TMP/output/00_metadata/sample_manifest.tsv"
printf 'cohort_id\tprimary_peak_caller\tprimary_peak_class\nC1\tmacs3\tbroad\n' \
    > "$TMP/output/00_metadata/cohort_manifest.tsv"

for sample in A B; do
    sample_root="$TMP/output/05_peaks/per_sample/$sample"
    mkdir -p "$sample_root/macs3"
    printf 'chr1\t10\t30\n' > "$sample_root/macs3/$sample.macs3.broad.bed"
    printf 'sample_key\tprimary_caller\tprimary_class\tstatus\treason\n%s\tmacs3\tbroad\tSUCCESS\t.\n' \
        "$sample" > "$sample_root/peakcall_metadata.tsv"
done

cat > "$TMP/config.sh" <<EOF
OUTPUT_DIR=$TMP/output
CONSENSUS_MIN_BIOLOGICAL_SAMPLES=2
ALLOW_SINGLE_SAMPLE_CONSENSUS=false
REQUIRE_ALL_ENABLED_TRACKS=false
MERGE_PARALLEL_JOBS=2
WRITE_COMMAND_LOG=false
EOF

PATH="$TMP/bin:$PATH" C2T_CONFIG="$TMP/config.sh" bash "$ROOT/scripts/consensus_batch.sh"
awk -F '\t' 'NR==2 {ok=($1=="C1" && $2=="SUCCESS" && $4==2 && $6==1)} END {exit !ok}' \
    "$TMP/output/05_peaks/consensus/consensus_status.tsv"

if PATH="$TMP/bin:$PATH" MOCK_CONSENSUS_FAIL=true C2T_CONFIG="$TMP/config.sh" \
    bash "$ROOT/scripts/consensus_batch.sh" >/dev/null 2>&1; then
    echo "ERROR: consensus batch did not propagate a worker failure" >&2
    exit 1
fi

echo "Consensus batch exit-status regression tests passed"
