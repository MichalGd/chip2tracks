#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
temporary="$(mktemp -d "$ROOT/.recovery-test.XXXXXX")"
portable_temp="$temporary"
if command -v cygpath >/dev/null 2>&1; then portable_temp="$(cygpath -m "$temporary")"; fi
trap 'rm -rf -- "$temporary"' EXIT
trap 'if [[ -f "$temporary/failed.log" ]]; then tail -n 30 "$temporary/failed.log"; fi' ERR
mkdir -p "$temporary/bin" "$temporary/output/00_metadata"
export TEST_CALLS="$temporary/trim_calls" TEST_AGGREGATE_FAIL="$temporary/fail" TEST_STARTUP_FAIL="$temporary/broken"
cat > "$temporary/bin/trim_galore" <<'TOOL'
#!/usr/bin/env bash
set -euo pipefail
[[ "${1:-}" != "--version" ]] || { echo 'Trim Galore test'; exit 0; }
inputs=()
paired=false
while (( $# )); do
    case "$1" in
        --output_dir) out="$2"; shift 2;;
        --cores|--length) shift 2;;
        --illumina) shift;;
        --paired) paired=true; shift;;
        *) inputs+=("$1"); shift;;
    esac
done
printf '%s\n' "${inputs[*]}" >> "$TEST_CALLS"
for i in "${!inputs[@]}"; do
    input="${inputs[$i]}"
    base="$(basename "$input" .fastq.gz)"
    if [[ "$paired" == true ]]; then
        cp "$input" "$out/${base}_val_$((i+1)).fq.gz"
    else
        cp "$input" "$out/${base}_trimmed.fq.gz"
    fi
done
TOOL
cat > "$temporary/bin/multiqc" <<'TOOL'
#!/usr/bin/env bash
set -euo pipefail
if [[ "${1:-}" == "--version" ]]; then
    [[ ! -f "$TEST_STARTUP_FAIL" ]] || exit 13
    echo 'MultiQC test'; exit 0
fi
[[ ! -f "$TEST_AGGREGATE_FAIL" ]] || exit 11
while (( $# )); do
    case "$1" in --outdir) out="$2"; shift 2;; *) shift;; esac
done
printf '<html>completed aggregation</html>\n' > "$out/multiqc_report.html"
TOOL
chmod +x "$temporary/bin/trim_galore" "$temporary/bin/multiqc"
export PATH="$temporary/bin:$PATH"
printf '@read1\nACGT\n+\nIIII\n' | gzip > "$temporary/A.fastq.gz"
cp "$temporary/A.fastq.gz" "$temporary/B.fastq.gz"
cp "$temporary/A.fastq.gz" "$temporary/B2.fastq.gz"
manifest="$temporary/output/00_metadata/sample_manifest.tsv"
printf '%s\n' $'sample_key\tsample_id\treplicate\tlayout\tgenome\tassay_profile\tfactor\tantibody_id\ttarget_class\tcondition\ttreatment\tcell_type\tis_control\tcontrol_type\tcontrol_id\tcontrol_key\tanalysis_duplicate_policy\tblacklist\tspikein_to_host_ratio\tspikein_stage\tspikein_lot\tbatch\tdonor\toutput_prefix\ttechnical_units\tfastq_1_list\tfastq_2_list\tcohort_id\tcohort_key\tprimary_peak_caller\tprimary_peak_class' > "$manifest"
for sample in A B; do
    layout=SE; assay=chipseq; r2=.
    if [[ "$sample" == B ]]; then layout=PE; assay=chipmentation; r2="$portable_temp/B2.fastq.gz"; fi
    printf '%s\t%s\t1\t%s\thg38\t%s\tF\tAB\tnarrow\tWT\tnone\tcells\tFALSE\tnone\t.\t.\tremove\t.\t.\t.\t.\tB1\tD1\t%s\t1\t%s\t%s\tC1\tkey\tmacs3\tnarrow\n' \
        "$sample.bioR1" "$sample" "$layout" "$assay" "$sample" "$portable_temp/$sample.fastq.gz" "$r2" >> "$manifest"
done
printf 'kind\tsha256\n' > "$temporary/output/00_metadata/reference_manifest.tsv"
config="$temporary/config.conf"
cat > "$config" <<EOF
OUTPUT_DIR=$portable_temp/output
RUN_ID=first
RUN_FASTQC=false
RUN_FASTQC_PER_TECHNICAL_UNIT=false
RUN_MULTIQC=true
TRIM_ADAPTERS=true
ADAPTER_PRESET=illumina
TRIMGALORE_EXTRA_ARGS=
ADAPTER_AUTO_DETECTION_REQUIRED=false
MIN_TRIMMED_LENGTH=20
THREADS_FASTQC=1
THREADS_TRIMGALORE=1
QC_SAMPLE_PARALLEL_JOBS=2
CHECKPOINT_PARALLEL_JOBS=1
WRITE_COMMAND_LOG=false
EOF
python3 "$ROOT/scripts/recovery.py" prepare --root "$ROOT" --config "$config" \
    --output "$temporary/output/00_metadata/recovery_manifest.json"
touch "$TEST_AGGREGATE_FAIL"
if C2T_CONFIG="$config" bash "$ROOT/scripts/preprocess_batch.sh" > "$temporary/failed.log" 2>&1; then
    echo 'ERROR: failed aggregation was accepted'; exit 1
fi
[[ "$(wc -l < "$TEST_CALLS")" -eq 2 ]]
test -s "$temporary/output/.checkpoints/preprocess_samples/A.bioR1.json"
test -s "$temporary/output/.checkpoints/preprocess_samples/B.bioR1.json"
test ! -e "$temporary/output/02_trimmed_fastq/stage_status.tsv"
sha256sum "$temporary/output/02_trimmed_fastq/"*.trimmed.fastq.gz > "$temporary/before.sha256"
rm "$TEST_AGGREGATE_FAIL"
C2T_CONFIG="$config" bash "$ROOT/scripts/preprocess_batch.sh" > "$temporary/resumed.log" 2>&1
[[ "$(wc -l < "$TEST_CALLS")" -eq 2 ]]
[[ "$(grep -c 'Reusing completed preprocessing library' "$temporary/resumed.log")" -eq 2 ]]
sha256sum -c "$temporary/before.sha256"
test -s "$temporary/output/01_fastq_qc/multiqc/multiqc_report.html"
# Corruption invalidates just the affected library rather than the entire batch.
printf 'broken gzip\n' > "$temporary/output/02_trimmed_fastq/A.bioR1.R1.trimmed.fastq.gz"
C2T_CONFIG="$config" bash "$ROOT/scripts/preprocess_batch.sh" > "$temporary/repaired.log" 2>&1
[[ "$(wc -l < "$TEST_CALLS")" -eq 3 ]]
sha256sum -c "$temporary/before.sha256"
# Startup failure must be detected before any additional library work.
touch "$TEST_STARTUP_FAIL"
if C2T_CONFIG="$config" bash "$ROOT/scripts/preprocess_batch.sh" > "$temporary/startup.log" 2>&1; then
    echo 'ERROR: broken MultiQC startup was accepted'; exit 1
fi
[[ "$(wc -l < "$TEST_CALLS")" -eq 3 ]]
grep -q 'MultiQC startup failed' "$temporary/startup.log"
echo 'Preprocessing aggregation-failure recovery regression passed'
