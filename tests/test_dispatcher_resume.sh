#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
temporary="$(mktemp -d "$ROOT/.dispatcher-test.XXXXXX")"
portable_temp="$temporary"
if command -v cygpath >/dev/null 2>&1; then portable_temp="$(cygpath -m "$temporary")"; fi
trap 'rm -rf -- "$temporary"' EXIT
trap 'for log in "$temporary/"*.log; do [[ ! -f "$log" ]] || tail -n 20 "$log"; done' ERR
workflow="$temporary/workflow"
mkdir -p "$workflow/config" "$temporary/output"
cp "$ROOT/chip2tracks.sh" "$ROOT/VERSION" "$workflow/"
cp "$ROOT/config/config.conf.template" "$workflow/config/"
cp -R "$ROOT/scripts" "$ROOT/common" "$workflow/"
# Use tiny deterministic stage commands to exercise the actual dispatcher,
# configuration validation, reference manifests, contracts and checkpoints.
cat > "$workflow/scripts/preflight.sh" <<'STAGE'
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
require_config
echo preflight >> "$OUTPUT_DIR/stub_calls"
printf 'status\nSUCCESS\n' > "$OUTPUT_DIR/00_metadata/preflight_status.tsv"
printf 'tool\tversion\nfixture\t1\n' > "$OUTPUT_DIR/00_metadata/software_versions.tsv"
printf 'stage\tstatus\nfixture\tPASS\n' > "$OUTPUT_DIR/00_metadata/resource_budget.tsv"
STAGE
cat > "$workflow/scripts/preprocess_batch.sh" <<'STAGE'
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
require_config
echo preprocess >> "$OUTPUT_DIR/stub_calls"
mkdir -p "$OUTPUT_DIR/02_trimmed_fastq" "$OUTPUT_DIR/01_fastq_qc"
raw="$(awk -F, 'NR==2 {print $2}' "$SAMPLESHEET")"
cp "$raw" "$OUTPUT_DIR/02_trimmed_fastq/S.bioR1.R1.merged.fastq.gz"
cp "$raw" "$OUTPUT_DIR/02_trimmed_fastq/S.bioR1.R1.trimmed.fastq.gz"
echo qc > "$OUTPUT_DIR/01_fastq_qc/fixture.txt"
printf 'sample_key\tevidence\nS.bioR1\tnot_trimmed\n' > "$OUTPUT_DIR/01_fastq_qc/adapter_detection.tsv"
# Keep the preprocessing-worker source boundary used by recovery.py.
merge_fastqs() { :; }
parallel_require_positive_integer() { :; }
STAGE
cat > "$workflow/scripts/align_batch.sh" <<'STAGE'
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
require_config
echo alignment >> "$OUTPUT_DIR/stub_calls"
mkdir -p "$OUTPUT_DIR/03_alignment/sorted" "$OUTPUT_DIR/03_alignment/spikein"
mkdir -p "$OUTPUT_DIR/03_alignment/metrics"
printf 'seed=%s\n' "$BOWTIE2_SEED" > "$OUTPUT_DIR/03_alignment/sorted/S.bioR1.host.sorted.bam"
printf 'index\n' > "$OUTPUT_DIR/03_alignment/sorted/S.bioR1.host.sorted.bam.bai"
echo metrics > "$OUTPUT_DIR/03_alignment/metrics/S.bioR1.alignment.tsv"
echo SUCCESS > "$OUTPUT_DIR/03_alignment/sorted/stage_status.tsv"
STAGE
cat > "$workflow/scripts/mark_filter_batch.sh" <<'STAGE'
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
require_config
echo filtering >> "$OUTPUT_DIR/stub_calls"
for name in analysis filtered marked metrics; do
    mkdir -p "$OUTPUT_DIR/03_alignment/$name"
done
cp "$OUTPUT_DIR/03_alignment/sorted/S.bioR1.host.sorted.bam" "$OUTPUT_DIR/03_alignment/marked/S.bioR1.host.marked.bam"
echo index > "$OUTPUT_DIR/03_alignment/marked/S.bioR1.host.marked.bam.bai"
echo metrics > "$OUTPUT_DIR/03_alignment/metrics/S.bioR1.duplicate_metrics.txt"
echo counts > "$OUTPUT_DIR/03_alignment/metrics/S.bioR1.filter_counts.tsv"
for branch in q0_dup-retained q0_dup-removed q30_dup-retained q30_dup-removed; do
    mkdir -p "$OUTPUT_DIR/03_alignment/filtered/$branch"
    quality="${branch%%_*}"; policy="${branch#*_}"
    cp "$OUTPUT_DIR/03_alignment/marked/S.bioR1.host.marked.bam" "$OUTPUT_DIR/03_alignment/filtered/$branch/S.bioR1.host.$quality.$policy.bam"
    echo index > "$OUTPUT_DIR/03_alignment/filtered/$branch/S.bioR1.host.$quality.$policy.bam.bai"
done
cp "$OUTPUT_DIR/03_alignment/marked/S.bioR1.host.marked.bam" "$OUTPUT_DIR/03_alignment/analysis/S.bioR1.host.analysis.bam"
echo index > "$OUTPUT_DIR/03_alignment/analysis/S.bioR1.host.analysis.bam.bai"
echo SUCCESS > "$OUTPUT_DIR/03_alignment/analysis/stage_status.tsv"
STAGE
printf 'reference\n' > "$temporary/ref"
printf 'index\n' > "$temporary/index.1.bt2"
printf '@read\nACGT\n+\nIIII\n' | gzip > "$temporary/raw.fastq.gz"
head -n 1 "$ROOT/config/samplesheet_template.csv" > "$temporary/samplesheet.csv"
printf 'S,%s,,SE,hg38,chipseq,F,AB,narrow,WT,none,cells,1,1,FALSE,none,,remove,,,,B1,D1,S\n' \
    "$portable_temp/raw.fastq.gz" >> "$temporary/samplesheet.csv"
config="$temporary/config.conf"
sed -E -e "s|^SAMPLESHEET=.*|SAMPLESHEET=$portable_temp/samplesheet.csv|" \
    -e "s|^OUTPUT_DIR=.*|OUTPUT_DIR=$portable_temp/output|" \
    -e 's|^RUN_MULTIQC=.*|RUN_MULTIQC=false|' \
    -e 's|^RUN_FASTQC=.*|RUN_FASTQC=false|' \
    -e 's|^TRIM_ADAPTERS=.*|TRIM_ADAPTERS=false|' \
    -e 's|^ALLOW_CONTROL_FREE_PEAKCALL=.*|ALLOW_CONTROL_FREE_PEAKCALL=true|' \
    -e "s|^INDEX_HG38=.*|INDEX_HG38=$portable_temp/index|" \
    -e "s@^(FASTA|CHROM_SIZES|CANONICAL_CONTIGS|GTF|BLACKLIST)_HG38=.*@\1_HG38=$portable_temp/ref@" \
    -e 's|^CCRE_BED_HG38=.*|CCRE_BED_HG38=|' \
    "$ROOT/config/config.conf.template" > "$config"
run_fixture() { bash "$workflow/chip2tracks.sh" --config "$config" "$@"; }
run_fixture --stop-after filtering > "$temporary/first.log" 2>&1 || { cat "$temporary/first.log"; exit 1; }
[[ "$(wc -l < "$temporary/output/stub_calls")" -eq 4 ]]
run_fixture --stop-after filtering > "$temporary/restart.log" 2>&1
[[ "$(wc -l < "$temporary/output/stub_calls")" -eq 4 ]]
grep -q '\[alignment\] SKIPPED' "$temporary/restart.log"
# A downstream-only source correction leaves all four earlier stages reusable.
echo '# downstream hotfix' >> "$workflow/scripts/differential_batch.sh"
run_fixture --stop-after filtering > "$temporary/hotfix.log" 2>&1
[[ "$(wc -l < "$temporary/output/stub_calls")" -eq 4 ]]
# Changed parameters cannot be hidden by --from-stage adoption.
sed -i 's/^BOWTIE2_SEED=.*/BOWTIE2_SEED=42/' "$config"
if run_fixture --from-stage filtering --stop-after filtering > "$temporary/unsafe.log" 2>&1; then
    echo 'ERROR: incompatible alignment was adopted'; exit 1
fi
grep -q 'cannot reuse alignment' "$temporary/unsafe.log"
run_fixture --stop-after filtering > "$temporary/changed.log" 2>&1
[[ "$(grep -c '^preprocess$' "$temporary/output/stub_calls")" -eq 1 ]]
[[ "$(grep -c '^alignment$' "$temporary/output/stub_calls")" -eq 2 ]]
[[ "$(grep -c '^filtering$' "$temporary/output/stub_calls")" -eq 2 ]]
grep -q 'seed=42' "$temporary/output/03_alignment/analysis/S.bioR1.host.analysis.bam"
# Rebuilding a corrupt upstream output must also rebuild downstream results.
echo corrupt > "$temporary/output/03_alignment/sorted/S.bioR1.host.sorted.bam"
run_fixture --stop-after filtering > "$temporary/corrupt.log" 2>&1
[[ "$(grep -c '^alignment$' "$temporary/output/stub_calls")" -eq 3 ]]
[[ "$(grep -c '^filtering$' "$temporary/output/stub_calls")" -eq 3 ]]
# Missing required preprocessing QC is not accepted by a valid FASTQ checkpoint.
rm "$temporary/output/01_fastq_qc/fixture.txt"
run_fixture --stop-after filtering > "$temporary/missing-qc.log" 2>&1
[[ "$(grep -c '^preprocess$' "$temporary/output/stub_calls")" -eq 2 ]]
echo 'Dispatcher restart, hotfix compatibility and invalidation regressions passed'
