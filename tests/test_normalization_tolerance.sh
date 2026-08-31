#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf -- "$TMP"' EXIT

mkdir -p "$TMP/bin" "$TMP/output/00_metadata"

cat > "$TMP/bin/Rscript" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
script="$1"; manifest="$2"; cohort="$3"; consensus="$4"; output="$5"; tables="$6"; policy="$7"
if [[ "$script" == *deseq2_enrichment_analysis.R ]]; then
    exit 0
fi
mkdir -p "$tables"
if [[ "$cohort" == "C_FAIL" ]]; then
    printf 'sample_key\tcohort_id\tpolicy\tconsensus_regions_with_any_count\tconsensus_count_sum\tstatus\n' \
        > "$tables/consensus_count_sums.tsv"
    printf 'FAIL_A.bioR1\tC_FAIL\tanalysis\t1\t0\tZERO\n' >> "$tables/consensus_count_sums.tsv"
    printf 'FAIL_B.bioR2\tC_FAIL\tanalysis\t1\t4\tNONZERO\n' >> "$tables/consensus_count_sums.tsv"
    echo 'zero consensus counts for samples: FAIL_A.bioR1' >&2
    exit 42
fi
printf 'sample_key\tcohort_id\tpolicy\tsignal_unit\tsize_factor\tdeseq2_consensus_scale\tconsensus_count_sum\tcohort_geometric_mean_column_sum\trobust_effective_library_size\tdeseq2_robust_cpm_scale\n' \
    > "$tables/normalization_factors.tsv"
printf 'OK_A.bioR1\tC_OK\tanalysis\tfragment\t1\t1\t10\t10\t10\t100000\n' >> "$tables/normalization_factors.tsv"
printf 'OK_B.bioR2\tC_OK\tanalysis\tfragment\t1\t1\t10\t10\t10\t100000\n' >> "$tables/normalization_factors.tsv"
printf 'region_id\tOK_A.bioR1\tOK_B.bioR2\nregion1\t10\t10\n' > "$tables/raw_counts.tsv.gz"
MOCK
chmod +x "$TMP/bin/Rscript"

header='sample_key\tis_control\tcohort_id'
printf '%b\n' "$header" > "$TMP/output/00_metadata/sample_manifest.tsv"
printf 'FAIL_A.bioR1\tFALSE\tC_FAIL\nFAIL_B.bioR2\tFALSE\tC_FAIL\nOK_A.bioR1\tFALSE\tC_OK\nOK_B.bioR2\tFALSE\tC_OK\n' \
    >> "$TMP/output/00_metadata/sample_manifest.tsv"

cohort_header='cohort_id\tcohort_key\tgenome\tassay_profile\tfactor\tantibody_id\tlayout\ttarget_class\tanalysis_duplicate_policy\tprimary_peak_caller\tprimary_peak_class\tn_biological_samples\tsample_keys\tconditions'
printf '%b\n' "$cohort_header" > "$TMP/output/00_metadata/cohort_manifest.tsv"
printf 'C_FAIL\tkey1\thg38\tchipseq\tF1\tA1\tPE\tbroad\tretain\tmacs3\tbroad\t2\tFAIL_A.bioR1,FAIL_B.bioR2\tA,B\n' \
    >> "$TMP/output/00_metadata/cohort_manifest.tsv"
printf 'C_OK\tkey2\thg38\tchipseq\tF2\tA2\tPE\tbroad\tretain\tmacs3\tbroad\t2\tOK_A.bioR1,OK_B.bioR2\tA,B\n' \
    >> "$TMP/output/00_metadata/cohort_manifest.tsv"

for cohort in C_FAIL C_OK; do
    directory="$TMP/output/05_peaks/consensus/$cohort/macs3/broad"
    mkdir -p "$directory"
    printf 'chr1\t0\t100\t%s.region1\t2\n' "$cohort" > "$directory/$cohort.consensus.bed"
done

cat > "$TMP/config.sh" <<EOF
OUTPUT_DIR=$TMP/output
REQUIRE_ALL_ENABLED_TRACKS=false
GENERATE_DESEQ2_CONSENSUS_TRACKS=false
RUN_DIFFBIND=false
RUN_DESEQ2_ENRICHMENT=true
GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS=false
GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS=false
GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS=false
RUN_CONTROL_SUBTRACTED_SENSITIVITY=false
RUN_TARGET_CONTROL_INTERACTION=false
DIFFERENTIAL_NORMALIZATION=internal
DIFFERENTIAL_MIN_REPLICATES_PER_CONDITION=2
DIFFERENTIAL_ALPHA=0.05
DIFFERENTIAL_MIN_ABS_LOG2FC=0
DIFFERENTIAL_BLOCK_COLUMNS=
DIFFERENTIAL_CONDITION_ORDER=
DIFFERENTIAL_REFERENCE_CONDITION=
WRITE_COMMAND_LOG=false
EOF

PATH="$TMP/bin:$PATH" C2T_CONFIG="$TMP/config.sh" bash "$ROOT/scripts/normalized_tracks_batch.sh"
grep -q '^COMPLETED_WITH_WARNINGS' "$TMP/output/04_tracks/stage_status.tsv"
grep -q 'zero consensus counts for samples: FAIL_A.bioR1' \
    "$TMP/output/logs/normalized_tracks/C_FAIL.analysis.factors.log"
test -s "$TMP/output/04_tracks/deseq2_consensus/C_FAIL/SKIPPED.json"
test -s "$TMP/output/04_tracks/deseq2_consensus/C_FAIL/tables/consensus_count_sums.tsv"
test -s "$TMP/output/04_tracks/deseq2_consensus/C_OK/tables/raw_counts.tsv.gz"
awk -F '\t' '$1=="C_FAIL" {found=($3=="SKIPPED")} END {exit !found}' \
    "$TMP/output/04_tracks/normalized_track_family_status.tsv"
awk -F '\t' '$1=="C_OK" {found=($3=="SUCCESS")} END {exit !found}' \
    "$TMP/output/04_tracks/normalized_track_family_status.tsv"

PATH="$TMP/bin:$PATH" C2T_CONFIG="$TMP/config.sh" bash "$ROOT/scripts/differential_batch.sh"
grep -q '^COMPLETED_WITH_WARNINGS' "$TMP/output/08_differential/stage_status.tsv"
test -s "$TMP/output/08_differential/C_FAIL/broad/SKIPPED.json"
test ! -e "$TMP/output/08_differential/C_FAIL/broad/FAILED.json"

disabled_output="$TMP/disabled_output"
disabled_cohort="C_DISABLED"
disabled_root="$disabled_output/08_differential/$disabled_cohort/broad"
mkdir -p "$disabled_output/00_metadata" "$disabled_root" \
    "$disabled_output/05_peaks/consensus/$disabled_cohort/epic2/broad"
printf 'sample_key\tis_control\tcohort_id\nS1.bioR1\tFALSE\t%s\n' "$disabled_cohort" \
    > "$disabled_output/00_metadata/sample_manifest.tsv"
printf '%b\n' "$cohort_header" > "$disabled_output/00_metadata/cohort_manifest.tsv"
printf '%s\tkey\tmm39\tchipmentation\tF3\tA3\tPE\tmixed\tremove\tepic2\tbroad\t1\tS1.bioR1\tbas\n' \
    "$disabled_cohort" >> "$disabled_output/00_metadata/cohort_manifest.tsv"
printf 'chr1\t10\t30\n' \
    > "$disabled_output/05_peaks/consensus/$disabled_cohort/epic2/broad/$disabled_cohort.consensus.bed"
cat > "$TMP/config.disabled.sh" <<EOF
OUTPUT_DIR=$disabled_output
RUN_DIFFBIND=false
RUN_DESEQ2_ENRICHMENT=false
RUN_TARGET_CONTROL_INTERACTION=false
RUN_CONTROL_SUBTRACTED_SENSITIVITY=false
DIFFERENTIAL_PARALLEL_JOBS=2
WRITE_COMMAND_LOG=false
EOF

# Reproduce the 0.2.1 false failure and verify that a fully disabled
# differential stage is a successful no-op that removes the stale marker.
printf '{"status":"FAILED","reason":"raw consensus counts unexpectedly unavailable"}\n' \
    > "$disabled_root/FAILED.json"
C2T_CONFIG="$TMP/config.disabled.sh" bash "$ROOT/scripts/differential_batch.sh"
awk -F '\t' 'NR==2 {ok=($1=="SUCCESS" && $2==0 && $3==0 && $4==2)} END {exit !ok}' \
    "$disabled_output/08_differential/stage_status.tsv"
test ! -e "$disabled_root/FAILED.json"
test ! -e "$disabled_root/SKIPPED.json"
test ! -e "$disabled_output/04_tracks/deseq2_consensus/$disabled_cohort/tables/raw_counts.tsv.gz"

# The no-op guard must not weaken validation when an analysis requiring those
# count tables is enabled.
sed 's/RUN_DESEQ2_ENRICHMENT=false/RUN_DESEQ2_ENRICHMENT=true/' \
    "$TMP/config.disabled.sh" > "$TMP/config.enabled-missing-counts.sh"
if C2T_CONFIG="$TMP/config.enabled-missing-counts.sh" \
        bash "$ROOT/scripts/differential_batch.sh" >/dev/null 2>&1; then
    echo "ERROR: enabled differential analysis accepted missing raw counts" >&2
    exit 1
fi
grep -q '^FAILED' "$disabled_output/08_differential/stage_status.tsv"
grep -q 'raw consensus counts unexpectedly unavailable' "$disabled_root/FAILED.json"

sed 's/REQUIRE_ALL_ENABLED_TRACKS=false/REQUIRE_ALL_ENABLED_TRACKS=true/' "$TMP/config.sh" > "$TMP/config.strict.sh"
if PATH="$TMP/bin:$PATH" C2T_CONFIG="$TMP/config.strict.sh" \
    bash "$ROOT/scripts/normalized_tracks_batch.sh" >/dev/null 2>&1; then
    echo "ERROR: strict normalized-track policy accepted a zero-count cohort" >&2
    exit 1
fi
grep -q '^FAILED' "$TMP/output/04_tracks/stage_status.tsv"

echo "Normalized-track cohort continuation and strict-policy regression tests passed"
