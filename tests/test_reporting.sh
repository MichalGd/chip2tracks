#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEMPORARY="$(mktemp -d)"
trap 'rm -rf -- "$TEMPORARY"' EXIT
OUTPUT="${TEMPORARY}/output"
mkdir -p \
    "$OUTPUT/00_metadata" "$OUTPUT/03_alignment/metrics" \
    "$OUTPUT/04_tracks" "$OUTPUT/05_peaks/per_sample" "$OUTPUT/05_peaks/consensus" \
    "$OUTPUT/06_qc/alignment_and_complexity" "$OUTPUT/06_qc/frip_and_peak_reproducibility" \
    "$OUTPUT/08_differential" "$TEMPORARY/bin"
mkdir -p "$OUTPUT/04_tracks/cpm"
mkdir -p \
    "$OUTPUT/08_differential/COHORT_A/broad/primary_target_only/deseq2_enrichment/comparisons/treated_vs_untreated" \
    "$OUTPUT/08_differential/COHORT_A/broad/primary_target_only/diffbind" \
    "$OUTPUT/08_differential/COHORT_A/broad/sensitivity_control_subtracted/diffbind" \
    "$OUTPUT/08_differential/COHORT_A/broad/sensitivity_target_control_interaction/deseq2"

printf 'sample_key\tcohort_id\nTARGET.bioR1\tCOHORT_A\n' > "$OUTPUT/00_metadata/sample_manifest.tsv"
printf 'cohort_id\tfactor\tantibody_id\tprimary_peak_class\tconditions\nCOHORT_A\tH3K27ac\tAB1\tbroad\tuntreated,treated\n' \
    > "$OUTPUT/00_metadata/cohort_manifest.tsv"
printf 'key\tvalue\nRUN_DESEQ2_ENRICHMENT\ttrue\nRUN_DIFFBIND\ttrue\nRUN_CONTROL_SUBTRACTED_SENSITIVITY\ttrue\nRUN_TARGET_CONTROL_INTERACTION\ttrue\nDIFFERENTIAL_ALPHA\t0.05\nDIFFERENTIAL_MIN_ABS_LOG2FC\t1\nDIFFERENTIAL_NORMALIZATION\tdeseq2\n' \
    > "$OUTPUT/00_metadata/resolved_config.tsv"
printf 'sample_key\tlayout\tsignal_unit\tanalysis_observations\nTARGET.bioR1\tPE\tfragment\t1234\n' \
    > "$OUTPUT/06_qc/alignment_and_complexity/observation_counts.tsv"
printf 'sample_key\tlayout\tgenome\talignment_records\tspikein_mode\nTARGET.bioR1\tPE\thg38\t2000\tnone\n' \
    > "$OUTPUT/03_alignment/metrics/TARGET.bioR1.alignment.tsv"
printf 'sample_key\tq0_dup_retained\tq0_dup_removed\tq30_dup_retained\tq30_dup_removed\tanalysis_policy\nTARGET.bioR1\t1800\t1700\t1300\t1234\tremove\n' \
    > "$OUTPUT/03_alignment/metrics/TARGET.bioR1.filter_counts.tsv"
printf 'sample_key\tcontrol_key\tprimary_caller\tprimary_class\tstatus\tprimary_peak_count\tcaller_warnings\treason\nTARGET.bioR1\tCTRL.bioR1\tmacs3\tbroad\tSUCCESS\t42\t.\t.\n' \
    > "$OUTPUT/05_peaks/per_sample/peakcall_status.tsv"
printf 'cohort_id\tstatus\ttotal_samples\tsuccessful_peak_samples\texcluded_samples\tregions\treason\nCOHORT_A\tSUCCESS\t1\t1\t0\t42\t.\n' \
    > "$OUTPUT/05_peaks/consensus/consensus_status.tsv"
printf 'cohort_id\tpolicy\tstatus\treason\tlog\nCOHORT_A\tanalysis\tSUCCESS\t.\t.\n' \
    > "$OUTPUT/04_tracks/normalized_track_family_status.tsv"
printf 'sample_key\tpolicy\tcpm_family\tnormalized_family_source\tbam_branch\tmapq_policy\tduplicates\tsignal_unit\ttotal_observations\tmapq0_observations\tmapq0_percent\tmapq_lt30_observations\tmapq_lt30_percent\txs_tagged_candidate_multimappers\txs_tagged_percent\nTARGET.bioR1\tpermissive\tcpm/permissive\tdeseq2_robust_cpm/permissive\tq0_dup-retained\t0\tretained\tfragment\t1800\t180\t10.0000\t360\t20.0000\t90\t5.0000\n' \
    > "$OUTPUT/04_tracks/cpm/mapping_composition.tsv"
printf 'status\tfailed_modules\tskipped_cohorts\nSUCCESS\t0\t0\n' > "$OUTPUT/08_differential/stage_status.tsv"
summary_header='comparison_id\tnumerator\treference\ttested\tsignificant\thigher_in_numerator\thigher_in_reference\tstatus\tall_results\tsignificant_results'
printf "%b\ntreated_vs_untreated\ttreated\tuntreated\t100\t12\t8\t4\tSUCCESS\tcomparisons/treated_vs_untreated/all.tsv.gz\tcomparisons/treated_vs_untreated/significant.tsv.gz\n" "$summary_header" \
    > "$OUTPUT/08_differential/COHORT_A/broad/primary_target_only/deseq2_enrichment/comparison_summary.tsv"
printf 'fixture\n' > "$OUTPUT/08_differential/COHORT_A/broad/primary_target_only/deseq2_enrichment/comparisons/treated_vs_untreated/all.tsv.gz"
printf 'fixture\n' > "$OUTPUT/08_differential/COHORT_A/broad/primary_target_only/deseq2_enrichment/comparisons/treated_vs_untreated/significant.tsv.gz"
printf "%b\ntreated_vs_untreated\ttreated\tuntreated\t100\t10\t7\t3\tSUCCESS\tcontrast_1_all.tsv\tcontrast_1_significant.tsv\n" "$summary_header" \
    > "$OUTPUT/08_differential/COHORT_A/broad/primary_target_only/diffbind/comparison_summary.tsv"
printf 'fixture\n' > "$OUTPUT/08_differential/COHORT_A/broad/primary_target_only/diffbind/contrast_1_all.tsv"
printf 'fixture\n' > "$OUTPUT/08_differential/COHORT_A/broad/primary_target_only/diffbind/contrast_1_significant.tsv"
printf '{"status":"SKIPPED","reason":"matched controls unavailable"}\n' \
    > "$OUTPUT/08_differential/COHORT_A/broad/sensitivity_control_subtracted/diffbind/SKIPPED.json"
printf "%b\ntreated_vs_untreated_target_control_interaction\ttreated\tuntreated\t100\t6\t5\t1\tSUCCESS\tinteraction_results_all.tsv.gz\tinteraction_results_significant.tsv.gz\n" "$summary_header" \
    > "$OUTPUT/08_differential/COHORT_A/broad/sensitivity_target_control_interaction/deseq2/comparison_summary.tsv"
printf 'fixture\n' > "$OUTPUT/08_differential/COHORT_A/broad/sensitivity_target_control_interaction/deseq2/interaction_results_all.tsv.gz"
printf 'fixture\n' > "$OUTPUT/08_differential/COHORT_A/broad/sensitivity_target_control_interaction/deseq2/interaction_results_significant.tsv.gz"
printf 'sample_key\tsignal_unit\ttotal\tin_consensus\tfrip\nTARGET.bioR1\tfragment\t1234\t321\t0.26012966\n' \
    > "$OUTPUT/06_qc/frip_and_peak_reproducibility/TARGET.bioR1.frip.tsv"

cat > "$TEMPORARY/bin/multiqc" <<'FAKE'
#!/usr/bin/env bash
set -euo pipefail
out=""; name=""
while (( $# )); do
    case "$1" in
        --outdir) out="$2"; shift 2 ;;
        --filename) name="$2"; shift 2 ;;
        --exclude|--ignore|--cl-config|--data-format) shift 2 ;;
        --export|--force) shift ;;
        *) shift ;;
    esac
done
[[ -n "$out" && -n "$name" ]]
mkdir -p "$out/${name}_data"
printf '<html><body>fake MultiQC report</body></html>\n' > "$out/${name}.html"
printf 'sample\tvalue\nTARGET.bioR1\t1\n' > "$out/${name}_data/multiqc_general_stats.txt"
echo "mqc_colour | Error converting color '55,126,184' to RGB: input #55,126,184 is not in #RRGGBB format" >&2
echo 'MultiQC complete'
FAKE
chmod +x "$TEMPORARY/bin/multiqc"

before="$(sha256sum "$OUTPUT/00_metadata/sample_manifest.tsv")"
PATH="$TEMPORARY/bin:$PATH" bash "$ROOT/utilities/regenerate_reports.sh" --output-dir "$OUTPUT"
after="$(sha256sum "$OUTPUT/00_metadata/sample_manifest.tsv")"
[[ "$before" == "$after" ]]

for required in \
    pipeline_report.html run_summary.tsv warning_summary.tsv coverage_mapping_composition.tsv \
    differential_occupancy_summary.tsv \
    chip2tracks_multiqc_report.html multiqc_status.tsv \
    multiqc_custom_content_manifest.tsv report_checksums.sha256; do
    [[ -s "$OUTPUT/10_reports/$required" ]]
done
[[ -d "$OUTPUT/10_reports/chip2tracks_multiqc_report_data" ]]
grep -q $'SUCCESS\t' "$OUTPUT/10_reports/multiqc_status.tsv"
grep -q 'chip2tracks_observations' "$OUTPUT/10_reports/multiqc_custom_content_manifest.tsv"
grep -q 'chip2tracks_coverage_mapping' "$OUTPUT/10_reports/multiqc_custom_content_manifest.tsv"
grep -q 'chip2tracks_comparisons' "$OUTPUT/10_reports/multiqc_custom_content_manifest.tsv"
grep -q 'Retained analysis observations' "$OUTPUT/10_reports/pipeline_report.html"
grep -q 'Mapping composition of coverage families' "$OUTPUT/10_reports/pipeline_report.html"
grep -q 'XS-tagged candidate multimappers' "$OUTPUT/10_reports/pipeline_report.html"
grep -q 'Differential occupancy analysis' "$OUTPUT/10_reports/pipeline_report.html"
grep -q 'sensitivity_control_subtracted' "$OUTPUT/10_reports/differential_occupancy_summary.tsv"
grep -q 'matched controls unavailable' "$OUTPUT/10_reports/differential_occupancy_summary.tsv"
grep -q 'sensitivity_target_control_interaction' "$OUTPUT/10_reports/differential_occupancy_summary.tsv"

echo "Unified reporting and completed-run regeneration test passed"
