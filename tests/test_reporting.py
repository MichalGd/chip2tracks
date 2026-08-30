#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ReportingTests(unittest.TestCase):
    def test_lightweight_and_custom_reports_recover_from_retained_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "output"
            files = {
                "00_metadata/sample_manifest.tsv": "sample_key\tcohort_id\nTARGET.bioR1\tCOHORT_A\n",
                "00_metadata/cohort_manifest.tsv": (
                    "cohort_id\tfactor\tantibody_id\tprimary_peak_class\tconditions\n"
                    "COHORT_A\tH3K27ac\tAB1\tbroad\tuntreated,treated\n"
                ),
                "00_metadata/resolved_config.tsv": (
                    "key\tvalue\nRUN_DESEQ2_ENRICHMENT\ttrue\nRUN_DIFFBIND\ttrue\n"
                    "RUN_CONTROL_SUBTRACTED_SENSITIVITY\ttrue\n"
                    "RUN_TARGET_CONTROL_INTERACTION\ttrue\nDIFFERENTIAL_ALPHA\t0.05\n"
                    "DIFFERENTIAL_MIN_ABS_LOG2FC\t1\nDIFFERENTIAL_NORMALIZATION\tdeseq2\n"
                ),
                "06_qc/alignment_and_complexity/observation_counts.tsv": (
                    "sample_key\tlayout\tsignal_unit\tanalysis_observations\n"
                    "TARGET.bioR1\tPE\tfragment\t1234\n"
                ),
                "05_peaks/per_sample/peakcall_status.tsv": (
                    "sample_key\tcontrol_key\tprimary_caller\tprimary_class\tstatus\t"
                    "primary_peak_count\tcaller_warnings\treason\n"
                    "TARGET.bioR1\tCTRL.bioR1\tmacs3\tbroad\tEMPTY\t0\tmacs3:broad=EMPTY\t"
                    "primary_caller_produced_no_peaks\n"
                ),
                "05_peaks/consensus/consensus_status.tsv": (
                    "cohort_id\tstatus\ttotal_samples\tsuccessful_peak_samples\texcluded_samples\tregions\treason\n"
                    "COHORT_A\tSKIPPED\t1\t0\t1\t0\tinsufficient successful samples\n"
                ),
                "04_tracks/normalized_track_family_status.tsv": (
                    "cohort_id\tpolicy\tstatus\treason\tlog\n"
                    "COHORT_A\tanalysis\tSKIPPED\tconsensus unavailable\t.\n"
                ),
                "04_tracks/cpm/mapping_composition.tsv": (
                    "sample_key\tpolicy\tcpm_family\tnormalized_family_source\tbam_branch\tmapq_policy\t"
                    "duplicates\tsignal_unit\ttotal_observations\tmapq0_observations\tmapq0_percent\t"
                    "mapq_lt30_observations\tmapq_lt30_percent\txs_tagged_candidate_multimappers\t"
                    "xs_tagged_percent\nTARGET.bioR1\tpermissive\tcpm/permissive\t"
                    "deseq2_robust_cpm/permissive\tq0_dup-retained\t0\tretained\tfragment\t1800\t"
                    "180\t10.0000\t360\t20.0000\t90\t5.0000\n"
                ),
                "04_tracks/cpm/mapping_composition_definitions.tsv": (
                    "metric\tdefinition\nxs_tagged_candidate_multimappers\tBowtie2 XS-tagged signal units.\n"
                ),
                "04_tracks/cpm/TARGET.bioR1.CPM.bw": "fixture\n",
                "06_qc/optional_qc_status.tsv": (
                    "sample_key\tmetric\tstatus\treason\n"
                    "TARGET.bioR1\tcross_correlation\tSKIPPED\tRUN_CROSS_CORRELATION=false\n"
                ),
                "06_qc/frip_and_peak_reproducibility/TARGET.bioR1.sample_primary_frip.tsv": (
                    "sample_key\tsignal_unit\ttotal\tin_sample_primary_peaks\tfrip\n"
                    "TARGET.bioR1\tfragment\t1234\t400\t0.32414911\n"
                ),
                "06_qc/frip_and_peak_reproducibility/TARGET.bioR1.frip.tsv": (
                    "sample_key\tsignal_unit\ttotal\tin_consensus\tfrip\n"
                    "TARGET.bioR1\tfragment\t1234\t321\t0.26012966\n"
                ),
                "06_qc/fragment_length_and_periodicity/TARGET.bioR1.fragment_lengths.tsv": (
                    "fragment_length\tcount\n100\t2\n150\t3\n200\t5\n"
                ),
                "06_qc/fragment_length_and_periodicity/TARGET.bioR1.phantompeak.tsv": (
                    "TARGET.tagAlign.gz\t1234\t150\t0.20\t36\t0.10\t500\t0.05\t1.40\t2.00\t2\n"
                ),
                "08_differential/stage_status.tsv": (
                    "status\tfailed_modules\tskipped_cohorts\nCOMPLETED_WITH_WARNINGS\t0\t1\n"
                ),
                "08_differential/COHORT_A/broad/primary_target_only/deseq2_enrichment/comparison_summary.tsv": (
                    "comparison_id\tnumerator\treference\ttested\tsignificant\thigher_in_numerator\t"
                    "higher_in_reference\tstatus\tall_results\tsignificant_results\n"
                    "treated_vs_untreated\ttreated\tuntreated\t100\t12\t8\t4\tSUCCESS\t"
                    "comparisons/treated_vs_untreated/all.tsv.gz\t"
                    "comparisons/treated_vs_untreated/significant.tsv.gz\n"
                ),
                "08_differential/COHORT_A/broad/primary_target_only/deseq2_enrichment/comparisons/treated_vs_untreated/all.tsv.gz": "fixture\n",
                "08_differential/COHORT_A/broad/primary_target_only/deseq2_enrichment/comparisons/treated_vs_untreated/significant.tsv.gz": "fixture\n",
                "08_differential/COHORT_A/broad/primary_target_only/diffbind/comparison_summary.tsv": (
                    "comparison_id\tnumerator\treference\ttested\tsignificant\thigher_in_numerator\t"
                    "higher_in_reference\tstatus\tall_results\tsignificant_results\n"
                    "treated_vs_untreated\ttreated\tuntreated\t100\t10\t7\t3\tSUCCESS\t"
                    "contrast_1_all.tsv\tcontrast_1_significant.tsv\n"
                ),
                "08_differential/COHORT_A/broad/primary_target_only/diffbind/contrast_1_all.tsv": "fixture\n",
                "08_differential/COHORT_A/broad/primary_target_only/diffbind/contrast_1_significant.tsv": "fixture\n",
                "08_differential/COHORT_A/broad/sensitivity_control_subtracted/diffbind/SKIPPED.json": (
                    '{"status":"SKIPPED","reason":"matched controls unavailable"}\n'
                ),
                "08_differential/COHORT_A/broad/sensitivity_target_control_interaction/deseq2/comparison_summary.tsv": (
                    "comparison_id\tnumerator\treference\ttested\tsignificant\thigher_in_numerator\t"
                    "higher_in_reference\tstatus\tall_results\tsignificant_results\n"
                    "treated_vs_untreated_target_control_interaction\ttreated\tuntreated\t100\t6\t5\t1\t"
                    "SUCCESS\tinteraction_results_all.tsv.gz\tinteraction_results_significant.tsv.gz\n"
                ),
                "08_differential/COHORT_A/broad/sensitivity_target_control_interaction/deseq2/interaction_results_all.tsv.gz": "fixture\n",
                "08_differential/COHORT_A/broad/sensitivity_target_control_interaction/deseq2/interaction_results_significant.tsv.gz": "fixture\n",
            }
            for relative, content in files.items():
                path = output / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(content, encoding="utf-8")

            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts/generate_report.py"), str(output)],
                text=True, capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            report_dir = output / "10_reports"
            for name in ("pipeline_report.html", "run_summary.tsv", "warning_summary.tsv",
                         "qc_module_summary.tsv", "coverage_mapping_composition.tsv",
                         "fragment_qc_summary.tsv", "cross_correlation_summary.tsv",
                         "track_inventory.tsv", "differential_occupancy_summary.tsv"):
                self.assertTrue((report_dir / name).is_file(), name)
                self.assertGreater((report_dir / name).stat().st_size, 0)
            warnings = (report_dir / "warning_summary.tsv").read_text(encoding="utf-8")
            self.assertIn("peakcalling:TARGET.bioR1:EMPTY", warnings)
            self.assertIn("normalization:COHORT_A:analysis:SKIPPED", warnings)
            html_report = (report_dir / "pipeline_report.html").read_text(encoding="utf-8")
            self.assertIn("Mapping composition of coverage families", html_report)
            self.assertIn("XS-tagged candidate multimappers", html_report)
            self.assertIn("QC modules and retained evidence", html_report)
            self.assertIn("FRiP against sample primary peaks", html_report)
            self.assertIn("Strand cross-correlation summary", html_report)
            self.assertIn("1.40", html_report)
            self.assertIn("Complete retained-track inventory", html_report)
            self.assertIn("TARGET.bioR1.CPM.bw", html_report)
            self.assertIn("Differential occupancy analysis", html_report)
            differential_summary = (report_dir / "differential_occupancy_summary.tsv").read_text(encoding="utf-8")
            self.assertIn("sensitivity_control_subtracted", differential_summary)
            self.assertIn("matched controls unavailable", differential_summary)
            self.assertIn("sensitivity_target_control_interaction", differential_summary)

            custom_dir = Path(temporary) / "custom"
            module = load_module(
                "prepare_multiqc_content",
                ROOT / "scripts/prepare_multiqc_content.py",
            )
            original_argv = sys.argv
            try:
                sys.argv = ["prepare_multiqc_content.py", str(output), str(custom_dir)]
                self.assertEqual(module.main(), 0)
            finally:
                sys.argv = original_argv
            manifest = (custom_dir / "custom_content_manifest.tsv").read_text(encoding="utf-8")
            self.assertIn("chip2tracks_observations_mqc.tsv", manifest)
            self.assertIn("chip2tracks_peakcalls_mqc.tsv", manifest)
            self.assertIn("chip2tracks_coverage_mapping_mqc.tsv", manifest)
            self.assertIn("chip2tracks_qc_modules_mqc.tsv", manifest)
            self.assertIn("chip2tracks_fragment_summary_mqc.tsv", manifest)
            self.assertIn("chip2tracks_cross_correlation_mqc.tsv", manifest)
            self.assertIn("chip2tracks_track_inventory_mqc.tsv", manifest)
            self.assertIn("chip2tracks_comparisons_mqc.tsv", manifest)
            custom_peakcalls = (custom_dir / "chip2tracks_peakcalls_mqc.tsv").read_text(encoding="utf-8")
            self.assertIn("# plot_type: table", custom_peakcalls)
            self.assertIn("TARGET.bioR1", custom_peakcalls)


if __name__ == "__main__":
    unittest.main()
