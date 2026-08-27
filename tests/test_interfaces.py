#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class InterfaceTests(unittest.TestCase):
    def test_chip_specific_interfaces_are_explicit(self) -> None:
        preprocess = (ROOT / "scripts/preprocess_batch.sh").read_text(encoding="utf-8")
        peakcall = (ROOT / "scripts/peakcall_batch.sh").read_text(encoding="utf-8")
        reproducibility = (ROOT / "scripts/reproducibility_batch.sh").read_text(encoding="utf-8")
        self.assertIn('ADAPTER_PRESET', preprocess)
        self.assertIn('adapter_detection.tsv', preprocess)
        self.assertIn('call_epic2', peakcall)
        self.assertIn('IDR_MACS3_PVALUE', peakcall)
        self.assertIn('.macs3.idr_peaks.narrowPeak', reproducibility)
        self.assertNotIn('SEACR_COMMAND', peakcall)

    def test_phantompeak_command_owns_its_r_runtime(self) -> None:
        qc = (ROOT / "scripts/qc_batch.sh").read_text(encoding="utf-8")
        self.assertIn('"$phantompeak_path" -c="$tagalign"', qc)
        self.assertNotIn('Rscript "$phantompeak_path"', qc)

        main_environment = (ROOT / "environment.yml").read_text(encoding="utf-8")
        self.assertNotIn("  - phantompeakqualtools", main_environment)
        self.assertNotIn("  - preseq", main_environment)
        self.assertTrue((ROOT / "environment.spp.yml").is_file())
        self.assertTrue((ROOT / "environment.preseq.yml").is_file())

        sidecar = (ROOT / "utilities/run_phantompeak_sidecar.sh").read_text(encoding="utf-8")
        launcher = (ROOT / "utilities/chip2tracks_shared_launcher.sh").read_text(encoding="utf-8")
        self.assertIn('exec "$RSCRIPT" "$RUN_SPP" "$@"', sidecar)
        self.assertIn('exec "$MAIN_ENV/bin/bash" "$WORKFLOW_ROOT/chip2tracks.sh" "$@"', launcher)

    def test_final_report_runs_unified_multiqc_and_validates_outputs(self) -> None:
        report = (ROOT / "scripts/report_batch.sh").read_text(encoding="utf-8")
        unified = (ROOT / "scripts/generate_multiqc_report.sh").read_text(encoding="utf-8")
        recovery = (ROOT / "utilities/regenerate_reports.sh").read_text(encoding="utf-8")
        self.assertIn('generate_multiqc_report.sh', report)
        self.assertIn('chip2tracks_multiqc_report.html', report)
        self.assertIn('--exclude deeptools', unified)
        self.assertIn('multiqc_custom_content_manifest.tsv', unified)
        self.assertIn('report_checksums.sha256', recovery)

    def test_filtering_uses_indexed_marked_bam_for_region_selection(self) -> None:
        script = (ROOT / "scripts/mark_filter_batch.sh").read_text(encoding="utf-8")
        self.assertNotIn('"$tmp/flags.bam"', script)
        self.assertIn('"$marked" "${contigs[@]}"', script)
        self.assertIn("Reusing validated marked BAM", script)
        self.assertNotIn("trap 'rm -rf \"$tmp\"' RETURN", script)
        self.assertIn("trap 'rm -rf -- \"$tmp\"' EXIT", script)

    def test_zero_consensus_counts_are_diagnostic_and_cohort_local(self) -> None:
        factors = (ROOT / "scripts/consensus_track_factors.R").read_text(encoding="utf-8")
        normalized = (ROOT / "scripts/normalized_tracks_batch.sh").read_text(encoding="utf-8")
        differential = (ROOT / "scripts/differential_batch.sh").read_text(encoding="utf-8")
        self.assertIn('"consensus_count_sums.tsv"', factors)
        self.assertIn("zero consensus counts for samples:", factors)
        self.assertIn("skip_or_fail_family", normalized)
        self.assertIn("normalized_track_family_status.tsv", normalized)
        self.assertIn("consensus normalization unavailable", differential)

    def test_annotation_uses_configured_genome_order(self) -> None:
        script = (ROOT / "scripts/annotate_browser.sh").read_text(encoding="utf-8")
        self.assertGreaterEqual(script.count('bedtools sort -faidx "$chrom_sizes"'), 2)
        self.assertIn('bedtools closest -a "$sorted_consensus" -b "$genes" -d -g "$chrom_sizes"', script)

    def test_config_template_is_complete_and_safe(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            config = directory / "config.conf"
            text = (ROOT / "config/config.conf.template").read_text(encoding="utf-8")
            text = text.replace("/absolute/path/to/samplesheet.csv", str(directory / "samples.csv"))
            text = text.replace("/absolute/path/to/results", str(directory / "results"))
            config.write_text(text, encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_config.py"), str(config),
                "--template", str(ROOT / "config/config.conf.template"),
                "--write-shell", str(directory / "resolved.conf"),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            resolved = (directory / "resolved.conf").read_text(encoding="utf-8")
            for key in (
                "GENERATE_CPM_TRACKS", "GENERATE_DESEQ2_CONSENSUS_TRACKS",
                "GENERATE_CPM_PERMISSIVE_TRACKS", "GENERATE_CPM_INTERMEDIATE_TRACKS",
                "GENERATE_CPM_STRINGENT_TRACKS",
                "GENERATE_DESEQ2_ROBUST_CPM_PERMISSIVE_TRACKS",
                "GENERATE_DESEQ2_ROBUST_CPM_INTERMEDIATE_TRACKS",
                "GENERATE_DESEQ2_ROBUST_CPM_STRINGENT_TRACKS",
            ):
                self.assertIn(f"{key}=true", resolved)
            self.assertIn("THREADS_FASTQC=10", resolved)
            self.assertIn("THREADS_TRIMGALORE=8", resolved)
            self.assertIn("PEAKCALL_FAILURE_POLICY=continue", resolved)
            self.assertIn("ASSAY_PROFILE=chipseq", resolved)
            self.assertIn("ADAPTER_PRESET=auto", resolved)
            self.assertIn("ALLOW_CONTROL_FREE_PEAKCALL=true", resolved)
            self.assertIn("SPIKEIN_MODE=none", resolved)
            self.assertIn("ALLOW_SHARED_CONTROLS=true", resolved)
            self.assertIn("SE_SIGNAL_MODE=extend", resolved)
            self.assertIn("SE_FRAGMENT_LENGTH=auto", resolved)
            self.assertIn("PEAK_CALLERS=macs3,epic2", resolved)
            self.assertIn("RUN_IDR=false", resolved)

    def test_config_rejects_umi_processing_in_v01(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            text = (ROOT / "config/config.conf.template").read_text(encoding="utf-8")
            text = text.replace("/absolute/path/to/samplesheet.csv", str(directory / "samples.csv"))
            text = text.replace("/absolute/path/to/results", str(directory / "results"))
            text = text.replace("UMI_BARCODE_TAG=", "UMI_BARCODE_TAG=RX")
            config = directory / "config.conf"
            config.write_text(text, encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_config.py"), str(config),
                "--template", str(ROOT / "config/config.conf.template"),
                "--write-shell", str(directory / "resolved.conf"),
            ], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("UMI_BARCODE_TAG is not supported", result.stderr)

    def test_config_rejects_unknown_peakcall_failure_policy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            text = (ROOT / "config/config.conf.template").read_text(encoding="utf-8")
            text = text.replace("/absolute/path/to/samplesheet.csv", str(directory / "samples.csv"))
            text = text.replace("/absolute/path/to/results", str(directory / "results"))
            text = text.replace("PEAKCALL_FAILURE_POLICY=continue", "PEAKCALL_FAILURE_POLICY=ignore")
            config = directory / "config.conf"
            config.write_text(text, encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_config.py"), str(config),
                "--template", str(ROOT / "config/config.conf.template"),
                "--write-shell", str(directory / "resolved.conf"),
            ], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("PEAKCALL_FAILURE_POLICY must be one of", result.stderr)

    def test_config_rejects_nonpositive_preprocessing_threads(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            text = (ROOT / "config/config.conf.template").read_text(encoding="utf-8")
            text = text.replace("/absolute/path/to/samplesheet.csv", str(directory / "samples.csv"))
            text = text.replace("/absolute/path/to/results", str(directory / "results"))
            text = text.replace("THREADS_TRIMGALORE=8", "THREADS_TRIMGALORE=0")
            config = directory / "config.conf"
            config.write_text(text, encoding="utf-8")
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_config.py"), str(config),
                "--template", str(ROOT / "config/config.conf.template"),
                "--write-shell", str(directory / "resolved.conf"),
            ], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("THREADS_TRIMGALORE must be a positive integer", result.stderr)

    def test_config_rejects_unknown_key(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            text = (ROOT / "config/config.conf.template").read_text(encoding="utf-8")
            text = text.replace("/absolute/path/to/samplesheet.csv", str(directory / "samples.csv"))
            text = text.replace("/absolute/path/to/results", str(directory / "results")) + "\nUNKNOWN_SETTING=yes\n"
            config = directory / "config.conf"; config.write_text(text, encoding="utf-8")
            result = subprocess.run([sys.executable, str(ROOT / "scripts/validate_config.py"), str(config),
                "--template", str(ROOT / "config/config.conf.template"), "--write-shell", str(directory / "resolved")],
                text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("unknown configuration key", result.stderr)

    def test_samplesheet_resolves_control_and_cohort(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "metadata"
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_samplesheet.py"),
                str(ROOT / "config/examples/chipseq_pe.csv"), "--assay-profile", "chipseq",
                "--spikein-mode", "none", "--spikein-reference-id", "", "--peak-callers", "macs3,epic2",
                "--primary-peak-caller", "auto", "--output-dir", str(out),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (out / "sample_manifest.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            target = next(row for row in rows if row["is_control"] == "FALSE")
            self.assertEqual(target["control_key"], "WT_Input_R1.bioR1")
            self.assertEqual(target["primary_peak_caller"], "macs3")
            self.assertNotEqual(target["cohort_id"], ".")

    def test_control_free_single_end_broad_uses_epic2(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            out = Path(temporary) / "metadata"
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/validate_samplesheet.py"),
                str(ROOT / "config/examples/chipseq_se.csv"), "--assay-profile", "chipseq",
                "--spikein-mode", "none", "--spikein-reference-id", "",
                "--peak-callers", "macs3,epic2", "--allow-control-free",
                "--primary-peak-caller", "auto", "--output-dir", str(out),
            ], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (out / "sample_manifest.tsv").open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["control_key"], ".")
            self.assertEqual(row["primary_peak_caller"], "epic2")

    def test_different_antibodies_require_separate_runs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary); source = ROOT / "config/examples/chipseq_pe.csv"
            with source.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            target = dict(rows[0]); second = dict(target)
            second.update(sample_id="WT_H3K27ac", factor="H3K27ac", antibody_id="AB_H3K27AC_01",
                          output_prefix="WT_H3K27ac_R1", fastq_1="/data/WT_H3K27ac_R1.fastq.gz",
                          fastq_2="/data/WT_H3K27ac_R2.fastq.gz")
            sheet = directory / "samples.csv"
            with sheet.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), lineterminator="\n")
                writer.writeheader(); writer.writerows([target, second, rows[1]])
            out = directory / "metadata"
            rejected = subprocess.run([sys.executable, str(ROOT / "scripts/validate_samplesheet.py"), str(sheet),
                "--assay-profile", "chipseq", "--spikein-mode", "none", "--spikein-reference-id", "", "--peak-callers", "macs3,epic2",
                "--primary-peak-caller", "auto", "--output-dir", str(directory / "rejected")],
                text=True, capture_output=True)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("assigned to multiple targets", rejected.stderr)
            result = subprocess.run([sys.executable, str(ROOT / "scripts/validate_samplesheet.py"), str(sheet),
                "--assay-profile", "chipseq", "--spikein-mode", "none", "--spikein-reference-id", "", "--peak-callers", "macs3,epic2",
                "--primary-peak-caller", "auto", "--allow-shared-controls", "--output-dir", str(out)],
                text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("one compatible target/peak universe", result.stderr)


if __name__ == "__main__":
    unittest.main()
