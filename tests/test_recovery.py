#!/usr/bin/env python3
"""Recovery contracts: invalidation, partial preprocessing and safe adoption."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import checkpoint
import recovery


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.output = self.directory / "output"
        self.metadata = self.output / "00_metadata"
        self.metadata.mkdir(parents=True)
        self.raw = self.directory / "raw.fastq.gz"
        with gzip.open(self.raw, "wb") as handle:
            handle.write(b"@read1\nACGT\n+\nIIII\n")
        self.sample = {"sample_key": "S.bioR1", "layout": "SE", "assay_profile": "chipseq",
                       "genome": "hg38", "analysis_duplicate_policy": "remove", "blacklist": "blacklist",
                       "fastq_1_list": str(self.raw), "fastq_2_list": ".", "technical_units": "1"}
        self.write_tsv("sample_manifest.tsv", [self.sample])
        self.reference = self.directory / "ref"
        self.reference.write_text("reference\n")
        self.write_tsv("reference_manifest.tsv", [{"kind": "bowtie2_index", "sha256": checkpoint.sha256(self.reference)}])
        self.config = self.directory / "resolved.conf"
        self.values = {"OUTPUT_DIR": str(self.output), "RUN_ID": "first", "RUN_FASTQC": "false",
                       "TRIM_ADAPTERS": "false", "RUN_MULTIQC": "false", "MIN_TRIMMED_LENGTH": "20",
                       "BOWTIE2_SEED": "0", "MIN_MAPQ": "30", "CHECKPOINT_PARALLEL_JOBS": "1"}
        self.manifest = self.metadata / "recovery_manifest.json"
        self.write_config()

    def write_tsv(self, name, rows):
        with (self.metadata / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def write_config(self):
        self.config.write_text("\n".join(f"{k}={shlex.quote(v)}" for k, v in self.values.items()) + "\n")

    def prepare(self, root=ROOT):
        self.write_config()
        with patch.object(recovery, "tool_identity", side_effect=lambda command: {"command": command}):
            recovery.prepare(root, self.config, self.manifest)
        return json.loads(self.manifest.read_text())

    def test_invocation_id_resources_and_cache_do_not_change_analysis_contracts(self):
        before = self.prepare()
        self.values.update(RUN_ID="second", THREADS_BOWTIE2="32", CHECKPOINT_PARALLEL_JOBS="2")
        after = self.prepare()
        for stage in recovery.SOURCES:
            if stage != "preflight":
                self.assertEqual(before["stages"][stage], after["stages"][stage], stage)
        self.assertEqual(before["preprocess_samples"], after["preprocess_samples"])

    def test_relocation_and_downstream_hotfix_preserve_upstream(self):
        before = self.prepare()
        copied = self.directory / "release"
        shutil.copytree(ROOT / "scripts", copied / "scripts", ignore=shutil.ignore_patterns("__pycache__"))
        shutil.copytree(ROOT / "common", copied / "common", ignore=shutil.ignore_patterns("__pycache__"))
        shutil.copyfile(ROOT / "VERSION", copied / "VERSION")
        cache = copied / "scripts/__pycache__"
        cache.mkdir()
        (cache / "ignored.pyc").write_bytes(b"cache")
        self.assertEqual(before, self.prepare(copied))
        with (copied / "scripts/differential_batch.sh").open("a") as handle:
            handle.write("\n# downstream hotfix\n")
        after = self.prepare(copied)
        for stage in ("preprocess", "alignment", "filtering", "cpm", "consensus", "qc"):
            self.assertEqual(before["stages"][stage], after["stages"][stage], stage)
        self.assertNotEqual(before["stages"]["differential"], after["stages"]["differential"])

    def test_raw_content_alignment_reference_and_mapq_changes(self):
        before = self.prepare()
        self.raw.write_bytes(self.raw.read_bytes() + b"changed content")
        after = self.prepare()
        self.assertNotEqual(before["preprocess_samples"], after["preprocess_samples"])
        self.values["BOWTIE2_SEED"] = "42"
        aligned = self.prepare()
        self.assertEqual(after["preprocess_samples"], aligned["preprocess_samples"])
        self.assertNotEqual(after["stages"]["alignment"], aligned["stages"]["alignment"])
        self.values["MIN_MAPQ"] = "40"
        filtered = self.prepare()
        self.assertEqual(aligned["stages"]["alignment"], filtered["stages"]["alignment"])
        self.assertNotEqual(aligned["stages"]["filtering"], filtered["stages"]["filtering"])
        self.write_tsv("reference_manifest.tsv", [{"kind": "bowtie2_index", "sha256": "changed-reference"}])
        reference = self.prepare()
        self.assertNotEqual(filtered["stages"]["alignment"], reference["stages"]["alignment"])
        self.assertEqual(filtered["stages"]["preprocess"], reference["stages"]["preprocess"])

    def test_multiqc_change_does_not_invalidate_finished_library(self):
        before = self.prepare()
        self.values["RUN_MULTIQC"] = "true"
        after = self.prepare()
        self.assertEqual(before["preprocess_samples"], after["preprocess_samples"])
        self.assertNotEqual(before["stages"]["preprocess"], after["stages"]["preprocess"])

    def test_aggregation_code_fix_preserves_library_completions(self):
        before = self.prepare()
        copied = self.directory / "release"
        shutil.copytree(ROOT / "scripts", copied / "scripts", ignore=shutil.ignore_patterns("__pycache__"))
        shutil.copytree(ROOT / "common", copied / "common", ignore=shutil.ignore_patterns("__pycache__"))
        shutil.copyfile(ROOT / "VERSION", copied / "VERSION")
        script = copied / "scripts/preprocess_batch.sh"
        script.write_text(script.read_text() + "\n# final aggregation correction\n")
        after = self.prepare(copied)
        self.assertEqual(before["preprocess_samples"], after["preprocess_samples"])
        self.assertNotEqual(before["stages"]["preprocess"], after["stages"]["preprocess"])

    def test_missing_index_or_changed_bam_fails_without_symlink_support(self):
        target = self.directory / "target.bam"
        index = self.directory / "target.bam.bai"
        target.write_bytes(b"bam-content")
        index.write_bytes(b"index")
        args = argparse.Namespace(checkpoint=self.directory / "checkpoint.json", stage="alignment",
                                  signature="inputs", outputs=[target, index], jobs=1)
        checkpoint.command_write(args)
        index.unlink()
        self.assertNotEqual(checkpoint.command_check(args), 0)
        index.write_bytes(b"index")
        target.write_bytes(b"changed-bam")
        self.assertNotEqual(checkpoint.command_check(args), 0)

    def test_rebuilt_upstream_invalidates_downstream_and_legacy_is_rejected(self):
        self.prepare()
        checkpoints = self.output / ".checkpoints"
        checkpoints.mkdir()
        parent = checkpoints / "preprocess.json"
        parent.write_text(json.dumps({"schema": 2, "signature": "same-inputs", "completed_utc": "first-generation"}))
        before = recovery.stage_signature(self.manifest, "alignment", checkpoints)
        parent.write_text(json.dumps({"schema": 2, "signature": "same-inputs", "completed_utc": "rebuilt-generation"}))
        self.assertNotEqual(before, recovery.stage_signature(self.manifest, "alignment", checkpoints))
        parent.write_text(json.dumps({"schema": 1, "signature": "legacy", "completed_utc": "old"}))
        with self.assertRaisesRegex(ValueError, "legacy checkpoint"):
            recovery.stage_signature(self.manifest, "alignment", checkpoints)

    def test_metadata_is_archived_before_replacement(self):
        old = self.metadata / "resolved_config.conf"
        old.write_text("RUN_ID=old\n")
        recovery.archive_metadata(self.output)
        old.write_text("RUN_ID=new\n")
        archives = list((self.metadata / "invocations").glob("*/resolved_config.conf"))
        self.assertEqual(len(archives), 1)
        self.assertEqual(archives[0].read_text(), "RUN_ID=old\n")

    def test_legacy_run_requires_explicit_rebuild_from_earliest_unproven_stage(self):
        checkpoints = self.output / ".checkpoints"
        checkpoints.mkdir()
        legacy = checkpoints / "preprocess.json"
        legacy.write_text(json.dumps({"schema": 1, "signature": "historical"}))
        original = legacy.read_bytes()
        with self.assertRaisesRegex(ValueError, "existing outputs and metadata were preserved"):
            recovery.guard_legacy(checkpoints, "")
        with self.assertRaisesRegex(ValueError, "legacy preprocess"):
            recovery.guard_legacy(checkpoints, "alignment")
        recovery.guard_legacy(checkpoints, "preprocess")
        self.assertEqual(original, legacy.read_bytes())

    def test_missing_index_changed_output_and_missing_analysis_link_are_rejected(self):
        target = self.directory / "target.bam"
        target.write_bytes(b"bam-content")
        index = self.directory / "target.bam.bai"
        index.write_bytes(b"index")
        link = self.directory / "analysis.bam"
        try:
            link.symlink_to(target)
        except OSError as exc:
            self.skipTest(f"host does not permit symlinks: {exc}")
        marker = self.directory / "checkpoint.json"
        args = argparse.Namespace(checkpoint=marker, stage="filtering", signature="inputs", outputs=[link, index], jobs=1)
        checkpoint.command_write(args)
        self.assertEqual(checkpoint.command_check(args), 0)
        link.unlink()
        self.assertNotEqual(checkpoint.command_check(args), 0)
        link.symlink_to(target)
        index.unlink()
        self.assertNotEqual(checkpoint.command_check(args), 0)
        index.write_bytes(b"index")
        target.write_bytes(b"changed-bam")
        self.assertNotEqual(checkpoint.command_check(args), 0)

    def test_fastqc_output_contract_requires_both_archive_and_html(self):
        self.values["RUN_FASTQC"] = "true"
        outputs = recovery.sample_outputs(self.values, self.sample)
        for path in outputs:
            path.parent.mkdir(parents=True, exist_ok=True)
            if path.suffix == ".zip":
                import zipfile
                with zipfile.ZipFile(path, "w") as archive:
                    archive.writestr("fastqc_data.txt", "test")
            else:
                path.write_bytes(b"content")
        recovery.validate_sample_outputs(outputs)
        html = next(p for p in outputs if p.suffix == ".html")
        html.unlink()
        with self.assertRaisesRegex(ValueError, "missing or empty"):
            recovery.validate_sample_outputs(outputs)

    def test_alignment_output_contract_rejects_index_that_was_never_created(self):
        bam = self.output / "03_alignment/sorted/S.bioR1.host.sorted.bam"
        bam.parent.mkdir(parents=True)
        bam.write_bytes(b"bam")
        with self.assertRaisesRegex(ValueError, "required BAM/index/QC artifact"):
            recovery.validate_stage(self.config, "alignment")
        Path(str(bam) + ".bai").write_bytes(b"index")
        (bam.parent / "stage_status.tsv").write_text("SUCCESS\n")
        metrics = self.output / "03_alignment/metrics/S.bioR1.alignment.tsv"
        metrics.parent.mkdir()
        metrics.write_text("metrics\n")
        recovery.validate_stage(self.config, "alignment")

    def test_later_outputs_in_shared_directories_do_not_invalidate_producers(self):
        self.test_alignment_output_contract_rejects_index_that_was_never_created()
        foreign = self.output / "03_alignment/spikein/filtered/derived.bam"
        foreign.parent.mkdir(parents=True)
        foreign.write_bytes(b"earlier spike-in result")
        own = recovery.validate_stage(self.config, "alignment")
        self.assertNotIn(foreign, own)
        marker = self.directory / "alignment.json"
        args = argparse.Namespace(checkpoint=marker, stage="alignment", signature="inputs", outputs=own, jobs=1)
        checkpoint.command_write(args)
        foreign.write_bytes(b"new spike-in result")
        self.assertEqual(checkpoint.command_check(args), 0)
        differential = self.output / "08_differential/cohort"
        differential.mkdir(parents=True)
        table = differential / "results.tsv"
        table.write_text("region_id\nregion1\n")
        annotated = differential / "results.annotated.tsv"
        annotated.write_text("old annotation\n")
        args = argparse.Namespace(checkpoint=self.directory / "differential.json", stage="differential",
                                  signature="inputs", outputs=[differential], jobs=1)
        checkpoint.command_write(args)
        annotated.write_text("new annotation\n")
        self.assertEqual(checkpoint.command_check(args), 0)


if __name__ == "__main__":
    unittest.main()
