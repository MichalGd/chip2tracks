#!/usr/bin/env python3
from __future__ import annotations

import csv
import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class UcscTrackTests(unittest.TestCase):
    def make_output(self, directory: Path) -> Path:
        output = directory / "output"
        relative_files = (
            "04_tracks/cpm/S1.CPM.bw",
            "04_tracks/cpm/permissive/S1.CPM.bw",
            "04_tracks/cpm/intermediate/S1.CPM.bw",
            "04_tracks/cpm/stringent/S1.CPM.bw",
            "04_tracks/deseq2_consensus/C1/S1.DESeq2Consensus.bw",
            "04_tracks/deseq2_robust_cpm/permissive/C1/S1.Robust.bw",
            "04_tracks/deseq2_robust_cpm/intermediate/C1/S1.Robust.bw",
            "04_tracks/deseq2_robust_cpm/stringent/C1/S1.Robust.bw",
            "04_tracks/spikein/C1/S1.SpikeIn.bw",
            "04_tracks/spikein_control/S1.dm6.CPM.bw",
            "04_tracks/control_normalized/S1.FE.bw",
            "04_tracks/future_family/S1.future.bw",
        )
        for relative in relative_files:
            path = output / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"nonempty bigwig fixture")
        return output

    def test_all_bigwig_families_use_valid_one_line_custom_tracks(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = self.make_output(Path(temporary))
            result = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts/generate_ucsc_tracks.py"),
                    str(output),
                    "--url-base",
                    "https://tracks.example.org/run7",
                    "--track-prefix",
                    "CHIP",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            ucsc = output / "09_browser/ucsc"
            combined = (ucsc / "all_bigwig_tracks.txt").read_text(encoding="utf-8")
            self.assertEqual(combined, (ucsc / "trackDb.txt").read_text(encoding="utf-8"))
            self.assertEqual(sum(line.startswith("track ") for line in combined.splitlines()), 12)
            self.assertIn("\n\n# family=cpm_permissive", combined)

            track_names: set[str] = set()
            for line in combined.splitlines():
                if not line or line.startswith("#"):
                    continue
                self.assertTrue(line.startswith("track "))
                tokens = shlex.split(line)
                attributes = dict(token.split("=", 1) for token in tokens[1:])
                self.assertEqual(attributes["type"], "bigWig")
                self.assertLessEqual(len(attributes["name"]), 15)
                self.assertLessEqual(len(attributes["description"]), 60)
                self.assertTrue(attributes["bigDataUrl"].startswith("https://tracks.example.org/run7/04_tracks/"))
                track_names.add(attributes["name"])
            self.assertEqual(len(track_names), 12)

            with (ucsc / "track_family_manifest.tsv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 12)
            self.assertTrue(all(row["url_mode"] == "remote" for row in rows))
            for row in rows:
                family_file = ucsc / row["descriptor_file"]
                self.assertTrue(family_file.is_file())
                self.assertEqual(sum(line.startswith("track ") for line in family_file.read_text(encoding="utf-8").splitlines()), 1)

    def test_empty_url_base_is_explicitly_local_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = self.make_output(Path(temporary))
            result = subprocess.run(
                [sys.executable, str(ROOT / "scripts/generate_ucsc_tracks.py"), str(output)],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            ucsc = output / "09_browser/ucsc"
            self.assertIn("local_path_not_ucsc_retrievable", (ucsc / "status.tsv").read_text(encoding="utf-8"))
            self.assertIn("public UCSC cannot retrieve", (ucsc / "README.txt").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
