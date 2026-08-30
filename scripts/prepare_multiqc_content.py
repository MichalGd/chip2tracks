#!/usr/bin/env python3
"""Prepare deterministic MultiQC custom content from chip2tracks outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import shutil
import sys
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def aggregate(
    pattern: str,
    root: Path,
    source_column: str | None = None,
    exclude_suffix: str = "",
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in sorted(root.glob(pattern)):
        if exclude_suffix and path.name.endswith(exclude_suffix):
            continue
        for row in read_tsv(path):
            if source_column:
                row[source_column] = path.relative_to(root).as_posix()
            rows.append(row)
    return rows


def safe_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "item"
    if len(cleaned) > 160:
        digest = hashlib.sha256(cleaned.encode()).hexdigest()[:12]
        cleaned = f"{cleaned[:147]}_{digest}"
    return cleaned


def write_custom_table(
    destination: Path,
    *,
    identifier: str,
    section_name: str,
    description: str,
    rows: list[dict[str, str]],
    columns: list[str] | None = None,
    key_columns: tuple[str, ...] = (),
) -> bool:
    if not rows:
        return False
    selected = columns or list(rows[0])
    if not selected:
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        handle.write(f"# id: {identifier}\n")
        handle.write(f"# section_name: {section_name}\n")
        handle.write(f"# description: {description}\n")
        handle.write("# format: tsv\n")
        handle.write("# plot_type: table\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["row_id", *selected])
        for index, row in enumerate(rows, 1):
            key = " | ".join(row.get(column, "") for column in key_columns).strip(" |")
            writer.writerow([key or f"row_{index}", *(row.get(column, "") for column in selected)])
    return True


def stage_images(root: Path, destination: Path) -> list[tuple[str, str]]:
    patterns = (
        "06_qc/controls/*.target_control_fingerprint.png",
        "06_qc/correlation_pca_fingerprint/*/spearman_heatmap.png",
        "06_qc/correlation_pca_fingerprint/*/pca.png",
        "06_qc/tss_signal_profile/*.descriptive_TSS_profile.png",
        "08_differential/**/pca.png",
        "08_differential/**/dispersion.png",
    )
    staged: list[tuple[str, str]] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for source in sorted(root.glob(pattern)):
            if not source.is_file() or source.stat().st_size == 0 or source in seen:
                continue
            seen.add(source)
            relative = source.relative_to(root).as_posix()
            name = f"chip2tracks_{safe_id(relative.removesuffix('.png'))}_mqc.png"
            shutil.copy2(source, destination / name)
            staged.append((relative, name))
    return staged


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("custom_dir", type=Path)
    args = parser.parse_args()
    root = args.output_dir.resolve()
    custom = args.custom_dir.resolve()
    try:
        if not (root / "00_metadata/sample_manifest.tsv").is_file():
            raise ValueError(f"sample manifest missing under output directory: {root}")
        custom.mkdir(parents=True, exist_ok=True)
        tables: list[tuple[str, str, str, list[dict[str, str]], list[str] | None, tuple[str, ...]]] = [
            ("run_summary", "Run summary", "Run-level chip2tracks result counts.",
             read_tsv(root / "10_reports/run_summary.tsv"), ["value"], ("metric",)),
            ("cohort_policy", "Cohort policy", "Researcher-selected cohort construction and hard compatibility boundaries.",
             read_tsv(root / "00_metadata/cohort_policy.tsv"), None, ("cohort_mode",)),
            ("cohort_membership", "Cohort membership", "Targets, controls, factor/antibody labels, and matched-control relationships used by the analysis.",
             read_tsv(root / "00_metadata/cohort_membership.tsv"), None, ("cohort_id", "sample_key", "role")),
            ("control_policy", "Control policy", "Whether control-free peak calling and shared controls were explicitly permitted.",
             read_tsv(root / "00_metadata/control_policy.tsv"), None, ("allow_control_free_peakcall",)),
            ("resource_budget", "Resource budget", "Maximum configured jobs × threads by long-running stage.",
             read_tsv(root / "00_metadata/resource_budget.tsv"), None, ("stage",)),
            ("software_versions", "Software versions", "Tool versions retained during preflight for provenance.",
             read_tsv(root / "00_metadata/software_versions.tsv"), None, ("tool",)),
            ("stage_timing", "Stage timing", "UTC start/end timestamps, elapsed time, and final status by workflow stage.",
             read_tsv(root / "00_metadata/stage_timing.tsv"), None, ("stage",)),
            ("warnings", "Workflow warnings and skips", "Recorded sample, cohort, and module warnings.",
             read_tsv(root / "10_reports/warning_summary.tsv"), ["severity", "item"], ("severity", "item")),
            ("observations", "Retained analysis observations", "Post-filter fragments for paired-end libraries or reads for single-end libraries.",
             read_tsv(root / "06_qc/alignment_and_complexity/observation_counts.tsv"),
             ["layout", "signal_unit", "analysis_observations"], ("sample_key",)),
            ("adapters", "Adapter trimming decisions", "Requested adapter preset and retained auto-detection evidence.",
             read_tsv(root / "01_fastq_qc/adapter_detection.tsv"), None, ("sample_key",)),
            ("complexity", "Library complexity", "Duplicate-retained q30 NRF and PCR bottleneck coefficients.",
             read_tsv(root / "06_qc/alignment_and_complexity/library_complexity.tsv"),
             ["layout", "total_observations", "distinct_observations", "NRF", "PBC1", "PBC2"], ("sample_key",)),
            ("alignment", "Host alignment records", "Host-alignment records before downstream analysis filtering.",
             aggregate("03_alignment/metrics/*.alignment.tsv", root),
             ["layout", "genome", "alignment_records", "spikein_mode"], ("sample_key",)),
            ("filtering", "Filtering and duplicate sensitivity", "Signal-unit counts across the four filtering policies.",
             aggregate("03_alignment/metrics/*.filter_counts.tsv", root),
             ["q0_dup_retained", "q0_dup_removed", "q30_dup_retained", "q30_dup_removed", "analysis_policy"], ("sample_key",)),
            ("peakcalls", "Per-sample peak calling", "Primary peak counts and caller warnings; empty and failed calls remain explicit.",
             read_tsv(root / "05_peaks/per_sample/peakcall_status.tsv"),
             ["primary_caller", "primary_class", "status", "primary_peak_count", "caller_warnings", "reason"], ("sample_key",)),
            ("consensus", "Consensus peak sets", "Successful contributions, exclusions, and retained consensus-region counts by cohort.",
             read_tsv(root / "05_peaks/consensus/consensus_status.tsv"),
             ["status", "total_samples", "successful_peak_samples", "excluded_samples", "regions", "reason"], ("cohort_id",)),
            ("reproducibility", "Replicate reproducibility", "Pairwise true-replicate IDR status for narrow MACS3 cohorts; broad cohorts are marked not applicable.",
             read_tsv(root / "05_peaks/reproducibility/status.tsv"),
             ["comparison", "method", "status", "passing_regions", "reason"], ("cohort_id", "comparison")),
            ("control_tracks", "Control-relative signal tracks", "MACS3 fold enrichment over matched-control/local lambda or local lambda alone.",
             read_tsv(root / "04_tracks/control_normalized/status.tsv"),
             ["background", "metric", "bigwig"], ("sample_key",)),
            ("cpm_tracks", "CPM coverage families", "Enabled analysis and fixed MAPQ/duplicate-policy CPM bedGraph/bigWig families.",
             read_tsv(root / "04_tracks/cpm/track_family_status.tsv"),
             ["enabled", "mapq", "duplicates", "path"], ("family",)),
            ("coverage_mapping", "Coverage-family mapping composition", "MAPQ 0, MAPQ <30, and Bowtie2 XS-tagged candidate-multimapper counts for the BAM policy underlying each CPM and normalized coverage family.",
             read_tsv(root / "04_tracks/cpm/mapping_composition.tsv"),
             ["policy", "cpm_family", "normalized_family_source", "bam_branch", "duplicates",
              "signal_unit", "total_observations", "mapq0_observations", "mapq0_percent",
              "mapq_lt30_observations", "mapq_lt30_percent",
              "xs_tagged_candidate_multimappers", "xs_tagged_percent"], ("sample_key", "policy")),
            ("mapping_definitions", "Mapping-composition definitions", "Interpretation limits for ambiguous and Bowtie2 XS-tagged signal-unit counts.",
             read_tsv(root / "04_tracks/cpm/mapping_composition_definitions.tsv"), None, ("metric",)),
            ("qc_modules", "QC module evidence", "Configured QC modules and non-empty retained artifact counts; retained evidence is not a biological pass/fail call.",
             read_tsv(root / "10_reports/qc_module_summary.tsv"), None, ("module",)),
            ("fragment_summary", "PE fragment-distribution summary", "Retained fragments in the plotted interval, mean/median/mode length, and percentage at or below 150 bp.",
             read_tsv(root / "10_reports/fragment_qc_summary.tsv"), None, ("sample_key",)),
            ("cross_correlation", "Strand cross-correlation summary", "phantompeakqualtools fragment-shift, phantom-peak, NSC, RSC, and quality-tag results; values are descriptive and target dependent.",
             read_tsv(root / "10_reports/cross_correlation_summary.tsv"), None, ("sample_key", "result_row")),
            ("track_inventory", "Retained track inventory", "Every retained bedGraph and bigWig, including family, format, size, and relative path.",
             read_tsv(root / "10_reports/track_inventory.tsv"), None, ("family", "path")),
            ("browser_track_families", "Browser track-family inventory", "UCSC descriptor and bigWig counts grouped by signal family.",
             read_tsv(root / "09_browser/ucsc/track_family_manifest.tsv"), None, ("family",)),
            ("normalization", "Normalized-track families", "Per-cohort normalization availability for each filtering policy.",
             read_tsv(root / "04_tracks/normalized_track_family_status.tsv"),
             ["status", "reason"], ("cohort_id", "policy")),
            ("differential", "Differential enrichment stage", "Run-level differential module failures and cohort skips.",
             read_tsv(root / "08_differential/stage_status.tsv"), None, ("status",)),
            ("frip", "Fraction of signal in consensus peaks", "Descriptive FRiP against each target cohort's consensus set.",
             aggregate("06_qc/frip_and_peak_reproducibility/*.frip.tsv", root,
                       exclude_suffix=".sample_primary_frip.tsv"),
             ["signal_unit", "total", "in_consensus", "frip"], ("sample_key",)),
            ("sample_frip", "Fraction of signal in sample primary peaks", "Descriptive FRiP against each sample's own primary caller output.",
             aggregate("06_qc/frip_and_peak_reproducibility/*.sample_primary_frip.tsv", root),
             ["signal_unit", "total", "in_sample_primary_peaks", "frip"], ("sample_key",)),
            ("optional_qc", "Optional QC status", "Explicit successes, skips and nonfatal failures for optional QC modules.",
             read_tsv(root / "06_qc/optional_qc_status.tsv"), ["metric", "status", "reason"], ("sample_key", "metric")),
            ("peak_annotation", "Peak feature annotation", "Feature-category counts and fractions for per-sample/caller and consensus peak sets.",
             read_tsv(root / "07_annotation/feature_summary/peak_feature_summary.tsv"),
             ["entity_type", "caller", "peak_class", "category", "count", "fraction", "percentage"],
             ("entity_id", "caller", "peak_class", "category")),
            ("spikein", "Spike-in calibration", "Host/spike observations, scale factors, and calibration status.",
             read_tsv(root / "06_qc/spikein/spikein_scaling.tsv"), None, ("sample_key",)),
            ("comparisons", "Differential occupancy summary", "All primary and sensitivity analysis variants, including disabled, skipped, failed, and completed comparisons.",
             read_tsv(root / "10_reports/differential_occupancy_summary.tsv"),
             ["factor", "peak_class", "analysis_family", "role", "method", "enabled",
              "comparison_id", "numerator", "reference", "tested_regions", "significant_regions",
              "higher_in_numerator", "higher_in_reference", "status", "reason", "normalization",
              "control_mode", "alpha", "min_abs_log2fc", "all_results", "significant_results"],
             ("cohort_id", "analysis_family", "method", "comparison_id")),
        ]
        manifest_rows: list[tuple[str, str, str]] = []
        for identifier, title, description, rows, columns, keys in tables:
            destination = custom / f"chip2tracks_{identifier}_mqc.tsv"
            if write_custom_table(destination, identifier=f"chip2tracks_{identifier}",
                                  section_name=title, description=description, rows=rows,
                                  columns=columns, key_columns=keys):
                manifest_rows.append(("table", identifier, destination.name))

        for source, destination in stage_images(root, custom):
            manifest_rows.append(("image", source, destination))

        with (custom / "custom_content_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["type", "source", "staged_name"])
            writer.writerows(manifest_rows)
        print(f"Prepared {len(manifest_rows)} MultiQC custom-content items in {custom}")
        return 0
    except (OSError, ValueError, csv.Error) as exc:
        print(f"MULTIQC CONTENT ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
