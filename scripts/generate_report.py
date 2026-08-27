#!/usr/bin/env python3
"""Generate the dependency-light HTML and TSV workflow report."""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
from datetime import datetime, timezone
from pathlib import Path


def table(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def html_table(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "<p>Not available.</p>"
    columns = list(rows[0])
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    def cell(row: dict[str, str], column: str) -> str:
        value = str(row.get(column, ""))
        escaped = html.escape(value)
        if column in {"all_results", "significant_results"} and value not in {"", "."}:
            escaped_href = html.escape(f"../{value}", quote=True)
            return f"<td><a href='{escaped_href}'>{escaped}</a></td>"
        return f"<td>{escaped}</td>"
    body = "".join("<tr>" + "".join(cell(row, column) for column in columns) + "</tr>" for row in rows)
    return f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


DIFFERENTIAL_COLUMNS = [
    "cohort_id", "factor", "antibody_id", "peak_class", "analysis_family", "role",
    "method", "enabled", "comparison_id", "numerator", "reference", "tested_regions",
    "significant_regions", "higher_in_numerator", "higher_in_reference", "status", "reason",
    "normalization", "control_mode", "alpha", "min_abs_log2fc", "all_results",
    "significant_results",
]


def config_values(root: Path) -> dict[str, str]:
    return {
        row.get("key", ""): row.get("value", "")
        for row in table(root / "00_metadata/resolved_config.tsv")
        if row.get("key")
    }


def enabled(config: dict[str, str], key: str, fallback: bool) -> bool:
    if key not in config:
        return fallback
    return config[key].lower() == "true"


def marker_status(directory: Path) -> tuple[str, str] | None:
    for filename, status in (("FAILED.json", "FAILED"), ("SKIPPED.json", "SKIPPED")):
        marker = directory / filename
        if not marker.is_file():
            continue
        try:
            payload = json.loads(marker.read_text(encoding="utf-8"))
            return status, str(payload.get("reason") or payload.get("module") or filename)
        except (OSError, ValueError, TypeError):
            return status, filename
    return None


def result_path(root: Path, module: Path, value: str) -> str:
    if not value or value in {".", "NA"}:
        return "."
    path = Path(value)
    if not path.is_absolute():
        path = path if path.is_file() else module / path
    if not path.is_file():
        return "."
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def result_counts(path: Path, direction_columns: tuple[str, ...] = ()) -> tuple[int, int, int]:
    if not path.is_file():
        return 0, 0, 0
    opener = gzip.open if path.suffix == ".gz" else open
    positive = negative = total = 0
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total += 1
            value = next((row.get(column, "") for column in direction_columns if row.get(column, "") not in {"", "NA"}), "")
            try:
                fold = float(value)
            except ValueError:
                continue
            positive += fold > 0
            negative += fold < 0
    return total, positive, negative


def differential_summary(root: Path, cohorts: list[dict[str, str]]) -> list[dict[str, str]]:
    config = config_values(root)
    differential_root = root / "08_differential"
    alpha = config.get("DIFFERENTIAL_ALPHA", ".")
    min_lfc = config.get("DIFFERENTIAL_MIN_ABS_LOG2FC", ".")
    run_deseq2 = enabled(config, "RUN_DESEQ2_ENRICHMENT", True)
    run_diffbind = enabled(config, "RUN_DIFFBIND", True)
    run_subtracted = run_diffbind and enabled(config, "RUN_CONTROL_SUBTRACTED_SENSITIVITY", True)
    run_interaction = enabled(config, "RUN_TARGET_CONTROL_INTERACTION", False)
    rows: list[dict[str, str]] = []

    for cohort in cohorts:
        cohort_id = cohort.get("cohort_id", ".")
        peak_class = cohort.get("primary_peak_class", ".")
        cohort_directory = differential_root / cohort_id
        if peak_class in {"", "."} and cohort_directory.is_dir():
            candidates = sorted(path.name for path in cohort_directory.iterdir() if path.is_dir())
            if candidates:
                peak_class = candidates[0]
        analysis_root = cohort_directory / peak_class
        variants = (
            ("primary_target_only", "PRIMARY", "DESeq2Enrichment", run_deseq2,
             analysis_root / "primary_target_only/deseq2_enrichment",
             config.get("DIFFERENTIAL_NORMALIZATION", "deseq2"), "target_only_counts", min_lfc),
            ("primary_target_only", "PRIMARY", "DiffBind_DESeq2", run_diffbind,
             analysis_root / "primary_target_only/diffbind",
             "DiffBind_DESeq2", "no_control_subtraction", "."),
            ("sensitivity_control_subtracted", "SENSITIVITY", "DiffBind_DESeq2", run_subtracted,
             analysis_root / "sensitivity_control_subtracted/diffbind",
             "DiffBind_DESeq2", "scaled_control_subtraction", "."),
            ("sensitivity_target_control_interaction", "SENSITIVITY", "DESeq2Interaction", run_interaction,
             analysis_root / "sensitivity_target_control_interaction/deseq2",
             "joint_DESeq2_poscounts", "target_control_interaction", "."),
        )
        root_marker = marker_status(analysis_root)

        for family, role, method, is_enabled, module, normalization, control_mode, variant_lfc in variants:
            base = {
                "cohort_id": cohort_id,
                "factor": cohort.get("factor", "."),
                "antibody_id": cohort.get("antibody_id", "."),
                "peak_class": peak_class,
                "analysis_family": family,
                "role": role,
                "method": method,
                "enabled": str(is_enabled).lower(),
                "comparison_id": ".",
                "numerator": ".",
                "reference": ".",
                "tested_regions": ".",
                "significant_regions": ".",
                "higher_in_numerator": ".",
                "higher_in_reference": ".",
                "status": "NOT_AVAILABLE",
                "reason": "enabled analysis produced no recognized result or status artifact",
                "normalization": normalization,
                "control_mode": control_mode,
                "alpha": alpha,
                "min_abs_log2fc": variant_lfc,
                "all_results": ".",
                "significant_results": ".",
            }
            if not is_enabled:
                rows.append({**base, "status": "DISABLED", "reason": "disabled in resolved config"})
                continue
            if root_marker:
                rows.append({**base, "status": root_marker[0], "reason": root_marker[1]})
                continue
            module_marker = marker_status(module) or marker_status(module.parent)
            if module_marker:
                rows.append({**base, "status": module_marker[0], "reason": module_marker[1]})
                continue

            summaries = table(module / "comparison_summary.tsv")
            if summaries:
                for summary_row in summaries:
                    comparison_id = summary_row.get("comparison_id", ".") or "."
                    all_value = summary_row.get("all_results", "")
                    significant_value = summary_row.get("significant_results", "")
                    if not all_value and method == "DESeq2Enrichment":
                        all_value = f"comparisons/{comparison_id}/all.tsv.gz"
                        significant_value = f"comparisons/{comparison_id}/significant.tsv.gz"
                    rows.append({
                        **base,
                        "comparison_id": comparison_id,
                        "numerator": summary_row.get("numerator", ".") or ".",
                        "reference": summary_row.get("reference", ".") or ".",
                        "tested_regions": summary_row.get("tested", ".") or ".",
                        "significant_regions": summary_row.get("significant", ".") or ".",
                        "higher_in_numerator": summary_row.get("higher_in_numerator", ".") or ".",
                        "higher_in_reference": summary_row.get("higher_in_reference", ".") or ".",
                        "status": summary_row.get("status", "SUCCESS") or "SUCCESS",
                        "reason": ".",
                        "all_results": result_path(root, module, all_value),
                        "significant_results": result_path(root, module, significant_value),
                    })
                continue

            legacy: list[tuple[str, Path, Path, str, str]] = []
            if method == "DESeq2Enrichment":
                for all_path in sorted(module.glob("comparisons/*/all.tsv.gz")):
                    comparison_id = all_path.parent.name
                    parts = comparison_id.split("_vs_", 1)
                    numerator, reference = parts if len(parts) == 2 else (".", ".")
                    legacy.append((comparison_id, all_path, all_path.with_name("significant.tsv.gz"), numerator, reference))
            elif method == "DiffBind_DESeq2":
                for all_path in sorted(module.glob("contrast_*_all.tsv")):
                    comparison_id = all_path.name.removesuffix("_all.tsv")
                    legacy.append((comparison_id, all_path, all_path.with_name(f"{comparison_id}_significant.tsv"), ".", "."))
            else:
                all_path = module / "interaction_results_all.tsv.gz"
                conditions = [value for value in cohort.get("conditions", "").split(",") if value]
                numerator = conditions[1] if len(conditions) == 2 else "."
                reference = conditions[0] if len(conditions) == 2 else "."
                comparison_id = f"{numerator}_vs_{reference}_target_control_interaction" if len(conditions) == 2 else "target_control_interaction"
                if all_path.is_file():
                    legacy.append((comparison_id, all_path, module / "interaction_results_significant.tsv.gz", numerator, reference))
            for comparison_id, all_path, significant_path, numerator, reference in legacy:
                tested, _, _ = result_counts(all_path)
                significant, higher_numerator, higher_reference = result_counts(
                    significant_path, ("log2FoldChange", "Fold")
                )
                rows.append({
                    **base,
                    "comparison_id": comparison_id,
                    "numerator": numerator,
                    "reference": reference,
                    "tested_regions": str(tested),
                    "significant_regions": str(significant),
                    "higher_in_numerator": str(higher_numerator) if numerator != "." else ".",
                    "higher_in_reference": str(higher_reference) if reference != "." else ".",
                    "status": "SUCCESS",
                    "reason": "legacy result summarized during report generation",
                    "all_results": result_path(root, module, str(all_path)),
                    "significant_results": result_path(root, module, str(significant_path)),
                })
            if not legacy:
                rows.append(base)

    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    root = args.output_dir
    report_dir = root / "10_reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    samples = table(root / "00_metadata/sample_manifest.tsv")
    cohorts = table(root / "00_metadata/cohort_manifest.tsv")
    counts = table(root / "06_qc/alignment_and_complexity/observation_counts.tsv")
    consensus = table(root / "05_peaks/consensus/consensus_status.tsv")
    peakcalls = table(root / "05_peaks/per_sample/peakcall_status.tsv")
    reproducibility = table(root / "05_peaks/reproducibility/status.tsv")
    adapters = table(root / "01_fastq_qc/adapter_detection.tsv")
    complexity = table(root / "06_qc/alignment_and_complexity/library_complexity.tsv")
    normalized = table(root / "04_tracks/normalized_track_family_status.tsv")
    mapping_composition = table(root / "04_tracks/cpm/mapping_composition.tsv")
    differential = table(root / "08_differential/stage_status.tsv")
    differential_results = differential_summary(root, cohorts)
    metagene = table(root / "06_qc/metagene/artifacts.tsv")
    warnings = []
    for failure in root.rglob("FAILED.json"):
        warnings.append({"severity": "ERROR", "item": str(failure.relative_to(root))})
    for skipped in root.rglob("SKIPPED.json"):
        warnings.append({"severity": "INFO", "item": str(skipped.relative_to(root))})
    for row in peakcalls:
        if row.get("status") != "SUCCESS" or row.get("caller_warnings") not in {"", "."}:
            warnings.append({
                "severity": "WARNING",
                "item": (f"peakcalling:{row.get('sample_key', '?')}:"
                         f"{row.get('status', '?')}:{row.get('caller_warnings', '.')}")
            })
    for row in consensus:
        if row.get("status") != "SUCCESS":
            warnings.append({"severity": "WARNING", "item": (
                f"consensus:{row.get('cohort_id', '?')}:{row.get('status', '?')}:"
                f"{row.get('reason', '.')}")})
    for row in reproducibility:
        if row.get("status") == "ERROR":
            warnings.append({"severity": "WARNING", "item": (
                f"reproducibility:{row.get('cohort_id', '?')}:{row.get('comparison', '?')}:"
                f"{row.get('reason', '.')}")})
    for row in adapters:
        if row.get("resolved_method") == "unresolved":
            warnings.append({"severity": "WARNING", "item": f"adapter_detection:{row.get('sample_key', '?')}:unresolved"})
    for row in normalized:
        if row.get("status") != "SUCCESS":
            warnings.append({"severity": "WARNING", "item": (
                f"normalization:{row.get('cohort_id', '?')}:{row.get('policy', '?')}:"
                f"{row.get('status', '?')}:{row.get('reason', '.')}")})
    for row in differential:
        if row.get("status") != "SUCCESS":
            warnings.append({"severity": "WARNING", "item": (
                f"differential:{row.get('status', '?')}:failed_modules="
                f"{row.get('failed_modules', '?')}:skipped_cohorts={row.get('skipped_cohorts', '0')}")})
    for row in differential_results:
        if row.get("status") in {"FAILED", "NOT_AVAILABLE"}:
            warnings.append({"severity": "WARNING", "item": (
                f"differential_variant:{row.get('cohort_id', '?')}:"
                f"{row.get('analysis_family', '?')}:{row.get('method', '?')}:"
                f"{row.get('status', '?')}:{row.get('reason', '.')}")})
    preseq_dir = root / "06_qc/alignment_and_complexity"
    for log in sorted((root / "logs/qc").glob("*.preseq.log")):
        sample = log.name.removesuffix(".preseq.log")
        if not (preseq_dir / f"{sample}.preseq.txt").is_file():
            warnings.append({"severity": "WARNING", "item": f"preseq:{sample}:FAILED_OR_EMPTY"})
    for log in sorted((root / "logs/qc").glob("*.phantompeak.log")):
        sample = log.name.removesuffix(".phantompeak.log")
        result = root / "06_qc/fragment_length_and_periodicity" / f"{sample}.phantompeak.tsv"
        if not result.is_file() or result.stat().st_size == 0:
            warnings.append({"severity": "WARNING", "item": f"cross_correlation:{sample}:FAILED_OR_EMPTY"})

    warnings.sort(key=lambda row: (row["severity"], row["item"]))
    write_tsv(report_dir / "warning_summary.tsv", warnings, ["severity", "item"])
    write_tsv(
        report_dir / "differential_occupancy_summary.tsv",
        differential_results,
        DIFFERENTIAL_COLUMNS,
    )
    mapping_columns = [
        "sample_key", "policy", "cpm_family", "normalized_family_source", "bam_branch",
        "mapq_policy", "duplicates", "signal_unit", "total_observations",
        "mapq0_observations", "mapq0_percent", "mapq_lt30_observations",
        "mapq_lt30_percent", "xs_tagged_candidate_multimappers", "xs_tagged_percent",
    ]
    write_tsv(report_dir / "coverage_mapping_composition.tsv", mapping_composition, mapping_columns)
    xs_percentages = []
    for row in mapping_composition:
        try:
            xs_percentages.append(float(row.get("xs_tagged_percent", "0")))
        except ValueError:
            pass
    differential_variants = {
        (row.get("cohort_id"), row.get("analysis_family"), row.get("method"))
        for row in differential_results if row.get("enabled") == "true"
    }
    incomplete_differential_variants = {
        (row.get("cohort_id"), row.get("analysis_family"), row.get("method"))
        for row in differential_results
        if row.get("enabled") == "true" and row.get("status") != "SUCCESS"
    }
    significant_region_calls = 0
    for row in differential_results:
        if row.get("status") != "SUCCESS":
            continue
        try:
            significant_region_calls += int(row.get("significant_regions", "0"))
        except ValueError:
            pass
    summary = [
        {"metric": "biological_libraries", "value": str(len(samples))},
        {"metric": "target_cohorts", "value": str(len(cohorts))},
        {"metric": "consensus_success", "value": str(sum(row.get("status") == "SUCCESS" for row in consensus))},
        {"metric": "consensus_skipped_or_failed", "value": str(sum(row.get("status") != "SUCCESS" for row in consensus))},
        {"metric": "idr_comparisons_success", "value": str(sum(row.get("method") == "idr" and row.get("status") == "SUCCESS" for row in reproducibility))},
        {"metric": "peakcall_samples_with_warnings", "value": str(sum(
            row.get("status") != "SUCCESS" or row.get("caller_warnings") not in {"", "."}
            for row in peakcalls
        ))},
        {"metric": "warnings_and_skips", "value": str(len(warnings))},
        {"metric": "normalized_track_families_skipped_or_failed", "value": str(sum(
            row.get("status") != "SUCCESS" for row in normalized
        ))},
        {"metric": "differential_stage_status", "value": differential[0].get("status", "NOT_AVAILABLE") if differential else "NOT_AVAILABLE"},
        {"metric": "differential_variants_enabled", "value": str(len(differential_variants))},
        {"metric": "differential_comparisons_success", "value": str(sum(
            row.get("status") == "SUCCESS" for row in differential_results
        ))},
        {"metric": "differential_significant_region_calls_across_comparisons", "value": str(significant_region_calls)},
        {"metric": "differential_variants_skipped_failed_or_missing", "value": str(len(incomplete_differential_variants))},
        {"metric": "metagene_plot_tasks", "value": str(len(metagene))},
        {"metric": "maximum_xs_tagged_candidate_multimapper_percent", "value": (
            f"{max(xs_percentages):.4f}" if xs_percentages else "NOT_AVAILABLE")},
    ]
    write_tsv(report_dir / "run_summary.tsv", summary, ["metric", "value"])
    document = f"""<!doctype html><html><head><meta charset='utf-8'><title>chip2tracks report</title>
<style>body{{font-family:sans-serif;margin:2rem}}table{{border-collapse:collapse;font-size:.85rem}}th,td{{border:1px solid #ccc;padding:.3rem}}th{{background:#eee}}code{{background:#f4f4f4}}</style></head><body>
<h1>chip2tracks report</h1><p>Generated {datetime.now(timezone.utc).isoformat()}.</p>
<p>QC thresholds are descriptive. Different antibody targets are never normalized together.</p>
<h2>Run summary</h2>{html_table(summary)}<h2>Cohorts</h2>{html_table(cohorts)}
<h2>Retained analysis observations</h2>{html_table(counts)}<h2>Consensus status</h2>{html_table(consensus)}
<h2>Mapping composition of coverage families</h2>
<p><code>MAPQ 0</code> and <code>MAPQ &lt; 30</code> report ambiguous/low-confidence signal units and are not universal multimapper definitions. <code>XS-tagged candidate multimappers</code> are coverage-driving representative alignments carrying Bowtie2's alternative-alignment-score tag. For paired-end data, one representative alignment is counted per fragment.</p>
{html_table(mapping_composition)}
<h2>Adapter decisions</h2>{html_table(adapters)}<h2>Library complexity</h2>{html_table(complexity)}
<h2>Replicate reproducibility</h2>{html_table(reproducibility)}
<h2>Per-sample peak-calling status</h2>{html_table(peakcalls)}
<h2>Normalized-track families</h2>{html_table(normalized)}
<h2>Differential occupancy analysis</h2>
<p>The target-only analyses are primary. Control-subtracted and target-control interaction models are sensitivity analyses and are labelled separately. Counts are comparison-specific and must not be summed as unique genomic regions across methods or contrasts.</p>
{html_table(differential_results)}
<h3>Differential stage status</h3>{html_table(differential)}
<h2>Metagene aggregate-signal outputs</h2>{html_table(metagene)}
<h2>Warnings and skips</h2>{html_table(warnings)}</body></html>"""
    temporary_html = report_dir / "pipeline_report.html.tmp"
    temporary_html.write_text(document, encoding="utf-8")
    temporary_html.replace(report_dir / "pipeline_report.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
