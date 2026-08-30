#!/usr/bin/env python3
"""Join consensus gene/cCRE reference annotations to differential result tables."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


EXTRA_FIELDS = [
    "primary_feature_category", "primary_feature_id", "primary_gene_name",
    "primary_gene_id", "nearest_gene_name", "nearest_gene_id",
    "nearest_tss_signed_distance", "distance_to_gene",
    "ccre_reference_overlaps",
]


def empty_annotation() -> dict[str, str]:
    return {key: "." for key in EXTRA_FIELDS}


def annotation_map(annotation_root: Path, feature_summary: Path, cohort: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = defaultdict(lambda: {
        **empty_annotation()
    })
    nearest = next(annotation_root.glob("*.nearest_gene.tsv"), None)
    if nearest:
        with nearest.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 12:
                    result[fields[3]].update({"nearest_gene_name": fields[8], "nearest_gene_id": fields[9], "distance_to_gene": fields[11]})
    ccre = next(annotation_root.glob("*.ccre_reference_overlaps.tsv"), None)
    if ccre:
        labels: dict[str, set[str]] = defaultdict(set)
        with ccre.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 10 and fields[8] not in {".", "-1"}:
                    labels[fields[3]].add(fields[8])
        for region, values in labels.items():
            result[region]["ccre_reference_overlaps"] = ",".join(sorted(values))
    assignments = feature_summary / "peak_feature_assignments.tsv.gz"
    if assignments.is_file():
        with gzip.open(assignments, "rt", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("entity_type") != "consensus" or row.get("entity_id") != cohort:
                    continue
                result[row["peak_id"]].update({
                    "primary_feature_category": row.get("primary_category", "."),
                    "primary_feature_id": row.get("primary_feature_id", "."),
                    "primary_gene_name": row.get("primary_gene_name", "."),
                    "primary_gene_id": row.get("primary_gene_id", "."),
                    "nearest_gene_name": row.get("nearest_gene_name", "."),
                    "nearest_gene_id": row.get("nearest_gene_id", "."),
                    "nearest_tss_signed_distance": row.get("nearest_tss_signed_distance", "."),
                })
    return result


def annotate(path: Path, annotations: dict[str, dict[str, str]]) -> None:
    compressed = path.name.endswith(".tsv.gz")
    suffix = ".tsv.gz" if compressed else ".tsv"
    destination = path.with_name(path.name.removesuffix(suffix) + ".annotated" + suffix)
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        if not reader.fieldnames or "region_id" not in reader.fieldnames:
            return
        with opener(destination, "wt", encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(target, fieldnames=reader.fieldnames + EXTRA_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                row.update(annotations.get(row["region_id"], empty_annotation()))
                writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(); parser.add_argument("output_dir", type=Path); args = parser.parse_args()
    differential = args.output_dir / "08_differential"
    annotation = args.output_dir / "07_annotation"
    if not differential.is_dir():
        return 0
    for cohort_dir in differential.iterdir():
        if not cohort_dir.is_dir():
            continue
        mapping = annotation_map(
            annotation / cohort_dir.name / "consensus",
            annotation / "feature_summary", cohort_dir.name,
        )
        for path in list(cohort_dir.rglob("*.tsv.gz")) + list(cohort_dir.rglob("*.tsv")):
            if ".annotated.tsv" not in path.name:
                annotate(path, mapping)
    return 0


if __name__ == "__main__": raise SystemExit(main())
