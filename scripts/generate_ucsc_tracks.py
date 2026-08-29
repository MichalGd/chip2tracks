#!/usr/bin/env python3
"""Generate UCSC custom-track descriptors for every retained bigWig family."""

from __future__ import annotations

import argparse
import hashlib
import re
import shlex
from collections import OrderedDict
from pathlib import Path
from urllib.parse import quote, urlparse


FAMILIES = OrderedDict(
    (
        ("cpm_analysis", "CPM analysis"),
        ("cpm_permissive", "CPM permissive MAPQ0 duplicates retained"),
        ("cpm_intermediate", "CPM intermediate MAPQ0 duplicates removed"),
        ("cpm_stringent", "CPM stringent MAPQ30 duplicates removed"),
        ("deseq2_consensus", "DESeq2 consensus normalized"),
        ("deseq2_robust_cpm_permissive", "DESeq2 robust CPM permissive"),
        ("deseq2_robust_cpm_intermediate", "DESeq2 robust CPM intermediate"),
        ("deseq2_robust_cpm_stringent", "DESeq2 robust CPM stringent"),
        ("spikein", "Spike-in normalized host"),
        ("spikein_control", "Spike-in control CPM"),
        ("control_normalized", "MACS3 control-relative fold enrichment"),
        ("other", "Other bigWig tracks"),
    )
)

COLORS = {
    "cpm_analysis": "0,0,0",
    "cpm_permissive": "160,160,160",
    "cpm_intermediate": "70,130,180",
    "cpm_stringent": "0,70,180",
    "deseq2_consensus": "190,0,190",
    "deseq2_robust_cpm_permissive": "186,120,186",
    "deseq2_robust_cpm_intermediate": "148,70,148",
    "deseq2_robust_cpm_stringent": "105,0,105",
    "spikein": "200,0,0",
    "spikein_control": "220,110,0",
    "control_normalized": "0,130,70",
    "other": "90,90,90",
}


def family_for(relative: Path) -> str:
    parts = relative.parts
    if len(parts) < 3 or parts[:2] != ("04_tracks", "cpm"):
        top = parts[1] if len(parts) > 1 and parts[0] == "04_tracks" else ""
        if top in {"deseq2_consensus", "spikein", "spikein_control", "control_normalized"}:
            return top
        if top == "deseq2_robust_cpm" and len(parts) > 2:
            candidate = f"deseq2_robust_cpm_{parts[2]}"
            return candidate if candidate in FAMILIES else "other"
        return "other"
    if len(parts) == 3:
        return "cpm_analysis"
    candidate = f"cpm_{parts[2]}"
    return candidate if candidate in FAMILIES else "other"


def quoted(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def track_name(prefix: str, relative: Path) -> str:
    safe_prefix = re.sub(r"[^A-Za-z0-9]", "", prefix).upper()[:3] or "C2T"
    digest = hashlib.sha256(relative.as_posix().encode()).hexdigest()[: 15 - len(safe_prefix)]
    return safe_prefix + digest


def description(family: str, path: Path) -> str:
    text = f"{FAMILIES[family]} | {path.stem}"
    text = re.sub(r"[\r\n\t]+", " ", text).replace('"', "'")
    return text[:60]


def bigwig_url(output_dir: Path, relative: Path, url_base: str) -> tuple[str, str]:
    if url_base:
        encoded = quote(relative.as_posix(), safe="/-._~")
        return f"{url_base.rstrip('/')}/{encoded}", "remote"
    return str((output_dir / relative).resolve()), "local_path_not_ucsc_retrievable"


def validate_track_line(line: str, remote_required: bool) -> None:
    if "\n" in line or "\r" in line:
        raise ValueError("UCSC track descriptor contains an embedded line break")
    tokens = shlex.split(line)
    if not tokens or tokens[0] != "track":
        raise ValueError(f"invalid UCSC track line: {line}")
    attributes: dict[str, str] = {}
    for token in tokens[1:]:
        if "=" not in token:
            raise ValueError(f"UCSC attribute lacks '=': {token}")
        key, value = token.split("=", 1)
        attributes[key] = value
    required = {"type", "name", "description", "bigDataUrl"}
    missing = required - attributes.keys()
    if missing:
        raise ValueError(f"UCSC track line lacks: {', '.join(sorted(missing))}")
    if attributes["type"] != "bigWig":
        raise ValueError("UCSC descriptor type must be bigWig")
    if len(attributes["name"]) > 15 or not re.fullmatch(r"[A-Za-z0-9]+", attributes["name"]):
        raise ValueError("UCSC name must be <=15 alphanumeric characters")
    if len(attributes["description"]) > 60:
        raise ValueError("UCSC description must be <=60 characters")
    if remote_required and urlparse(attributes["bigDataUrl"]).scheme not in {"http", "https", "ftp"}:
        raise ValueError("remote bigDataUrl must use HTTP, HTTPS, or FTP")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--url-base", nargs="?", default="", const="")
    parser.add_argument("--track-prefix", default="CHIP")
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    tracks_root = output_dir / "04_tracks"
    ucsc_root = output_dir / "09_browser/ucsc"
    if not tracks_root.is_dir():
        parser.error(f"track directory is missing: {tracks_root}")
    if args.url_base and urlparse(args.url_base).scheme not in {"http", "https", "ftp"}:
        parser.error("--url-base must use HTTP, HTTPS, or FTP")
    ucsc_root.mkdir(parents=True, exist_ok=True)
    for index, family in enumerate(FAMILIES, 1):
        (ucsc_root / f"{index:02d}_{family}.txt").unlink(missing_ok=True)

    grouped: dict[str, list[tuple[Path, str]]] = {family: [] for family in FAMILIES}
    names: set[str] = set()
    url_mode = "remote" if args.url_base else "local_path_not_ucsc_retrievable"
    for bigwig in sorted(tracks_root.rglob("*.bw")):
        if not bigwig.is_file() or bigwig.stat().st_size == 0:
            continue
        relative = bigwig.relative_to(output_dir)
        family = family_for(relative)
        name = track_name(args.track_prefix, relative)
        if name in names:
            raise ValueError(f"duplicate generated UCSC track name: {name}")
        names.add(name)
        url, current_mode = bigwig_url(output_dir, relative, args.url_base)
        if current_mode != url_mode:
            raise ValueError("inconsistent UCSC URL modes")
        line = " ".join(
            (
                "track",
                "type=bigWig",
                f"name={quoted(name)}",
                f"description={quoted(description(family, bigwig))}",
                "visibility=full",
                "autoScale=on",
                f"color={COLORS[family]}",
                f"bigDataUrl={quoted(url)}",
            )
        )
        validate_track_line(line, remote_required=bool(args.url_base))
        grouped[family].append((relative, line))

    combined_sections: list[str] = []
    manifest_rows: list[tuple[str, str, int, str]] = []
    for index, (family, label) in enumerate(FAMILIES.items(), 1):
        entries = grouped[family]
        if not entries:
            continue
        filename = f"{index:02d}_{family}.txt"
        content = f"# family={family}; {label}\n" + "\n".join(line for _, line in entries) + "\n"
        (ucsc_root / filename).write_text(content, encoding="utf-8", newline="\n")
        combined_sections.append(content.rstrip("\n"))
        manifest_rows.append((family, filename, len(entries), url_mode))

    combined = "\n\n".join(combined_sections)
    if combined:
        combined += "\n"
    else:
        combined = "# No non-empty bigWig files were found under 04_tracks.\n"
    (ucsc_root / "all_bigwig_tracks.txt").write_text(combined, encoding="utf-8", newline="\n")
    # Preserve the historical path, but make its contents valid custom-track
    # text rather than multi-line trackDb hub stanzas.
    (ucsc_root / "trackDb.txt").write_text(combined, encoding="utf-8", newline="\n")

    with (ucsc_root / "track_family_manifest.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("family\tdescriptor_file\ttrack_count\turl_mode\n")
        for row in manifest_rows:
            handle.write("\t".join(map(str, row)) + "\n")

    readme = [
        "UCSC custom-track descriptors generated by chip2tracks.",
        "Each non-comment track descriptor occupies exactly one physical line.",
        "all_bigwig_tracks.txt and trackDb.txt contain all families separated by",
        "a comment plus a blank line. UCSC ignores comments and empty lines.",
        f"URL mode: {url_mode}",
    ]
    if not args.url_base:
        readme.extend(
            (
                "WARNING: UCSC_BIGDATA_URL_BASE was empty. Descriptor syntax is valid,",
                "but public UCSC cannot retrieve server-local filesystem paths. Set an",
                "HTTP/HTTPS/FTP base with byte-range support and regenerate these files.",
            )
        )
    (ucsc_root / "README.txt").write_text("\n".join(readme) + "\n", encoding="utf-8", newline="\n")
    with (ucsc_root / "status.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("status\ttotal_tracks\tfamilies\turl_mode\n")
        handle.write(f"SUCCESS\t{len(names)}\t{len(manifest_rows)}\t{url_mode}\n")
    print(f"Generated {len(names)} UCSC bigWig descriptors across {len(manifest_rows)} families")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
