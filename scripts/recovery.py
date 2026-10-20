#!/usr/bin/env python3
"""Stage input contracts and validated preprocessing/marking recovery.

The run signature remains provenance; these contracts decide scientific reuse.
Unknown configuration settings are conservatively included in every contract.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import zipfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from checkpoint import sha256

SOURCES = {
    "preflight": "preflight.sh",
    "preprocess": "preprocess_batch.sh",
    "alignment": "align_batch.sh",
    "filtering": "mark_filter_batch.sh",
    "cpm": "coverage_batch.sh",
    "peakcalling": "peakcall_batch.sh",
    "consensus": "reproducibility_batch.sh consensus_batch.sh build_consensus.py",
    "spikein": "spikein_batch.sh spikein_qc.py",
    "normalized_tracks": "normalized_tracks_batch.sh consensus_track_factors.R",
    "metagene": "metagene_batch.sh build_metagene_track_manifest.py",
    "qc": "qc_batch.sh prepare_tss_bed.py",
    "differential": "differential_batch.sh deseq2_enrichment_analysis.R diffbind_analysis.R control_interaction_analysis.R",
    "annotation": "annotate_browser.sh summarize_peak_annotations.py generate_ucsc_tracks.py annotate_differential_results.py",
    "report": "report_batch.sh generate_report.py generate_multiqc_report.sh prepare_multiqc_content.py",
    "cleanup": "cleanup.sh",
    "finalize": "finalize.py",
}
PARENTS = {
    "preflight": (), "preprocess": (), "alignment": ("preprocess",),
    "filtering": ("alignment",), "cpm": ("filtering",),
    "peakcalling": ("filtering", "cpm"), "consensus": ("peakcalling",),
    "spikein": ("alignment", "filtering", "cpm"),
    "normalized_tracks": ("filtering", "consensus", "spikein", "cpm"),
    "metagene": ("cpm", "spikein"),
    "qc": ("filtering", "peakcalling", "consensus", "cpm", "metagene", "spikein"),
    "differential": ("normalized_tracks", "consensus", "filtering", "spikein"),
    "annotation": ("peakcalling", "consensus", "differential", "cpm", "normalized_tracks", "spikein"),
    "report": ("preflight", "qc", "annotation", "normalized_tracks", "differential"),
    "cleanup": ("report",), "finalize": ("cleanup",),
}
TOOLS = {
    "preflight": (), "preprocess": ("fastqc", "trim_galore", "cutadapt", "multiqc"),
    "alignment": ("bowtie2", "samtools"),
    "filtering": ("samtools", "bedtools", "PICARD_COMMAND"),
    "cpm": ("bamCoverage", "bedtools", "bedGraphToBigWig", "MACS3_COMMAND"),
    "peakcalling": ("MACS3_COMMAND", "EPIC2_COMMAND", "bedtools", "bedGraphToBigWig"),
    "consensus": ("IDR_COMMAND",),
    "spikein": ("PICARD_COMMAND", "samtools", "bedtools", "bamCoverage", "bedGraphToBigWig"),
    "normalized_tracks": ("Rscript", "bamCoverage", "bedtools", "bedGraphToBigWig"),
    "metagene": ("computeMatrix", "plotProfile", "plotHeatmap"),
    "qc": ("samtools", "bedtools", "preseq", "PHANTOMPEAK_COMMAND", "plotFingerprint", "multiBamSummary", "plotCorrelation", "plotPCA", "computeMatrix", "plotProfile"),
    "differential": ("Rscript",), "annotation": ("bedtools",),
    "report": ("multiqc",), "cleanup": (), "finalize": (),
}
REFS = {
    "preprocess": set(), "alignment": {"bowtie2_index", "composite_bowtie2_index", "chrom_sizes", "spikein_allowed_contigs"},
    "filtering": {"canonical_contigs", "blacklist"}, "cpm": {"chrom_sizes"},
    "peakcalling": {"chrom_sizes"}, "consensus": set(),
    "spikein": {"spikein_chrom_sizes", "spikein_allowed_contigs", "spikein_blacklist"},
    "normalized_tracks": {"chrom_sizes"},
    "qc": {"gtf", "tss_bed", "chrom_sizes"},
    "annotation": {"gtf", "ccre_bed", "chrom_sizes"},
}
SAMPLE_FIELDS = {
    "preprocess": ("sample_key", "layout", "assay_profile", "fastq_1_list", "fastq_2_list", "technical_units"),
    "alignment": ("sample_key", "layout", "genome"),
    "filtering": ("sample_key", "layout", "genome", "analysis_duplicate_policy", "blacklist"),
    "cpm": ("sample_key", "layout", "genome", "analysis_duplicate_policy"),
}


def digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read_config(path: Path) -> dict[str, str]:
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("#"):
            key, value = line.split("=", 1)
            tokens = shlex.split(value)
            values[key] = tokens[0] if tokens else ""
    return values


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def operational(key: str) -> bool:
    return (key in {"RUN_ID", "SAMPLESHEET", "TOTAL_CPU_BUDGET", "RESOURCE_CHECK_MODE",
                    "WRITE_COMMAND_LOG", "WRITE_STRUCTURED_LOG", "WRITE_CONSOLE_LOG",
                    "WRITE_FILE_CHECKSUMS", "CHECKSUM_PARALLEL_JOBS", "CHECKPOINT_PARALLEL_JOBS"}
            or key.startswith("THREADS_") or key.endswith("_PARALLEL_JOBS"))


def tool_identity(command: str) -> dict[str, str]:
    path = shutil.which(command)
    if not path:
        return {"command": command, "executable": "unavailable"}
    resolved = Path(path).resolve()
    result = {"command": command, "executable": str(resolved), "sha256": sha256(resolved)}
    # Capture package versions as well as launchers (Conda scripts can retain
    # identical bytes across a package upgrade). Never execute analysis commands.
    if Path(command).name in {"fastqc", "trim_galore", "cutadapt", "multiqc", "bowtie2", "samtools", "bedtools", "bamCoverage", "macs3", "epic2", "Rscript"}:
        try:
            probe = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=30)
            result["version"] = (probe.stdout + probe.stderr).strip()
            result["version_exit"] = str(probe.returncode)
        except (OSError, subprocess.TimeoutExpired) as exc:
            result["version"] = type(exc).__name__
    return result


def source_files(root: Path, stage: str) -> list[Path]:
    files = [root / "scripts" / name for name in SOURCES[stage].split()]
    files += [root / "scripts/lib/common.sh", root / "scripts/lib/parallel_jobs.sh"]
    if stage == "metagene":
        files += [p for p in (root / "common/metagene").rglob("*")
                  if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc"]
    return files


def prepare(root: Path, config_path: Path, output: Path) -> None:
    config = read_config(config_path)
    metadata = Path(config["OUTPUT_DIR"]) / "00_metadata"
    samples = read_tsv(metadata / "sample_manifest.tsv")
    references = read_tsv(metadata / "reference_manifest.tsv")
    texts = {stage: "\n".join((root / "scripts" / name).read_text(encoding="utf-8") for name in SOURCES[stage].split())
             for stage in SOURCES}
    used = {stage: set(re.findall(r"\b[A-Z][A-Z0-9_]+\b", text)) for stage, text in texts.items()}
    # Common coverage helper consumes these options even when they do not
    # appear in the calling shell script. Dynamic reference keys are handled below.
    for stage in ("cpm", "spikein", "normalized_tracks"):
        used[stage].update({"SE_SIGNAL_MODE", "SE_FRAGMENT_LENGTH", "SE_FRAGMENT_LENGTH_FALLBACK"})
    known = set().union(*used.values())
    reference_key = re.compile(r"^(INDEX|FASTA|CHROM_SIZES|CANONICAL_CONTIGS|GTF|BLACKLIST|TSS_BED|CCRE_BED|EFFECTIVE_GENOME_SIZE)_")
    known.update({"COHORT_MODE", "ALLOW_SHARED_CONTROLS", "ALLOW_CONTROL_FREE_PEAKCALL", "ALLOW_MIXED_LAYOUTS", "ALLOW_MIXED_GENOMES"})
    unknown = {k for k in config if k not in known and not operational(k) and not reference_key.match(k)}
    raw_paths = sorted({Path(p) for sample in samples for field in ("fastq_1_list", "fastq_2_list")
                        for p in sample.get(field, "").split(";") if p and p != "."}, key=str)
    with ThreadPoolExecutor(max_workers=int(config.get("CHECKPOINT_PARALLEL_JOBS", "1"))) as executor:
        raw_hashes = dict(zip(raw_paths, executor.map(sha256, raw_paths), strict=True))
    raw = {}
    for sample in samples:
        paths = [Path(p) for field in ("fastq_1_list", "fastq_2_list")
                 for p in sample.get(field, "").split(";") if p and p != "."]
        raw[sample["sample_key"]] = [{"path": str(p.resolve()), "sha256": raw_hashes[p]} for p in paths]
    tools: dict[str, dict[str, str]] = {}
    contracts = {}
    for stage in SOURCES:
        keys = set(config) if stage in {"preflight", "report", "finalize"} else used[stage] | unknown
        # Reference configuration keys must be included for stages using the
        # resource even if they are accessed via reference_value dynamically.
        for key in config:
            if reference_key.match(key):
                if stage not in {"preprocess", "alignment", "filtering", "cpm"}:
                    keys.add(key)
                if stage == "alignment" and key.startswith(("INDEX_", "CHROM_SIZES_")):
                    keys.add(key)
                if stage == "filtering" and key.startswith(("CANONICAL_CONTIGS_", "BLACKLIST_")):
                    keys.add(key)
                if stage == "cpm" and key.startswith(("CHROM_SIZES_", "EFFECTIVE_GENOME_SIZE_")):
                    keys.add(key)
        selected = {k: v for k, v in config.items() if k in keys and not operational(k)}
        if stage == "preflight":
            selected.update({k: v for k, v in config.items() if k != "RUN_ID" and operational(k)})
        stage_tools = {}
        for name in TOOLS[stage]:
            command = config.get(name, name)
            if command not in tools:
                tools[command] = tool_identity(command)
            stage_tools[name] = tools[command]
        fields = SAMPLE_FIELDS.get(stage)
        rows = [{k: row.get(k, "") for k in fields} for row in samples] if fields else samples
        contract = {
            "contract_schema": 1, "config": selected, "samples": rows,
            "sources": {p.relative_to(root).as_posix(): sha256(p) for p in source_files(root, stage)},
            "references": [r for r in references if stage not in REFS or r["kind"] in REFS[stage]],
            "tools": stage_tools,
        }
        if stage == "preprocess":
            contract["raw_inputs"] = raw
        if stage in {"report", "finalize"}:
            contract["workflow_version"] = (root / "VERSION").read_text().strip()
        contracts[stage] = contract
    # Exclude schema implementation bytes from stage contracts: changes to the
    # fingerprint algorithm require an explicit contract_schema bump instead.
    sample_contract = dict(contracts["preprocess"])
    sample_contract["config"] = {k: v for k, v in sample_contract["config"].items()
                                 if k not in {"RUN_MULTIQC", "MULTIQC_EXPORT_PLOTS"}}
    sample_contract["tools"] = {k: v for k, v in sample_contract["tools"].items() if k != "multiqc"}
    worker_text = (root / "scripts/preprocess_batch.sh").read_text(encoding="utf-8")
    worker_text = worker_text[worker_text.index("merge_fastqs() {"):worker_text.index("parallel_require_positive_integer")]
    sample_contract["sources"] = {k: v for k, v in sample_contract["sources"].items()
                                  if k != "scripts/preprocess_batch.sh"}
    sample_contract["sources"]["preprocess_worker"] = digest(worker_text)
    per_sample = {}
    for row in sample_contract["samples"]:
        key = row["sample_key"]
        per_sample[key] = digest({**sample_contract, "samples": [row], "raw_inputs": raw[key]})
    atomic_json(output, {"stages": contracts, "preprocess_samples": per_sample})


def stage_signature(manifest: Path, stage: str, checkpoints: Path) -> str:
    contract = json.loads(manifest.read_text())["stages"][stage]
    upstream = {}
    for parent in PARENTS[stage]:
        payload = json.loads((checkpoints / f"{parent}.json").read_text())
        if payload.get("schema") != 2:
            raise ValueError(f"{parent}: legacy checkpoint has no compatible input provenance; restart at or before {parent}")
        upstream[parent] = {"signature": payload["signature"], "generation": payload["completed_utc"]}
    return digest({"contract": contract, "upstream": upstream})


def fastqc_paths(fastq: Path, directory: Path) -> list[Path]:
    name = re.sub(r"\.(fastq|fq)(\.gz)?$", "", fastq.name)
    return [directory / f"{name}_fastqc.{suffix}" for suffix in ("zip", "html")]


def sample_outputs(config: dict[str, str], row: dict[str, str]) -> list[Path]:
    root = Path(config["OUTPUT_DIR"])
    key = row["sample_key"]
    work = root / "02_trimmed_fastq"
    outputs = []
    for mate in (1, 2) if row["layout"] == "PE" else (1,):
        raw = work / f"{key}.R{mate}.merged.fastq.gz"
        trimmed = work / f"{key}.R{mate}.trimmed.fastq.gz"
        outputs += [raw, trimmed]
        if config.get("RUN_FASTQC") == "true":
            outputs += fastqc_paths(raw, root / "01_fastq_qc/raw")
            outputs += fastqc_paths(trimmed, root / "01_fastq_qc/trimmed")
            if config.get("RUN_FASTQC_PER_TECHNICAL_UNIT") == "true":
                for unit, path in enumerate(row[f"fastq_{mate}_list"].split(";"), 1):
                    outputs += fastqc_paths(Path(path), root / f"01_fastq_qc/raw_units/{key}/techR{unit}")
    if config.get("TRIM_ADAPTERS") == "true":
        outputs.append(root / f"01_fastq_qc/adapters/{key}.tsv")
    return outputs


def validate_sample_outputs(paths: list[Path]) -> None:
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"required preprocessing output missing or empty: {path}")
        if path.suffix == ".zip":
            with zipfile.ZipFile(path) as archive:
                if archive.testzip() is not None or not archive.namelist():
                    raise ValueError(f"invalid FastQC archive: {path}")


def validate_stage(config_path: Path, stage: str) -> list[Path]:
    config = read_config(config_path)
    root = Path(config["OUTPUT_DIR"])
    samples = read_tsv(root / "00_metadata/sample_manifest.tsv")
    required = []

    def bam(path: Path, picard: bool = False) -> None:
        required.append(path)
        index = Path(str(path) + ".bai")
        if picard and not index.is_file():
            index = path.with_suffix(".bai")
        required.append(index)

    for row in samples:
        key = row["sample_key"]
        if stage == "preprocess":
            validate_sample_outputs(sample_outputs(config, row))
        elif stage == "alignment":
            bam(root / f"03_alignment/sorted/{key}.host.sorted.bam")
            required.append(root / f"03_alignment/metrics/{key}.alignment.tsv")
            if config.get("SPIKEIN_MODE", "none") != "none":
                bam(root / f"03_alignment/spikein/composite/{key}.composite.sorted.bam")
                bam(root / f"03_alignment/spikein/spike/{key}.{config['SPIKEIN_MODE']}.bam")
        elif stage == "filtering":
            bam(root / f"03_alignment/marked/{key}.host.marked.bam", picard=True)
            required += [root / f"03_alignment/metrics/{key}.duplicate_metrics.txt",
                         root / f"03_alignment/metrics/{key}.filter_counts.tsv"]
            for branch in ("q0_dup-retained", "q0_dup-removed", "q30_dup-retained", "q30_dup-removed"):
                quality, policy = branch.split("_", 1)
                bam(root / f"03_alignment/filtered/{branch}/{key}.host.{quality}.{policy}.bam")
            analysis = root / f"03_alignment/analysis/{key}.host.analysis.bam"
            bam(analysis)
            if analysis.is_symlink():
                policy = "retained" if row["analysis_duplicate_policy"] == "retain" else "removed"
                expected = root / f"03_alignment/filtered/q30_dup-{policy}/{key}.host.q30.dup-{policy}.bam"
                if analysis.resolve() != expected.resolve():
                    raise ValueError(f"analysis BAM link does not match configured duplicate policy: {analysis}")
    if stage == "alignment":
        required.append(root / "03_alignment/sorted/stage_status.tsv")
    elif stage == "filtering":
        required.append(root / "03_alignment/analysis/stage_status.tsv")
    elif stage == "annotation":
        required += sorted((root / "08_differential").rglob("*.annotated.tsv*"))
    if stage == "preprocess":
        required.append(root / "01_fastq_qc/adapter_detection.tsv")
        if config.get("RUN_MULTIQC") == "true":
            required.append(root / "01_fastq_qc/multiqc/multiqc_report.html")
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"{stage}: required BAM/index/QC artifact missing or empty: {path}")
    return required


def archive_metadata(output: Path) -> None:
    metadata = output / "00_metadata"
    names = ["resolved_config.conf", "resolved_config.tsv", "sanitized_samplesheet.csv",
             "sample_manifest.tsv", "cohort_manifest.tsv", "reference_manifest.tsv",
             "run_manifest.tsv", "recovery_manifest.json", "software_versions.tsv"]
    existing = [metadata / name for name in names if (metadata / name).is_file()]
    if existing:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        archive = metadata / "invocations" / f"{stamp}-{os.getpid()}"
        archive.mkdir(parents=True)
        for path in existing:
            shutil.copy2(path, archive / path.name)


def guard_legacy(checkpoints: Path, from_stage: str) -> None:
    """Do not silently overwrite unproven legacy runs during an ordinary retry."""
    names = list(SOURCES)
    from_stage = "consensus" if from_stage == "reproducibility" else from_stage
    requested = names.index(from_stage) if from_stage else len(names)
    for stage in names:
        if stage in {"preflight", "cleanup", "finalize"}:
            continue  # Operational audits/finalization can safely be repeated.
        path = checkpoints / f"{stage}.json"
        if path.is_file():
            try:
                payload = json.loads(path.read_text())
            except ValueError:
                continue  # Normal checkpoint validation diagnoses corruption.
            if payload.get("schema") != 2 and requested > names.index(stage):
                raise ValueError(f"legacy {stage} checkpoint lacks input provenance; existing outputs and metadata were preserved. "
                                 f"Revalidate the original run, or explicitly rebuild with --from-stage {stage} (or earlier). "
                                 "Do not rewrite checkpoint signatures.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    archive = sub.add_parser("archive")
    archive.add_argument("output", type=Path)
    guard = sub.add_parser("guard-legacy")
    guard.add_argument("--checkpoints", type=Path, required=True)
    guard.add_argument("--from-stage", default="")
    plan = sub.add_parser("prepare")
    plan.add_argument("--root", type=Path, required=True)
    plan.add_argument("--config", type=Path, required=True)
    plan.add_argument("--output", type=Path, required=True)
    stage = sub.add_parser("stage")
    stage.add_argument("--manifest", type=Path, required=True)
    stage.add_argument("--stage", choices=SOURCES, required=True)
    stage.add_argument("--checkpoints", type=Path, required=True)
    validate = sub.add_parser("validate-stage")
    validate.add_argument("--config", type=Path, required=True)
    validate.add_argument("--stage", choices=SOURCES, required=True)
    validate.add_argument("--print-paths", action="store_true")
    sample = sub.add_parser("sample")
    sample.add_argument("--config", type=Path, required=True)
    sample.add_argument("--sample-key", required=True)
    sample.add_argument("--validate", action="store_true")
    marked = sub.add_parser("marked")
    marked.add_argument("--config", type=Path, required=True)
    marked.add_argument("--bam", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "archive":
            archive_metadata(args.output)
        elif args.command == "guard-legacy":
            guard_legacy(args.checkpoints, args.from_stage)
        elif args.command == "prepare":
            prepare(args.root.resolve(), args.config, args.output)
        elif args.command == "stage":
            print(stage_signature(args.manifest, args.stage, args.checkpoints))
        elif args.command == "validate-stage":
            paths = validate_stage(args.config, args.stage)
            if args.print_paths:
                for path in paths:
                    print(path)
        elif args.command == "marked":
            config = read_config(args.config)
            print(digest({"source": sha256(args.bam), "picard": tool_identity(config["PICARD_COMMAND"]),
                          "settings": "MarkDuplicates REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"}))
        else:
            config = read_config(args.config)
            metadata = Path(config["OUTPUT_DIR"]) / "00_metadata"
            row = next(r for r in read_tsv(metadata / "sample_manifest.tsv") if r["sample_key"] == args.sample_key)
            paths = sample_outputs(config, row)
            if args.validate:
                validate_sample_outputs(paths)
                for path in paths:
                    print(path)
            else:
                manifest = json.loads((metadata / "recovery_manifest.json").read_text())
                print(manifest["preprocess_samples"][args.sample_key])
        return 0
    except (OSError, ValueError, KeyError, StopIteration, zipfile.BadZipFile) as exc:
        print(f"RECOVERY ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
