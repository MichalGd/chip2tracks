# Outputs and recovery

Results use numbered directories from `00_metadata` through `10_reports`.
Sample files use `<output_prefix>.bioR<replicate>`. Cohort IDs contain a readable
prefix and an eight-character SHA-256 suffix; the complete key is retained in
`cohort_manifest.tsv`.

Per-caller results are recorded in
`05_peaks/per_sample/<sample_key>/caller_status.tsv`; the stable run-wide table
is `05_peaks/per_sample/peakcall_status.tsv`. A continuation-mode exclusion from
consensus is recorded under each cohort as `excluded_peak_samples.tsv`. These
records distinguish a biological zero-peak result (`EMPTY`) from a caller
execution failure (`ERROR`) and must be reviewed with the ordinary QC outputs.

## Output map

```text
<OUTPUT_DIR>/
├── 00_metadata/                    run, sample, cohort and reference manifests
├── 01_fastq_qc/                    raw/trimmed FastQC and MultiQC
├── 02_trimmed_fastq/               merged/trimmed FASTQs; retained by default
├── 03_alignment/                   sorted, marked, filtered and analysis BAMs
├── 04_tracks/                      CPM, DESeq2-derived and spike-in tracks
├── 05_peaks/                       per-sample and consensus peaks
├── 06_qc/                          alignment, fragment, FRiP, control and metagene QC
├── 07_annotation/                  consensus and differential annotations
├── 08_differential/                primary and control-aware sensitivity results
├── 09_browser/                     UCSC and IGV assets
├── 10_reports/                     unified MultiQC and lightweight HTML reports
├── logs/                           stage/tool logs
└── .checkpoints/                   signature-and-output JSON checkpoints
```

The priority coverage tree is:

```text
04_tracks/
├── cpm/                            analysis-policy *.CPM.{bw,bedGraph}
│   ├── permissive/                 MAPQ 0, duplicates retained
│   ├── intermediate/               MAPQ 0, duplicates removed
│   └── stringent/                  MAPQ 30, duplicates removed
└── deseq2_robust_cpm/
    └── stringent/<cohort>/         stringent robust-CPM tracks and tables
```

The complete family matrix, scaling formulas, file-format behavior, and track
selection guidance are in
[Genomic coverage tracks and normalization](12_tracks_and_normalization.md).

The root README gives a more detailed
[stage and output overview](../README.md#stages-and-main-outputs), and
[Pipeline stages](07_pipeline_stages.md) maps every stage to its declared
checkpoint output.

`04_tracks/cpm/mapping_composition.tsv` and its companion definitions table
describe the mapping composition of every coverage policy. Each row reports
the total signal observations, MAPQ 0 and MAPQ <30 observations, and Bowtie2
`XS`-tagged candidate multimappers. A bedGraph and bigWig generated from the
same policy share these counts because both represent the same filtered BAM
signal.

## Checkpoints

Each stage checkpoint is JSON containing the complete run signature and hashes
of declared outputs. The signature covers the sanitized CSV, resolved config,
workflow version, scripts, reusable `common/` modules, and reference manifest.
A changed or missing output invalidates the checkpoint.

`00_metadata/stage_timing.tsv` records UTC start/end timestamps and elapsed
seconds for completed, skipped, reused, and failed stages in the current
invocation. `CHECKPOINT_PARALLEL_JOBS` controls hash-validation concurrency;
`CHECKSUM_PARALLEL_JOBS` controls optional terminal checksum generation. The
timing table is intentionally excluded from `final_checksums.sha256` because
its final row is written after checksum generation finishes.

Run `--plan` to validate the metadata model and inspect
`00_metadata/planned_stages.tsv` without checking files/tools. Run
`--preflight-only` for the complete input, tool, and reference audit.

Reuse validated outputs before a stage, then force that stage and all later
stages:

```bash
bash chip2tracks.sh --config /path/to/config.conf --from-stage qc
```

This is an explicit recovery override. Before the requested stage, each existing
checkpoint is validated against its stored output sizes and SHA-256 hashes and
its signature adoption is recorded in the checkpoint JSON. If any earlier
checkpoint or output is missing or changed, recovery stops instead of silently
rerunning or accepting it. The named stage and all later stages always rerun.
Select a starting stage at or before the earliest output that could be affected
by the configuration or code change.

If a stage command completes every output but fails during terminal bookkeeping
before its checkpoint is written, preserve the outputs and validate them with
the relevant native tools before recovery. A replacement checkpoint may be
written with `scripts/checkpoint.py write` only after expected file counts,
indexes, links, and content integrity have all been confirmed and the incident
has been recorded. Deploy the code correction, then use `--from-stage` on the
next stage so the replacement checkpoint is explicitly adopted into the new
workflow signature. Never use this procedure for partial or merely
size-checked outputs.

Stop after a stage without editing the workflow:

```bash
bash chip2tracks.sh --config /path/to/config.conf --stop-after reproducibility
```

## Stable reporting interfaces

Metagene outputs are retained under `06_qc/metagene`. Each sample/gene-set/mode
directory contains the deepTools matrix, exported profile values, sorted BED,
PNG/PDF profile, PNG/PDF heatmap, and task metadata. The run-wide
`artifacts.tsv` is the stable reporting interface.

The final interactive report is
`10_reports/chip2tracks_multiqc_report.html`; its parsed data directory,
exported plots, log, custom-content manifest, and status table are retained
beside it. `10_reports/pipeline_report.html` remains a dependency-light summary,
with `run_summary.tsv` and `warning_summary.tsv` as stable tabular companions.
Both reports include coverage-family mapping composition, and
`10_reports/coverage_mapping_composition.tsv` is its stable final-report table.
The counts describe the BAM policy underlying both the bedGraph and bigWig,
not information recovered from either coverage file after conversion.
`10_reports/differential_occupancy_summary.tsv` similarly provides the stable
all-variant differential table, including completed comparisons and explicit
disabled, skipped, failed, or missing variants.
Machine-readable TSVs in `00_metadata`, `04_tracks`, `05_peaks/consensus`,
`06_qc`, and `08_differential` remain authoritative and should be preferred over
parsing HTML or filenames. Browser definitions are written under `09_browser/`.
See [Quality control](11_quality_control.md),
[Differential occupancy](14_differential_occupancy.md), and
[Annotation](15_annotation.md) for the stable scientific outputs in those
directories.

### UCSC custom-track descriptors

`09_browser/ucsc/` contains:

- numbered family `.txt` files, one for each generated bigWig family;
- `all_bigwig_tracks.txt`, containing every family separated by an ignored
  `# family=...` comment and a blank line;
- `trackDb.txt`, an identical backward-compatible custom-track-text filename;
- `track_family_manifest.tsv`, `status.tsv`, and `README.txt`.

Every descriptor occupies exactly one physical line and begins with `track
type=bigWig`. Required `name`, `description`, and `bigDataUrl` attributes are
validated before the annotation stage succeeds. Empty and `#` lines are valid
UCSC separators and are ignored by its parser. Despite its historical name,
this `trackDb.txt` is a custom-track submission file, not a multi-file track-hub
`trackDb` stanza file.

Public UCSC servers cannot read `/home/...` paths. Set
`UCSC_BIGDATA_URL_BASE` to the HTTP/HTTPS/FTP URL corresponding to the analysis
output root, and ensure the web server supports byte-range requests. With an
empty base, the files remain syntax-valid but `status.tsv` marks them
`local_path_not_ucsc_retrievable`.

Regenerate browser descriptors for a completed run without rerunning analysis:

```bash
PYTHON_COMMAND=/opt/miniconda/envs/chip2tracks-0.1.0/bin/python3 \
  bash /opt/bioinformatics/workflows/chip2tracks/current/utilities/regenerate_ucsc_tracks.sh \
  --output-dir /absolute/path/to/chip2tracks \
  --url-base https://tracks.example.org/path/to/analysis
```

The optional URL override avoids changing the completed run's historical
resolved configuration. Omit it only when that configuration already contains
the correct public base URL.

Format and hosting behavior follow the official
[UCSC custom-track specification](https://genome.ucsc.edu/goldenPath/help/customTrack.html)
and [bigWig documentation](https://genome.ucsc.edu/goldenPath/help/bigWig.html).

`RUN_MULTIQC=true` enables both the preprocessing FastQC aggregation and this
final unified report. The final scan includes retained FastQC, Trim Galore,
Bowtie2, Picard, samtools, and preseq outputs plus ChIP-specific custom-content
tables and selected QC plots. Native deepTools parsing is excluded because it
is unreliable for some MultiQC 1.35 tables; original deepTools files are kept,
and selected plots plus the metagene summary are supplied as custom content.

`MULTIQC_EXPORT_PLOTS=false` is recommended for routine server runs. The HTML
report remains complete, but MultiQC avoids exporting every plot to extra
PNG/SVG/PDF files. Enable it only when those static assets are required.

Reports can be recovered from a completed output directory without rerunning
alignment, filtering, peak calling, normalization, or differential models:

```bash
bash /path/to/chip2tracks/utilities/regenerate_reports.sh \
  --output-dir /absolute/path/to/chip2tracks_results
```

Regeneration writes `10_reports/report_checksums.sha256`. It intentionally does
not rewrite the historical `00_metadata/final_checksums.sha256` or workflow
checkpoints. If automatic cleanup already removed an intermediate, MultiQC can
only report from the retained logs and QC artifacts; normal final outputs are
otherwise sufficient for the ChIP-specific summary sections.

Normalization continuation is recorded in
`04_tracks/normalized_track_family_status.tsv`. A failed factor calculation has
a dedicated log under `logs/normalized_tracks/`, a `SKIPPED.json` marker in the
affected family, and `tables/consensus_count_sums.tsv` when counting completed.
With `REQUIRE_ALL_ENABLED_TRACKS=false`, unaffected cohorts and later stages
continue; strict mode records `FAILED.json` and stops.

## Cleanup and retention

Cleanup begins only after the final report exists, targets only explicit
children of the resolved output root, and records every deletion in
`00_metadata/cleanup_manifest.tsv`.

`ENABLE_AUTOMATIC_CLEANUP=false` is the default. When cleanup is explicitly
enabled, the `KEEP_*` switches control individual intermediate families and the
workflow invalidates checkpoints whose declared outputs were deleted. A later
full rerun rebuilds those stages; an incompatible `--from-stage` request fails
early rather than reusing an incomplete checkpoint.
Final tracks, peaks, QC, differential results, annotations, browser assets,
reports, metadata, logs, and checksums are not cleanup targets under the default
policy.
