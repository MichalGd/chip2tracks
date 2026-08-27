# Implementation and release plan

## Decisions implemented in v0.1

| Decision | Implementation |
|---|---|
| Repository name | `chip2tracks` |
| License | MIT |
| Source baseline | File-level SHA-256 manifest plus immutable import-tag recommendation |
| Assays | Separate `chipseq` and `chipmentation` run profiles |
| Adapter metadata | Auto-detect common adapters by default, retain evidence, permit explicit/custom overrides |
| Controls | Optional; absence is explicit in metadata and output labels |
| UMIs | Rejected in v0.1; planned extension point |
| Spike-in | Drosophila cells, competitive host+dm6 alignment |
| Broad domains | epic2 primary, MACS3 broad companion, local pilot required |
| Reproducibility | `ATACseq2tracks`-style support consensus by default; optional supplementary IDR |
| Engine | Bash retained until real-data validation passes |

## Reuse disposition

| Category | Components |
|---|---|
| Reuse unchanged | safe config parser, checkpoint signatures, technical-replicate merge, Bowtie2 execution, BAM policy branches, CPM/spike-in/DESeq2 track engine, metagene helpers, differential/annotation/report framework |
| Reuse with modification | samplesheet validator, preflight, preprocessing, peak calling, consensus stage, QC, MultiQC/report labels, browser session, environments and tests |
| Replace | CUT assay enums; SEACR primary path; CUT-specific SE shift/extension defaults; chromosome-length sum as MACS genome size; ATAC-derived QC contract |
| Create new | adapter evidence table, epic2 broad path, effective-genome-size reference key, true-replicate IDR, ChIP library complexity/cross-correlation/correlation QC, server audit, GitHub metadata and source-baseline provenance |

## Phases, files, dependencies, and acceptance criteria

### Phase 0 — trustworthy baseline

Files: `provenance/cutnrun2tracks_v0.2.8.sha256`,
`provenance/SOURCE_BASELINE.md`, `LICENSE`, `CITATION.cff`.

Acceptance: all inherited files verify against the manifest; a future Git import
tag points to the untouched snapshot and is never rewritten.

### Phase 1 — assay-aware metadata and preflight

Files: `chip2tracks.sh`, `config/config.conf.template`,
`config/samplesheet_template.csv`, `config/examples/*.csv`,
`scripts/validate_config.py`, `scripts/validate_samplesheet.py`,
`scripts/preflight.sh`.

Acceptance: plan mode passes independently for both profiles, PE/SE layouts and
controlled/control-free targets; invalid profiles, UMI settings, missing dm6
metadata, absent effective genome size, and incompatible controls fail clearly.

### Phase 2 — preprocessing, alignment, and spike-in

Files: `scripts/preprocess_batch.sh`, inherited `scripts/align_batch.sh`,
`scripts/mark_filter_batch.sh`, `scripts/spikein_batch.sh`,
`utilities/prepare_composite_reference.sh`.

Acceptance: common Illumina and Nextera synthetic libraries resolve or record a
useful unresolved warning; explicit adapters are honored; host and dm6 counts
are mutually exclusive under competitive alignment; q30 retained/removed BAMs
are indexed and restartable.

### Phase 3 — peaks, tracks, and reproducibility

Files: `scripts/peakcall_batch.sh`, `scripts/reproducibility_batch.sh`,
`scripts/build_consensus.py`, `scripts/coverage_batch.sh`,
`scripts/normalized_tracks_batch.sh`.

Acceptance: PE MACS3 is invoked with BAMPE; SE uses modeling unless fixed mode is
explicit; narrow, broad, and mixed targets route correctly; control-free runs
remain labeled; MACS3 fold-enrichment bigWigs are valid; two-replicate narrow
cohorts produce biological-support consensus; `RUN_IDR=true` additionally
produces IDR status and passing BED; broad cohorts produce support consensus and
signal correlation.

### Phase 4 — ChIP QC and reporting

Files: `scripts/qc_batch.sh`, `scripts/prepare_multiqc_content.py`,
`scripts/generate_multiqc_report.sh`, `scripts/generate_report.py`,
`scripts/report_batch.sh`, `docs/04_outputs_and_recovery.md`.

Acceptance: flagstat/stats, NRF/PBC, preseq, fragment length, NSC/RSC, FRiP,
target/control fingerprints, correlations, PCA, peak status, IDR/consensus, and
spike-in warnings are discoverable from `10_reports`; missing optional outputs
produce warnings rather than false success.

### Phase 5 — server pilot and dependency lock

Files: `utilities/audit_server_environment.sh`, `environment.yml`,
`environment.lock.yml`, `environment.epic2.yml`, `environment.idr.yml`,
`docs/08_server_audit.md`.

Acceptance: the audit is reviewed; reused executable versions and reference
paths are recorded; a Linux-specific explicit lock is resolved; no server-wide
environment is modified; small ChIP-seq narrow, ChIP-seq broad, and
ChIPmentation datasets pass end to end.

### Phase 6 — GitHub release

Files: `.github/workflows/ci.yml`, `.github/ISSUE_TEMPLATE/*`,
`CHANGELOG.md`, `CONTRIBUTING.md`, `CITATION.cff`, release archive/checksums.

Acceptance: CI passes syntax/unit/plan tests; the release archive is generated
from a clean commit; version/tag/CITATION agree; checksums and pilot matrix are
attached; v0.1 remains marked prevalidation until real-data acceptance is
signed off.

## Real-data pilot matrix

At minimum: PE narrow ChIP-seq with two biological replicates and matched input;
PE or SE broad ChIP-seq; PE ChIPmentation; one control-free target; one dm6
spike-in low-count warning; and one adapter override. Where permitted, use a
small public reference dataset in CI/integration testing and compare peak count,
FRiP, NSC/RSC, replicate correlation and representative track values to frozen
tolerance ranges rather than exact byte identity.

## Deferred or unresolved

- Full pooled-replicate and pseudoreplicate IDR is deferred to v0.2.
- UMI extraction/deduplication is deferred until UMI location and whitelist
  semantics can be represented in the sample sheet.
- ChIPmentation kit/adapter name is optional for common automatically detected
  adapters, but custom or ambiguous adapters still require laboratory input.
- Drosophila-to-host cell ratio and lot must be recorded; scale-factor behavior
  must be validated against the lab protocol and a no-biological-change pilot.
- epic2 is the best-aligned v0.1 choice for diffuse domains, not a universal
  ground truth. Mark-specific pilot comparison against MACS3 broad remains an
  acceptance condition.
- Nextflow/Snakemake migration is considered only after the Bash release passes
  the real-data matrix, so engine migration cannot obscure scientific changes.
