# Changelog

## Unreleased

- Corrected Git executable modes for the workflow entrypoint and all tracked
  shell utilities/tests, and made CI reject non-executable shell archives. This
  ensures `git archive` produces a directly runnable shared release.
- Fixed duplicate-policy-aware signal counting across filtering metrics, CPM
  scale factors, coverage mapping-composition summaries, QC/FRiP denominators,
  paired-end fragment QC, and spike-in host/control ratios. Duplicate-retained
  branches now include duplicate-flagged observations, while duplicate-removed
  branches continue to exclude them.
- Added paired-end and single-end regression coverage for signal-count masks and
  strengthened filtering/coverage tests so retained and removed branches must
  produce distinct counts when duplicates are present.
- Corrected the public browser-track prefix default from inherited `CUT` to
  `CHIP`.
- Made the shared launcher resolve `current` to an immutable release at process
  start, preventing a later atomic promotion from changing code underneath an
  active run.
- Isolated phantompeakqualtools/SPP and preseq from the modern main environment
  and allowed the configured cross-correlation command to select its own R
  runtime, matching the validated all-user server deployment.
- Expanded the task-oriented documentation with dedicated guides for reference
  reuse, blacklist and alignment filtering, ChIP-specific QC, genomic coverage
  normalization, biological/technical replicates, controls and design,
  differential occupancy, and annotation.
- Added a track-family matrix covering CPM, DESeq2 consensus, robust CPM,
  Drosophila spike-in, and MACS3 fold-enrichment outputs, including bedGraph vs
  bigWig behavior and multimapper/ambiguity reporting.
- Documented exact implemented behavior and explicit v0.1 boundaries, including
  the lighter annotation layer relative to ATACseq2tracks. No workflow behavior
  or configuration default changed.

## chip2tracks 0.1.0 - 2026-08-25

- Forked the checksum-recorded local `cutnrun2tracks` 0.2.8 source into an
  independent MIT-licensed repository.
- Added explicit ChIP-seq and ChIPmentation profiles, adapter auto-detection
  evidence/overrides, optional matched controls, and v0.1 UMI rejection.
- Added effective genome sizes, MACS3 ChIP modes, epic2 broad domains, MACS3
  fold-enrichment tracks, true-replicate IDR, NRF/PBC, duplicate-retained
  preseq/cross-correlation, and replicate correlation/PCA.
- Made the `ATACseq2tracks`-style biological-support consensus the default
  reproducibility method; pairwise true-replicate IDR is optional and disabled
  by default.
- Configured competitive host+dm6 cell spike-in support and added a read-only
  server environment/reference audit.
- Added priority analysis, permissive (MAPQ 0/duplicates retained), intermediate
  (MAPQ 0/duplicates removed), and stringent (MAPQ 30/duplicates removed) CPM
  bigWig/bedGraph families plus stringent consensus-derived robust CPM.
- Added per-policy coverage mapping-composition tables and final-report/MultiQC
  sections for MAPQ 0, MAPQ <30, and Bowtie2 `XS`-tagged candidate multimappers.
- Added a consolidated differential-occupancy report covering every primary
  and sensitivity variant, including comparison counts, result paths, and
  explicit disabled/skipped/failed states in HTML, TSV, and MultiQC.
- Made uncompressed bedGraph and bigWig retention independently configurable,
  added MACS3-estimated SE coverage extension with a recorded fallback, and
  made spike-in optional/off in the public template.
- Added real spike-in parallel execution, lightweight CPU-budget preflight,
  one-universe-per-samplesheet validation, dead-option removal, and
  cleanup-aware checkpoint invalidation.

The entries below are inherited `cutnrun2tracks` history and are preserved for
source provenance; their version numbers are not `chip2tracks` releases.

## 0.2.8 - 2026-08-25

- Treat MultiQC colour-conversion diagnostics as non-fatal when MultiQC exits
  successfully and creates the required report and data directory.
- Retain strict failure handling for non-zero MultiQC exits, module crashes,
  validation errors, and missing report artifacts.
- Add a reporting regression test covering the non-fatal colour diagnostic seen
  in completed human and mouse CUT&Tag report recovery.

## 0.2.7 - 2026-08-25

- Added a final unified MultiQC report over the complete retained workflow
  output, matching the reporting role used by ATACseq2tracks while preserving
  CUT-specific sample, cohort, peak, normalization, and differential semantics.
- Added deterministic MultiQC custom-content tables for retained observations,
  alignment/filtering counts, FRiP, peak calls, consensus sets, normalized-track
  families, spike-in calibration, differential status, and comparison summaries.
- Added selected target/control fingerprints, descriptive TSS profiles, and
  differential PCA/dispersion images as MultiQC custom content while excluding
  the fragile native deepTools parser; authoritative source files are unchanged.
- Made the report stage validate its lightweight HTML/TSV outputs and its final
  MultiQC HTML before cleanup can proceed.
- Added `utilities/regenerate_reports.sh` to recover the lightweight and unified
  reports from completed analyses without rerunning upstream stages. A separate
  report checksum inventory avoids rewriting historical workflow checksums.
- Expanded the lightweight report to include normalization, differential, and
  preseq warnings, and added end-to-end reporting regression coverage.

## 0.2.6 - 2026-08-25

- Fixed nearest-gene annotation failures caused by inconsistent chromosome
  ordering between consensus and GTF-derived BED inputs.
- Genome-sorted both inputs with the configured chromosome-sizes file and
  passed that same order to `bedtools closest -g`, making annotation independent
  of shell locale and suitable for hg38, mm39, and compatible custom genomes.
- Added an executable annotation-order regression test.

## 0.2.5 - 2026-08-25

- Made consensus-normalization failures cohort-local when
  `REQUIRE_ALL_ENABLED_TRACKS=false`, so a zero-count or otherwise
  non-normalizable cohort is recorded and skipped without terminating
  unaffected cohorts, QC, metagene, annotation, or reporting.
- Kept normalization scientifically strict: samples with zero consensus counts
  are not silently removed and no scaling factors are fabricated. Strict
  fail-fast behavior remains available with `REQUIRE_ALL_ENABLED_TRACKS=true`.
- Added per-sample `consensus_count_sums.tsv`, per-family status tables, and
  dedicated factor-calculation logs to make normalization exclusions auditable.
- Made differential analysis recognize an upstream normalization skip as an
  expected cohort skip while continuing to treat unexpectedly missing count
  tables as failures.
- Added continuation/strict-mode regression coverage with one failing and one
  successful cohort.

## 0.2.4 - 2026-08-24

- Added configurable `PEAKCALL_FAILURE_POLICY=continue|fail`; continuation mode
  records per-caller `SUCCESS`, `EMPTY`, or `ERROR` results without allowing one
  failed sample or auxiliary caller to terminate all other samples.
- Made consensus construction exclude only failed/empty primary peak
  contributions, proceed when enough successful biological samples remain, and
  record every exclusion and reason in machine-readable TSV files.
- Preserved all BAM, coverage, QC, metagene, and downstream consensus-counting
  eligibility for a peak-call-excluded sample; peak failure is reported and is
  not treated as an automatic biological QC exclusion.
- Added completed-with-warnings reporting and regression coverage for both
  continuation and strict peak-calling policies.

## 0.2.3 - 2026-08-24

- Scoped each filtering temporary directory and its cleanup trap to a
  subshell, preventing the trap from firing again after its local `tmp`
  variable had expired when a parallel worker returned.
- Added an executable filtering-stage regression fixture that exercises all
  four BAM branches and fails if temporary cleanup leaks into the worker.
- Documented recovery of fully validated filtering outputs without repeating
  completed preprocessing, alignment, or filtering work.

## 0.2.2 - 2026-08-24

- Fixed canonical-contig filtering by applying region selection directly to
  the indexed duplicate-marked BAM instead of an unindexed temporary BAM.
- Reused validated marked BAMs on filtering recovery, avoiding unnecessary
  repeated Picard duplicate marking.
- Made `--from-stage` validate and adopt unchanged earlier-stage checkpoints
  across a workflow signature update, while still forcing the named and later
  stages.
- Added regression coverage for canonical filtering and checkpoint adoption.

## 0.2.1 - 2026-08-24

- Added validated `THREADS_FASTQC` and `THREADS_TRIMGALORE` settings and passed
  them to FastQC and Trim Galore, matching the proven ATACseq2tracks resource
  model for high-capacity servers.
- Made the shared bounded worker pool work-conserving by refilling a slot when
  any worker finishes instead of waiting for submission order.
- Fixed `set -u` failures caused by referencing variables within the same
  `local` declaration in the workflow dispatcher and preprocessing worker.
- Added regression coverage for preprocessing thread validation and parallel
  pool scheduling.

## Unreleased - documentation

- Expanded the GitHub README with a Mermaid workflow diagram, installation,
  stage, output, validation, recovery, and limitations guidance.
- Added a task-oriented documentation index and an exact pipeline-stage guide.
- Clarified output discovery and checkpoint/cleanup recovery without changing
  workflow behavior or configuration defaults.
- Added a server installation runbook covering isolated environment cloning,
  pinned SEACR, shared host-reference reuse, derived canonical/dm6 coordinate
  files, validation, pilot deployment, and guarded promotion.
- Added an explicit all-user deployment model with immutable runtime
  permissions, an optional `/usr/local/bin` launcher, parent-path auditing, and
  acceptance tests executed as an ordinary server user.

## 0.2.0 - 2026-08-22

- Added a reusable, headless deepTools metagene module for TSS, TES, and scaled
  gene-body profiles with matching PNG/PDF heatmaps.
- Added deterministic GTF-to-BED12 protein-coding reference preparation,
  blacklist/contig/length/overlap filters, and auditable exclusion tables.
- Added HPA broadly expressed reference-set and Ensembl BioMart mouse-ortholog
  builders.
- Added manifest validation, explicit CPM/spike-in track selection, reporting
  tables, resource bounds, and signature-aware workflow integration.

## 0.1.0 - 2026-08-22

- Initial independent CUT&RUN/CUT&Tag workflow implementation.
- Added strict samplesheet/config validation and target-specific cohort IDs.
- Added PE fragment-aware alignment, filtering, coverage, SEACR, and MACS3 paths.
- Added five independent coverage families and generalized spike-in calibration.
- Added target-only differential analysis plus separately labelled control-aware
  sensitivity analyses.
