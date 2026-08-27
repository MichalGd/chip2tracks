# Replicates, controls, and experimental design

[Documentation index](README.md) | [Differential occupancy](14_differential_occupancy.md)

The samplesheet distinguishes sequencing units, biological libraries, controls,
and conditions. Correct metadata are essential because the workflow uses these
fields to decide what may be merged, what remains independent, and which samples
share a peak and differential-analysis universe.

## Technical and biological replicates

One CSV row is one sequencing unit. Rows with the same `sample_id` and
`replicate` and different `tech_replicate` values are technical units of one
biological library. They may differ only in `fastq_1`, `fastq_2`, and
`tech_replicate`; every biological and analysis field must match.

The FASTQs are concatenated before trimming, FastQC, alignment, and all
downstream work. Technical units therefore increase depth but do not increase
the biological replicate count or the degrees of freedom in differential
analysis.

Biological replicates use different `replicate` values and remain separate for:

- trimming and alignment outputs;
- duplicate marking and filtering;
- coverage tracks and peak calling;
- consensus support;
- QC and replicate correlation; and
- raw counting and differential occupancy.

Do not merge biological BAMs before statistical analysis. A pooled browser
track can be useful for presentation, but v0.1 does not create pooled-condition
tracks and such a file must not be presented as an independent replicate.

## One samplesheet, one target universe

Version 0.1 permits one compatible target/antibody and peak universe per
samplesheet. Target rows must agree on genome, assay profile, factor, antibody,
layout, target class, analysis duplicate policy, primary caller, and primary
peak class. Spike mode/reference/stage/lot also participate in the cohort
identity when calibration is enabled.

Conditions, treatments, donors, batches, and biological replicates belong
together when they form one valid comparison. Use separate runs for different
antibodies, factors, narrow/broad target policies, layouts, or genomes. This
prevents unrelated regions from being combined into one consensus and count
matrix.

## Matched and shared controls

Input, IgG, and mock controls are supported and optional. A target's
`control_id` links to the control `sample_id`.

Control resolution follows this order:

1. use exactly one control with the requested `sample_id` and the same
   biological replicate number;
2. if none exists and `ALLOW_SHARED_CONTROLS=true`, use exactly one row whose
   condition is `shared`; and
3. otherwise fail metadata validation rather than guessing a replicate.

Target and selected control must match genome, assay profile, layout, treatment,
cell type, batch, and spike-in ratio/stage/lot. A replicate-matched control must
also match condition; a `shared` control is deliberately condition-independent.
Shared controls should be used only when the same experimental input is
scientifically appropriate for every linked target.

Controls contribute to MACS3/epic2 background modeling, target/control QC, and
optional sensitivity analyses. They do not become biological target replicates
and are not columns in the primary target-only count model.

## Consensus and reproducibility

Technical units are already collapsed when per-library peaks are called. The
default consensus then retains genomic intervals supported by peaks from at
least `CONSENSUS_MIN_BIOLOGICAL_SAMPLES` distinct target sample keys (two by
default). The consensus is the fixed region universe for normalized tracks,
FRiP, annotation, and differential analysis.

`ALLOW_SINGLE_SAMPLE_CONSENSUS=false` prevents one target from being promoted
to a reproducible cohort by default. Enabling it can support a technical pilot,
but it does not establish biological reproducibility.

`RUN_IDR=false` is the default. When enabled, pairwise true-biological-replicate
IDR is a supplementary result for narrow MACS3 cohorts. It does not replace the
consensus, and v0.1 does not implement pooled or pseudoreplicate IDR. Broad
domains use biological-support consensus and signal concordance rather than
IDR. See the [IDR framework paper](https://doi.org/10.1214/11-AOAS466) when
designing a full IDR analysis.

The following files report reproducibility explicitly:

- `05_peaks/consensus/consensus_status.tsv`;
- `05_peaks/consensus/<cohort>/excluded_peak_samples.tsv`;
- `05_peaks/consensus/<cohort>/.../*.consensus.bed`;
- `05_peaks/consensus/<cohort>/.../consensus_support.tsv.gz`; and
- `05_peaks/reproducibility/status.tsv`.

## Conditions, batch, and donor

The primary DESeq2 enrichment model requires at least
`DIFFERENTIAL_MIN_REPLICATES_PER_CONDITION` target libraries in at least two
conditions. The default minimum is two; three or more biological replicates per
condition are preferable when feasible, especially for heterogeneous samples.

Use `DIFFERENTIAL_BLOCK_COLUMNS=batch`, `donor`, or a comma-separated
combination when those covariates are part of a valid design. The model is
constructed as:

```text
~ block_1 + block_2 + ... + condition
```

Missing, constant, or rank-deficient block designs fail rather than silently
dropping a covariate. A batch or donor completely confounded with condition
cannot be corrected computationally; the experimental design must provide the
necessary crossing or pairing.

DiffBind v0.1 does not map arbitrary block columns. If blocks are requested,
the DiffBind variant records `SKIPPED`, while the block-aware target-only DESeq2
analysis remains primary.

## Contrast direction

Without an explicit order, condition levels follow first appearance in the
sanitized samplesheet. Set `DIFFERENTIAL_CONDITION_ORDER` to make the complete
order reproducible, or `DIFFERENTIAL_REFERENCE_CONDITION` to place one eligible
condition first. The workflow exports every pair among eligible conditions and
records `numerator`, `reference`, and direction-specific significant counts.

Never infer direction only from a folder name. Confirm the `numerator` and
`reference` fields in `comparison_summary.tsv` or
`10_reports/differential_occupancy_summary.tsv` before interpreting the sign of
a fold change.

## Pre-analysis checklist

- Technical units represent the same biological library and differ only in
  sequencing files/lane identity.
- Biological replicates are genuinely independent.
- Every condition intended for testing has sufficient biological replication.
- Target and control metadata describe the same experiment and spike protocol.
- Batch/donor variables are not confounded with condition.
- All samples use the same genome, blacklist, target definition, and signal
  unit within the run.
- Narrow versus broad target class and caller policy are biologically justified.
- QC and within-condition concordance are reviewed before differential results.
