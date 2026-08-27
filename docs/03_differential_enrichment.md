# Controls and differential enrichment

Controls have three separate uses:

1. optional matched background during MACS3/epic2 peak calling;
2. target-control QC and fingerprint plots;
3. optional, separately labelled sensitivity models.

The primary analysis counts raw target fragments over the condition-unbiased
cohort consensus. Controls are not biological replicates and are not subtracted
from DESeq2 input. DESeq2Enrichment supports validated block columns and all
eligible pairwise condition contrasts.

DiffBind receives `bamReads`, `bamControl`, and `ControlID`. Its primary run
sets `bSubControl=FALSE`; the sensitivity run explicitly uses scaled control
subtraction. Version 0.1 skips DiffBind when arbitrary blocking columns are
requested, leaving the block-aware DESeq2Enrichment analysis primary.

The optional interaction model tests whether the condition change in target
signal exceeds the condition change in matched control signal. Version 0.1
requires exactly two conditions and one-to-one condition-matched controls. It
uses joint DESeq2 normalization and must be pilot-validated.

Outputs are separated into `primary_target_only`,
`sensitivity_control_subtracted`, and
`sensitivity_target_control_interaction`; only the first is primary.

For exact eligibility, design formulas, normalization, contrast direction,
result trees, and interpretation, see
[Differential occupancy analysis](14_differential_occupancy.md). Replicate and
shared-control metadata rules are in
[Replicates, controls, and experimental design](13_replicates_and_design.md).

## Final differential summary

`10_reports/differential_occupancy_summary.tsv` is the stable run-wide summary
for all four analysis variants: primary target-only DESeq2Enrichment, primary
target-only DiffBind, control-subtracted DiffBind sensitivity, and the
target-control interaction sensitivity model. Every cohort/variant is retained
as one or more rows, including variants that are `DISABLED`, `SKIPPED`,
`FAILED`, or `NOT_AVAILABLE`.

Completed comparison rows report the condition direction, tested and
significant consensus regions, direction-specific significant counts,
normalization/control mode, thresholds, and paths to the complete and
significant result tables. The same table is embedded in the lightweight HTML
and final MultiQC report. Significant counts are comparison-specific calls and
must not be added across methods or contrasts as though they were unique
genomic regions.
