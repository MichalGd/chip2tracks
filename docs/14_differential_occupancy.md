# Differential occupancy analysis

[Documentation index](README.md) | [Replicates and design](13_replicates_and_design.md)

Differential occupancy asks whether raw ChIP fragment/read counts over a fixed
consensus region universe differ between biological conditions. It does not
test bigWig or bedGraph values. The workflow separates primary target-only
models from control-aware sensitivity analyses so their biological meanings are
not conflated.

## Common input and eligibility

For each target cohort, raw integer counts are generated over the
biological-support consensus from the analysis-policy BAMs. All target samples
help define the consensus and count matrix. Statistical modeling retains only
conditions with at least `DIFFERENTIAL_MIN_REPLICATES_PER_CONDITION` biological
samples; at least two eligible conditions are required.

The relevant configuration is:

```text
RUN_DESEQ2_ENRICHMENT=true
RUN_DIFFBIND=true
DIFFERENTIAL_MIN_REPLICATES_PER_CONDITION=2
DIFFERENTIAL_ALPHA=0.05
DIFFERENTIAL_MIN_ABS_LOG2FC=0
DIFFERENTIAL_BLOCK_COLUMNS=
DIFFERENTIAL_REFERENCE_CONDITION=
DIFFERENTIAL_CONDITION_ORDER=
DIFFERENTIAL_NORMALIZATION=deseq2
RUN_CONTROL_SUBTRACTED_SENSITIVITY=true
RUN_TARGET_CONTROL_INTERACTION=false
```

A consensus-normalization family may be generated even when no differential
contrast is eligible. Conversely, differential analysis cannot run if consensus
counting was skipped or failed. Status files distinguish these cases.

## Analysis variants

| Analysis family | Method | Role | Control treatment | Important limitation |
|---|---|---|---|---|
| `primary_target_only` | DESeq2Enrichment | primary | controls excluded from count model | assumes chosen normalization is valid |
| `primary_target_only` | DiffBind/DESeq2 | primary peer | complete control set attached when available; `bSubControl=FALSE` | arbitrary block columns unsupported in v0.1 |
| `sensitivity_control_subtracted` | DiffBind/DESeq2 | sensitivity | scaled control subtraction | requires matched controls for every target |
| `sensitivity_target_control_interaction` | DESeq2 interaction | sensitivity | tests target condition effect minus control condition effect | exactly two conditions and one-to-one matched controls in v0.1 |

Only the target-only analyses are primary. Agreement across variants can
increase confidence, but disagreement is diagnostic rather than a reason to
select whichever method gives more significant regions.

## Primary DESeq2Enrichment

The model is `~ condition` or `~ blocks + condition`. It supports any number of
eligible conditions and exports every pairwise contrast. With
`DIFFERENTIAL_NORMALIZATION=deseq2`, size factors use DESeq2 `poscounts`, which
is suitable for sparse peak-count matrices while retaining the usual relative
normalization assumption.

With `DIFFERENTIAL_NORMALIZATION=spikein`, size factors are derived from
retained spike observations divided by the declared spike-to-host ratio and
then centered geometrically. This mode requires valid spike-in results for the
cohort. It is intended for experiments where the external reference was added
consistently and a global target shift is plausible.

Each contrast applies adjusted P-value `DIFFERENTIAL_ALPHA` and absolute
log2-fold-change `DIFFERENTIAL_MIN_ABS_LOG2FC`. The complete result retains
independent-filtered `padj=NA` rows; only finite rows passing both configured
criteria enter the significant table. Zero significant regions is a successful
result, not a workflow failure.

Principal outputs are:

```text
08_differential/<cohort>/<peak_class>/primary_target_only/deseq2_enrichment/
|-- raw_counts.tsv.gz
|-- normalized_counts.tsv.gz
|-- size_factors.tsv
|-- comparison_summary.tsv
|-- pca.tsv
|-- pca.png
|-- dispersion.png
|-- deseq2_object.rds
|-- session_info.txt
`-- comparisons/<numerator>_vs_<reference>/
    |-- all.tsv.gz
    `-- significant.tsv.gz
```

See the [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
for statistical assumptions and result-column definitions.

## Primary DiffBind peer analysis

DiffBind receives target BAMs, the same fixed consensus, and matched control
BAMs when every eligible target has one. If the cohort has incomplete control
coverage, the primary sheet omits controls for the complete cohort. The primary
run uses `bSubControl=FALSE`, so target counts remain the tested signal. The
workflow disables DiffBind's additional blacklist because upstream BAMs and
consensus regions are already blacklist-filtered; it enables greylisting when a
complete matched-control set is present.

DiffBind performs DESeq2 contrasts and writes an R object, complete/significant
TSV files, a comparison summary, session information, and the generated
DiffBind samplesheet under:

```text
08_differential/<cohort>/<peak_class>/primary_target_only/diffbind/
```

This peer analysis is not expected to return identical statistics to the
direct DESeq2 implementation because DiffBind manages its own counting and
normalization workflow. See the
[DiffBind Bioconductor documentation](https://bioconductor.org/packages/release/bioc/html/DiffBind.html).

## Control-aware sensitivity analyses

The control-subtracted DiffBind variant sets `bSubControl=TRUE` and
`bScaleControl=TRUE`. It is skipped when any eligible target lacks a matched
control. Subtraction can reduce background-driven differences but also adds
control noise and changes the estimand, so it remains a sensitivity output.

The interaction variant jointly counts targets and one-to-one matched controls
and fits:

```text
~ library_type + condition + library_type:condition
```

Its coefficient represents the condition effect in target signal minus the
condition effect in control signal. Version 0.1 requires exactly two conditions,
sufficient target replication, unique matched controls, and condition-matched
control rows. It uses joint DESeq2 `poscounts` normalization and must be
pilot-validated before being promoted to a primary result.

## Annotation and final summary

The later annotation stage creates annotated copies of gzipped DESeq2 and
interaction result tables when nearest-gene/cCRE annotation is enabled. Original
statistical tables remain unchanged. DiffBind's uncompressed TSV results are not
annotated automatically in v0.1. See [Annotation](15_annotation.md) for the
implemented columns and limitations.

`10_reports/differential_occupancy_summary.tsv` is the stable run-wide index. It
contains every enabled or disabled primary and sensitivity variant and retains
`SUCCESS`, `DISABLED`, `SKIPPED`, `FAILED`, and `NOT_AVAILABLE` states. Completed
comparison rows record:

- cohort, factor, and peak class;
- analysis family, role, method, and control mode;
- comparison ID, numerator, and reference;
- tested and significant region counts;
- significant regions higher in either direction;
- normalization and thresholds; and
- paths to all and significant result tables.

The same information appears in the lightweight report and final MultiQC
report. Counts are method- and contrast-specific and must not be summed as if
they were unique genomic loci.

## Interpretation checklist

1. Confirm contrast direction, design, size factors, and session information.
2. Review raw depth, complexity, FRiP, peak-call status, consensus support,
   within-condition correlation, PCA, and batch/donor balance.
3. Inspect MA-like effect patterns and dispersion behavior, not only the number
   of significant regions.
4. Treat `padj=NA` as untested after independent filtering, not significant or
   non-significant evidence.
5. Interpret target-only results as primary and control-aware variants as
   sensitivity analyses unless the project prespecified another model.
6. Do not interpret nearest-gene or cCRE overlap as a proven regulatory target.
