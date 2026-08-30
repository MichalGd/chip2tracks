# Peak and differential-result annotation

[Documentation index](README.md) | [Differential occupancy](14_differential_occupancy.md)

Annotation is descriptive and runs after peak calling/statistical testing. It
does not alter peaks, counts, P-values, adjusted P-values or fold changes.

## Implemented layers

The lightweight cohort view is retained:

- nearest whole-gene span for each consensus interval;
- all cCRE reference overlaps for each consensus interval.

The comprehensive feature summary additionally processes every successful
per-sample/caller peak set and, optionally, every cohort consensus. It assigns:

- promoter, enhancer, exon, intron, gene-end, other-regulatory, intergenic or
  unclassified context;
- all overlapping feature records;
- one deterministic primary feature using the configured precedence;
- primary gene identifiers/names where applicable;
- strand-aware signed distance to the nearest TSS;
- cCRE-derived enhancer context when a compatible cCRE reference is available.

Configuration:

```text
RUN_SIMPLE_PEAK_ANNOTATION=true
RUN_CCRE_ANNOTATION=true
RUN_FEATURE_ANNOTATION_SUMMARY=true
PEAK_ANNOTATION_PROMOTER_UPSTREAM=2000
PEAK_ANNOTATION_PROMOTER_DOWNSTREAM=500
PEAK_ANNOTATION_GENE_END_WINDOW=2000
PEAK_ANNOTATION_FEATURE_PRECEDENCE=promoter,enhancer,exon,intron,gene_end,other_regulatory,intergenic,unclassified
PEAK_ANNOTATION_PLOT_FORMATS=png,pdf,svg
PEAK_ANNOTATION_INCLUDE_CONSENSUS=true
```

GTF and chromosome-size references are required per genome. cCRE references are
required only when cCRE annotation is enabled. Gzip-compressed GTF/cCRE inputs
are accepted.

## Outputs

`07_annotation/feature_summary/` contains:

- `peak_feature_assignments.tsv.gz`: one deterministic assignment per peak;
- `peak_feature_all_overlaps.tsv.gz`: exhaustive peak/feature overlaps;
- `peak_feature_counts.tsv`, `peak_feature_fractions.tsv` and
  `peak_feature_bp_coverage.tsv`;
- `peak_feature_summary.tsv` and `peak_feature_colors.tsv`;
- `peak_annotation_status.tsv`;
- composition plots in the selected formats.

Cohort-local lightweight files remain under
`07_annotation/<cohort>/consensus/`.

Differential result tables containing `region_id` receive adjacent annotated
copies. Both gzip-compressed and uncompressed TSVs are supported, including
DiffBind outputs. Comprehensive primary category/gene/TSS fields and lightweight
nearest-gene/cCRE fields are propagated when available; original statistical
tables remain unchanged.

## Interpretation limits

A nearest or overlapping gene is not a validated regulatory target. cCRE overlap
is reference context, not evidence that an element is active in the assayed
cell type. Category results depend on the exact GTF/cCRE releases, promoter
windows and precedence ordering, all of which must be retained in provenance.

Motif enrichment and enhancer-to-gene assignment remain out of scope.
