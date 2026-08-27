# Peak and differential-result annotation

[Documentation index](README.md) | [Differential occupancy](14_differential_occupancy.md)

The annotation stage adds two lightweight reference-overlap layers to each
cohort consensus: the nearest GTF `gene` interval and overlapping cCRE reference
records. Annotation is descriptive, runs after statistical testing, and never
changes peak calls, counts, P-values, adjusted P-values, or fold changes.

## Configuration and required references

```text
RUN_SIMPLE_PEAK_ANNOTATION=true
RUN_CCRE_ANNOTATION=true
GTF_HG38=/path/to/annotation.gtf
CCRE_BED_HG38=/path/to/hg38.ccre.bed.gz
```

Equivalent assembly-specific keys are required for other genomes. With
`RUN_CCRE_ANNOTATION=true`, preflight requires a readable nonempty cCRE BED.
Set it to `false` for nearest-gene-only annotation. Setting
`RUN_SIMPLE_PEAK_ANNOTATION=false` disables nearest-gene annotation but does not
implicitly disable the separately configured cCRE overlap.

`RUN_MOTIF_ENRICHMENT` must remain `false`: motif enrichment is not implemented
in v0.1, and preflight fails if it is enabled.

## Nearest-gene implementation

For each cohort, the workflow:

1. extracts GTF records whose feature type is exactly `gene`;
2. converts them to BED with chromosome, gene span, `gene_name`, `gene_id`, and
   strand;
3. sorts both genes and consensus intervals using the configured chromosome-size
   order; and
4. runs `bedtools closest -d` against whole gene spans.

The output is:

```text
07_annotation/<cohort>/consensus/<cohort>.nearest_gene.tsv
```

The first five columns are the consensus interval and support count, followed
by the six GTF-derived gene BED columns and the unsigned distance reported by
`bedtools closest`. Distance is zero for an overlap. This is nearest **gene
span**, not nearest TSS, and it is not strand-oriented regulatory distance.

The result depends on the exact GTF release and requires `gene` features with
usable `gene_id`; `gene_name` is retained when present. A nearest gene is a
proximity annotation, not evidence that the ChIP-enriched region regulates or
is regulated by that gene.

## cCRE overlap implementation

When enabled, the workflow runs `bedtools intersect -wao` between the cohort
consensus and `CCRE_BED_<GENOME>`:

```text
07_annotation/<cohort>/consensus/<cohort>.ccre_reference_overlaps.tsv
```

The table preserves every consensus/cCRE overlap, all supplied cCRE BED fields,
and the overlap length. Non-overlapping consensus regions remain visible with
the normal `-wao` placeholder fields.

For annotated differential tables, v0.1 collects unique values from cCRE BED
column 4 into `ccre_reference_overlaps`. Therefore column 4 should contain a
stable element identifier or the value the project wants propagated. The
workflow does not parse the full ENCODE cCRE class vocabulary or choose a
primary cCRE when several overlap.

cCRE reference overlap is cell-type-agnostic reference context. It is not proof
that an element is active in the assayed cells, and it does not establish an
enhancer-to-gene relationship. Record the source release and assembly in the
reference manifest. For current ENCODE registry concepts and downloads, see
[SCREEN](https://screen.wenglab.org/) and the
[ENCODE candidate regulatory element resources](https://www.encodeproject.org/).

## Propagation to differential tables

After consensus annotation, `scripts/annotate_differential_results.py` finds
gzipped TSV result tables under `08_differential/<cohort>/` that contain a
`region_id` column and writes an adjacent `*.annotated.tsv.gz` copy. It adds:

| Column | Meaning |
|---|---|
| `nearest_gene_name` | nearest whole-gene interval's name, when present |
| `nearest_gene_id` | nearest whole-gene interval's GTF identifier |
| `distance_to_gene` | unsigned `bedtools closest -d` distance; zero on overlap |
| `ccre_reference_overlaps` | unique overlapping cCRE BED column-4 values |

The original table remains the statistical source of truth. In v0.1, direct
DESeq2Enrichment and target/control interaction outputs are gzipped TSVs and
receive annotated copies. DiffBind writes uncompressed TSVs and does not receive
automatic annotation.

## What is not implemented

Unlike the richer ATACseq2tracks annotation module, `chip2tracks` v0.1 does not
currently provide:

- promoter, exon, intron, or distal-intergenic classification;
- strand-aware nearest-TSS distance;
- a deterministic primary cCRE and parsed cCRE class fields;
- annotation summaries per differential contrast;
- automatic annotation of DiffBind TSVs;
- motif enrichment; or
- enhancer-to-gene assignment.

These are explicit limitations, not implied outputs. The current lightweight
annotation is sufficient for interval lookup and preliminary review, while a
future release can reuse the tested ATACseq2tracks annotation design after
ChIP-specific validation.
