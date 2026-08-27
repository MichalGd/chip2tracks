# References, blacklist, and alignment filtering

[Documentation index](README.md) | [Track normalization](12_tracks_and_normalization.md)

This page describes the reference contract and the exact observations retained
after alignment. These choices affect every downstream BAM, track, peak set, QC
metric, and differential result.

## Reference contract

Each genome used in a run needs the following `config.conf` values:

| Key | Role | Main consumers |
|---|---|---|
| `INDEX_<GENOME>` | Bowtie2 index prefix | host alignment |
| `FASTA_<GENOME>` | genome sequence | reference validation and provenance |
| `CHROM_SIZES_<GENOME>` | chromosome order and lengths | sorted bedGraph, bigWig, epic2, annotation |
| `CANONICAL_CONTIGS_<GENOME>` | allowed sequence names | all filtered host BAM branches |
| `GTF_<GENOME>` | gene features | generated TSS BED and nearest-gene annotation |
| `BLACKLIST_<GENOME>` | assembly-matched artifact regions | reference validation and run setup |
| `EFFECTIVE_GENOME_SIZE_<GENOME>` | mappable genome size | MACS3 and epic2 |
| `TSS_BED_<GENOME>` | optional prebuilt one-base TSS BED | descriptive TSS profiles |
| `CCRE_BED_<GENOME>` | optional cCRE intervals | cCRE overlap annotation |

The template points to the hg38 and mm39 resources already shared by
`cutnrun2tracks` and `ATACseq2tracks` on the target server. Reuse those files
when their assembly, chromosome naming, and provenance are appropriate. Do not
rebuild or download a second copy merely to give it a `chip2tracks` path.

`EFFECTIVE_GENOME_SIZE` is the estimated mappable genome size, not the sum of
the chromosome-size file. Changing the read length, mapping policy, or genome
assembly can change the defensible value. See the
[deepTools effective-genome-size guidance](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html)
and record the selected value with the run.

The preflight requires readable nonempty core references and checks the Bowtie2
index. It also writes `00_metadata/reference_manifest.tsv` with paths, sizes,
and checksums. This makes the exact installed files auditable; it does not by
itself certify that a file is the biologically correct release.

## Blacklist selection

Genome blacklists identify regions that frequently produce anomalously high or
unstructured signal across many sequencing assays. They are assembly-specific
artifact masks, not lists of biologically uninteresting loci. The rationale and
construction are described by
[Amemiya, Kundaje, and Boyle (2019)](https://doi.org/10.1038/s41598-019-45839-z),
and maintained lists are available from the
[Boyle Lab blacklist repository](https://github.com/Boyle-Lab/Blacklist).

The samplesheet `blacklist` column is required on every row and is the path
actually used to filter that biological library. Use the same validated,
assembly-matched blacklist for every target and control in one samplesheet.
The corresponding `BLACKLIST_<GENOME>` config value records the run-level
reference. A custom per-row path is technically possible, but mixing blacklist
definitions within a comparison changes the observable universe and is not a
valid default analysis design.

For an unlisted custom genome:

1. use matching names in the Bowtie2 index, FASTA, chromosome sizes, canonical
   contigs, GTF, blacklist, TSS, and cCRE files;
2. coordinate-sort BED references in the order of `CHROM_SIZES_<GENOME>`;
3. provide a valid, nonempty, curated blacklist BED; if none exists, treat the
   custom genome as not production-ready rather than inventing artifact regions;
   and
4. validate biological provenance before production use.

Do not substitute an hg19, mm10, or lifted blacklist without documenting and
validating the coordinate conversion.

## Exact host filtering order

Picard first marks duplicates without removing them. From that common marked
BAM, `scripts/mark_filter_batch.sh` builds four branches:

| Branch | Default MAPQ | Duplicate records | Main role |
|---|---:|---|---|
| `q0_dup-retained` | 0 | retained | permissive coverage |
| `q0_dup-removed` | 0 | removed | intermediate coverage |
| `q30_dup-retained` | 30 | retained | duplicate-sensitive analysis option and ChIP QC |
| `q30_dup-removed` | 30 | removed | stringent analysis and coverage |

For each branch the workflow:

1. retains mapped primary alignments and removes secondary, supplementary, and
   QC-failed records;
2. requires proper pairing for PE libraries;
3. applies the branch MAPQ threshold;
4. retains only configured canonical contigs, excluding `chrM`, `MT`, or `M`
   when `REMOVE_MITO=true`;
5. removes alignments overlapping the samplesheet blacklist with
   `bedtools intersect -v`; and
6. repairs PE mate information, removes orphaned pairs, coordinate-sorts,
   validates, and indexes the result.

For PE data, a fragment loses proper-pair eligibility if blacklist filtering
removes one mate, so a one-mate blacklist overlap does not survive as an orphan
fragment. For SE data, the overlapping read itself is removed.

The analysis BAM is a symlink to the q30 duplicate-retained or q30
duplicate-removed branch selected by the samplesheet
`analysis_duplicate_policy`. The public examples use duplicate removal. Peak
callers receive this already filtered analysis BAM and are told to keep the
observations supplied to them, avoiding a second duplicate-removal policy.
DiffBind also disables its additional blacklist step because both the BAMs and
consensus regions have already passed the workflow blacklist; when every
eligible target has a matched control, DiffBind greylisting remains enabled.

## Auditing filter impact

The following files make filtering visible rather than implicit:

- `03_alignment/metrics/<sample_key>.duplicate_metrics.txt`: Picard duplicate
  metrics from the marked BAM;
- `03_alignment/metrics/<sample_key>.filter_counts.tsv`: signal-unit counts in
  all four branches;
- `04_tracks/cpm/mapping_composition.tsv`: MAPQ 0, MAPQ less than 30, and
  Bowtie2 `XS` alternative-alignment evidence for every coverage policy; and
- `00_metadata/reference_manifest.tsv`: blacklist and reference identities.

`XS` evidence and MAPQ categories are useful ambiguity measures, but neither is
an aligner-independent exact definition of a multimapper. See
[Tracks and normalization](12_tracks_and_normalization.md) for how these BAM
branches map to bedGraph and bigWig families.

## Spike-in filtering

With `SPIKEIN_MODE=dm6`, spike reads come from competitive host-plus-dm6
alignment and use their own `SPIKEIN_MIN_MAPQ`, `SPIKEIN_DUPLICATE_POLICY`,
allowed-contig list, and optional `SPIKEIN_BLACKLIST`. The host analysis still
uses the host blacklist above. An absent fly blacklist does not mean that the
host blacklist should be applied to dm6 coordinates.

## Failure behavior

Preflight stops for a missing core reference, an unreadable samplesheet
blacklist, a missing required cCRE file, or an invalid enabled spike reference.
Filtering stops if a branch BAM is invalid or contains zero signal units unless
`ALLOW_EMPTY_FILTERED_BAM=true`. The latter is a diagnostic escape hatch, not a
recommended production setting.
