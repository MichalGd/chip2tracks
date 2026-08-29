# Genomic coverage tracks and normalization

[Documentation index](README.md) | [Quality control](11_quality_control.md)

`chip2tracks` produces several bedGraph and bigWig families because no single
normalization answers every ChIP-seq or ChIPmentation question. The file suffix
names the scaling method; the directory and metadata name the filtering policy.
Do not relabel one family as another.

## Signal unit and file format

PE coverage represents observed properly paired fragments counted once. SE
coverage represents reads; by default they are extended to the MACS3 `predictd`
fragment-length estimate, with a recorded configurable fallback. Set
`SE_SIGNAL_MODE=read` for unextended SE reads or provide an explicit
`SE_FRAGMENT_LENGTH` when required.

All coverage is generated at `TRACK_BIN_SIZE` resolution with deepTools
`bamCoverage` and then chromosome-sorted. BigWig and bedGraph made from the same
stem contain the same binned values:

- bigWig is compressed, indexed, and preferred for UCSC, IGV, and routine
  visualization;
- bedGraph is uncompressed text and is useful for inspection or downstream
  tools that require text intervals, but consumes substantially more storage.

`GENERATE_COVERAGE_BIGWIGS` and `GENERATE_COVERAGE_BEDGRAPHS` independently
control retention. Both are `true` in the public template.

## Track-family matrix

| Family | Default retained observations | Scale applied to coverage | Intended comparison |
|---|---|---|---|
| analysis CPM | MAPQ 30; sample `analysis_duplicate_policy` | `1,000,000 / retained signal units` | routine per-library visualization |
| permissive CPM | MAPQ 0; duplicates retained | `1,000,000 / retained signal units` | mapping-plus-duplicate sensitivity |
| intermediate CPM | MAPQ 0; duplicates removed | `1,000,000 / retained signal units` | mapping ambiguity sensitivity after deduplication |
| stringent CPM | MAPQ 30; duplicates removed | `1,000,000 / retained signal units` | conservative filtering comparison |
| DESeq2 consensus | analysis-policy BAM | reciprocal DESeq2 consensus-count size factor | relative cohort signal on the analysis policy |
| DESeq2 robust CPM | selected permissive/intermediate/stringent BAM | `1,000,000 / robust effective library size` | cohort-normalized filtering sensitivity |
| spike-in host | analysis-policy host BAM | `SPIKEIN_SCALE_TARGET * declared ratio / retained spike observations` | external-reference calibration |
| spike-control CPM | filtered spike BAM | `1,000,000 / retained spike observations` | fly coverage-shape and calibration QC |
| MACS3 fold enrichment | MACS3 treatment pileup and lambda | treatment pileup divided by control/local lambda | enrichment over modeled background |

The defaults above are recorded rather than inferred. If MAPQ, duplicate, bin,
or SE settings are changed, the per-track metadata under `04_tracks/` and the
resolved configuration are authoritative.

## CPM families

The priority CPM outputs are:

```text
04_tracks/cpm/<sample_key>.CPM.{bw,bedGraph}
04_tracks/cpm/permissive/<sample_key>.CPM.{bw,bedGraph}
04_tracks/cpm/intermediate/<sample_key>.CPM.{bw,bedGraph}
04_tracks/cpm/stringent/<sample_key>.CPM.{bw,bedGraph}
```

For each family, `C` is unscaled bin coverage and `L` is the number of retained
PE fragments or SE reads in that exact BAM policy. The output is `C * 10^6 / L`.
Consequently every CPM family has its own denominator; the permissive and
stringent tracks are not made by applying an extra filter to already normalized
values.

Use CPM for within-sample visualization and approximate comparisons when global
occupancy and composition are not expected to shift substantially. CPM is not
an external calibration and does not correct antibody efficiency, IP yield, or
global biological changes.

## Consensus-derived DESeq2 families

For each cohort, the workflow counts raw integer fragments/reads over the fixed
biological-support consensus. DESeq2 `poscounts` size factors are then estimated
from that matrix.

`04_tracks/deseq2_consensus/<cohort>/` applies `1 / size_factor` to the
analysis-policy BAM. These values are relative to the fitted cohort and have no
absolute per-million unit.

The robust-CPM families are under:

```text
04_tracks/deseq2_robust_cpm/permissive/<cohort>/
04_tracks/deseq2_robust_cpm/intermediate/<cohort>/
04_tracks/deseq2_robust_cpm/stringent/<cohort>/
```

For one policy, let `s_j` be sample `j`'s DESeq2 size factor and let `G` be the
geometric mean of cohort consensus-count column sums. The effective library size
is `s_j * G`, and the robust-CPM scale is `1,000,000 / (s_j * G)`. Each policy
gets its own raw count matrix and factors over the same consensus intervals.

These tracks are useful for comparing relative patterns within that cohort and
for checking sensitivity to MAPQ/duplicate policy. They are not absolute
cross-study calibration. A cohort with fewer than two targets, no usable
consensus, or a zero-count sample cannot receive valid factors; the status and
reason are recorded in `04_tracks/normalized_track_family_status.tsv`. No sample
is silently dropped to force normalization.

## Drosophila spike-in tracks

With `SPIKEIN_MODE=dm6`, competitive alignment assigns reads to namespaced host
or fly sequences. For a sample with retained fly count `D` and declared
`spikein_to_host_ratio` `R`, the host scale is:

```text
SPIKEIN_SCALE_TARGET * R / D
```

Host tracks are written below `04_tracks/spikein/<cohort>/`; fly CPM controls
are written below `04_tracks/spikein_control/`. This reproduces the inherited
`cutnrun2tracks` ratio-aware scaling policy. Interpret the tracks together with
`06_qc/spikein/spikein_scaling.tsv`; low fly counts, extreme fractions, lot
changes, or inconsistent cell addition can invalidate calibration.

Use spike-in-scaled host tracks when the experimental protocol supports global
occupancy comparisons. They are not interchangeable with CPM or DESeq2 tracks
and do not remove unrelated batch effects.

## MACS3 fold-enrichment tracks

When `MACS3_GENERATE_SIGNAL_TRACKS=true`, MACS3 creates treatment pileup and
control-lambda bedGraphs, and `macs3 bdgcmp -m FE` writes
`04_tracks/control_normalized/<sample_key>.FE.{bw,bedGraph}`. With a matched
control, lambda includes that control and MACS3 local background; without one,
it is a local-lambda-only background. The status table labels which case was
used.

Fold enrichment is useful for viewing signal relative to modeled background.
It is not a library-size normalization and should not be used as the count input
for differential occupancy. RPGC and log-likelihood/log-fold track families are
not implemented in v0.1.

## Mapping composition and multimappers

`04_tracks/cpm/mapping_composition.tsv` reports the observations underlying
each CPM policy and its corresponding normalized family:

- total signal observations;
- MAPQ 0 count and percent;
- MAPQ less than 30 count and percent; and
- count and percent with a Bowtie2 `XS` alternative-alignment score.

For PE libraries, one representative mate is counted per proper fragment; for
SE, one retained alignment is counted. MAPQ is an aligner confidence score, and
`XS` indicates a reported alternative score, so neither field alone is labeled
an exact universal multimapper count. A bedGraph and bigWig from the same policy
share one composition row because conversion changes format, not observations.

The same table is copied to
`10_reports/coverage_mapping_composition.tsv` and included in the lightweight
and MultiQC final reports.

Every retained `.bw` from the families above receives a UCSC one-line custom
track descriptor under `09_browser/ucsc/`. Family-specific files and the
grouped `all_bigwig_tracks.txt` are generated from the actual output tree, so
disabled or skipped track families do not create stale descriptors.

## Choosing a track

| Goal | Preferred starting family |
|---|---|
| inspect an individual library | analysis CPM |
| determine whether ambiguous alignments or duplicates drive a feature | compare permissive, intermediate, and stringent CPM |
| compare relative signal among cohort samples | DESeq2 consensus or the matching robust-CPM family |
| assess a plausible global occupancy shift with validated fly-cell addition | spike-in host tracks |
| inspect enrichment over input/IgG/local background | MACS3 fold enrichment |
| statistical differential occupancy | raw consensus counts, never a bedGraph or bigWig |

For exact deepTools coverage semantics, see the
[bamCoverage documentation](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html).
