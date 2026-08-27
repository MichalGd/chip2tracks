# Methods and assay-specific choices

## Shared processing

Technical FASTQ units are concatenated only within one biological library.
FastQC is run before and after Trim Galore. Trim Galore's auto mode examines up
to the first million reads for common adapter families; `illumina`, `nextera`,
and custom overrides are supported. Auto-detection does not eliminate the value
of recording the kit and oligos: uncommon, partial, or protocol-specific
adapters still require an explicit sequence.

Reads align with Bowtie2 to the selected host reference, or to a namespaced
host+dm6 competitive index when spike-in is enabled. BAM branches preserve
permissive/stringent MAPQ and duplicate-retained/removed variants. Secondary,
supplementary, QC-failed, noncanonical, blacklisted and optionally mitochondrial
alignments are removed according to configuration.

The configured `EFFECTIVE_GENOME_SIZE_<ASSEMBLY>` is supplied to MACS3. The
total length of the chromosome sizes file is not substituted for mappable
effective genome size.

## ChIP-seq profile

Paired-end MACS3 uses `BAMPE`, which makes MACS3 use observed fragments.
Single-end MACS3 uses its ChIP fragment model by default; fixed `--shift` and
`--extsize` are available only as an explicit validated override. Narrow marks
use MACS3 q-value calls. Broad marks use epic2 as the primary diffuse-domain
caller and MACS3 `--broad` as a companion.

The default analysis branch removes duplicates, while duplicate-retained q30
BAMs feed NRF/PBC, preseq and phantompeakqualtools. This avoids circularly
measuring library complexity after duplicate removal. MACS3 and epic2 are told
to retain the reads they receive because duplicate policy has already been
applied upstream; callers must not silently apply a second policy.

## ChIPmentation profile

ChIPmentation shares the ChIP statistical model; it is not treated as CUT&Tag.
The profile requires adapter trimming in v0.1 and defaults to read-based common
adapter-family detection; `nextera` or custom sequences remain explicit
overrides. After trimming, alignment, filtering, control modeling, peak-class
routing, and replicate logic follow ChIP principles rather than ATAC/Tn5 cut-site
shifting.

The profile defaults to adapter auto-detection because both standard Illumina
and Nextera/Tn5-derived sequences may occur locally. No blanket Tn5 insertion
shift is applied: peak calling models immunoprecipitated chromatin fragments,
not open-chromatin insertion sites. PE fragment distributions, overrepresented
sequences, residual adapters and short-insert rates receive additional review.

## Controls and tracks

Matched input, IgG, or mock controls improve background modeling but are
optional. When present, they must match genome, profile, layout, treatment,
cell type, batch and spike-in metadata. When absent, MACS3 local-lambda and
epic2 control-free modes are used and the limitation is explicitly recorded.
This workflow flexibility does not redefine experimental standards: ENCODE
expects a corresponding input control with matching run type/read length and
replicate structure for released ChIP-seq experiments.

Tracks include CPM coverage, spike-in-scaled coverage, DESeq2
consensus-normalized coverage, and MACS3 fold enrichment over its local/control
lambda where available. CPM and spike-in scaling answer different questions
and are never labeled interchangeably.

The priority CPM families use explicit BAM policies:

| Family | MAPQ | Duplicates | Location |
|---|---:|---|---|
| analysis | configured | configured per biological library | `04_tracks/cpm/` |
| permissive | 0 | retained | `04_tracks/cpm/permissive/` |
| intermediate | 0 | removed | `04_tracks/cpm/intermediate/` |
| stringent | 30 | removed | `04_tracks/cpm/stringent/` |

For each policy, the denominator is the number of retained PE fragments or SE
reads in that policy BAM, and the scale is `1,000,000 / denominator`. The
stringent consensus-derived robust-CPM family is written below
`04_tracks/deseq2_robust_cpm/stringent/`. BigWig and uncompressed bedGraph are
independently enabled in `config.conf`.

The workflow also records mapping composition for the BAM policy underlying
each CPM and normalized coverage family. It counts one representative
alignment per proper PE fragment, matching the fragment signal used for
coverage, or one alignment per SE read. Reported fields include total
observations, MAPQ 0 observations, MAPQ <30 observations, and alignments with a
Bowtie2 `XS` tag. Bowtie2 uses `XS` for the best-scoring alternative alignment
when more than one alignment was found, so this is the most direct retained
evidence of candidate multimapping available in the BAM. MAPQ 0 and MAPQ <30
are reported separately as useful ambiguity/filter-impact measures; they are
not labeled as exact multimapper counts because MAPQ is an aligner-specific
confidence score rather than a universal classification. For PE libraries,
the `XS` count refers to the representative mate used by the fragment-counting
policy, not both SAM records.

bedGraph and bigWig are derived from the same filtered observations and
therefore share one mapping-composition row per policy. Normalization changes
signal weights, not read membership: the analysis policy maps to the DESeq2
consensus family, while permissive, intermediate, and stringent policies map
to the corresponding robust-CPM families. Generation status for optional
normalized tracks remains authoritative in
`04_tracks/normalized_track_family_status.tsv`.

PE coverage uses observed paired fragments. For SE, `SE_SIGNAL_MODE=extend` and
`SE_FRAGMENT_LENGTH=auto` run MACS3 `predictd` once per sample. A failed model
uses `SE_FRAGMENT_LENGTH_FALLBACK` (200 bp by default) with a warning and a
metadata record; an explicit positive length or unextended `read` mode remains
available.

## Replicates and QC

The primary reproducible peak set follows the `ATACseq2tracks` support rule:
technical replicates are collapsed first, and an interval is retained when it
is covered by peaks from at least `CONSENSUS_MIN_BIOLOGICAL_SAMPLES` distinct
biological sample keys (two by default). This consensus drives quantification,
normalized tracks, FRiP, annotation, and differential analysis.

IDR is disabled by default and never replaces the consensus set. When
`RUN_IDR=true`, narrow MACS3 replicate pairs receive additional relaxed peak
calls and pairwise IDR using ranked narrowPeak signal. Broad peaks continue to
use biological-support consensus and signal correlation. Pooled and
pseudoreplicate IDR are planned for v0.2.

QC includes FastQC/MultiQC, mapping/flag statistics, duplicate metrics,
NRF/PBC1/PBC2, preseq complexity extrapolation, PE fragment length,
phantompeakqualtools NSC/RSC, FRiP, target/control fingerprints, replicate
Spearman correlation and PCA. Thresholds are warnings rather than automatic
biological truth; broad marks, low-input libraries, and ChIPmentation may need
mark- and protocol-specific interpretation.

## Primary references

- [ATACseq2tracks consensus definition](https://github.com/MichalGd/ATACseq2tracks#consensus-peaks)
- [Trim Galore documentation](https://github.com/FelixKrueger/TrimGalore)
- [Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [SAMtools documentation](https://www.htslib.org/doc/)
- [MACS3 callpeak documentation](https://macs3-project.github.io/MACS/docs/callpeak.html)
- [MACS3 predictd documentation](https://macs3-project.github.io/MACS/docs/predictd.html)
- [epic2 paper](https://doi.org/10.1093/bioinformatics/btz232)
- [IDR framework paper](https://doi.org/10.1214/11-AOAS466)
- [ENCODE ChIP-seq guidelines](https://www.encodeproject.org/chip-seq/transcription_factor/)
- [phantompeakqualtools](https://github.com/kundajelab/phantompeakqualtools)
- [deepTools effective genome sizes](https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html)
