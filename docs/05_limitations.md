# Limitations and pilot decisions

- Adapter auto-detection covers common families but cannot reconstruct unknown
  custom oligos or replace laboratory records. Residual-adapter FastQC warnings
  must be reviewed, especially for short-insert ChIPmentation libraries.
- Controls are optional, but control-free peaks and local-lambda fold enrichment
  are not equivalent to matched-input normalization.
- epic2 is purpose-built for diffuse domains but is not universally superior.
  Pilot each broad mark against the retained MACS3 broad companion.
- IDR is optional and disabled by default. When enabled, v0.1 pairwise IDR does
  not implement pooled or pseudoreplicate rescue/self-consistency ratios. Those
  are required before claiming full ENCODE pipeline equivalence.
- Broad domains use overlap consensus and signal correlation, not narrow-peak
  IDR. Thresholds need mark-specific validation.
- UMI-aware libraries are rejected in v0.1.
- Spike-in calibration assumes Drosophila cells were added consistently before
  the experiment and that competitive host+dm6 alignment is valid. Low counts,
  extreme fractions, lot changes, or biology that affects both species can
  invalidate scale factors.
- `EFFECTIVE_GENOME_SIZE_*` must be curated for each reference/read setup; the
  workflow does not infer it from chromosome-length totals.
- Mixed assay profiles are not allowed in one run. Mixed layouts/genomes remain
  opt-in because they complicate cohort comparisons.
- Differential models and normalized consensus tracks skip rather than silently
  discard zero-count samples. Sparse or discordant cohorts require review.
- The v0.1 annotation layer reports nearest whole-gene spans and raw cCRE
  overlaps. It does not yet reproduce ATACseq2tracks promoter/exon/intron,
  nearest-TSS, primary-cCRE, contrast-summary, or DiffBind annotation features;
  motif enrichment is not implemented.
- The inherited Bash engine is retained for v0.1. A workflow-engine migration is
  deferred until real ChIP-seq and ChIPmentation pilots pass.
