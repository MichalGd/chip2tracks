# Limitations and pilot decisions

- Adapter auto-detection covers common families but cannot reconstruct unknown
  custom oligos or replace laboratory records. Residual-adapter FastQC warnings
  must be reviewed, especially for short-insert ChIPmentation libraries.
- Controls are required by default, but the explicit control-free opt-in remains;
  control-free peaks and local-lambda fold enrichment
  are not equivalent to matched-input normalization.
- epic2 is purpose-built for diffuse domains but is not universally superior.
  Pilot each broad mark against the retained MACS3 broad companion.
- IDR is optional and disabled by default. When enabled, pairwise IDR does
  not implement pooled or pseudoreplicate rescue/self-consistency ratios. Those
  are required before claiming full ENCODE pipeline equivalence.
- Broad domains use overlap consensus and signal correlation, not narrow-peak
  IDR. Thresholds need mark-specific validation.
- UMI-aware libraries are rejected in this release.
- Spike-in calibration assumes Drosophila cells were added consistently before
  the experiment and that competitive host+dm6 alignment is valid. Low counts,
  extreme fractions, lot changes, or biology that affects both species can
  invalidate scale factors.
- `EFFECTIVE_GENOME_SIZE_*` must be curated for each reference/read setup; the
  workflow does not infer it from chromosome-length totals.
- Mixed assay profiles can share a run but remain separate cohorts. Mixed
  layouts/genomes remain opt-in and are also separated by hard compatibility.
- `global-compatible` cohorting is researcher-controlled and can combine
  different factors/antibodies into one peak universe. It does not establish
  biological equivalence and can yield misleading models when target identity
  is confounded with condition.
- Differential models and normalized consensus tracks skip rather than silently
  discard zero-count samples. Sparse or discordant cohorts require review.
- Feature annotation now covers per-sample/caller and consensus peak sets, but
  gene/cCRE assignments remain reference- and rule-dependent; motif enrichment
  and enhancer-to-gene inference are not implemented.
- The inherited Bash engine is retained. A workflow-engine migration is
  deferred until real ChIP-seq and ChIPmentation pilots pass.
