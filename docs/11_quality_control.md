# Quality control and interpretation

[Documentation index](README.md) | [Replicates and design](13_replicates_and_design.md)

`chip2tracks` reports QC at the biological-library level after technical FASTQ
units have been merged. It deliberately does not convert one universal set of
thresholds into automatic sample exclusion: acceptable enrichment, duplicate
rate, fragment distribution, and cross-correlation depend on the target,
antibody, cell number, assay profile, and whether enrichment is narrow or
broad. The final decision must consider the metrics together.

## QC inventory

| Question | Metric or artifact | Input used | Main output |
|---|---|---|---|
| Are raw and trimmed reads technically sound? | FastQC and Trim Galore reports, adapter evidence | merged biological-library FASTQs | `01_fastq_qc/` |
| How many alignments remain? | Bowtie2 summary, `samtools flagstat`, `samtools stats` | raw alignment and analysis BAM | `03_alignment/metrics/`, `06_qc/alignment_and_complexity/` |
| What did filtering remove? | q0/q30 and duplicate-retained/removed counts | four filtered BAM branches | `03_alignment/metrics/*.filter_counts.tsv` |
| Are ambiguous alignments contributing to tracks? | MAPQ 0, MAPQ less than 30, Bowtie2 `XS` evidence | each coverage-policy BAM | `04_tracks/cpm/mapping_composition.tsv` |
| Is the library complex? | Picard duplicate metrics, NRF, PBC1, PBC2, preseq | q30 duplicate-retained BAM | `03_alignment/metrics/`, `06_qc/alignment_and_complexity/` |
| Is fragment structure plausible? | PE insert-length histogram | analysis BAM | `06_qc/fragment_length_and_periodicity/` |
| Is ChIP enrichment detectable? | NSC/RSC cross-correlation | q30 duplicate-retained tagAlign | `06_qc/fragment_length_and_periodicity/*.phantompeak.tsv` |
| Is signal concentrated in reproducible peaks? | FRiP against cohort consensus | target analysis BAM | `06_qc/frip_and_peak_reproducibility/*.frip.tsv` |
| Does target separate from background? | target/control fingerprint | matched analysis BAMs | `06_qc/controls/` |
| Do biological replicates agree? | genome-bin Spearman heatmap and PCA | all target analysis BAMs | `06_qc/correlation_pca_fingerprint/` |
| Is there aggregate TSS-proximal signal? | descriptive TSS profile | analysis CPM bigWig | `06_qc/tss_signal_profile/` |
| Is fly calibration credible? | host/fly observations, fractions, scale factors, warnings | competitive-alignment branches | `06_qc/spikein/` |

Raw and trimmed FastQC run after technical units have been concatenated into one
biological library. This catches run-level adapter and quality problems but does
not preserve separate post-merge FastQC for each lane. Investigate the original
sequencing-run QC when a merged report suggests lane-specific failure.

## Library complexity

Complexity metrics use the duplicate-retained q30 BAM. Calculating them after
duplicate removal would erase the multiplicity that they are meant to measure.

- **NRF** is distinct genomic signal units divided by total signal units.
- **PBC1** is locations observed once divided by distinct observed locations.
- **PBC2** is locations observed once divided by locations observed twice.
- **preseq** extrapolates expected distinct observations with additional
  sequencing.

For PE libraries, one properly paired fragment is one signal unit; for SE, one
retained alignment is one signal unit. Picard metrics and the workflow's
coordinate-multiplicity metrics answer related but not identical questions, so
small numerical differences are expected. Interpret unusually low complexity
with library depth, cell input, target breadth, and duplicate sensitivity.

## Cross-correlation and fragments

When `RUN_CROSS_CORRELATION=true`, phantompeakqualtools runs on duplicate-retained
q30 alignments and produces strand cross-correlation results and a PDF. NSC and
RSC are enrichment diagnostics, not universal pass/fail rules. Broad histone
domains and low-input ChIPmentation can behave differently from sharp
transcription-factor ChIP. See the
[phantompeakqualtools documentation](https://github.com/kundajelab/phantompeakqualtools)
and the applicable
[ENCODE ChIP-seq guidance](https://www.encodeproject.org/chip-seq/transcription_factor/)
when defining project-specific acceptance criteria.

PE fragment histograms come from observed insert lengths in the analysis BAM.
They are useful for detecting short inserts, unexpected long-fragment tails,
and ChIPmentation library anomalies. The workflow does not apply an ATAC-style
Tn5 cut-site shift or require nucleosomal periodicity.

## FRiP and peak reproducibility

FRiP is calculated for each target against the run's biological-support
consensus, not against that sample's own peaks. This makes replicates comparable
within the same universe, but its value depends on the consensus construction,
peak class, caller, and filtering policy. Broad marks can have a lower or less
directly comparable FRiP than focused marks.

Review together:

- per-sample peak-call status and peak count;
- the number of successful samples contributing to consensus;
- excluded peak contributions in `excluded_peak_samples.tsv`;
- consensus region count and support;
- FRiP; and
- replicate correlation/PCA.

A zero-peak result may be recorded as `EMPTY` and allowed to continue under the
configured peak-call failure policy. It is not automatically relabeled as a
failed library, but it must be explained before biological interpretation.

## Replicate correlation, PCA, and controls

`multiBamSummary bins` summarizes target analysis BAMs in genome bins;
`plotCorrelation` reports Spearman correlation and `plotPCA` reports the main
between-library axes. These plots are descriptive. With multiple conditions,
good data may separate by condition; disagreement should therefore be judged
within condition and alongside batch/donor metadata.

When a matched control exists, `plotFingerprint` compares target and control.
Control-free samples have no fingerprint plot. Controls are not biological
replicates and are not included as target columns in the primary differential
model.

## TSS profile is not an ATAC TSS-enrichment score

The TSS module plots average CPM signal around strand-aware TSS positions from
`TSS_BED_<GENOME>` or a BED generated from GTF `gene` records. It is explicitly
descriptive and does not calculate the standard ATAC-seq TSS-enrichment ratio.
For many ChIP targets, absence of a TSS-centered pattern is biologically
expected and must not be treated as automatic failure.

## Spike-in QC

Spike-in is disabled by default. When enabled, the workflow records retained
host and fly observations, fly fraction, declared spike-to-host ratio, host
scale, fly CPM scale, reference median scale, status, failure reasons, and
warnings. In v0.1 the field named `cohort_median_host_scale` is calculated over
all valid libraries in the run; use one samplesheet per compatible spike batch.
The important configurable checks are:

- `SPIKEIN_MIN_OBSERVATIONS_FAIL`;
- `SPIKEIN_MIN_OBSERVATIONS_WARN`;
- `SPIKEIN_WARN_LOW_FRACTION` and `SPIKEIN_WARN_HIGH_FRACTION`; and
- a greater-than-tenfold scale-factor deviation from the cohort median.

`ALLOW_FAILED_SPIKEIN=false` makes a failed calibration fatal. Passing these
checks confirms computational count sufficiency; it cannot prove that fly cells
were added consistently or that the spike behaves independently of the biology.

## Final report and decision record

The unified report is `10_reports/chip2tracks_multiqc_report.html`. The lighter
`10_reports/pipeline_report.html` and the machine-readable tables remain useful
when MultiQC is unavailable. Important stable tables include:

- `10_reports/run_summary.tsv`;
- `10_reports/warning_summary.tsv`;
- `10_reports/coverage_mapping_composition.tsv`; and
- `10_reports/differential_occupancy_summary.tsv`.

Before accepting a cohort, record the project-specific reason for accepting or
excluding each biological library. The workflow reports evidence but does not
silently remove a sample from differential analysis based on a single QC
threshold.
