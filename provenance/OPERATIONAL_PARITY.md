# Operational parity provenance

The 0.2.0 parity work was reviewed against:

- chip2tracks `02c586fa2c5c265e480a4eb626a5a3a47ff344da` (0.1.0);
- cutnrun2tracks `5c55e8fdbb2d386a90451f112a44baac95b60c7e` (0.3.1);
- `method_review/chip2tracks_vs_cutnrun2tracks_review_30aug2026.md` in the
  parent workspace.

Ported concepts include the CLI/stage vocabulary, multi-cohort metadata model,
persistent telemetry, bounded long-stage workers, technical-unit FastQC,
cohort-local replicate QC, and comprehensive peak-feature annotation. ChIP
caller, control, coverage, adapter, IDR and QC behavior was retained and adapted
rather than replaced with CUT&RUN/CUT&Tag-specific assumptions.

Deliberate non-ports:

- SEACR and ataqv are CUT/ATAC-specific and are not general ChIP requirements;
- full pooled/self-pseudoreplicate ENCODE IDR is not claimed or implemented
  without a separately specified and validated design;
- exact site GTF/cCRE/blacklist releases cannot be inferred from placeholder
  paths. Each run's reference manifest records resolved paths, sizes and hashes;
  deployment must additionally record source release, URL and transformation.
