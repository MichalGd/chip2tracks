# Migration from 0.1 to 0.2

1. Remove `ASSAY_PROFILE` from new configs. Existing 0.1 configs may still
   contain it; the validator accepts and ignores it while assay values are read
   from samplesheet rows.
2. Remove the `blacklist` samplesheet column and its value from every row.
   Keep `BLACKLIST_<GENOME>` in config.
3. Set `COHORT_MODE=automatic` unless a shared cross-antibody peak universe is
   an explicit scientific decision. Use `global-compatible` only after review.
4. Decide whether control-free or shared-control behavior is justified. Both
   now default to `false`.
5. Review new job limits, annotation settings, technical-unit FastQC and logging
   settings. Old configs receive conservative backward-compatible values during
   validation.
6. Use `consensus` in `--from-stage` and `--stop-after`. The old
   `reproducibility` argument remains a temporary alias.
7. Run `--plan` into a new output directory, inspect `cohort_manifest.tsv`,
   `cohort_membership.tsv`, `cohort_policy.tsv` and `control_policy.tsv`,
   then run full preflight.

Do not reuse a 0.1 sanitized samplesheet with the old 25-column header. For a
content-validated recovery of existing outputs, retain the original output
directory and use the updated raw 24-column samplesheet/config with
`--from-stage`; earlier stages are adopted only if hashes and declared outputs
remain valid.
