# Troubleshooting and safe recovery

## Five-minute triage

1. Read the external scheduler/nohup log and `logs/chip2tracks.console.log`.
2. Check the final rows of `00_metadata/workflow_events.tsv` and
   `stage_timing.tsv`.
3. Inspect `00_metadata/command_events.tsv` for a nonzero exit code.
4. Open the affected stage's `stage_status.tsv`, per-sample/cohort status and
   log directory.
5. Re-run `--preflight-only` after environment/reference changes.

Do not delete checkpoints to force recovery. `--from-stage NAME` validates and
adopts all earlier content, then reruns the named and later stages. Common safe
restart points are `preprocess`, `alignment`, `filtering`, `peakcalling`,
`consensus`, `normalized_tracks`, `qc`, `differential`, `annotation`
and `report`.

```bash
chip2tracks --config /project/config/config.conf --from-stage consensus
```

The deprecated name `reproducibility` is accepted as a temporary alias for
`consensus`.

A failed reuse can mean changed inputs/settings/code, missing output artifacts,
or an old checkpoint without producing-input provenance. The diagnostic names
the earliest incompatible stage. Restart at or before that stage; do not edit
checkpoint JSON. An ordinary legacy retry stops before replacing the historical
metadata. Automatic cleanup invalidates affected checkpoints deliberately.

If final preprocessing MultiQC failed after library work completed, correct the
executable/environment and retry with the same config and no `--from-stage`.
Validated library records are then reused and aggregation is retried. Check
`00_metadata/multiqc_preflight.log` for the selected executable and startup
errors. An explicit `--from-stage preprocess` intentionally forces library work.

## Frequent causes

- Metadata failure: compare the CSV header byte-for-byte with the 24-column
  template; inspect control linkage and cohort policy.
- Missing blacklist/reference: define the assembly-specific config key and
  verify chromosome naming across all references.
- Chipmentation trimming failure: ensure `TRIM_ADAPTERS=true` and select a
  protocol-compatible adapter preset.
- Normalization skipped: inspect consensus status, excluded peak samples and
  normalization-factor logs.
- Differential skipped: verify two or more biological replicates per condition
  and a successful cohort consensus/count matrix.
- Annotation failure: verify GTF/cCRE compression, contig compatibility and the
  matplotlib dependency.
- CPU overcommit: use `resource_budget.tsv` and reduce jobs or per-tool
  threads; R differential jobs also require memory headroom.

## Support bundle

Provide resolved config/manifest tables (after checking paths for sensitivity),
workflow/command/stage events, resource budget, software versions, relevant
stage status tables and logs, workflow version/commit, and the exact failing
command. Do not include FASTQs or BAMs unless explicitly requested.
