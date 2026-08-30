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

A failed reuse means an earlier checkpoint or one of its hashed outputs is
missing or changed. Restart from that earlier stage; do not edit checkpoint
JSON. Automatic cleanup invalidates affected checkpoints deliberately.

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
