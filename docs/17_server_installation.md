# Shared-server installation and promotion

Install immutable workflow releases and point one stable launcher at
`/opt/bioinformatics/workflows/chip2tracks/current`. The supplied launcher
resolves that target once at process start and prepends versioned main and
sidecar environments; it does not activate Conda or modify user shells.

Recommended layout:

```text
/opt/bioinformatics/workflows/chip2tracks/releases/0.2.0/
/opt/bioinformatics/workflows/chip2tracks/current -> releases/0.2.0
/opt/bioinformatics/envs/chip2tracks/0.2.0/
/usr/local/bin/chip2tracks
```

Installation checklist:

1. Copy a clean release into the immutable release directory and record its Git
   commit and archive checksum.
2. Build the pinned main environment from `environment.lock.yml`.
3. Build/verify enabled epic2, IDR, preseq and phantompeak sidecars from their
   manifests; keep legacy runtimes isolated.
4. Install `utilities/chip2tracks_shared_launcher.sh` as
   `/usr/local/bin/chip2tracks` and the SPP sidecar launcher as
   `/usr/local/bin/run_spp.R`.
5. Point `current` at the candidate only after validation.
6. Verify `chip2tracks --version`, `--plan`, `--preflight-only`, a small
   representative run, stop/restart from `consensus`, checksums, logs and final
   reports without manual environment activation.

Promote by atomically changing the `current` symlink. Existing jobs remain
pinned because the launcher resolved the release before execution. Roll back by
repointing `current` to the last validated immutable release; never mutate an
already promoted release in place.

For a dedicated 140-CPU/500-GB node, start conservatively and keep each
`jobs × threads` request near or below 120–128 CPUs. For example, eight
alignment jobs with 8 Bowtie2 + 4 samtools threads request 96 CPUs. Differential
jobs are memory-heavy; begin with 2–4 even when many CPUs are available.
