# Shared-server installation and promotion

Install immutable workflow releases and point one stable launcher at
`/opt/bioinformatics/workflows/chip2tracks/current`. The supplied launcher
resolves that target once at process start and prepends versioned main and
sidecar environments; it does not activate Conda or modify user shells.

Recorded biolserv layout (verify the active target on the server before each
update; these paths do not establish which release is currently active):

```text
/opt/bioinformatics/workflows/chip2tracks/releases/<version>-<commit>/
/opt/bioinformatics/workflows/chip2tracks/current -> releases/<version>-<commit>
/opt/miniconda/envs/chip2tracks-0.1.0/
/usr/local/bin/chip2tracks
```

The environment name is independent of the workflow version. The supplied
launcher defaults to the paths above and supports `CHIP2TRACKS_ROOT` and
`CHIP2TRACKS_MAIN_ENV` overrides. The recovery correction adds no dependencies:
reuse the validated main environment and installed sidecar launchers rather
than recreating or upgrading them for this code-only update. Never replace an
existing release directory in place.

Installation checklist:

1. Copy a clean release into the immutable release directory and record its Git
   commit and archive checksum.
2. For a fresh installation, build the pinned main environment from
   `environment.lock.yml`. For a code-only update, verify and reuse the existing
   environment instead.
3. Build/verify enabled epic2, IDR, preseq and phantompeak sidecars from their
   manifests; keep legacy runtimes isolated.
4. Install `utilities/chip2tracks_shared_launcher.sh` as
   `/usr/local/bin/chip2tracks` and the SPP sidecar launcher as
   `/usr/local/bin/run_spp.R`.
5. Point `current` at the candidate only after validation.
6. Verify `chip2tracks --version`, `--plan`, `--preflight-only`, a small
   representative run, stop/restart from `consensus`, checksums, logs and final
   reports without manual environment activation.

For a code-only recovery update, first run `bash tests/run_tests.sh` from the
candidate on Linux, with the existing main environment's `bin` followed by
`/usr/local/bin:/usr/bin:/bin` on `PATH`. The regression suite uses temporary
synthetic fixtures, not existing analysis results. Resolve any Linux-specific
failures before promotion. Candidate `--plan`/`--preflight-only` checks must use
a separate test configuration/output directory: they write metadata and are not
read-only checks of an existing run.

Promote by atomically changing the `current` symlink. Existing jobs remain
pinned because the launcher resolved the release before execution. Roll back by
repointing `current` to the last validated immutable release; never mutate an
already promoted release in place.

### Existing runs and recovery-checkpoint compatibility

The recovery correction introduces schema-2 scientific checkpoint provenance.
Older schema-1 checkpoints cannot safely prove that their producing inputs and
settings match the new invocation. Keep each historical run's original release
and metadata available. Do not delete checkpoints, rewrite their signatures, or
copy new fingerprints into old markers to make a run appear compatible.

Use the new release for new runs. If an old run needs continuation, pin its
original compatible release explicitly, or assess a rebuild from the earliest
unproven stage separately; switching `current` is not a checkpoint migration.
See [Outputs and recovery](04_outputs_and_recovery.md#checkpoints).

Before promoting, review active jobs, record the previous absolute `current`
target, and retain that release for rollback. A versioned candidate should record
its source commit and archive checksum; publish the same validated source to
GitHub, not an independently edited server copy. Uncommitted local corrections
are not installed by `git pull` on the server.

For a dedicated 140-CPU/500-GB node, start conservatively and keep each
`jobs × threads` request near or below 120–128 CPUs. For example, eight
alignment jobs with 8 Bowtie2 + 4 samtools threads request 96 CPUs. Differential
jobs are memory-heavy; begin with 2–4 even when many CPUs are available.
