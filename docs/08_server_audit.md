# Server audit before dependency resolution

Run the audit on the Linux analysis server before creating the final lockfile.
It is read-only except for writing its report. It does not install, update, or
remove anything.

```bash
cd /path/to/chip2tracks
bash utilities/audit_server_environment.sh \
  --output "$PWD/chip2tracks_server_audit.txt" \
  --reference-root /opt/bioinformatics/references \
  --host-index /path/to/hg38/bowtie2/hg38 \
  --composite-index /path/to/hg38_dm6/bowtie2/hg38_dm6
```

If the reference prefixes are not yet known, omit those options on the first
run. Send back `chip2tracks_server_audit.txt`; it contains paths and versions,
but no FASTQ contents. Review it for site-specific secrets before sharing.

The audit determines which dependencies can be reused unchanged, which need an
isolated add-on environment, and whether the existing host+dm6 index is truly a
competitive composite index. Do not finalize `environment.lock.yml` before this
report is reviewed.

The repository isolates current Bioconda epic2, phantompeakqualtools/SPP,
preseq, and optional IDR in `environment.epic2.yml`, `environment.spp.yml`,
`environment.preseq.yml`, and `environment.idr.yml`. Their legacy Python, R, or
zlib constraints can conflict with the modern main environment. Prefer already
installed compatible commands; otherwise create only the required sidecars.

The epic2 sidecar pins `setuptools=80.9.0` because epic2 0.0.54 imports the
legacy `pkg_resources` module. After creating or updating that sidecar, run:

```bash
EPIC2_COMMAND=/opt/miniconda/envs/chip2tracks-epic2-0.0.54/bin/epic2 \
  bash utilities/smoke_test_epic2.sh
```

This creates temporary paired-end BAMs, exercises epic2's `--guess-bampe`
path, and deletes the fixtures. Normal workflow preflight also executes
`epic2 --help`; an import-broken launcher is a preflight failure.

`PHANTOMPEAK_COMMAND` is executed as a command rather than forced through the
main environment's `Rscript`. A site launcher may therefore expose a small
`run_spp.R` wrapper that executes both the script and R from the SPP sidecar.
Likewise, an existing shared preseq executable can be added to the launcher's
`PATH`; no per-user installation is needed. Point `EPIC2_COMMAND` and optional
`IDR_COMMAND` to absolute executables when they are not exposed by that launcher.

The repository supplies `utilities/chip2tracks_shared_launcher.sh` and
`utilities/run_phantompeak_sidecar.sh` for an all-user `/opt` installation.
Install copies under `/usr/local/bin` so their executable modes and targets are
controlled by the administrator; do not add per-user Conda activation to shell
startup files.
