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

The repository isolates current Bioconda epic2 and optional IDR in
`environment.epic2.yml` and `environment.idr.yml` because their Python 3.10 and
3.9 constraints conflict with the main Python 3.11 environment. Prefer already
installed compatible commands; otherwise create only the required sidecars and
point `EPIC2_COMMAND`/`IDR_COMMAND` to their absolute executables.
