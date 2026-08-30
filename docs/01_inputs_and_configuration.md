# Inputs and configuration

The samplesheet header must exactly match `config/samplesheet_template.csv`.
One row is a sequencing unit; identical `sample_id + replicate` rows with
different `tech_replicate` values are merged before trimming. Biological
replicates remain independent.

Each row declares `assay_profile=chipseq` or `chipmentation`. Both profiles may
coexist in a run, but are always assigned to separate compatible cohorts. Targets may
reference a control `sample_id` through `control_id`; an empty value is valid
when `ALLOW_CONTROL_FREE_PEAKCALL=true`. The validator
first requires the same replicate. A condition=`shared` control is eligible
only when `ALLOW_SHARED_CONTROLS=true` and all other biological context matches.
There is no replicate-1 fallback.

One samplesheet may define multiple isolated peak universes. With the recommended
`COHORT_MODE=automatic`, factor and antibody are part of the cohort identity.
With explicit `COHORT_MODE=global-compatible`, they are ignored so a researcher
can define a shared consensus across different antibodies; hard technical and
analysis compatibility boundaries are still enforced and recorded.
See [Replicates, controls, and experimental design](13_replicates_and_design.md)
for the merge, shared-control, consensus, block, and contrast rules.

Configuration is plain `KEY=VALUE`, not a sourced user shell script. Unknown
keys, duplicate keys, shell expansions, invalid booleans, and incomplete spike
metadata fail before tools run. The resolved configuration is written to
`00_metadata/resolved_config.conf` and `resolved_config.tsv`.

Run `--plan` for metadata/cohort validation without checking FASTQ/reference
existence or installed tools. Run `--preflight-only` for the full environment
and reference audit.

`ADAPTER_PRESET=auto` delegates common adapter-family detection to Trim Galore
and records evidence per sample. Explicit `illumina`, `nextera`, and custom
sequences are supported. Auto mode cannot reliably identify every custom kit,
so unresolved detection requires read-level review and an override.

When `SPIKEIN_MODE=dm6`, every row records a positive spike-to-host ratio,
addition stage, and lot. Reference configuration points to a competitive
host+dm6 index; post hoc independent alignment is not treated as equivalent.
Spike-in is optional and `SPIKEIN_MODE=none` is the public default.

Blacklist paths are resolved centrally from `BLACKLIST_<GENOME>` config keys
and written into the generated sample manifest; they are not samplesheet columns.
Reference compatibility and the exact blacklist/canonical/MAPQ/duplicate
filtering order are documented in
[References, blacklist, and filtering](10_references_blacklist_and_filtering.md).

The optional metagene stage uses a separate tab-delimited gene-set manifest so
annotation releases and HPA-derived subsets can be replaced without modifying
workflow code. Filtering thresholds belong to the offline reference-build
command and are recorded in its reference manifest; plotting-window and
rendering settings belong to `config.conf`.
