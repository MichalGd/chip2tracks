# Configuration

Copy `config.conf.template` outside the repository and replace every path. It
is a non-executable `KEY=VALUE` file: unknown keys, shell expansion, unsafe
values, and missing required options are rejected. Set `SAMPLESHEET` in this
file; the normal launch command needs only `--config`.

Assay profile is a samplesheet field, not a run-wide config value. A run may
contain `chipseq` and `chipmentation` libraries; they always remain in
separate compatible cohorts. If any chipmentation library is present,
`TRIM_ADAPTERS=true` is required. The adapter preset remains run-wide, so
mixed-profile runs must use a preset that is suitable for every library.

## Cohort policy

`COHORT_MODE=automatic` is the recommended default. It creates independent
peak universes for each compatible factor and antibody. Multiple factors,
antibodies, conditions and replicates can therefore be processed in one run
without mixing their downstream count matrices.

`COHORT_MODE=global-compatible` is an explicit researcher override. It ignores
factor and antibody when grouping targets, allowing a shared consensus peak
universe across different antibodies. It still separates genomes, assay
profiles, PE/SE layouts, target classes, duplicate policies, primary callers,
primary peak classes and spike-in policies. Generated `cohort_policy.tsv` and
`cohort_membership.tsv` record this choice, and mixed factor/antibody fields in
`cohort_manifest.tsv` are reported as `MULTIPLE`. This option defines a
technical region universe; it does not assert biological equivalence.

Matched controls are the template default. Set
`ALLOW_CONTROL_FREE_PEAKCALL=true` only for a justified control-free design.
Shared controls likewise require `ALLOW_SHARED_CONTROLS=true` and still must
pass metadata compatibility checks.

## Samplesheet and references

The exact CSV header has 24 columns. Use `config/examples/` as schema examples.
PE requires both FASTQs; SE requires an empty `fastq_2`. Each row is a
technical unit. Rows sharing `sample_id + replicate` are merged only when all
biological metadata agree.

Blacklist paths do not belong in the samplesheet. Every assembly needs
`INDEX_<GENOME>`, `FASTA_<GENOME>`, `CHROM_SIZES_<GENOME>`,
`CANONICAL_CONTIGS_<GENOME>`, `GTF_<GENOME>`,
`BLACKLIST_<GENOME>` and `EFFECTIVE_GENOME_SIZE_<GENOME>`. The workflow
resolves the configured blacklist into the generated sample manifest.

`ADAPTER_PRESET=auto` records Trim Galore evidence. Use `illumina`, `nextera`
or `custom` when the protocol is known. `RUN_FASTQC_PER_TECHNICAL_UNIT=true`
adds pre-merge unit QC while preserving merged-raw and trimmed FastQC.

Spike-in is off by default. Enabling it requires the competitive reference
settings and all three spike-in samplesheet fields.

## Resources and observability

Separate bounded job limits cover preprocessing/QC, alignment, tracks,
peak-calling, consensus/IDR, normalized tracks, differential analysis and
annotation. `TOTAL_CPU_BUDGET=auto` uses visible CPUs. Preflight writes
`resource_budget.tsv`; `RESOURCE_CHECK_MODE=fail` makes overcommit fatal.

`WRITE_CONSOLE_LOG`, `WRITE_COMMAND_LOG` and `WRITE_STRUCTURED_LOG` preserve
console output, command start/end/exit records and workflow events across
resumes. Stage timing includes the run ID and is append-only.

Further guidance:

- [Inputs and configuration](../docs/01_inputs_and_configuration.md)
- [Quality control](../docs/11_quality_control.md)
- [Replicates, controls, and design](../docs/13_replicates_and_design.md)
- [Annotation](../docs/15_annotation.md)
