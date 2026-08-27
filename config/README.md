# Configuration

Copy `config.conf.template` outside the repository and replace every path. It is
a non-executable `KEY=VALUE` file: unknown keys, shell expansion, unsafe values,
and missing required options are rejected.

Choose exactly one run-wide profile:

```text
ASSAY_PROFILE=chipseq
ASSAY_PROFILE=chipmentation
```

`ADAPTER_PRESET=auto` is preferred when the local kit is unknown. The workflow
records Trim Galore detection evidence. Use `illumina`, `nextera`, or `custom`
when known; custom mode requires `CUSTOM_ADAPTER_R1` and optionally R2. Adapter
auto-detection is a convenience, not a substitute for laboratory provenance.

Matched controls are optional by default. A target with no control uses an
empty `control_id`. Set `ALLOW_CONTROL_FREE_PEAKCALL=false` if a project
requires controls. Compatible
shared controls are allowed by default and still undergo metadata validation.

The supplied public template disables spike-in. To enable the laboratory's dm6
cell spike-in protocol, set `SPIKEIN_MODE=dm6`, replace the composite-index and
dm6 paths, and populate the three spike-in sample-sheet columns.

Every selected assembly needs `INDEX`, `FASTA`, `CHROM_SIZES`,
`CANONICAL_CONTIGS`, `GTF`, `BLACKLIST`, and `EFFECTIVE_GENOME_SIZE` reference
keys. Effective genome size is not the same as the chromosome-length sum.

The CSV header is exact and contains 25 columns. Use the files under
`config/examples/` as schema examples. PE requires both FASTQs; SE requires an
empty `fastq_2`. Each row is a technical unit. Rows sharing
`sample_id + replicate` are merged only when all biological metadata agree.

Relevant biological fields include assay, factor, antibody, narrow/broad/mixed
class, condition, treatment, cell type, biological/technical replicate,
control type/linkage, duplicate policy, blacklist, spike-in ratio/stage/lot,
batch, donor and a safe unique output prefix.

One samplesheet is one compatible target/antibody and peak universe. Put
different factors, antibodies, layouts, target classes, or analysis policies in
separate runs. Conditions and replicates remain together for consensus and
differential analysis.

The template enables bigWig and uncompressed bedGraph output for the analysis,
MAPQ-0 duplicate-retained, MAPQ-0 duplicate-removed, and MAPQ-30
duplicate-removed CPM families. The corresponding `GENERATE_*` switches are
ordinary `config.conf` values. Stringent robust-CPM is also enabled.

`TOTAL_CPU_BUDGET=auto` uses the CPUs visible to the process. Preflight records
per-stage maximum CPU requests and warns on overcommit by default; use
`RESOURCE_CHECK_MODE=fail` for strict enforcement.
