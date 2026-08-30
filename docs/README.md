# chip2tracks documentation

This index groups the documentation by user task. Start with the root
[README](../README.md) for the workflow map, installation, quick start, stage
summary, and principal output tree.

## Prepare a run

| Topic | Document |
|---|---|
| Samplesheet rows, technical and biological replicates, matched controls, configuration safety, and plan/preflight validation | [Inputs and configuration](01_inputs_and_configuration.md) |
| Complete key list and example samplesheets | [`config/README.md`](../config/README.md) and [`config/config.conf.template`](../config/config.conf.template) |
| Exact stage order, enable/disable behavior, checkpoints, and declared outputs | [Pipeline stages](07_pipeline_stages.md) |
| Read-only inventory before deciding what server software and references can be reused | [Server audit](08_server_audit.md) |
| Genome-reference contract, reuse of server resources, blacklist choice, and exact BAM filtering order | [References, blacklist, and filtering](10_references_blacklist_and_filtering.md) |

## Understand the analysis

| Topic | Document |
|---|---|
| PE/SE signal units, filtering branches, CPM and DESeq2 track formulas, cohort isolation, and caller behavior | [Methods and normalization](02_methods.md) |
| The distinct roles of IgG/input/mock controls, primary raw-count models, DiffBind, and sensitivity analyses | [Controls and differential enrichment](03_differential_enrichment.md) |
| TSS/TES/gene-body aggregate plots, BED12 reference preparation, HPA subsets, manifests, and standalone reuse | [Metagene aggregate-signal module](06_metagene.md) |
| FastQC, alignment, complexity, fragment, cross-correlation, FRiP, control, replicate, TSS, and spike-in QC | [Quality control and interpretation](11_quality_control.md) |
| CPM, DESeq2 consensus, robust-CPM, spike-in, fold-enrichment, bedGraph/bigWig, and multimapper reporting | [Genomic coverage tracks and normalization](12_tracks_and_normalization.md) |
| Biological versus technical replicates, shared controls, consensus, IDR, blocks, and contrast direction | [Replicates, controls, and experimental design](13_replicates_and_design.md) |
| Inputs, models, outputs, statuses, and interpretation for all four differential variants | [Differential occupancy analysis](14_differential_occupancy.md) |
| Comprehensive feature/cCRE annotations, composition summaries, differential propagation, and interpretation limits | [Peak and differential-result annotation](15_annotation.md) |

## Operate and validate the workflow

| Topic | Document |
|---|---|
| Output organization, checkpoint recovery, partial reruns, and guarded cleanup | [Outputs and recovery](04_outputs_and_recovery.md) |
| Five-minute triage, stage diagnostics, safe restart points, and support bundle | [Troubleshooting and safe recovery](16_troubleshooting.md) |
| Immutable releases, environments, launchers, promotion, rollback, and large-node profile | [Shared-server installation and promotion](17_server_installation.md) |
| Samplesheet/config changes, compatibility defaults, and recovery from 0.1 | [Migration from 0.1 to 0.2](18_migration_0.1_to_0.2.md) |
| Scientific limitations and real-data pilot decisions | [Limitations and pilot decisions](05_limitations.md) |
| Synthetic coverage and outstanding Linux/real-data fixtures | [Test matrix](06_test_matrix.md) |

## Project information

- [Release history](../CHANGELOG.md)
- [Contribution requirements](../CONTRIBUTING.md)
- [Metagene shared-module interface](../common/metagene/README.md)
- [Implementation and release plan](09_implementation_and_release.md)

Documentation describes `chip2tracks` 0.2.0 unless a page explicitly says otherwise.
The executable behavior is defined by `chip2tracks.sh`, the scripts under
`scripts/`, and the validated configuration template.
