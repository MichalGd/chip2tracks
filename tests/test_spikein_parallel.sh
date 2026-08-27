#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf -- "$TMP"' EXIT
mkdir -p "$TMP/bin" "$TMP/output/00_metadata" "$TMP/output/03_alignment/analysis" \
    "$TMP/output/03_alignment/spikein/spike" "$TMP/output/03_alignment/metrics"
printf 'chr1\t1000\n' > "$TMP/host.chrom.sizes"
printf 'dm6_chr2L\t1000\n' > "$TMP/spike.chrom.sizes"

for key in S1 S2; do
    printf 'host\n' > "$TMP/output/03_alignment/analysis/$key.host.analysis.bam"
    printf 'spike\n' > "$TMP/output/03_alignment/spikein/spike/$key.dm6.bam"
done

cat > "$TMP/bin/picard" <<'EOF'
#!/usr/bin/env bash
for arg in "$@"; do case "$arg" in I=*) input="${arg#I=}" ;; O=*) output="${arg#O=}" ;; M=*) metrics="${arg#M=}" ;; esac; done
cp "$input" "$output"; printf 'metrics\n' > "$metrics"; printf 'index\n' > "${output}.bai"
EOF
cat > "$TMP/bin/samtools" <<'EOF'
#!/usr/bin/env bash
command="$1"; shift
case "$command" in
  view)
    if [[ " $* " == *" -c "* ]]; then echo 100; exit 0; fi
    while (( $# )); do if [[ "$1" == "-o" ]]; then output="$2"; shift 2; else input="$1"; shift; fi; done
    cp "$input" "$output"
    ;;
  index) output="${@: -1}"; printf 'index\n' > "${output}.bai" ;;
  quickcheck) exit 0 ;;
  *) exit 0 ;;
esac
EOF
cat > "$TMP/bin/bamCoverage" <<'EOF'
#!/usr/bin/env bash
while (( $# )); do if [[ "$1" == "--outFileName" ]]; then output="$2"; shift 2; else shift; fi; done
printf 'chr1\t0\t10\t1\n' > "$output"
EOF
cat > "$TMP/bin/bedtools" <<'EOF'
#!/usr/bin/env bash
while (( $# )); do if [[ "$1" == "-i" || "$1" == "-abam" ]]; then input="$2"; shift 2; else shift; fi; done
cat "$input"
EOF
cat > "$TMP/bin/bedGraphToBigWig" <<'EOF'
#!/usr/bin/env bash
cp "$1" "$3"
EOF
cat > "$TMP/bin/python3" <<'EOF'
#!/usr/bin/env bash
counts="$2"; output="$3"
printf 'sample_key\tlayout\tcohort_id\tspikein_stage\tspikein_lot\tspikein_to_host_ratio\thost_observations\tspike_observations\tspike_fraction\thost_scale_factor\tspike_cpm_scale_factor\tcohort_median_host_scale\tstatus\tfailure_reasons\twarnings\n' > "$output"
tail -n +2 "$counts" | while IFS=$'\t' read -r key layout cohort stage lot ratio host spike; do
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t0.5\t500\t10000\t500\tPASS\t.\t.\n' \
        "$key" "$layout" "$cohort" "$stage" "$lot" "$ratio" "$host" "$spike"
done >> "$output"
EOF
chmod +x "$TMP/bin/"*

cat > "$TMP/output/00_metadata/sample_manifest.tsv" <<'EOF'
sample_key	sample_id	replicate	layout	genome	assay_profile	factor	antibody_id	target_class	condition	treatment	cell_type	is_control	control_type	control_id	control_key	analysis_duplicate_policy	blacklist	spikein_to_host_ratio	spikein_stage	spikein_lot	batch	donor	output_prefix	technical_units	fastq_1_list	fastq_2_list	cohort_id	cohort_key	primary_peak_caller	primary_peak_class
S1	S1	1	PE	hg38	chipseq	CTCF	AB	narrow	WT	none	CELL	FALSE	none	.	.	remove	blacklist	0.05	cells	LOT	B1	D1	S1	1	fq1	fq2	C1	K1	macs3	narrow
S2	S2	1	PE	hg38	chipseq	CTCF	AB	narrow	KO	none	CELL	FALSE	none	.	.	remove	blacklist	0.05	cells	LOT	B1	D2	S2	1	fq1	fq2	C2	K2	macs3	narrow
EOF
cat > "$TMP/config.sh" <<EOF
OUTPUT_DIR='$TMP/output'
SPIKEIN_MODE=dm6
SPIKEIN_DUPLICATE_POLICY=remove
SPIKEIN_MIN_MAPQ=30
SPIKEIN_BLACKLIST=
SPIKEIN_SCALE_TARGET=1000000
SPIKEIN_MIN_OBSERVATIONS_FAIL=1
SPIKEIN_MIN_OBSERVATIONS_WARN=1
SPIKEIN_WARN_LOW_FRACTION=0.001
SPIKEIN_WARN_HIGH_FRACTION=0.2
ALLOW_FAILED_SPIKEIN=false
PICARD_COMMAND=picard
THREADS_SAMTOOLS=1
SPIKEIN_PARALLEL_JOBS=2
THREADS_BAMCOVERAGE=1
TRACK_BIN_SIZE=10
GENERATE_COVERAGE_BIGWIGS=true
GENERATE_COVERAGE_BEDGRAPHS=true
GENERATE_SPIKEIN_CONTROL_TRACKS=true
SPIKEIN_CHROM_SIZES='$TMP/spike.chrom.sizes'
CHROM_SIZES_HG38='$TMP/host.chrom.sizes'
SE_SIGNAL_MODE=extend
SE_FRAGMENT_LENGTH=auto
SE_FRAGMENT_LENGTH_FALLBACK=200
WRITE_COMMAND_LOG=true
EOF

PATH="$TMP/bin:$PATH" C2T_CONFIG="$TMP/config.sh" bash "$ROOT/scripts/spikein_batch.sh"
for key in S1 S2; do
    cohort="C${key#S}"
    test -s "$TMP/output/04_tracks/spikein/$cohort/$key.SpikeInScaled.dm6.bw"
    test -s "$TMP/output/04_tracks/spikein/$cohort/$key.SpikeInScaled.dm6.bedGraph"
    test -s "$TMP/output/04_tracks/spikein_control/$key.dm6.CPM.bw"
    test -s "$TMP/output/04_tracks/spikein_control/$key.dm6.CPM.bedGraph"
done
test "$(tail -n +2 "$TMP/output/06_qc/spikein/spikein_counts.tsv" | wc -l)" -eq 2
echo "Parallel spike-in stage regression test passed"
