#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
if [[ -z "${KEEP_TEST_TMP:-}" ]]; then trap 'rm -rf -- "$TMP"' EXIT; else echo "Keeping test directory: $TMP"; fi
mkdir -p "$TMP/bin" "$TMP/output/00_metadata" "$TMP/output/03_alignment/analysis" \
    "$TMP/output/03_alignment/filtered/q0_dup-retained" \
    "$TMP/output/03_alignment/filtered/q0_dup-removed" \
    "$TMP/output/03_alignment/filtered/q30_dup-removed"

printf 'chr1\t1000\n' > "$TMP/chrom.sizes"
for bam in \
    "$TMP/output/03_alignment/analysis/S1.host.analysis.bam" \
    "$TMP/output/03_alignment/filtered/q0_dup-retained/S1.host.q0.dup-retained.bam" \
    "$TMP/output/03_alignment/filtered/q0_dup-removed/S1.host.q0.dup-removed.bam" \
    "$TMP/output/03_alignment/filtered/q30_dup-removed/S1.host.q30.dup-removed.bam"; do
    printf 'bam\n' > "$bam"
done

cat > "$TMP/bin/samtools" <<'EOF'
#!/usr/bin/env bash
if [[ "$*" == *q0.dup-retained.bam* ]]; then count=200; mapq0=20; low=40; xs=30
elif [[ "$*" == *q0.dup-removed.bam* ]]; then count=160; mapq0=16; low=32; xs=24
elif [[ "$*" == *q30.dup-removed.bam* ]]; then count=100; mapq0=0; low=0; xs=5
else count=100; mapq0=0; low=0; xs=5
fi
if [[ " $* " == *" -c "* ]]; then echo "$count"; exit 0; fi
for ((i=1; i<=count; i++)); do
    if (( i <= mapq0 )); then mapq=0; elif (( i <= low )); then mapq=10; else mapq=30; fi
    if (( i <= xs )); then tag=$'\tXS:i:20'; else tag=''; fi
    printf 'read%s\t66\tchr1\t1\t%s\t10M\t=\t1\t10\tAAAAAAAAAA\tFFFFFFFFFF%s\n' "$i" "$mapq" "$tag"
done
EOF
cat > "$TMP/bin/bamCoverage" <<'EOF'
#!/usr/bin/env bash
while (( $# )); do
    if [[ "$1" == "--outFileName" ]]; then output="$2"; shift 2; else shift; fi
done
printf 'chr1\t0\t10\t1\nchr1\t10\t20\t2\n' > "$output"
EOF
cat > "$TMP/bin/bedtools" <<'EOF'
#!/usr/bin/env bash
while (( $# )); do
    if [[ "$1" == "-i" ]]; then input="$2"; shift 2; else shift; fi
done
cat "$input"
EOF
cat > "$TMP/bin/bedGraphToBigWig" <<'EOF'
#!/usr/bin/env bash
cp "$1" "$3"
EOF
chmod +x "$TMP/bin/"*

cat > "$TMP/output/00_metadata/sample_manifest.tsv" <<'EOF'
sample_key	sample_id	replicate	layout	genome	rest
S1	S1	1	PE	hg38	.
EOF
cat > "$TMP/config.sh" <<EOF
OUTPUT_DIR='$TMP/output'
CHROM_SIZES_HG38='$TMP/chrom.sizes'
GENERATE_CPM_TRACKS=true
GENERATE_CPM_PERMISSIVE_TRACKS=true
GENERATE_CPM_INTERMEDIATE_TRACKS=true
GENERATE_CPM_STRINGENT_TRACKS=true
GENERATE_COVERAGE_BIGWIGS=true
GENERATE_COVERAGE_BEDGRAPHS=true
TRACK_BIN_SIZE=10
THREADS_BAMCOVERAGE=1
TRACK_PARALLEL_JOBS=1
BAMCOVERAGE_COMMON_ARGS=
SE_SIGNAL_MODE=extend
SE_FRAGMENT_LENGTH=auto
SE_FRAGMENT_LENGTH_FALLBACK=200
MACS3_COMMAND=macs3
EOF

PATH="$TMP/bin:$PATH" C2T_CONFIG="$TMP/config.sh" \
    bash "$ROOT/scripts/coverage_batch.sh"

for relative in \
    S1.CPM.bw S1.CPM.bedGraph \
    permissive/S1.CPM.bw permissive/S1.CPM.bedGraph \
    intermediate/S1.CPM.bw intermediate/S1.CPM.bedGraph \
    stringent/S1.CPM.bw stringent/S1.CPM.bedGraph; do
    test -s "$TMP/output/04_tracks/cpm/$relative"
done
if find "$TMP/output/04_tracks/cpm" -name '*.bedGraph.gz' -print -quit | grep -q .; then
    echo "ERROR: CPM stage emitted legacy compressed bedGraph output" >&2
    exit 1
fi
grep -q $'permissive\tq0_dup-retained\t0\tretained\tfragment\t200\t5000' \
    "$TMP/output/04_tracks/cpm/permissive/S1.normalization_metadata.tsv"
grep -q $'stringent\tq30_dup-removed\t30\tremoved\tfragment\t100\t10000' \
    "$TMP/output/04_tracks/cpm/stringent/S1.normalization_metadata.tsv"
grep -q $'S1\tpermissive\tcpm/permissive\tdeseq2_robust_cpm/permissive\tq0_dup-retained\t0\tretained\tfragment\t200\t20\t10.0000\t40\t20.0000\t30\t15.0000' \
    "$TMP/output/04_tracks/cpm/mapping_composition.tsv"
grep -q $'S1\tstringent\tcpm/stringent\tdeseq2_robust_cpm/stringent\tq30_dup-removed\t30\tremoved\tfragment\t100\t0\t0.0000\t0\t0.0000\t5\t5.0000' \
    "$TMP/output/04_tracks/cpm/mapping_composition.tsv"
test -s "$TMP/output/04_tracks/cpm/mapping_composition_definitions.tsv"
echo "Coverage-policy track regression test passed"
