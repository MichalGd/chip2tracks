#!/usr/bin/env bash
set -euo pipefail

EPIC2_COMMAND="${EPIC2_COMMAND:-epic2}"
SAMTOOLS_COMMAND="${SAMTOOLS_COMMAND:-samtools}"

usage() {
    cat <<'USAGE'
Usage: smoke_test_epic2.sh

Run a tiny paired-end epic2 call using EPIC2_COMMAND (default: epic2) and
SAMTOOLS_COMMAND (default: samtools). No user data are read or modified.
USAGE
}

case "${1:-}" in
    -h|--help) usage; exit 0 ;;
    "") ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
esac

command -v "$EPIC2_COMMAND" >/dev/null 2>&1 || {
    echo "ERROR: epic2 command not found: $EPIC2_COMMAND" >&2
    exit 2
}
command -v "$SAMTOOLS_COMMAND" >/dev/null 2>&1 || {
    echo "ERROR: samtools command not found: $SAMTOOLS_COMMAND" >&2
    exit 2
}

temporary="$(mktemp -d "${TMPDIR:-/tmp}/chip2tracks-epic2-smoke.XXXXXX")"
trap 'rm -rf -- "$temporary"' EXIT
printf 'chrSmoke\t10000\n' > "$temporary/chrom.sizes"

write_sam() {
    local destination="$1" offset="$2" index start mate
    {
        printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrSmoke\tLN:10000\n'
        for index in $(seq 1 20); do
            start=$((offset + index * 100))
            mate=$((start + 100))
            printf 'pair%s\t99\tchrSmoke\t%s\t60\t50M\t=\t%s\t150\t%s\t%s\n' \
                "$index" "$start" "$mate" \
                'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA' \
                'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
            printf 'pair%s\t147\tchrSmoke\t%s\t60\t50M\t=\t%s\t-150\t%s\t%s\n' \
                "$index" "$mate" "$start" \
                'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT' \
                'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
        done
    } > "$destination"
}

write_sam "$temporary/treatment.sam" 100
write_sam "$temporary/control.sam" 5000
"$SAMTOOLS_COMMAND" view -b -o "$temporary/treatment.bam" "$temporary/treatment.sam"
"$SAMTOOLS_COMMAND" view -b -o "$temporary/control.bam" "$temporary/control.sam"
"$SAMTOOLS_COMMAND" index "$temporary/treatment.bam"
"$SAMTOOLS_COMMAND" index "$temporary/control.bam"

"$EPIC2_COMMAND" \
    --treatment "$temporary/treatment.bam" \
    --control "$temporary/control.bam" \
    --chromsizes "$temporary/chrom.sizes" \
    --effective-genome-fraction 1 \
    --bin-size 200 \
    --gaps-allowed 3 \
    --fragment-size 200 \
    --false-discovery-rate-cutoff 0.05 \
    --keep-duplicates --mapq 0 --guess-bampe \
    --output "$temporary/epic2.tsv"

echo "EPIC2_PAIRED_END_SMOKE_TEST_PASSED"
