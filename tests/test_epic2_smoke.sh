#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf -- "$TMP"' EXIT
mkdir -p "$TMP/bin"

cat > "$TMP/bin/samtools" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
[[ "$1" == "view" ]]
shift
output=""; input=""
while (( $# )); do
    case "$1" in
        -b) shift ;;
        -o) output="$2"; shift 2 ;;
        *) input="$1"; shift ;;
    esac
done
[[ -n "$output" && -s "$input" ]]
cp "$input" "$output"
MOCK

cat > "$TMP/bin/epic2" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
treatment=""; control=""; chromsizes=""; output=""
guess=false
while (( $# )); do
    case "$1" in
        --treatment) treatment="$2"; shift 2 ;;
        --control) control="$2"; shift 2 ;;
        --chromsizes) chromsizes="$2"; shift 2 ;;
        --output) output="$2"; shift 2 ;;
        --guess-bampe) guess=true; shift ;;
        --keep-duplicates) shift ;;
        --effective-genome-fraction|--bin-size|--gaps-allowed|--fragment-size|--false-discovery-rate-cutoff|--mapq) shift 2 ;;
        *) echo "unexpected epic2 argument: $1" >&2; exit 2 ;;
    esac
done
[[ -s "$treatment" && -s "$control" && -s "$chromsizes" && "$guess" == "true" && -n "$output" ]]
: > "$output"
MOCK
chmod +x "$TMP/bin/samtools" "$TMP/bin/epic2"

result="$(EPIC2_COMMAND="$TMP/bin/epic2" SAMTOOLS_COMMAND="$TMP/bin/samtools" \
    bash "$ROOT/utilities/smoke_test_epic2.sh")"
grep -q '^EPIC2_PAIRED_END_SMOKE_TEST_PASSED$' <<< "$result"
echo "epic2 paired-end smoke-test utility regression test passed"
