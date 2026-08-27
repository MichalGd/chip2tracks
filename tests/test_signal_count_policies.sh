#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
temporary="$(mktemp -d)"
trap 'rm -rf -- "$temporary"' EXIT
mkdir -p "$temporary/bin"

cat > "$temporary/bin/samtools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
[[ "${1:-}" == "view" ]]
shift
exclude=""
while (( $# )); do
    if [[ "$1" == "-F" ]]; then exclude="$2"; shift 2; else shift; fi
done
case "$exclude" in
    2828) echo 120 ;;
    3852) echo 100 ;;
    2820) echo 240 ;;
    3844) echo 200 ;;
    *) echo "unexpected exclusion mask: $exclude" >&2; exit 95 ;;
esac
EOF
chmod +x "$temporary/bin/samtools"

PATH="$temporary/bin:$PATH"
# shellcheck source=../scripts/lib/common.sh
source "$ROOT/scripts/lib/common.sh"

[[ "$(signal_count example.bam PE retain)" == 120 ]]
[[ "$(signal_count example.bam PE removed)" == 100 ]]
[[ "$(signal_count example.bam SE retained)" == 240 ]]
[[ "$(signal_count example.bam SE remove)" == 200 ]]
if (signal_count example.bam PE ambiguous >/dev/null 2>&1); then
    echo "ERROR: invalid duplicate policy was accepted" >&2
    exit 1
fi

echo "Signal-count duplicate-policy regression test passed"
