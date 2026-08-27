#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
temporary="$(mktemp -d)"
trap 'rm -rf -- "$temporary"' EXIT
release="$temporary/releases/release-one"
environment="$temporary/environment"
mkdir -p "$release" "$environment/bin"

cat > "$release/chip2tracks.sh" <<'EOF'
#!/usr/bin/env bash
exit 0
EOF
cat > "$environment/bin/bash" <<'EOF'
#!/bin/sh
printf '%s\n' "$1"
EOF
cat > "$environment/bin/python3" <<'EOF'
#!/bin/sh
exit 0
EOF
chmod +x "$release/chip2tracks.sh" "$environment/bin/bash" "$environment/bin/python3"
ln -s "$release" "$temporary/current"
if [[ ! -L "$temporary/current" ]]; then
    echo "Shared-launcher test skipped: this platform cannot create symbolic links"
    exit 0
fi

resolved="$(CHIP2TRACKS_ROOT="$temporary/current" CHIP2TRACKS_MAIN_ENV="$environment" \
    bash "$ROOT/utilities/chip2tracks_shared_launcher.sh")"
[[ "$resolved" == "$release/chip2tracks.sh" ]] || {
    echo "ERROR: shared launcher did not pin the immutable release: $resolved" >&2
    exit 1
}

echo "Shared-launcher immutable-release regression test passed"
