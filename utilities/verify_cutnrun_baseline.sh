#!/usr/bin/env bash
set -euo pipefail
source_dir="${1:?Usage: verify_cutnrun_baseline.sh /path/to/cutnrun2tracks}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
manifest="$(cd "$script_dir/../provenance" && pwd)/cutnrun2tracks_v0.2.8.sha256"
[[ -d "$source_dir" && -s "$manifest" ]] || { echo "ERROR: source directory or manifest missing" >&2; exit 1; }
(
    cd "$source_dir"
    sha256sum --check "$manifest"
)
