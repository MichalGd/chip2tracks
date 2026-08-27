#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mapfile -t scripts < <(find "$ROOT" -type f -name '*.sh' -print | sort)
(( ${#scripts[@]} > 0 ))
for script in "${scripts[@]}"; do
    bash -n "$script"
    [[ -x "$script" ]] || {
        echo "ERROR: shell entrypoint is not executable: $script" >&2
        exit 1
    }
done
echo "Bash syntax and executable modes OK: ${#scripts[@]} files"
