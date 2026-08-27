#!/usr/bin/env bash
set -euo pipefail

WORKFLOW_ROOT="${CHIP2TRACKS_ROOT:-/opt/bioinformatics/workflows/chip2tracks/current}"
MAIN_ENV="${CHIP2TRACKS_MAIN_ENV:-/opt/miniconda/envs/chip2tracks-0.1.0}"

[[ -x "$WORKFLOW_ROOT/chip2tracks.sh" ]] || {
    echo "ERROR: chip2tracks release is not executable: $WORKFLOW_ROOT/chip2tracks.sh" >&2
    exit 127
}
[[ -x "$MAIN_ENV/bin/bash" && -x "$MAIN_ENV/bin/python3" ]] || {
    echo "ERROR: chip2tracks main environment is incomplete: $MAIN_ENV" >&2
    exit 127
}

# Do not activate Conda or alter the caller's shell. The workflow receives its
# versioned main tools, followed by system-wide sidecar launchers.
export PATH="$MAIN_ENV/bin:/usr/local/bin:/usr/bin:/bin"
exec "$MAIN_ENV/bin/bash" "$WORKFLOW_ROOT/chip2tracks.sh" "$@"
