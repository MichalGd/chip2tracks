#!/usr/bin/env bash
set -euo pipefail

SPP_ENV="${CHIP2TRACKS_SPP_ENV:-/opt/miniconda/envs/chip2tracks-spp-1.2.2}"
RSCRIPT="${SPP_ENV}/bin/Rscript"
RUN_SPP="${SPP_ENV}/bin/run_spp.R"

[[ -x "$RSCRIPT" ]] || {
    echo "ERROR: SPP-sidecar Rscript is not executable: $RSCRIPT" >&2
    exit 127
}
[[ -f "$RUN_SPP" ]] || {
    echo "ERROR: phantompeakqualtools script is missing: $RUN_SPP" >&2
    exit 127
}

exec "$RSCRIPT" "$RUN_SPP" "$@"
