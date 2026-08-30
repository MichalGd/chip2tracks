#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bash "$ROOT/tests/check_bash_syntax.sh"
bash "$ROOT/tests/test_parallel_jobs.sh"
bash "$ROOT/tests/test_consensus_batch.sh"
bash "$ROOT/tests/test_signal_count_policies.sh"
bash "$ROOT/tests/test_shared_launcher.sh"
bash "$ROOT/tests/test_filter_cleanup.sh"
bash "$ROOT/tests/test_cleanup_checkpoints.sh"
bash "$ROOT/tests/test_coverage_policies.sh"
bash "$ROOT/tests/test_spikein_parallel.sh"
bash "$ROOT/tests/test_peakcall_tolerance.sh"
bash "$ROOT/tests/test_epic2_smoke.sh"
bash "$ROOT/tests/test_normalization_tolerance.sh"
bash "$ROOT/tests/test_annotation_sorting.sh"
bash "$ROOT/tests/test_reporting.sh"
python3 -m unittest discover -s "$ROOT/tests" -p 'test_*.py' -v
echo "All chip2tracks tests passed"
