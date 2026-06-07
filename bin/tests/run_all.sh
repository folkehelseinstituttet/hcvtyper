#!/bin/bash

# run_all.sh - Single-command runner for the HCVTyper R regression guard (D-01)
# Runs every bin/tests/test_*.R via Rscript and exits non-zero if ANY fails.
# This is BOTH the local developer command and the command the CI job (Plan 03)
# invokes. New test_*.R files are picked up automatically by the glob (D-04:
# test_denovo_confirm.R folds in unchanged).
#
# Usage (from any cwd):
#   bash bin/tests/run_all.sh

# Resolve this script's own directory so it runs from any cwd.
here="$(cd "$(dirname "$0")" && pwd)"

# NOTE: deliberately NOT using `set -e`. We want EVERY test to run even if an
# early one fails, then propagate failure via the explicit accumulator below.
status=0

for t in "$here"/test_*.R; do
  echo "=================================================================="
  echo ">>> Running $(basename "$t")"
  echo "=================================================================="
  # Invoke as `Rscript "$t"` (real path) so each test's --file= self-location
  # resolves bin/ correctly. Any non-zero exit flips the accumulator.
  Rscript "$t" || status=1
done

echo "=================================================================="
if [ "$status" -eq 0 ]; then
  echo ">>> ALL R TESTS PASSED"
else
  echo ">>> ONE OR MORE R TESTS FAILED"
fi
echo "=================================================================="

exit $status
