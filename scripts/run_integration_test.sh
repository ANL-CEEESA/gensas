#!/usr/bin/env bash
# Integration test: run power flow on d_70k_070.mat and verify:
#   1. The solver converges to alpha=1
#   2. "Computation time:" is printed (clean exit after computation)
#   3. Output .mat file is created and non-empty
#
# Usage: bash scripts/run_integration_test.sh
#   Must be run from the repository root.
#
# Note: the binary is invoked directly (not via "bazel run") to avoid
# the Bazel darwin sandbox restricting write access to arbitrary paths.

set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# The app prepends GetCurrentWorkingDir() to the -o argument, so we pass
# a bare filename and let the app place it in the CWD (REPO_DIR).
OUTPUT_NAME="gensas_integration_test_output.mat"
OUTPUT_MAT="${REPO_DIR}/${OUTPUT_NAME}"

cd "${REPO_DIR}"
rm -f "${OUTPUT_MAT}"

echo "=== GenSAS power flow integration test ==="
echo "Input : resources/psat_mat/d_70k_070.mat"
echo "Output: ${OUTPUT_MAT}"
echo ""

# Build first (exits non-zero on build failure, aborting the test).
bazel build //app:app

echo ""

# Run the binary directly — bypasses the Bazel sandbox so the output
# .mat file can be written to an arbitrary filesystem path.
# Pass a relative output name; the app prepends CWD, yielding REPO_DIR/OUTPUT_NAME.
set +e
APP_OUTPUT=$("${REPO_DIR}/bazel-bin/app/app" \
  -p \
  -f "${REPO_DIR}/resources/psat_mat/d_70k_070.mat" \
  -s 0.5 -l 28 -d 1e-5 \
  -o "${OUTPUT_NAME}" \
  2>&1)
APP_EXIT=$?
set -e

echo "${APP_OUTPUT}"
echo ""

PASS=1

# ── Check 1: process exited cleanly ────────────────────────────────────────
if [ "${APP_EXIT}" -ne 0 ]; then
  echo "FAIL: app exited with code ${APP_EXIT} (crash or error)"
  PASS=0
else
  echo "PASS: clean exit (exit code 0)"
fi

# ── Check 2: solver converged to alpha=1 ────────────────────────────────────
if echo "${APP_OUTPUT}" | grep -q "Step=1,"; then
  echo "PASS: solver converged to alpha=1"
else
  echo "FAIL: 'Step=1,' not found — solver did not reach alpha=1"
  PASS=0
fi

# ── Check 3: timing line present (confirms no crash after computation) ───────
if echo "${APP_OUTPUT}" | grep -q "Computation time:"; then
  echo "PASS: 'Computation time:' found"
else
  echo "FAIL: 'Computation time:' not found — likely crashed before clean exit"
  PASS=0
fi

# ── Check 4: output .mat file exists and is non-empty ───────────────────────
if [ -f "${OUTPUT_MAT}" ]; then
  MAT_SIZE=$(wc -c < "${OUTPUT_MAT}")
  if [ "${MAT_SIZE}" -gt 0 ]; then
    echo "PASS: output .mat created (${MAT_SIZE} bytes)"
  else
    echo "FAIL: output .mat file is empty"
    PASS=0
  fi
else
  echo "FAIL: output .mat file was not created at ${OUTPUT_MAT}"
  PASS=0
fi

# ── Summary ─────────────────────────────────────────────────────────────────
echo ""
if [ "${PASS}" -eq 1 ]; then
  echo "=== All checks passed ==="
  exit 0
else
  echo "=== One or more checks FAILED ==="
  exit 1
fi
