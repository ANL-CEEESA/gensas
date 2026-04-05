#!/usr/bin/env bash
# Run clang-format on all C++ source files in the repository.
#
# Usage:
#   bash scripts/format.sh           # format files in-place
#   bash scripts/format.sh --check   # check only, exit non-zero if any file
#                                    # needs formatting (used by CI)
#
# Requirements:
#   clang-format must be on PATH.
#   Install (both macOS and Linux):
#     pip install clang-format==22.1.2
#
# Override the binary with:
#   CLANG_FORMAT=/path/to/clang-format bash scripts/format.sh

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

CLANG_FORMAT="${CLANG_FORMAT:-clang-format}"

CHECK_ONLY=0

# ── Parse arguments ──────────────────────────────────────────────────────────
for arg in "$@"; do
  case "${arg}" in
    --check) CHECK_ONLY=1 ;;
    *) echo "Unknown argument: ${arg}"; exit 1 ;;
  esac
done

# ── Verify clang-format is available ─────────────────────────────────────────
if ! command -v "${CLANG_FORMAT}" &>/dev/null; then
  echo "error: '${CLANG_FORMAT}' not found."
  echo "  Install: pip install clang-format==22.1.2"
  exit 1
fi

CF_VERSION=$("${CLANG_FORMAT}" --version)
echo "Using: ${CF_VERSION}"
echo ""

# ── Collect source files ──────────────────────────────────────────────────────
mapfile -t FILES < <(
  find "${REPO_DIR}" \
    -not \( -path "${REPO_DIR}/bazel-*" -prune \) \
    \( -name "*.cpp" -o -name "*.h" \) \
  | sort
)
echo "Found ${#FILES[@]} file(s) to process."
echo ""

# ── Format or check ───────────────────────────────────────────────────────────
if [ "${CHECK_ONLY}" -eq 1 ]; then
  FAILED=0
  for f in "${FILES[@]}"; do
    if ! "${CLANG_FORMAT}" --dry-run --Werror "${f}" 2>/dev/null; then
      echo "  needs formatting: ${f#${REPO_DIR}/}"
      FAILED=1
    fi
  done
  echo ""
  if [ "${FAILED}" -eq 1 ]; then
    echo "FAIL: one or more files need formatting."
    echo "      Run  bash scripts/format.sh  to fix."
    exit 1
  else
    echo "PASS: all files are correctly formatted."
  fi
else
  for f in "${FILES[@]}"; do
    "${CLANG_FORMAT}" -i "${f}"
  done
  echo "Done. All files formatted in-place."
fi
