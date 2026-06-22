#!/usr/bin/env bash
#
# build_and_run.sh - build the linux/amd64 image and run the MSAMDD smoke test.
#
# Run this from anywhere; it figures out the repo root itself. Prerequisite: you
# have placed the CPLEX Linux x86-64 installer at
#     MSAMDD_MWE/docker/cplex_linux.bin
# See docker/README.md for how to obtain it.
#
# Results (alignment + logs) are written back to MSAMDD_MWE/results/ on the host.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # .../MSAMDD_MWE/docker
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"               # .../MSAMDD-Ling
IMAGE="msamdd-ling:amd64"

cd "${REPO_ROOT}"

if [ ! -f "${SCRIPT_DIR}/cplex_linux.bin" ]; then
  echo "ERROR: missing CPLEX installer at MSAMDD_MWE/docker/cplex_linux.bin" >&2
  echo "       Download the Linux x86-64 CPLEX .bin from IBM, move it there, and" >&2
  echo "       rename it to cplex_linux.bin. See MSAMDD_MWE/docker/README.md." >&2
  exit 1
fi

echo "== [1/2] Building image ${IMAGE} (linux/amd64; emulated on Apple Silicon) =="
echo "   This is slow the first time (CPLEX install + compile under emulation)."
docker build --platform linux/amd64 \
  -f MSAMDD_MWE/docker/Dockerfile \
  -t "${IMAGE}" \
  .

echo ""
echo "== [2/2] Running smoke test in a container =="
mkdir -p "${REPO_ROOT}/MSAMDD_MWE/results"
docker run --rm --platform linux/amd64 \
  -v "${REPO_ROOT}/MSAMDD_MWE/results:/work/MSAMDD_MWE/results" \
  "${IMAGE}"

echo ""
echo "Done. Logs and the alignment are in MSAMDD_MWE/results/ on your Mac."
