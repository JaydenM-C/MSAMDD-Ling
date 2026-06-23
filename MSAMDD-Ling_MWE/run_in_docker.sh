#!/usr/bin/env bash
#
# run_in_docker.sh - run the Phase 4 validation inside the linux/amd64 image.
#
# Why Docker: msa_aff (CPLEX) and the bundled muscle are x86-64 Linux binaries
# and do not run natively on Apple Silicon.
#
# This REUSES the image already built in Phase 0 (msamdd-ling:amd64) and
# bind-mounts your live repo at /work, so the current scripts and source are
# always used WITHOUT rebuilding the image. The image only supplies the
# toolchain + CPLEX; everything in /work comes from your working tree. That
# avoids the slow, network-dependent apt/CPLEX build on every run.
#
# The image is built only if it is missing (or FORCE_BUILD=1). Building needs the
# CPLEX installer at MSAMDD_MWE/docker/cplex_linux.bin and a working network
# inside Docker; running a stage needs neither.
#
# Prerequisites:
#   * Docker Desktop, "Use Rosetta for x86_64/amd64 emulation" enabled.
#   * The SeqAn submodule initialised:  git submodule update --init
#   * Image present from Phase 0, OR cplex_linux.bin in place to build it once.
#
# Usage (from the repository root):
#   bash MSAMDD-Ling_MWE/run_in_docker.sh                 # all stages
#   bash MSAMDD-Ling_MWE/run_in_docker.sh incumbent       # one stage
#   bash MSAMDD-Ling_MWE/run_in_docker.sh optimum
#   bash MSAMDD-Ling_MWE/run_in_docker.sh linguistic
#   FORCE_BUILD=1 bash MSAMDD-Ling_MWE/run_in_docker.sh   # rebuild the image first
#   TIME_LIMIT=1800 ENTRIES="Data/bali1/1tgxA.fa Data/bali1/1aab.fa" \
#     bash MSAMDD-Ling_MWE/run_in_docker.sh optimum

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="msamdd-ling:amd64"
STAGE="${1:-all}"

cd "${REPO_ROOT}"

if ! command -v docker >/dev/null 2>&1; then
  echo "ERROR: docker not found. Install Docker Desktop and enable Rosetta emulation." >&2
  exit 1
fi
if [ ! -f "includes/seqan/include/seqan/version.h" ]; then
  echo "ERROR: SeqAn submodule not initialised (includes/seqan is empty)." >&2
  echo "       Run:  git submodule update --init" >&2
  exit 1
fi

# --- ensure the image exists; build only if missing or forced -----------------
if [ "${FORCE_BUILD:-0}" != "1" ] && docker image inspect "${IMAGE}" >/dev/null 2>&1; then
  echo "== Reusing existing image ${IMAGE} (set FORCE_BUILD=1 to rebuild). =="
else
  if [ ! -f "MSAMDD_MWE/docker/cplex_linux.bin" ]; then
    echo "ERROR: image ${IMAGE} not found and no CPLEX installer to build it." >&2
    echo "       Put the Linux x86-64 CPLEX .bin at MSAMDD_MWE/docker/cplex_linux.bin" >&2
    echo "       (see MSAMDD_MWE/docker/README.md), then re-run." >&2
    exit 1
  fi
  echo "== Building ${IMAGE} (one-time; needs network inside Docker) =="
  if ! docker build --platform linux/amd64 \
        -f MSAMDD_MWE/docker/Dockerfile -t "${IMAGE}" .; then
    echo "" >&2
    echo "Image build failed. The log above usually shows an apt/DNS error such as" >&2
    echo "  'Could not resolve archive.ubuntu.com' -- a transient Docker network issue" >&2
    echo "under amd64 emulation, not a problem with this repo. Try, in order:" >&2
    echo "  1. Re-run this command (apt layers resume from cache)." >&2
    echo "  2. Restart Docker Desktop; disable any VPN; then re-run." >&2
    echo "  3. Docker Desktop > Settings > Resources > Network, or add a DNS server" >&2
    echo "     (e.g. 8.8.8.8) in Docker's daemon settings, then re-run." >&2
    echo "Once the image exists, later runs reuse it and need no network." >&2
    exit 1
  fi
fi

# --- run the stage with the live repo bind-mounted ----------------------------
# The container builds msa_aff / seqan_warmstart fresh from the mounted source
# (Linux x86-64), so a couple of Linux build artifacts (obj/, msa_aff,
# seqan_warmstart) will appear in your tree. They are git-ignored and cannot run
# on macOS; `make clean` clears them if you want.
mkdir -p "${REPO_ROOT}/MSAMDD-Ling_MWE/results"
echo "== Running Phase 4 validation (stage: ${STAGE}) =="
docker run --rm --platform linux/amd64 \
  -e STAGE -e TIME_LIMIT -e ENTRIES -e OP -e EP -e SUBMAT \
  -v "${REPO_ROOT}:/work" \
  -w /work \
  "${IMAGE}" \
  bash MSAMDD-Ling_MWE/validate_all.sh "${STAGE}"

echo ""
echo "Done. Logs and alignments are in MSAMDD-Ling_MWE/results/ on your Mac."
