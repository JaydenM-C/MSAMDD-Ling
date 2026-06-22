#!/usr/bin/env bash
#
# smoke_test.sh - minimal working example / toolchain smoke test for MSAMDD-Ling.
#
# Goal: confirm the bundled exact affine MSA solver (msa_aff) builds and produces
# a genuine multiple sequence alignment on the tiny real dataset shipped with the
# repo, BEFORE any code changes (e.g. the planned SeqAn warm-start swap).
#
# Steps:
#   1. Build ./msa_aff via the Makefile if missing (FORCE_BUILD=1 forces rebuild).
#   2. Probe the MUSCLE warm-start binary STANDALONE: run ./muscle directly with
#      the exact args the solver uses, and report whether it actually executes
#      and produces a valid alignment on this host.
#   3. Run the exact affine solver on Data/bali1/1tgxA.fa with a bounded time
#      limit, and report what MUSCLE did DURING the solve (from the run log).
#   4. Validate that the solver output is a *real* alignment.
#
# PASS  => the toolchain works: solver builds, runs, emits a valid alignment.
#          (MUSCLE working is reported but is NOT required for PASS - its failure
#          is non-fatal; CPLEX still solves the MSA exactly from scratch.)
# FAIL  => something in build/run/output is broken; diagnostics are printed.
#
# Why MUSCLE may not run here: the bundled `muscle` is a Linux x86-64 ELF. A Linux
# ELF cannot execute on macOS at all (Rosetta translates x86-64 *macOS* binaries,
# not Linux ones), so on bare macOS it will fail with an exec-format error. It only
# runs inside a Linux x86-64 environment (container/VM) or if replaced by a native
# build. The solver tolerates that failure: it prints
#     "MUSCLE did not execute! Using previous solution if available"
# and solves from scratch. This script reports MUSCLE's true status either way.
#
# Config (env vars, all optional):
#   TIME_LIMIT   solver -time in seconds          (default 600)
#   INPUT_REL    input .fa relative to repo root   (default Data/bali1/1tgxA.fa)
#   OP           opening gap penalty               (default 12,   msa_aff default)
#   EP           extension gap penalty             (default 2.22, msa_aff default)
#   BASISDIR     CPLEX install root for detection  (default /opt/ibm/ILOG)
#   FORCE_BUILD  set to 1 to `make clean` first    (default unset)
#
# Usage:
#   bash MSAMDD_MWE/smoke_test.sh
#   TIME_LIMIT=120 FORCE_BUILD=1 bash MSAMDD_MWE/smoke_test.sh

set -euo pipefail

# --- locate paths -----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

TIME_LIMIT="${TIME_LIMIT:-600}"
INPUT_REL="${INPUT_REL:-Data/bali1/1tgxA.fa}"
OP="${OP:-12}"
EP="${EP:-2.22}"
BASISDIR="${BASISDIR:-/opt/ibm/ILOG}"
SUBMAT_REL="Data/blosum.ncbi"

INPUT_ABS="${REPO_ROOT}/${INPUT_REL}"
RESULTS_DIR="${SCRIPT_DIR}/results"
OUT_FILE="${RESULTS_DIR}/$(basename "${INPUT_REL%.*}")_aligned.fa"
RUN_LOG="${RESULTS_DIR}/last_run.log"
MUSCLE_OUT="${RESULTS_DIR}/muscle_probe.fa"
MUSCLE_LOG="${RESULTS_DIR}/muscle_probe.log"
VALIDATOR="${SCRIPT_DIR}/validate_alignment.py"

mkdir -p "${RESULTS_DIR}"

# Gap params as the solver passes them to MUSCLE (see src_aff/utility.cpp):
#   -gapopen = -(op + ep) ,  -gapextend = -ep
GAPOPEN="$(awk -v o="${OP}" -v e="${EP}" 'BEGIN{printf "%.6g", -(o+e)}')"
GAPEXTEND="$(awk -v e="${EP}" 'BEGIN{printf "%.6g", -e}')"

# --- helpers ----------------------------------------------------------------
say()  { printf '%s\n' "$*"; }
rule() { printf '%s\n' "------------------------------------------------------------"; }
fail() { say ""; rule; say "RESULT: FAIL  -  $*"; rule; exit 1; }

# Prefer a timeout wrapper if present (GNU coreutils; gtimeout on macOS via brew).
TIMEOUT=""
if command -v timeout  >/dev/null 2>&1; then TIMEOUT="timeout 60"; fi
if command -v gtimeout >/dev/null 2>&1; then TIMEOUT="gtimeout 60"; fi

rule
say "MSAMDD-Ling smoke test (exact affine solver)"
rule
say "repo root   : ${REPO_ROOT}"
say "input       : ${INPUT_REL}"
say "output      : ${OUT_FILE#${REPO_ROOT}/}"
say "time limit  : ${TIME_LIMIT}s   gap penalties: op=${OP} ep=${EP}"
say ""

# --- environment diagnostics ------------------------------------------------
say "[env] host        : $(uname -s) $(uname -m)"
if command -v g++ >/dev/null 2>&1; then
  say "[env] g++         : $(g++ --version | head -1)"
else
  fail "g++ not found on PATH - cannot build the solver."
fi
if [ -d "${BASISDIR}" ]; then
  CPLEX_HIT="$(find "${BASISDIR}" -maxdepth 1 -type d -name 'CPLEX_Studio*' 2>/dev/null | sort -V | tail -1 || true)"
  if [ -n "${CPLEX_HIT}" ]; then
    say "[env] CPLEX        : ${CPLEX_HIT}"
  else
    say "[env] CPLEX        : none under ${BASISDIR} (build will fail without it)"
  fi
else
  say "[env] CPLEX        : ${BASISDIR} absent (build will fail without it)"
fi
[ -f "${INPUT_ABS}" ] || fail "input file not found: ${INPUT_ABS}"
[ -f "${VALIDATOR}" ] || fail "validator not found: ${VALIDATOR}"
command -v python3 >/dev/null 2>&1 || fail "python3 not found - needed for validation."
say ""

# --- step 1: build ----------------------------------------------------------
say "[1/4] build msa_aff"
cd "${REPO_ROOT}"
if [ "${FORCE_BUILD:-}" = "1" ]; then
  say "      FORCE_BUILD=1 -> make clean"
  make clean >/dev/null 2>&1 || true
fi
if [ -x "${REPO_ROOT}/msa_aff" ] && [ -z "${FORCE_BUILD:-}" ]; then
  say "      ./msa_aff already present - skipping build (FORCE_BUILD=1 to rebuild)."
else
  say "      running: make msa_aff"
  if ! make msa_aff >"${RESULTS_DIR}/build.log" 2>&1; then
    say ""
    say "      --- build.log (tail) ---"
    tail -n 25 "${RESULTS_DIR}/build.log" | sed 's/^/      /'
    say "      ------------------------"
    fail "build failed. On Apple Silicon / macOS the Makefile is hard-wired for
       x86-64 Linux + CPLEX at ${BASISDIR} (SYSTEM=x86-64_linux, -no-pie). You
       likely need a Linux x86-64 environment with CPLEX, or to retarget the
       Makefile. Full log: ${RESULTS_DIR}/build.log"
  fi
  say "      build OK -> ${REPO_ROOT}/msa_aff"
fi
[ -x "${REPO_ROOT}/msa_aff" ] || fail "msa_aff binary missing after build step."
say ""

# --- step 2: probe MUSCLE standalone ---------------------------------------
# Definitively answer "is MUSCLE actually running?" by invoking ./muscle exactly
# as the solver does, to our own output file, and inspecting the result. The
# solver deletes its warm-start file at the end, so a standalone probe is the only
# way to inspect what MUSCLE produced.
say "[2/4] probe MUSCLE warm-start binary (standalone)"
MUSCLE_STATUS="unknown"
if [ ! -f "${REPO_ROOT}/muscle" ]; then
  say "      no ./muscle binary present - warm-start unavailable."
  MUSCLE_STATUS="absent"
else
  say "      ./muscle       : $(file "${REPO_ROOT}/muscle" 2>/dev/null | cut -d: -f2- | sed 's/^ //')"
  chmod +x "${REPO_ROOT}/muscle" 2>/dev/null || true
  rm -f "${MUSCLE_OUT}"
  say "      invoking: ./muscle -in ./${INPUT_REL} -out <results>/$(basename "${MUSCLE_OUT}") \\"
  say "                -matrix ./${SUBMAT_REL} -gapopen ${GAPOPEN} -gapextend ${GAPEXTEND} -center 0.0 -quiet"
  set +e
  # shellcheck disable=SC2086
  ${TIMEOUT} "${REPO_ROOT}/muscle" -in "./${INPUT_REL}" -out "${MUSCLE_OUT}" \
    -matrix "./${SUBMAT_REL}" -gapopen "${GAPOPEN}" -gapextend "${GAPEXTEND}" \
    -center 0.0 -quiet </dev/null >"${MUSCLE_LOG}" 2>&1
  MUSCLE_RC=$?
  set -e
  if [ "${MUSCLE_RC}" -eq 126 ] || grep -qiE "exec format error|cannot execute|bad CPU type|no such file" "${MUSCLE_LOG}"; then
    MUSCLE_STATUS="cannot-execute"
    say "      -> MUSCLE CANNOT EXECUTE on this host (rc=${MUSCLE_RC}; exec-format/arch mismatch)."
    say "         This is the expected case for the Linux x86-64 ELF on macOS."
  elif [ "${MUSCLE_RC}" -ne 0 ]; then
    MUSCLE_STATUS="ran-error"
    say "      -> MUSCLE ran but exited non-zero (rc=${MUSCLE_RC}). Log tail:"
    tail -n 8 "${MUSCLE_LOG}" | sed 's/^/         /'
  elif [ ! -s "${MUSCLE_OUT}" ]; then
    MUSCLE_STATUS="no-output"
    say "      -> MUSCLE exited 0 but produced no output file. Suspicious."
  elif python3 "${VALIDATOR}" "${INPUT_ABS}" "${MUSCLE_OUT}" >/dev/null 2>&1; then
    MUSCLE_STATUS="working"
    say "      -> MUSCLE IS RUNNING: produced a valid warm-start alignment."
    say "         (saved for inspection: ${MUSCLE_OUT#${REPO_ROOT}/})"
  else
    MUSCLE_STATUS="bad-output"
    say "      -> MUSCLE exited 0 but its output is not a valid alignment:"
    python3 "${VALIDATOR}" "${INPUT_ABS}" "${MUSCLE_OUT}" 2>&1 | sed 's/^/         /' || true
  fi
fi
say ""

# --- step 3: run the exact affine solver -----------------------------------
say "[3/4] run exact affine solver"
say "      ./msa_aff -in ./${INPUT_REL} -out <results>/$(basename "${OUT_FILE}") -time ${TIME_LIMIT} -op ${OP} -ep ${EP}"
rm -f "${OUT_FILE}"
# Solver calls cin.get() on a load error and would hang; feed it /dev/null.
set +e
./msa_aff -in "./${INPUT_REL}" -out "${OUT_FILE}" -time "${TIME_LIMIT}" \
  -op "${OP}" -ep "${EP}" </dev/null >"${RUN_LOG}" 2>&1
RUN_RC=$?
set -e

if [ "${RUN_RC}" -ne 0 ]; then
  say "      --- last_run.log (tail) ---"
  tail -n 20 "${RUN_LOG}" | sed 's/^/      /'
  say "      ---------------------------"
  fail "solver exited with code ${RUN_RC}. Full log: ${RUN_LOG}"
fi
if ! grep -q "END" "${RUN_LOG}"; then
  say "      warning: solver did not print END; output validity will be checked anyway."
fi

# Report what MUSCLE did DURING this solve, straight from the run log.
say "      MUSCLE during the solve:"
if grep -q "MUSCLE did not execute" "${RUN_LOG}"; then
  say "        - solver reported: \"MUSCLE did not execute\" (warm-start skipped, non-fatal)."
else
  say "        - no \"MUSCLE did not execute\" message: warm-start binary exited 0."
fi
if grep -q "No solution file exists" "${RUN_LOG}"; then
  say "        - no warm-start file was produced (Load_sol found nothing)."
fi
HVAL="$(grep -E "heuristic solution value:" "${RUN_LOG}" | tail -1 | sed 's/.*value: *//' || true)"
if [ -n "${HVAL}" ]; then
  # -IloInfinity (~ -1e+20) means MUSCLE contributed no usable warm start.
  IS_INF="$(awk -v v="${HVAL}" 'BEGIN{print (v+0 <= -1e18) ? 1 : 0}')"
  if [ "${IS_INF}" = "1" ]; then
    say "        - heuristic (warm-start) value: ${HVAL}  => MUSCLE gave NO usable solution."
  else
    say "        - heuristic (warm-start) value: ${HVAL}  => MUSCLE warm-start WAS used."
  fi
fi
say "      run complete (rc=0). Log: ${RUN_LOG#${REPO_ROOT}/}"
say ""

# --- step 4: validate solver output ----------------------------------------
say "[4/4] validate solver output is a real alignment"
[ -s "${OUT_FILE}" ] || fail "output file empty or missing: ${OUT_FILE}"
if ! python3 "${VALIDATOR}" "${INPUT_ABS}" "${OUT_FILE}"; then
  fail "solver output is not a valid alignment (see problems above). File: ${OUT_FILE}"
fi

say ""
rule
say "RESULT: PASS  -  toolchain works; valid alignment written to:"
say "  ${OUT_FILE#${REPO_ROOT}/}"
say "MUSCLE warm-start status (standalone probe): ${MUSCLE_STATUS}"
rule
