#!/usr/bin/env bash
#
# validate_all.sh - Phase 4 validation driver. Runs INSIDE the linux/amd64
# container (muscle + CPLEX are x86-64 Linux), from the repository root. Builds
# the binaries once, then runs the selected validation stage(s).
#
# Stages (pass as the first argument, or via STAGE=...; default "all"):
#   incumbent   SeqAn-vs-MUSCLE warm-start incumbent quality (no CPLEX solve)
#   optimum     end-to-end: same proven optimum from either warm start (CPLEX)
#   linguistic  end-to-end exact MSA over a phoneme alphabet + custom matrix
#   all         all three, in that order
#
# Env knobs are forwarded to the stage scripts: TIME_LIMIT, ENTRIES, OP, EP,
# SUBMAT (see each script's header).

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS="${SCRIPT_DIR}/results"
cd "${REPO_ROOT}"

STAGE="${1:-${STAGE:-all}}"
mkdir -p "${RESULTS}"
say()  { printf '%s\n' "$*"; }
rule() { printf '%s\n' "####################################################################"; }

rule
say "Phase 4 validation  |  stage: ${STAGE}  |  host: $(uname -s) $(uname -m)"
rule

# --- build once; stages run with SKIP_BUILD=1 -------------------------------
say "[build] rm -f msa_aff seqan_warmstart && make msa_aff seqan_warmstart"
rm -f msa_aff seqan_warmstart
if ! make msa_aff seqan_warmstart >"${RESULTS}/build.log" 2>&1; then
  say "  BUILD FAILED. tail of ${RESULTS}/build.log:"
  tail -n 30 "${RESULTS}/build.log" | sed 's/^/    /'
  say "  (Need CPLEX for msa_aff and the initialised includes/seqan submodule for"
  say "   seqan_warmstart. Run inside the amd64 image per MSAMDD-Ling_MWE/README.md.)"
  exit 1
fi
export SKIP_BUILD=1
say "  build OK: $(file msa_aff | cut -d, -f1-2 | sed 's/.*: //'), seqan_warmstart present."
say ""

rc_inc=0; rc_opt=0; rc_ling=0
run_incumbent() {
  rule; say "STAGE incumbent: warm-start incumbent quality (SeqAn vs MUSCLE)"; rule
  python3 "${SCRIPT_DIR}/balibase/incumbent_bench.py"; rc_inc=$?
}
run_optimum() {
  rule; say "STAGE optimum: same proven optimum from either warm start"; rule
  bash "${SCRIPT_DIR}/balibase/optimum_compare.sh"; rc_opt=$?
}
run_linguistic() {
  rule; say "STAGE linguistic: exact MSA over a phoneme alphabet + custom matrix"; rule
  bash "${SCRIPT_DIR}/linguistic/linguistic_mwe.sh"; rc_ling=$?
}

case "${STAGE}" in
  incumbent)  run_incumbent ;;
  optimum)    run_optimum ;;
  linguistic) run_linguistic ;;
  all)        run_incumbent; say ""; run_optimum; say ""; run_linguistic ;;
  *) say "unknown stage '${STAGE}' (use: incumbent | optimum | linguistic | all)"; exit 2 ;;
esac

say ""
rule
say "Phase 4 validation summary"
[ "${STAGE}" = "all" ] || [ "${STAGE}" = "incumbent" ]  && say "  incumbent  : exit ${rc_inc}"
[ "${STAGE}" = "all" ] || [ "${STAGE}" = "optimum" ]    && say "  optimum    : exit ${rc_opt}"
[ "${STAGE}" = "all" ] || [ "${STAGE}" = "linguistic" ] && say "  linguistic : exit ${rc_ling}"
say "  artifacts  : MSAMDD-Ling_MWE/results/"
rule
# Non-zero overall exit if any run stage failed.
exit $(( rc_inc | rc_opt | rc_ling ))
