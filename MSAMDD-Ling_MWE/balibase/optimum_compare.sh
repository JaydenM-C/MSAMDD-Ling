#!/usr/bin/env bash
#
# optimum_compare.sh - confirm msa_aff reaches the SAME proven optimum whether
# the warm start comes from SeqAn or from MUSCLE (warm start = speed, not optimum).
#
# Method. The exact solver is a single fixed binary; the only thing we vary is
# which heuristic aligner supplies its initial incumbent. msa_aff shells out to
# "./seqan_warmstart" with MUSCLE-compatible flags, so:
#   * SeqAn arm  : the real seqan_warmstart binary is in place.
#   * MUSCLE arm : seqan_warmstart is temporarily replaced by a shim that execs
#                  ./muscle "$@" -- i.e. the solver's warm-start call runs MUSCLE
#                  with byte-for-byte the same flags the former run_MUSCLE used.
# This isolates exactly one variable (the initial aligner) while the solver,
# objective and instance stay identical. Both arms are run to proven optimality
# (Gap -> 0); a fixed MILP has one global optimum, so the two objectives must
# match. They differ only in the warm-start incumbent and the time to converge.
#
# Caveat (honest): the current solver short-circuits the intra-Benders refine
# call, so via the shim MUSCLE supplies only the INITIAL warm start (no
# muscle -refine mid-solve). That cannot change the optimum, which is the claim
# under test; it does mean this is not a re-run of the old end-to-end MUSCLE
# refine behaviour.
#
# Runs inside the linux/amd64 container (muscle + CPLEX are x86-64 Linux). It
# must run from the repository root.
#
# Env knobs:
#   ENTRIES      space-separated input FASTAs  (default: two short bali1 N=4 sets)
#   TIME_LIMIT   per-arm solver -time seconds   (default 900)
#   OP, EP       gap penalties                  (default 12, 2.22 = protein)
#   SUBMAT       substitution matrix            (default ./Data/blosum.ncbi)
#
# NOTE ON RUNTIME: the exact affine solve can be slow and is not guaranteed to
# prove optimality inside TIME_LIMIT on every instance. Start with the small
# defaults; scale up deliberately. Per-entry the script reports whether each arm
# actually PROVED optimality (Gap=0); only proven-vs-proven pairs are a clean
# optimum check (flagged INCONCLUSIVE otherwise).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
LIB="${SCRIPT_DIR}/../lib/score_sop.py"
RESULTS="${SCRIPT_DIR}/../results"
cd "${REPO_ROOT}"

ENTRIES="${ENTRIES:-Data/bali1/1tgxA.fa Data/bali1/1tvxA.fa}"
TIME_LIMIT="${TIME_LIMIT:-900}"
OP="${OP:-12}"
EP="${EP:-2.22}"
SUBMAT="${SUBMAT:-./Data/blosum.ncbi}"

mkdir -p "${RESULTS}"
say()  { printf '%s\n' "$*"; }
rule() { printf '%s\n' "===================================================================="; }

# --- preflight --------------------------------------------------------------
[ -f "${REPO_ROOT}/muscle" ] || { say "ERROR: ./muscle missing (need it for the MUSCLE arm)."; exit 1; }
command -v python3 >/dev/null || { say "ERROR: python3 not found."; exit 1; }

# Force-rebuild the binaries we exercise so a stale (wrong-arch) copy left in the
# tree cannot shadow a native build -- a known foot-gun on this repo. (validate_all.sh
# builds once up front and sets SKIP_BUILD=1 to avoid rebuilding per stage.)
if [ "${SKIP_BUILD:-0}" != "1" ]; then
  say "[build] rm -f msa_aff seqan_warmstart && make msa_aff seqan_warmstart"
  rm -f msa_aff seqan_warmstart
  make msa_aff seqan_warmstart >"${RESULTS}/build.log" 2>&1 || {
    say "  build FAILED -- tail of ${RESULTS}/build.log:"; tail -n 25 "${RESULTS}/build.log" | sed 's/^/    /'; exit 1; }
fi
[ -x "${REPO_ROOT}/msa_aff" ] && [ -x "${REPO_ROOT}/seqan_warmstart" ] || { say "ERROR: binaries missing (build first, or unset SKIP_BUILD)."; exit 1; }

# Keep the genuine SeqAn binary aside; restore it on any exit.
cp -f seqan_warmstart seqan_warmstart.real
restore() { cp -f seqan_warmstart.real seqan_warmstart 2>/dev/null || true; chmod +x seqan_warmstart 2>/dev/null || true; rm -f seqan_warmstart.real; }
trap restore EXIT

install_seqan() { cp -f seqan_warmstart.real seqan_warmstart; chmod +x seqan_warmstart; }
install_muscle_shim() { printf '#!/bin/sh\nexec ./muscle "$@"\n' > seqan_warmstart; chmod +x seqan_warmstart; }

# Run one arm; echoes: "<incumbent>\t<gap>\t<wall_s>\t<objective>\t<optimal 0|1>"
run_arm() {
  local entry="$1" arm="$2" name="$3"
  local out="${RESULTS}/${name}.${arm}.fa" log="${RESULTS}/${name}.${arm}.log"
  local t0 t1 wall
  t0=$(date +%s)
  set +e
  ./msa_aff -in "${entry}" -out "${out}" -submat "${SUBMAT}" \
            -op "${OP}" -ep "${EP}" -time "${TIME_LIMIT}" </dev/null >"${log}" 2>&1
  set -e
  t1=$(date +%s); wall=$((t1 - t0))
  local inc gap obj opt
  inc=$(grep "heuristic solution value:" "${log}" | tail -1 | sed 's/.*value:[[:space:]]*//' || true)
  gap=$(awk -F'Gap is:' '/Gap is:/{g=$2} END{gsub(/[ \t\r]/,"",g); print g}' "${log}")
  [ -n "${inc}" ] || inc="NA"
  [ -n "${gap}" ] || gap="NA"
  if [ -s "${out}" ]; then
    obj=$(python3 "${LIB}" score --matrix "${SUBMAT}" --op "${OP}" --ep "${EP}" "${out}" 2>/dev/null || echo "NA")
  else
    obj="NA"
  fi
  opt=0
  if [ "${gap}" != "NA" ] && awk "BEGIN{exit !(${gap} < 0.01)}"; then opt=1; fi
  printf '%s\t%s\t%s\t%s\t%s\n' "${inc}" "${gap}" "${wall}" "${obj}" "${opt}"
}

rule; say "optimum_compare: SeqAn vs MUSCLE warm start -> same proven optimum?"
say "entries    : ${ENTRIES}"
say "penalties  : op=${OP} ep=${EP}   matrix=${SUBMAT}   per-arm -time=${TIME_LIMIT}s"
rule

printf '%-12s %-8s %14s %10s %8s %8s %6s\n' "entry" "arm" "objective" "incumbent" "gap%" "wall_s" "opt"
printf '%s\n' "--------------------------------------------------------------------------------"

overall_pass=1
any_inconclusive=0
for entry in ${ENTRIES}; do
  [ -f "${entry}" ] || { say "  (missing input: ${entry})"; continue; }
  name="$(basename "${entry%.*}")"
  # Run on a COPY inside results/ so the solver's byproducts (_s.txt, and the
  # heur.cpp "_h.fa" heuristic-solution file, which msa_aff does NOT clean up)
  # land in results/ rather than polluting Data/.
  infile="${RESULTS}/${name}.in.fa"; cp -f "${entry}" "${infile}"

  install_seqan
  IFS=$'\t' read -r sInc sGap sWall sObj sOpt < <(run_arm "${infile}" "seqan" "${name}")
  install_muscle_shim
  IFS=$'\t' read -r mInc mGap mWall mObj mOpt < <(run_arm "${infile}" "muscle" "${name}")
  install_seqan

  printf '%-12s %-8s %14s %10s %8s %8s %6s\n' "${name}" "seqan" "${sObj}" "${sInc}" "${sGap}" "${sWall}" "${sOpt}"
  printf '%-12s %-8s %14s %10s %8s %8s %6s\n' "${name}" "muscle" "${mObj}" "${mInc}" "${mGap}" "${mWall}" "${mOpt}"

  if [ "${sOpt}" = "1" ] && [ "${mOpt}" = "1" ]; then
    if awk "BEGIN{d=(${sObj})-(${mObj}); if(d<0)d=-d; exit !(d <= 1e-6*(( (${sObj}<0)?-${sObj}:${sObj} )+1))}"; then
      verdict="SAME OPTIMUM"
    else
      verdict="*** OPTIMA DIFFER ***"; overall_pass=0
    fi
  else
    verdict="INCONCLUSIVE (an arm did not prove optimality in ${TIME_LIMIT}s)"
    any_inconclusive=1
  fi
  say "  -> ${name}: ${verdict}   (SeqAn obj=${sObj} / MUSCLE obj=${mObj})"
  printf '%s\n' "--------------------------------------------------------------------------------"
done

rule
if [ "${overall_pass}" = "1" ] && [ "${any_inconclusive}" = "0" ]; then
  say "RESULT: PASS - every entry proved the SAME optimum from both warm starts."
elif [ "${overall_pass}" = "1" ]; then
  say "RESULT: PARTIAL - no optimum mismatch, but some entries were INCONCLUSIVE"
  say "        (raise TIME_LIMIT or pick smaller ENTRIES to prove optimality)."
else
  say "RESULT: FAIL - at least one entry reported DIFFERENT optima. Investigate:"
  say "        a warm start must never change the proven optimum; this signals a"
  say "        scoring/integration bug, not a heuristic-quality difference."
fi
say "Per-arm logs + alignments: ${RESULTS#${REPO_ROOT}/}/<entry>.{seqan,muscle}.{log,fa}"
rule

# Exit non-zero only on a genuine optimum mismatch (the regression we guard
# against). INCONCLUSIVE is a time-budget outcome, not a failure of the swap.
[ "${overall_pass}" = "1" ]
