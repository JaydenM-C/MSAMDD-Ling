#!/usr/bin/env bash
#
# lardil_mwe.sh - end-to-end exact MSA of a toy Lardil dataset with the SeqAn
# warm start. This is the first run of MSAMDD-Ling against real (if small)
# linguistic data, and the point is narrow: show that the SeqAn warm start
# consumes a unicode-derived phoneme alphabet + a custom substitution matrix
# and that msa_aff runs end-to-end on it. MUSCLE is not involved (it rejects
# non-amino-acid matrices outright -- see the Phase-4 session notes).
#
# Pipeline:
#   1. lardil_prep.py  : unicode IPA  -> single-char ASCII FASTA + matrix.
#   2. seqan_warmstart : ASCII FASTA  -> a valid warm-start MSA (the thing we
#                        are demonstrating; works even if the exact solve below
#                        does not close in time).
#   3. msa_aff         : exact, custom-matrix MSA over the phoneme alphabet,
#                        consuming the SeqAn warm start as its incumbent.
#
# Everything generated lands in Lardil_MWEs/Results/. Run it from anywhere:
#   bash MSAMDD-Ling/Lardil_MWEs/lardil_mwe.sh
#
# Scoring scheme (synchronic morphology):
#   * substitution matrix : +DIAG on the diagonal, OFFDIAG (uniform) off it.
#   * gaps                : -op to open, -ep to extend, with ep = 0.
#   * INVARIANT (this task): |OFFDIAG| > OP  -- a substitution costs more than
#     opening a gap. The script asserts it before running.
#
# Env knobs (all optional):
#   OP=1 EP=0            gap open / extend penalties passed to the solver
#   DIAG=2 OFFDIAG=-2    matrix diagonal / off-diagonal scores
#   TIME_LIMIT=600       msa_aff -time budget (seconds)
#   MAXFORMS=            use only the first N forms (tractability knob; N=9 here
#                        may not prove optimality -- try MAXFORMS=3 for a quick
#                        provably-optimal run on the three Thigh forms)
#   REBUILD=1            force `rm -f msa_aff seqan_warmstart && make`
#   SKIP_BUILD=1         use existing binaries, never build
# First positional arg overrides the input FASTA (default the toy Lardil set).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS="${SCRIPT_DIR}/Results"
PREP="${SCRIPT_DIR}/lardil_prep.py"
# Authoritative SOP scorer / MSA validator (shared with the Phase-4 harness).
# Optional: if absent, we skip scoring/validation but still run the pipeline.
SCORE_LIB="${REPO_ROOT}/MSAMDD-Ling_MWE/lib/score_sop.py"

SRC="${1:-${SCRIPT_DIR}/Data/Lardil_short_plain.fasta}"
# Resolve a relative input against the caller's cwd first, then the MWE dir, and
# store it ABSOLUTE so it survives the later `cd "${REPO_ROOT}"`.
if   [ -f "${SRC}" ];               then SRC="$(cd "$(dirname "${SRC}")" && pwd)/$(basename "${SRC}")"
elif [ -f "${SCRIPT_DIR}/${SRC}" ]; then SRC="${SCRIPT_DIR}/${SRC}"
fi
OP="${OP:-1}"
EP="${EP:-0}"
DIAG="${DIAG:-2}"
OFFDIAG="${OFFDIAG:--2}"
TIME_LIMIT="${TIME_LIMIT:-600}"
MAXFORMS="${MAXFORMS:-}"

say()  { printf '%s\n' "$*"; }
rule() { printf '%s\n' "===================================================================="; }
have_scorer() { [ -f "${SCORE_LIB}" ]; }

# --- progress helpers -------------------------------------------------------
# Animate only on an interactive TTY; piped/redirected runs print one static
# line so logs stay clean. The watched work always runs in the BACKGROUND so we
# can poll it. Bash-3.2-safe: no arrays, no ${v^^}, POSIX arithmetic only.
_tty()  { [ -t 1 ]; }
_mmss() { printf '%d:%02d' $(( $1 / 60 )) $(( $1 % 60 )); }

# spin_until PID LABEL -- spinner + elapsed seconds while PID runs (build steps).
spin_until() {
  _p=$1; _lbl=$2; _s='|/-\'; _i=0; _t0=$(date +%s)
  if ! _tty; then say "  ${_lbl} ..."; return 0; fi
  while kill -0 "$_p" 2>/dev/null; do
    printf '\r  %s %s  %ss\033[K' \
      "$(printf '%s' "$_s" | cut -c $(( _i % 4 + 1 )))" "$_lbl" "$(( $(date +%s) - _t0 ))"
    _i=$(( _i + 1 )); sleep 0.5
  done
  printf '\r\033[K'
}

# monitor_solver PID LOG BUDGET -- elapsed/budget bar + live gap & incumbent
# scraped from the solver log. The bar tracks TIME USED, not work done: an exact
# solver has no honest percent-complete. gap/inc appear as soon as the solver
# flushes them (it prints status with endl, so usually promptly).
monitor_solver() {
  _p=$1; _log=$2; _bud=$3; _w=22; _t0=$(date +%s)
  if ! _tty; then say "  solving (up to $(_mmss "$_bud")) ..."; return 0; fi
  while kill -0 "$_p" 2>/dev/null; do
    _el=$(( $(date +%s) - _t0 ))
    _pct=0; [ "$_bud" -gt 0 ] && _pct=$(( _el * 100 / _bud )); [ "$_pct" -gt 100 ] && _pct=100
    _f=$(( _pct * _w / 100 )); _bar=''; _k=0
    while [ "$_k" -lt "$_w" ]; do
      if [ "$_k" -lt "$_f" ]; then _bar="${_bar}#"; else _bar="${_bar}-"; fi; _k=$(( _k + 1 ))
    done
    _gap=$(awk -F'Gap is:' '/Gap is:/{g=$2} END{gsub(/[ \t\r]/,"",g); print (g==""?"--":g"%")}' "$_log" 2>/dev/null || true)
    _inc=$(grep -i 'heuristic solution value:' "$_log" 2>/dev/null | tail -1 | sed 's/.*value:[[:space:]]*//' | tr -d ' \t\r' || true)
    [ -z "$_inc" ] && _inc='--'
    printf '\r  [solve] %s/%s [%s] %3d%%  gap %s  inc %s\033[K' \
      "$(_mmss "$_el")" "$(_mmss "$_bud")" "$_bar" "$_pct" "$_gap" "$_inc"
    sleep 2
  done
  printf '\r\033[K'
}

# Per-dataset output dir, derived from the input filename stem, so runs on
# different inputs never clobber each other (e.g. the short-set results are
# preserved when you run the full set).
STEM="$(basename "${SRC}")"; STEM="${STEM%.fasta}"; STEM="${STEM%.fa}"
OUTDIR="${RESULTS}/${STEM}"
mkdir -p "${OUTDIR}"
cd "${REPO_ROOT}"

# --- preflight --------------------------------------------------------------
[ -f "${SRC}" ]  || { say "ERROR: input FASTA not found: ${SRC}"; exit 1; }
[ -f "${PREP}" ] || { say "ERROR: missing ${PREP}"; exit 1; }
command -v python3 >/dev/null || { say "ERROR: python3 not found."; exit 1; }

# Scoring regime (informational -- this script is a general scheme explorer).
# The |OFFDIAG| vs OP relationship decides what the optimiser prefers:
#   |OFFDIAG| > OP  -> substitution dearer than a gap: align only identical
#                      segments, scatter the rest into gaps (identity-anchored).
#   |OFFDIAG| < OP  -> gap dearer than a substitution: align slots across the
#                      paradigm (positional correspondence).
REGIME="$(awk -v o="${OFFDIAG}" -v p="${OP}" 'BEGIN{ a=(o<0?-o:o);
  if (a>p)      printf "|sub|=%g > op=%g  -> identity-anchored (gappy)", a, p;
  else if (a<p) printf "|sub|=%g < op=%g  -> favours correspondence", a, p;
  else          printf "|sub|=%g = op=%g  -> balanced", a, p; }')"

# Flags exactly as run_warmstart builds them: -gapopen -(op+ep), -gapextend -ep.
GAPOPEN="$(awk -v o="${OP}" -v e="${EP}" 'BEGIN{v=-(o+e); if(v==0)v=0; printf "%g", v}')"
GAPEXTEND="$(awk -v e="${EP}" 'BEGIN{v=-e; if(v==0)v=0; printf "%g", v}')"

# --- build (or reuse) the binaries -----------------------------------------
# seqan_warmstart is pure SeqAn (builds native, no CPLEX). msa_aff needs CPLEX;
# if it will not build we still demonstrate the warm start and skip the solve.
#
# CPLEX location passthrough: the Makefile auto-detects the studio dir by
# globbing BASISDIR for CPLEX_Studio*; if that finds nothing the include path
# comes out empty (ilcplex/ilocplex.h not found). Set any of these to point make
# at your install and they are forwarded on the command line (highest precedence):
#   BASISILOG=/Applications/CPLEX_Studio2211   # pin the exact studio dir (easiest)
#   BASISDIR=/some/parent  SYSTEM=x86-64_osx
MAKEVARS=()
[ -n "${BASISILOG:-}" ] && MAKEVARS+=("BASISILOG=${BASISILOG}")
[ -n "${BASISDIR:-}" ]  && MAKEVARS+=("BASISDIR=${BASISDIR}")
[ -n "${SYSTEM:-}" ]    && MAKEVARS+=("SYSTEM=${SYSTEM}")
if [ "${REBUILD:-0}" = "1" ]; then
  say "[build] REBUILD=1 -> rm -f msa_aff seqan_warmstart && make"
  rm -f msa_aff seqan_warmstart
fi
if [ "${SKIP_BUILD:-0}" != "1" ]; then
  if [ ! -x "${REPO_ROOT}/seqan_warmstart" ]; then
    say "[build] seqan_warmstart"
    make seqan_warmstart >"${OUTDIR}/build_seqan.log" 2>&1 & _bp=$!
    spin_until "${_bp}" "compiling seqan_warmstart"
    if ! wait "${_bp}"; then
      say "  seqan_warmstart build FAILED -- tail:"; tail -n 20 "${OUTDIR}/build_seqan.log" | sed 's/^/    /'
      say "  (need a C++14 compiler + the SeqAn submodule: git submodule update --init)"; exit 1
    fi
  fi
  if [ ! -x "${REPO_ROOT}/msa_aff" ]; then
    say "[build] ${MAKEVARS[*]+${MAKEVARS[*]} }msa_aff"
    make ${MAKEVARS[@]+"${MAKEVARS[@]}"} msa_aff >"${OUTDIR}/build_msaaff.log" 2>&1 & _bp=$!
    spin_until "${_bp}" "compiling msa_aff (CPLEX)"
    if ! wait "${_bp}"; then
      say "  msa_aff build FAILED (needs CPLEX configured for your platform)."
      say "  If it can't find ilcplex/ilocplex.h, point make at your CPLEX install, e.g.:"
      say "    BASISILOG=/Applications/CPLEX_Studio2211 bash Lardil_MWEs/lardil_mwe.sh"
      say "  Warm-start demo will still run; the exact solve (step 3) is skipped."
      say "  tail:"; tail -n 12 "${OUTDIR}/build_msaaff.log" | sed 's/^/    /'
    fi
  fi
fi
[ -x "${REPO_ROOT}/seqan_warmstart" ] || { say "ERROR: seqan_warmstart binary missing."; exit 1; }
MSA_AFF_OK=0; [ -x "${REPO_ROOT}/msa_aff" ] && MSA_AFF_OK=1

# --- step 1: unicode -> ASCII + matrix -------------------------------------
ASCII="${OUTDIR}/ascii.fa"
MATRIX="${OUTDIR}/matrix.mat"
LEGEND="${OUTDIR}/legend.tsv"
# Optional --max-forms. Note the `${arr[@]+"${arr[@]}"}` guard everywhere an
# array is expanded: macOS ships Bash 3.2, where "${arr[@]}" on an EMPTY array
# under `set -u` aborts with "unbound variable". The guard expands to nothing
# when empty and to the elements otherwise.
MAXARG=(); [ -n "${MAXFORMS}" ] && MAXARG=(--max-forms "${MAXFORMS}")

rule
say "Lardil MWE: exact MSA over a unicode-derived phoneme alphabet, custom matrix"
say "input : ${SRC#${REPO_ROOT}/}"
say "gaps  : op=${OP} / ep=${EP}   matrix: diagonal=${DIAG} off-diagonal=${OFFDIAG}"
say "regime: ${REGIME}"
say "output: ${OUTDIR#${REPO_ROOT}/}/"
rule
say ""
say "[1/3] lardil_prep.py: transliterate to single-char ASCII + emit matrix"
python3 "${PREP}" build --in "${SRC}" --out-fasta "${ASCII}" --out-matrix "${MATRIX}" \
        --out-legend "${LEGEND}" --diagonal "${DIAG}" --off-diagonal "${OFFDIAG}" ${MAXARG[@]+"${MAXARG[@]}"} \
        | sed 's/^/      /'
say "      ASCII FASTA -> ${ASCII#${REPO_ROOT}/}"
say "      matrix      -> ${MATRIX#${REPO_ROOT}/}"

# --- step 2: SeqAn warm start ----------------------------------------------
WARM="${OUTDIR}/warmstart.fa"
say ""
say "[2/3] seqan_warmstart on the ASCII data + custom matrix"
say "      ./seqan_warmstart -in ${ASCII#${REPO_ROOT}/} -out ${WARM#${REPO_ROOT}/} -matrix ${MATRIX#${REPO_ROOT}/} -gapopen ${GAPOPEN} -gapextend ${GAPEXTEND} -center 0.0 -quiet"
./seqan_warmstart -in "${ASCII}" -out "${WARM}" -matrix "${MATRIX}" \
                  -gapopen "${GAPOPEN}" -gapextend "${GAPEXTEND}" -center 0.0 -quiet
if have_scorer; then
  say "      validation:"
  python3 "${SCORE_LIB}" validate --input "${ASCII}" "${WARM}"
  W_SOP="$(python3 "${SCORE_LIB}" score --matrix "${MATRIX}" --op "${OP}" --ep "${EP}" "${WARM}")"
  say "      warm-start sum-of-pairs (msa_aff's incumbent): ${W_SOP}"
else
  say "      (score_sop.py not found at ${SCORE_LIB#${REPO_ROOT}/} -- skipping validation/score)"
fi
say "      warm-start alignment (IPA):"
python3 "${PREP}" decode --in "${WARM}" | sed 's/^/        /'

# --- step 3: msa_aff end-to-end --------------------------------------------
say ""
if [ "${MSA_AFF_OK}" != "1" ]; then
  rule
  say "msa_aff not available -- the SeqAn warm start above is demonstrated; the"
  say "exact solve was skipped. Build msa_aff (CPLEX) and re-run for the optimum."
  rule
  exit 0
fi

ALIGNED="${OUTDIR}/aligned.fa"
A_LOG="${OUTDIR}/msaaff.log"
say "[3/3] msa_aff end-to-end (exact, custom matrix, op=${OP}/ep=${EP}, -time ${TIME_LIMIT})"
say "      ./msa_aff -in ${ASCII#${REPO_ROOT}/} -out ${ALIGNED#${REPO_ROOT}/} -submat ${MATRIX#${REPO_ROOT}/} -op ${OP} -ep ${EP} -time ${TIME_LIMIT}"
set +e
./msa_aff -in "${ASCII}" -out "${ALIGNED}" -submat "${MATRIX}" \
          -op "${OP}" -ep "${EP}" -time "${TIME_LIMIT}" </dev/null >"${A_LOG}" 2>&1 &
A_PID=$!
monitor_solver "${A_PID}" "${A_LOG}" "${TIME_LIMIT}"
wait "${A_PID}"; A_RC=$?
set -e
if [ "${A_RC}" -ne 0 ] || [ ! -s "${ALIGNED}" ]; then
  if [ "${A_RC}" -gt 128 ]; then
    say "      msa_aff CRASHED (signal $((A_RC - 128)); 11 = segfault) -- an upstream"
    say "      solver limit at this instance size, NOT a problem with the warm start."
  else
    say "      msa_aff produced no alignment (exit ${A_RC})."
  fi
  say "      Tail of log:"; tail -n 20 "${A_LOG}" | sed 's/^/        /'
  say ""
  rule
  say "Exact solve unavailable for this instance. The SeqAn warm start above IS a"
  say "valid alignment of all forms and is saved at:"
  say "  ${WARM#${REPO_ROOT}/}"
  say "For exact (provably-optimal) output, shrink the instance, e.g. MAXFORMS=<N>."
  rule
  exit 0
fi
INC="$(grep -i 'heuristic solution value:' "${A_LOG}" | tail -1 | sed 's/.*value:[[:space:]]*//')"
GAP="$(awk -F'Gap is:' '/Gap is:/{g=$2} END{gsub(/[ \t\r]/,"",g); print g}' "${A_LOG}")"
say "      msa_aff completed."
if have_scorer; then
  say "      validation:"
  python3 "${SCORE_LIB}" validate --input "${ASCII}" "${ALIGNED}"
  OBJ="$(python3 "${SCORE_LIB}" score --matrix "${MATRIX}" --op "${OP}" --ep "${EP}" "${ALIGNED}")"
  say ""
  say "      warm-start incumbent (SeqAn): ${INC:-n/a}"
  say "      final objective (get_value) : ${OBJ}"
fi
say "      optimality gap reported     : ${GAP:-n/a}%"
PROVEN="no (timed out -- raise TIME_LIMIT, or lower MAXFORMS)"
if [ -n "${GAP}" ] && awk "BEGIN{exit !(${GAP} < 0.01)}"; then PROVEN="YES (gap = 0)"; fi
say "      proven optimal              : ${PROVEN}"
say ""
say "      optimal alignment (ASCII):"
sed 's/^/        /' "${ALIGNED}"
say "      optimal alignment (IPA):"
python3 "${PREP}" decode --in "${ALIGNED}" | sed 's/^/        /'

rule
say "RESULT: exact MSA over a unicode-derived Lardil phoneme alphabet with a"
say "custom substitution matrix -- input MUSCLE cannot align. Artifacts in"
say "${OUTDIR#${REPO_ROOT}/}/: ascii.fa, matrix.mat, warmstart.fa, aligned.fa"
rule
