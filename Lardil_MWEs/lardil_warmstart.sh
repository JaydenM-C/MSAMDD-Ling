#!/usr/bin/env bash
#
# lardil_warmstart.sh - SeqAn warm start ONLY (approximate MSA), no exact solve.
#
# The full pipeline (lardil_mwe.sh) ends in msa_aff, the exact CPLEX solver,
# which segfaults inside its Benders cut-separation above 12 forms -- so the full
# 33-form Lardil paradigm cannot be solved as one instance today. This script
# stops one step earlier and emits just the SeqAn warm start: a valid, heuristic
# alignment of ALL forms. It never builds or runs msa_aff and needs no CPLEX, so
# the >12-form solver crash structurally cannot occur here.
#
# Pipeline (steps 1-2 of the MWE; step 3 deliberately omitted):
#   1. lardil_prep.py  : unicode IPA -> single-char ASCII FASTA + matrix.
#   2. seqan_warmstart : ASCII FASTA -> a valid progressive-MSA alignment of all
#                        forms (the deliverable), then decoded back to IPA.
#
# This is the same warm start msa_aff would consume as its incumbent; the
# downstream solver re-scores it with its own objective, so the alignment only
# needs to be *valid* (equal-length rows that de-gap to the inputs), which the
# bundled validator checks. It is NOT provably optimal -- that is what the exact
# solve buys, at the cost of the form ceiling.
#
# Output lands in Results/<dataset>_warmstart/ -- a separate namespace from the
# MWE's Results/<dataset>/ so a full-set (33-form) warm start never clobbers, or
# is confused with, a capped (<=12-form) exact result sitting next to it.
#
# Run from anywhere (paths resolve relative to the script):
#   bash Lardil_MWEs/lardil_warmstart.sh                              # full 33-form set
#   bash Lardil_MWEs/lardil_warmstart.sh Lardil_MWEs/Data/Lardil_short_plain.fasta
#
# Scoring scheme (mirrors lardil_mwe.sh; synchronic morphology default):
#   * matrix : +DIAG on the diagonal, OFFDIAG (uniform) off it.
#   * gaps   : -op to open, -ep to extend (passed to SeqAn as -gapopen/-gapextend
#              exactly as run_warmstart builds them: -(op+ep) and -ep).
#   * |OFFDIAG| > OP (default -2 vs 1) -> identity-anchored: align only identical
#     segments, gap the rest apart. Set OFFDIAG=-1 OP=2 for positional/diachronic.
#
# Env knobs (all optional):
#   OP=1 EP=0            gap open / extend penalties
#   DIAG=2 OFFDIAG=-2    matrix diagonal / off-diagonal scores
#   MAXFORMS=            use only the first N forms (default: all)
#   RESULTS_DIR=         output root (default: Lardil_MWEs/Results)
#   REBUILD=1            force `rm -f seqan_warmstart && make seqan_warmstart`
#   SKIP_BUILD=1         use the existing binary, never build
#   CXX=g++-13           override the compiler passed to make
# First positional arg overrides the input FASTA (default the full Lardil set).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RESULTS_DIR="${RESULTS_DIR:-${SCRIPT_DIR}/Results}"
PREP="${SCRIPT_DIR}/lardil_prep.py"
# Authoritative SOP scorer / MSA validator (shared with the Phase-4 harness).
# Optional: if absent we skip scoring/validation but still produce the alignment.
SCORE_LIB="${REPO_ROOT}/MSAMDD-Ling_MWE/lib/score_sop.py"

# Default input is the FULL 33-form set -- the whole point of this script is the
# approximate alignment of the complete paradigm that the exact solver can't take.
SRC="${1:-${SCRIPT_DIR}/Data/Lardil_plain.fasta}"
# Resolve a relative input against the caller's cwd first, then the MWE dir, and
# store it ABSOLUTE so it survives the later `cd "${REPO_ROOT}"`.
if   [ -f "${SRC}" ];               then SRC="$(cd "$(dirname "${SRC}")" && pwd)/$(basename "${SRC}")"
elif [ -f "${SCRIPT_DIR}/${SRC}" ]; then SRC="${SCRIPT_DIR}/${SRC}"
fi

OP="${OP:-1}"
EP="${EP:-0}"
DIAG="${DIAG:-2}"
OFFDIAG="${OFFDIAG:--2}"
MAXFORMS="${MAXFORMS:-}"

say()  { printf '%s\n' "$*"; }
rule() { printf '%s\n' "===================================================================="; }
have_scorer() { [ -f "${SCORE_LIB}" ]; }

# Per-dataset output dir, namespaced with a _warmstart suffix so it never shares a
# directory with the MWE's exact-solve artifacts for the same input.
STEM="$(basename "${SRC}")"; STEM="${STEM%.fasta}"; STEM="${STEM%.fa}"
OUTDIR="${RESULTS_DIR}/${STEM}_warmstart"
mkdir -p "${OUTDIR}"
cd "${REPO_ROOT}"

# --- preflight --------------------------------------------------------------
[ -f "${SRC}" ]  || { say "ERROR: input FASTA not found: ${SRC}"; exit 1; }
[ -f "${PREP}" ] || { say "ERROR: missing ${PREP}"; exit 1; }
command -v python3 >/dev/null || { say "ERROR: python3 not found."; exit 1; }

# Gap flags exactly as run_warmstart builds them: -gapopen -(op+ep), -gapextend -ep.
GAPOPEN="$(awk -v o="${OP}" -v e="${EP}" 'BEGIN{v=-(o+e); if(v==0)v=0; printf "%g", v}')"
GAPEXTEND="$(awk -v e="${EP}" 'BEGIN{v=-e; if(v==0)v=0; printf "%g", v}')"

# --- build (or reuse) seqan_warmstart --------------------------------------
# Pure SeqAn: a C++14 compiler + the vendored SeqAn v2 headers, no CPLEX. We only
# ever `make seqan_warmstart` -- msa_aff is never a target here, so a missing or
# broken CPLEX install is irrelevant and there is no exact solver to crash.
if [ "${REBUILD:-0}" = "1" ]; then
  say "[build] REBUILD=1 -> rm -f seqan_warmstart && make seqan_warmstart"
  rm -f seqan_warmstart
fi
if [ "${SKIP_BUILD:-0}" != "1" ] && [ ! -x "${REPO_ROOT}/seqan_warmstart" ]; then
  say "[build] make seqan_warmstart (vendored SeqAn v2; no CPLEX)"
  if ! make seqan_warmstart >"${OUTDIR}/build_seqan.log" 2>&1; then
    say "  build FAILED -- tail:"; tail -n 20 "${OUTDIR}/build_seqan.log" | sed 's/^/    /'
    say "  (need a C++14 compiler + the SeqAn submodule: git submodule update --init)"
    exit 1
  fi
fi
[ -x "${REPO_ROOT}/seqan_warmstart" ] || { say "ERROR: seqan_warmstart binary missing."; exit 1; }

# Informational: what the |OFFDIAG| vs OP relationship makes the aligner prefer.
REGIME="$(awk -v o="${OFFDIAG}" -v p="${OP}" 'BEGIN{ a=(o<0?-o:o);
  if (a>p)      printf "|sub|=%g > op=%g  -> identity-anchored (gappy)", a, p;
  else if (a<p) printf "|sub|=%g < op=%g  -> favours correspondence", a, p;
  else          printf "|sub|=%g = op=%g  -> balanced", a, p; }')"
NFORMS="$(grep -c '^>' "${SRC}")"

rule
say "Lardil SeqAn warm start -- APPROXIMATE MSA ONLY (no exact solve, no CPLEX)"
say "input : ${SRC#${REPO_ROOT}/}  (${NFORMS} forms${MAXFORMS:+, using first ${MAXFORMS}})"
say "gaps  : op=${OP} / ep=${EP}   matrix: diagonal=${DIAG} off-diagonal=${OFFDIAG}"
say "regime: ${REGIME}"
say "output: ${OUTDIR#${REPO_ROOT}/}/"
rule

# --- step 1: unicode -> ASCII + matrix -------------------------------------
ASCII="${OUTDIR}/ascii.fa"
MATRIX="${OUTDIR}/matrix.mat"
LEGEND="${OUTDIR}/legend.tsv"
# `${arr[@]+"${arr[@]}"}` guard: macOS ships Bash 3.2, where "${arr[@]}" on an
# EMPTY array under `set -u` aborts. The guard expands to nothing when empty.
MAXARG=(); [ -n "${MAXFORMS}" ] && MAXARG=(--max-forms "${MAXFORMS}")

say ""
say "[1/2] lardil_prep.py: transliterate to single-char ASCII + emit matrix"
python3 "${PREP}" build --in "${SRC}" --out-fasta "${ASCII}" --out-matrix "${MATRIX}" \
        --out-legend "${LEGEND}" --diagonal "${DIAG}" --off-diagonal "${OFFDIAG}" ${MAXARG[@]+"${MAXARG[@]}"} \
        | sed 's/^/      /'
say "      ASCII FASTA -> ${ASCII#${REPO_ROOT}/}"
say "      matrix      -> ${MATRIX#${REPO_ROOT}/}"

# --- step 2: SeqAn warm start ----------------------------------------------
WARM="${OUTDIR}/warmstart.fa"
WARM_IPA="${OUTDIR}/warmstart_ipa.fa"
say ""
say "[2/2] seqan_warmstart on the ASCII data + custom matrix"
say "      ./seqan_warmstart -in ${ASCII#${REPO_ROOT}/} -out ${WARM#${REPO_ROOT}/} -matrix ${MATRIX#${REPO_ROOT}/} -gapopen ${GAPOPEN} -gapextend ${GAPEXTEND} -center 0.0 -quiet"
./seqan_warmstart -in "${ASCII}" -out "${WARM}" -matrix "${MATRIX}" \
                  -gapopen "${GAPOPEN}" -gapextend "${GAPEXTEND}" -center 0.0 -quiet

if have_scorer; then
  say "      validation:"
  python3 "${SCORE_LIB}" validate --input "${ASCII}" "${WARM}"
  W_SOP="$(python3 "${SCORE_LIB}" score --matrix "${MATRIX}" --op "${OP}" --ep "${EP}" "${WARM}")"
  say "      warm-start sum-of-pairs (would be msa_aff's incumbent): ${W_SOP}"
else
  say "      (score_sop.py not found at ${SCORE_LIB#${REPO_ROOT}/} -- skipping validation/score)"
fi

# Decode to readable IPA, saved alongside the ASCII alignment and printed.
python3 "${PREP}" decode --in "${WARM}" --out "${WARM_IPA}"
say ""
say "      warm-start alignment (IPA) -> ${WARM_IPA#${REPO_ROOT}/}:"
sed 's/^/        /' "${WARM_IPA}"

# If a capped exact result is sitting in the sibling MWE dir, say so plainly so
# the two are never confused. We never touch that file.
MWE_ALN="${RESULTS_DIR}/${STEM}/aligned.fa"
if [ -f "${MWE_ALN}" ]; then
  N_WARM="$(grep -c '^>' "${WARM}")"
  N_ALN="$(grep -c '^>' "${MWE_ALN}" 2>/dev/null || echo 0)"
  if [ "${N_WARM}" != "${N_ALN}" ]; then
    say ""
    say "note: an exact-solve alignment exists at ${MWE_ALN#${REPO_ROOT}/} (${N_ALN} forms);"
    say "      this warm start has ${N_WARM}. They are different instances -- the exact"
    say "      solver is capped at 12 forms. This script leaves that file untouched."
  fi
fi

say ""
rule
say "DONE -- approximate (SeqAn) MSA of all ${MAXFORMS:-${NFORMS}} forms:"
say "  ASCII : ${WARM#${REPO_ROOT}/}"
say "  IPA   : ${WARM_IPA#${REPO_ROOT}/}"
say "Produced without the exact solve, so the >12-form msa_aff crash cannot occur."
say "It is a valid heuristic alignment, not a proven optimum -- for that, run"
say "lardil_mwe.sh with MAXFORMS<=12."
rule
