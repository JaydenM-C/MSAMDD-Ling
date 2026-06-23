#!/usr/bin/env bash
#
# linguistic_mwe.sh - end-to-end exact MSA on NON-protein data: a small
# morphological paradigm over a 13-symbol phoneme alphabet with a custom
# substitution matrix. This is the capability the project exists to add: MUSCLE
# cannot take an arbitrary alphabet + custom matrix, so the former warm start
# refuses this input; the SeqAn warm start consumes it and msa_aff proves an
# optimum.
#
# Three steps:
#   1. MUSCLE on the SAME input + custom matrix, with the exact flags msa_aff
#      would pass. We capture whatever MUSCLE actually does and report it (it is
#      built for the 20 amino acids and a protein matrix; a phoneme matrix with
#      non-amino-acid labels is outside what it accepts).
#   2. seqan_warmstart on the same input + matrix -> a valid warm-start MSA.
#   3. msa_aff end-to-end (-op 1 -ep 0, the custom matrix) -> a proven-optimal
#      alignment over the phoneme alphabet.
#
# Runs inside the linux/amd64 container, from the repository root.
#
# Env knobs: TIME_LIMIT (solver -time, default 600), OP (default 1), EP (default 0).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
LIB="${SCRIPT_DIR}/../lib/score_sop.py"
RESULTS="${SCRIPT_DIR}/../results"
cd "${REPO_ROOT}"

LING="MSAMDD-Ling_MWE/linguistic"
FASTA="${LING}/paradigm.fa"
MATRIX="${LING}/phon.mat"
TIME_LIMIT="${TIME_LIMIT:-600}"
OP="${OP:-1}"
EP="${EP:-0}"
# Flags exactly as run_warmstart builds them: -gapopen -(op+ep), -gapextend -ep.
GAPOPEN="$(awk -v o="${OP}" -v e="${EP}" 'BEGIN{printf "%g", -(o+e)}')"
GAPEXTEND="$(awk -v e="${EP}" 'BEGIN{printf "%g", -e}')"

mkdir -p "${RESULTS}"
say()  { printf '%s\n' "$*"; }
rule() { printf '%s\n' "===================================================================="; }

[ -f "${FASTA}" ]  || { say "ERROR: missing ${FASTA}"; exit 1; }
[ -f "${MATRIX}" ] || { say "ERROR: missing ${MATRIX}"; exit 1; }
[ -f "${REPO_ROOT}/muscle" ] || { say "ERROR: ./muscle missing."; exit 1; }
command -v python3 >/dev/null || { say "ERROR: python3 not found."; exit 1; }

# validate_all.sh builds once up front and sets SKIP_BUILD=1.
if [ "${SKIP_BUILD:-0}" != "1" ]; then
  say "[build] rm -f msa_aff seqan_warmstart && make msa_aff seqan_warmstart"
  rm -f msa_aff seqan_warmstart
  make msa_aff seqan_warmstart >"${RESULTS}/ling_build.log" 2>&1 || {
    say "  build FAILED -- tail:"; tail -n 25 "${RESULTS}/ling_build.log" | sed 's/^/    /'; exit 1; }
fi
[ -x "${REPO_ROOT}/msa_aff" ] && [ -x "${REPO_ROOT}/seqan_warmstart" ] || { say "ERROR: binaries missing (build first, or unset SKIP_BUILD)."; exit 1; }

ALPHA="$(head -1 "${MATRIX}" | wc -w | tr -d ' ')"
rule
say "Linguistic MWE: ${ALPHA}-symbol phoneme alphabet, custom matrix, gaps op=${OP}/ep=${EP}"
say "input : ${FASTA}"
say "matrix: ${MATRIX}  (positive diagonal, negative off-diagonal; integer => lossless)"
rule

# --- step 1: MUSCLE on the phoneme input + custom matrix --------------------
say ""
say "[1/3] MUSCLE on the phoneme input + custom matrix"
say "      ./muscle -in ${FASTA} -matrix ${MATRIX} -gapopen ${GAPOPEN} -gapextend ${GAPEXTEND} -center 0.0 -quiet"
chmod +x ./muscle 2>/dev/null || true
M_OUT="${RESULTS}/muscle_phoneme.fa"; M_LOG="${RESULTS}/muscle_phoneme.log"
rm -f "${M_OUT}"
set +e
./muscle -in "${FASTA}" -out "${M_OUT}" -matrix "${MATRIX}" \
         -gapopen "${GAPOPEN}" -gapextend "${GAPEXTEND}" -center 0.0 -quiet </dev/null >"${M_LOG}" 2>&1
M_RC=$?
set -e
if [ "${M_RC}" -ne 0 ]; then
  say "      -> MUSCLE REFUSED the input (exit ${M_RC}). Message:"
  grep -iE "invalid|matrix|residue|character|assert|error|cannot|fatal" "${M_LOG}" | head -4 | sed 's/^/         /' \
    || tail -4 "${M_LOG}" | sed 's/^/         /'
elif [ ! -s "${M_OUT}" ]; then
  say "      -> MUSCLE exited 0 but produced NO alignment (it could not use this input)."
elif python3 "${LIB}" validate --input "${FASTA}" "${M_OUT}" >/dev/null 2>&1; then
  say "      -> MUSCLE ran and emitted an alignment, BUT it cannot apply a phoneme"
  say "         matrix: it treats these symbols as amino acids and scores them with"
  say "         protein chemistry. The residues survive only because they happen to"
  say "         be letters; the linguist's substitution costs are ignored."
else
  say "      -> MUSCLE exited 0 but its output is NOT a valid MSA of the input"
  say "         (it mangled the non-protein symbols). Detail:"
  python3 "${LIB}" validate --input "${FASTA}" "${M_OUT}" 2>&1 | grep -A6 PROBLEMS | sed 's/^/         /' || true
fi
say "      (full MUSCLE log: ${M_LOG#${REPO_ROOT}/})"

# --- step 2: seqan_warmstart on the same input + matrix --------------------
say ""
say "[2/3] seqan_warmstart on the same input + custom matrix"
S_OUT="${RESULTS}/seqan_phoneme.fa"
./seqan_warmstart -in "${FASTA}" -out "${S_OUT}" -matrix "${MATRIX}" \
                  -gapopen "${GAPOPEN}" -gapextend "${GAPEXTEND}" -center 0.0 -quiet
say "      -> seqan_warmstart produced an alignment. Validation:"
python3 "${LIB}" validate --input "${FASTA}" "${S_OUT}"
S_SOP="$(python3 "${LIB}" score --matrix "${MATRIX}" --op "${OP}" --ep "${EP}" "${S_OUT}")"
say "      warm-start sum-of-pairs (what msa_aff would load as its incumbent): ${S_SOP}"

# --- step 3: msa_aff end-to-end over the phoneme alphabet ------------------
say ""
say "[3/3] msa_aff end-to-end (exact, custom matrix, op=${OP}/ep=${EP}, -time ${TIME_LIMIT})"
cp -f "${FASTA}" "${RESULTS}/ling_paradigm.fa"
A_OUT="${RESULTS}/ling_paradigm_aligned.fa"; A_LOG="${RESULTS}/ling_paradigm.log"
set +e
./msa_aff -in "${RESULTS}/ling_paradigm.fa" -out "${A_OUT}" -submat "${MATRIX}" \
          -op "${OP}" -ep "${EP}" -time "${TIME_LIMIT}" </dev/null >"${A_LOG}" 2>&1
A_RC=$?
set -e
if [ "${A_RC}" -ne 0 ] || [ ! -s "${A_OUT}" ]; then
  say "      msa_aff FAILED (exit ${A_RC}). Tail of log:"; tail -n 20 "${A_LOG}" | sed 's/^/        /'; exit 1
fi
INC="$(grep 'heuristic solution value:' "${A_LOG}" | tail -1 | sed 's/.*value:[[:space:]]*//')"
GAP="$(awk -F'Gap is:' '/Gap is:/{g=$2} END{gsub(/[ \t\r]/,"",g); print g}' "${A_LOG}")"
OBJ="$(python3 "${LIB}" score --matrix "${MATRIX}" --op "${OP}" --ep "${EP}" "${A_OUT}")"
say "      -> msa_aff completed. Output validation:"
python3 "${LIB}" validate --input "${FASTA}" "${A_OUT}"
say ""
say "      warm-start incumbent (SeqAn) : ${INC}"
say "      final objective (get_value)  : ${OBJ}"
say "      optimality gap reported      : ${GAP}%"
PROVEN="no (timed out -- raise TIME_LIMIT)"
if [ -n "${GAP}" ] && awk "BEGIN{exit !(${GAP} < 0.01)}"; then PROVEN="YES (gap = 0)"; fi
say "      proven optimal               : ${PROVEN}"
say ""
say "      optimal alignment:"
sed 's/^/        /' "${A_OUT}"

rule
say "RESULT: exact, provably-optimal MSA over a ${ALPHA}-symbol phoneme alphabet"
say "with a custom substitution matrix -- input MUSCLE cannot align. See step 1"
say "above for exactly how MUSCLE handled (or refused) the same input."
say "Artifacts in ${RESULTS#${REPO_ROOT}/}/: muscle_phoneme.*, seqan_phoneme.fa, ling_paradigm_aligned.fa"
rule
