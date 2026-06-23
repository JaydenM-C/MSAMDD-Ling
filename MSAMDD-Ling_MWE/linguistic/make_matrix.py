#!/usr/bin/env python3
"""Generate the custom phoneme substitution matrix (phon.mat).

Scoring convention for synchronic morphology (no natural-class structure for
now): a positive diagonal (identical phonemes match) and a uniform negative
off-diagonal (any substitution is penalised equally). Pair this with Erich's
gap scheme -- open 1, extend 0 -- passed to the solver as `-op 1 -ep 0`. All
values are integers, so they round-trip losslessly through seqan_warmstart's
integer Score (only the protein defaults 12 / 2.22 ever round).

The alphabet is a 13-symbol single-character phoneme inventory. The matrix
declares the full inventory whether or not a given paradigm uses every symbol;
the format mirrors Data/blosum.ncbi so both Load_subsmat (the solver) and
seqan_warmstart's loader parse it identically.

    python3 make_matrix.py > phon.mat
"""
import sys

# stops / nasals / fricative / liquids, then vowels
ALPHABET = list("ptkmnsrlaeiou")
DIAGONAL = 2     # identical phoneme
OFF_DIAG = -1    # any substitution

def main():
    out = sys.stdout
    out.write("  " + "  ".join(ALPHABET) + "\n")
    for a in ALPHABET:
        scores = [DIAGONAL if a == b else OFF_DIAG for b in ALPHABET]
        out.write(a + " " + " ".join(f"{v:>2}" for v in scores) + "\n")

if __name__ == "__main__":
    main()
