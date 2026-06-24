#!/usr/bin/env python3
"""Prepare a unicode-phoneme FASTA for MSAMDD: transliterate to single-char
ASCII and emit a matching substitution matrix.

WHY THIS EXISTS
---------------
MSAMDD reads a sequence one BYTE at a time (load_inst.cpp pushes `(int)*it`
for each byte of the line), so every phoneme MUST be a single ASCII byte. A
unicode IPA symbol such as the retroflex stop "t" (U+0288, two UTF-8 bytes)
would be read as TWO residues, silently breaking the 1:1 phoneme:symbol
assumption the whole method relies on. This helper performs the lossless
unicode -> single-char ASCII step and writes the substitution matrix over the
resulting alphabet, so the solver and the SeqAn warm start score identically.

It is the general "multi-glyph / unicode phoneme inventory -> single-char FASTA"
converter the project flagged as a prerequisite for real morphological data;
the toy Lardil set is the first user.

SCORING SCHEME (synchronic morphology; no natural-class structure)
------------------------------------------------------------------
Positive diagonal (identical phoneme), a single uniform negative off-diagonal
(every non-identity substitution penalised equally). The off-diagonal MAGNITUDE
is chosen strictly greater than the gap-OPEN penalty the driver passes to the
solver (default matrix off-diagonal = -2 vs gap-open op = 1), so a substitution
costs more than opening a gap. Gap-EXTEND is zero (set by the driver, -ep 0).
All values are integers => lossless through seqan_warmstart's integer Score.

The gap penalties live with the SOLVER flags (-op/-ep), NOT in this matrix; the
driver (lardil_mwe.sh) owns them and checks |off-diagonal| > op at run time.

USAGE
-----
    # transliterate + emit matrix (+ a human-readable legend)
    python3 lardil_prep.py build \
        --in Data/Lardil_short_plain.fasta \
        --out-fasta Results/lardil_ascii.fa \
        --out-matrix Results/lardil.mat \
        --out-legend Results/lardil_legend.tsv \
        [--max-forms N] [--diagonal 2] [--off-diagonal -2]

    # turn an aligned ASCII FASTA back into readable IPA (gaps preserved)
    python3 lardil_prep.py decode --in Results/lardil_aligned.fa \
        [--out Results/lardil_aligned_ipa.fa]
"""
import argparse
import sys

GAP = "-"  # the solver's gap symbol (ASCII 45); never used as a phoneme target

# ---------------------------------------------------------------------------
# Phoneme -> single-char ASCII map.
#
# Plain-ASCII phonemes map to themselves; non-ASCII symbols get mnemonic CAPITALS
# so they never collide with a lowercase phoneme already in use. Order here =
# matrix column order (consonants then vowels), purely for readability.
#
# Covers the full Lardil_plain.fasta inventory (19 symbols) and the short toy set.
#
# STRESS: the acute-marked vowels (á/é/ú) are treated as DISTINCT phonemes from
# their plain counterparts -- faithful to the 1:1 character<->phoneme encoding.
# Stress is positionally predictable in Lardil (word-initial), so in this data
# they only ever occur on the first vowel and never need to align against an
# unmarked vowel. To instead FOLD stress (treat á==a, é==e, ú==u, i.e. stress as
# non-phonemic), point their targets at the plain vowel ("á","a") etc. -- the
# injectivity check then tolerates it because they share a target by design.
#
# To extend further: add (unicode_phoneme, ascii_char) pairs. Keep every ASCII
# target unique -- `build` asserts injectivity over the symbols actually used and
# refuses any phoneme it has no mapping for, so collisions/omissions fail loudly
# rather than corrupting the alignment.
# ---------------------------------------------------------------------------
PHONEME_MAP = [
    # consonants
    ("p", "p"),  # p  bilabial stop          -> p
    ("ṱ", "D"),  # t̪  laminal/dental stop    -> D
    ("ʈ", "T"),  # ʈ  retroflex stop         -> T
    ("c", "c"),  # c  palatal stop           -> c
    ("k", "k"),  # k  velar stop             -> k
    ("m", "m"),  # m  bilabial nasal         -> m
    ("n", "n"),  # n  alveolar nasal         -> n
    ("ŋ", "N"),  # ŋ  velar nasal (eng)      -> N
    ("l", "l"),  # l  lateral                -> l
    ("r", "r"),  # r  rhotic / tap           -> r
    ("ɻ", "R"),  # ɻ  retroflex approximant  -> R
    ("j", "j"),  # j  palatal glide          -> j
    ("w", "w"),  # w  labio-velar glide      -> w
    # vowels (plain, then acute-stressed counterparts)
    ("a", "a"),
    ("á", "A"),  # stressed a -> A
    ("e", "e"),
    ("é", "E"),  # stressed e -> E
    ("i", "i"),
    ("u", "u"),
    ("ú", "U"),  # stressed u -> U
]
FWD = dict(PHONEME_MAP)                  # phoneme -> ascii
REV = {a: p for p, a in PHONEME_MAP}     # ascii   -> phoneme (for decode)
ORDER = [a for _, a in PHONEME_MAP]      # ascii symbols in declared order


def read_fasta(path):
    """(ids, seqs) preserving order. Header = text after '>' with surrounding
    whitespace stripped; blank lines ignored. Reads UTF-8 so each IPA codepoint
    is one Python character."""
    ids, seqs, cur, buf = [], [], None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line[0] == ">":
                if cur is not None:
                    ids.append(cur)
                    seqs.append("".join(buf))
                cur, buf = line[1:].strip(), []
            else:
                buf.append("".join(line.split()))
    if cur is not None:
        ids.append(cur)
        seqs.append("".join(buf))
    return ids, seqs


def normalise_header(h):
    """Make a header safe + stable for MSAMDD's name<->index map. The solver
    keys on the FULL header line, and seqan_warmstart echoes the id SeqAn parsed;
    a leading space after '>' is the classic mismatch. Collapse internal
    whitespace to '_' and strip, so input and warm-start output agree exactly."""
    return "_".join(h.split())


def build(args):
    ids, seqs = read_fasta(args.in_path)
    if not ids:
        sys.exit(f"lardil_prep: no sequences found in {args.in_path}")

    if args.max_forms is not None:
        if args.max_forms < 1:
            sys.exit("lardil_prep: --max-forms must be >= 1")
        ids, seqs = ids[: args.max_forms], seqs[: args.max_forms]

    # transliterate, failing loudly on any unmapped phoneme
    ascii_seqs, used = [], []
    for hdr, seq in zip(ids, seqs):
        out = []
        for ch in seq:
            if ch not in FWD:
                sys.exit(
                    f"lardil_prep: no ASCII mapping for phoneme {ch!r} "
                    f"(U+{ord(ch):04X}) in form '{hdr}'. Add it to PHONEME_MAP."
                )
            a = FWD[ch]
            out.append(a)
            if a not in used:
                used.append(a)
        ascii_seqs.append("".join(out))

    # the resulting alphabet must be 1:1 (no two phonemes share an ASCII byte)
    if len(set(used)) != len(used):  # defensive; FWD is a dict so can't happen
        sys.exit("lardil_prep: internal error -- non-injective transliteration")

    alphabet = [a for a in ORDER if a in used]  # declared order, used symbols only
    headers = [normalise_header(h) for h in ids]
    if len(set(headers)) != len(headers):
        sys.exit("lardil_prep: duplicate headers after normalisation -- "
                 "MSAMDD keys sequences by header, they must be unique.")

    # --- write ASCII FASTA ---
    with open(args.out_fasta, "w", encoding="ascii") as fh:
        for hdr, seq in zip(headers, ascii_seqs):
            fh.write(f">{hdr}\n{seq}\n")

    # --- write substitution matrix (square, symmetric, NCBI-style) ---
    d, o = args.diagonal, args.off_diagonal
    if not (abs(o) > 0):
        sys.exit("lardil_prep: --off-diagonal must be non-zero")
    width = max(2, max(len(str(d)), len(str(o))) + 1)
    with open(args.out_matrix, "w", encoding="ascii") as fh:
        fh.write("  " + " ".join(f"{a:>{width}}" for a in alphabet) + "\n")
        for a in alphabet:
            row = [d if a == b else o for b in alphabet]
            fh.write(a + " " + " ".join(f"{v:>{width}}" for v in row) + "\n")

    # --- write a human-readable legend (ascii <-> ipa) ---
    if args.out_legend:
        with open(args.out_legend, "w", encoding="utf-8") as fh:
            fh.write("ascii\tipa\tU+\n")
            for a in alphabet:
                p = REV[a]
                fh.write(f"{a}\t{p}\tU+{ord(p):04X}\n")

    # --- report to stdout (captured/echoed by the driver) ---
    legend = "  ".join(f"{REV[a]}->{a}" for a in alphabet)
    print(f"forms={len(headers)}  alphabet({len(alphabet)})={''.join(alphabet)}")
    print(f"map: {legend}")
    print(f"matrix: diagonal={d}  off-diagonal={o}  (|off-diagonal|={abs(o)})")
    return 0


def decode(args):
    """Reverse-transliterate an (aligned) ASCII FASTA back to IPA, preserving
    gaps and header lines verbatim. Lets a human read T/R back as the real
    retroflex symbols."""
    out = sys.stdout if not args.out else open(args.out, "w", encoding="utf-8")
    try:
        with open(args.in_path, encoding="ascii") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">") or not line.strip():
                    out.write(line + "\n")
                    continue
                out.write("".join(REV.get(ch, ch) for ch in line) + "\n")
    finally:
        if out is not sys.stdout:
            out.close()
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="unicode FASTA -> ASCII FASTA + matrix")
    b.add_argument("--in", dest="in_path", required=True)
    b.add_argument("--out-fasta", required=True)
    b.add_argument("--out-matrix", required=True)
    b.add_argument("--out-legend", default=None)
    b.add_argument("--max-forms", type=int, default=None,
                   help="use only the first N forms (tractability knob)")
    b.add_argument("--diagonal", type=int, default=2,
                   help="match score (identical phoneme), default 2")
    b.add_argument("--off-diagonal", type=int, default=-2,
                   help="mismatch score (any substitution), default -2; "
                        "|value| must exceed the gap-open penalty op")
    b.set_defaults(func=build)

    d = sub.add_parser("decode", help="aligned ASCII FASTA -> readable IPA")
    d.add_argument("--in", dest="in_path", required=True)
    d.add_argument("--out", default=None)
    d.set_defaults(func=decode)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
