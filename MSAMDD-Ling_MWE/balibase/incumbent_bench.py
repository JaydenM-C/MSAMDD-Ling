#!/usr/bin/env python3
"""Head-to-head warm-start incumbent quality: SeqAn vs MUSCLE (no CPLEX needed).

Runs ./seqan_warmstart and ./muscle standalone on each input with the identical
flags msa_aff passes, then scores both alignments with MSA-MDD's own sum-of-pairs
objective (lib/score_sop.py). The score is exactly the incumbent the solver would
load from each tool, so "higher = better warm start". Width / gap% expose whether
one tool over-inserts gaps. This answers "are the SeqAn warm-start incumbents in
the same ballpark as MUSCLE's?" without running the (slow) exact solve.

Run from the repository root, where ./muscle, ./seqan_warmstart and ./Data live:

    python3 MSAMDD-Ling_MWE/balibase/incumbent_bench.py            # default handful
    python3 MSAMDD-Ling_MWE/balibase/incumbent_bench.py Data/bali4/BB1100*.fa

Knobs via env: OPN_PEN (default 12), EXT_PEN (default 2.22), MATRIX
(default ./Data/blosum.ncbi).
"""
import glob
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "lib"))
from score_sop import load_matrix, read_fasta, get_value, gap_stats  # noqa: E402

OPN = float(os.environ.get("OPN_PEN", "12"))
EXT = float(os.environ.get("EXT_PEN", "2.22"))
MATRIX = os.environ.get("MATRIX", "./Data/blosum.ncbi")
MUSCLE = "./muscle"
SEQAN = "./seqan_warmstart"

# A small default handful spanning short, easy instances (where the exact solve
# is also tractable) up to wider ones. Override by passing FASTA paths / globs.
DEFAULT = [
    "Data/bali1/1tgxA.fa", "Data/bali1/1aab.fa", "Data/bali1/1idy.fa",
    "Data/bali4/BB11001.fa", "Data/bali4/BB11025.fa", "Data/bali4/BB11021.fa",
]


def run_aligner(binary, infile, outfile):
    cmd = [binary, "-in", infile, "-out", outfile, "-matrix", MATRIX,
           "-gapopen", f"{-(OPN + EXT):.6f}", "-gapextend", f"{-EXT:.6f}",
           "-center", "0.0", "-quiet"]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True)
    except OSError as e:
        return -1, f"cannot exec {binary}: {e}"
    return p.returncode, (p.stderr or p.stdout).strip()


def main():
    for label, path in ((MUSCLE, MUSCLE), (SEQAN, SEQAN), ("matrix", MATRIX)):
        if not os.path.exists(path):
            sys.exit(f"missing {label}: {path}\n"
                     "Run from the repo root; `make seqan_warmstart` if the "
                     "binary is absent, and ensure ./muscle is present.")
    os.chmod(MUSCLE, 0o755)

    args = sys.argv[1:]
    files = []
    for a in (args or DEFAULT):
        files += sorted(glob.glob(a)) if any(c in a for c in "*?[") else [a]
    files = [f for f in files if os.path.exists(f)]
    if not files:
        sys.exit("no input files matched.")

    sub, _ = load_matrix(MATRIX)
    hdr = (f"{'file':<14}{'N':>3} {'SeqAn SOP':>12}{'MUSCLE SOP':>12}  win  "
           f"{'Swid':>5}{'Mwid':>5} {'Sgap%':>6}{'Mgap%':>6}")
    print(hdr)
    print("-" * len(hdr))

    sw = mw = ties = n_ok = 0
    dsop = 0.0
    with tempfile.TemporaryDirectory() as td:
        for f in files:
            name = os.path.splitext(os.path.basename(f))[0]
            _, inseqs = read_fasta(f)
            if len(inseqs) < 2:
                print(f"{name:<14}{len(inseqs):>3}  (skipped: <2 sequences)")
                continue
            so, mo = os.path.join(td, "s.fa"), os.path.join(td, "m.fa")
            rcS, errS = run_aligner(SEQAN, f, so)
            rcM, errM = run_aligner(MUSCLE, f, mo)
            if rcS == -1:
                sys.exit(f"\nseqan_warmstart will not execute ({errS}).\n"
                         "Rebuild it for THIS machine's arch, forcing past any "
                         "stale binary:\n    rm -f seqan_warmstart && make seqan_warmstart")
            if rcM == -1:
                sys.exit(f"\nmuscle will not execute ({errM}).\n"
                         "It is an x86-64 Linux ELF; run inside the amd64 container.")
            if rcS or not os.path.exists(so):
                print(f"{name:<14}{len(inseqs):>3}  SeqAn FAILED: {errS[:40]}")
                continue
            if rcM or not os.path.exists(mo):
                print(f"{name:<14}{len(inseqs):>3}  MUSCLE FAILED: {errM[:40]}")
                continue
            _, Sseq = read_fasta(so)
            _, Mseq = read_fasta(mo)
            Ss, Ms = get_value(Sseq, sub, OPN, EXT), get_value(Mseq, sub, OPN, EXT)
            Sw, _, Sg = gap_stats(Sseq)
            Mw, _, Mg = gap_stats(Mseq)
            win = "S" if Ss > Ms else ("M" if Ms > Ss else "=")
            sw += win == "S"
            mw += win == "M"
            ties += win == "="
            dsop += Ss - Ms
            n_ok += 1
            print(f"{name:<14}{len(inseqs):>3} {Ss:>12.2f}{Ms:>12.2f}   {win}   "
                  f"{Sw:>5}{Mw:>5} {Sg:>5.1f} {Mg:>5.1f}")

    if n_ok:
        print("-" * len(hdr))
        print(f"files scored: {n_ok} | SeqAn-higher: {sw}  MUSCLE-higher: {mw}  "
              f"ties: {ties}")
        print(f"mean SOP delta (SeqAn - MUSCLE): {dsop / n_ok:+.2f}  "
              "(positive => SeqAn gives the better incumbent)")
        print("\nReminder: the warm start sets only the starting incumbent (=> solver"
              "\nspeed). The exact solver proves the SAME optimum from either tool;"
              "\nuse optimum_compare.sh to confirm that end-to-end.")


if __name__ == "__main__":
    main()
