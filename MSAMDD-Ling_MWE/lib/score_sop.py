#!/usr/bin/env python3
"""Score an alignment with MSA-MDD's own sum-of-pairs objective.

This is a faithful Python port of the affine solver's scoring path:

  * get_value()    -> src/src_aff/utility.cpp
  * Load_subsmat() -> src/src_aff/load_inst.cpp  (matrix + subsref index map)

Because the solver re-scores every warm start with get_value(), the number this
module returns for an aligned FASTA is exactly the incumbent msa_aff would
receive from that alignment (higher = better warm start). Keeping one
authoritative port here means the benchmark and the end-to-end checks cannot
drift apart in how they score.

Importable (load_matrix / read_fasta / get_value / validate) and runnable:

    # print the sum-of-pairs value of an alignment (for shell capture)
    python3 score_sop.py score --matrix phon.mat --op 1 --ep 0 aligned.fa

    # check an alignment is a genuine MSA of its input (any alphabet)
    python3 score_sop.py validate --input paradigm.fa aligned.fa
"""
import argparse
import sys

GAP = "-"  # ASCII 45, the solver's gap symbol


def load_matrix(path):
    """Parse an NCBI-format matrix exactly as Load_subsmat / subsref do.

    Row 0 holds the column-header characters (first char of each token); each
    later row is a label char followed by integer scores. subsref maps a header
    char to its 0-based column index; an UNKNOWN character maps to index 0,
    because the solver looks it up in a std::map<int,int> whose operator[]
    default-constructs a missing key to 0. Returns (sub, alphabet).
    """
    cols, mat = [], []
    with open(path) as fh:
        for line in fh:
            toks = line.split()
            if not toks:
                continue
            if not cols:
                cols = [t[0] for t in toks]
            else:
                mat.append([int(v) for v in toks[1:]])
    if not cols or not mat:
        sys.exit(f"score_sop: empty or invalid matrix file: {path}")
    idx = {c: i for i, c in enumerate(cols)}

    def sub(a, b):
        return mat[idx.get(a, 0)][idx.get(b, 0)]

    return sub, set(cols)


def read_fasta(path):
    """Return (ids, seqs) preserving order. ids keep the full header text after
    '>'; sequences are whitespace-stripped and case-preserved (the solver is
    case-sensitive over raw chars)."""
    ids, seqs, cur, buf = [], [], None, []
    with open(path) as fh:
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


def get_value(seqs, sub, opn, ext):
    """Faithful port of src_aff/utility.cpp get_value: sum of pairs with an
    affine gap of length L costing (opn + ext*L). A column where both rows are
    gaps is skipped and does NOT break a running gap; a non-gap column flushes
    any pending gap run; a trailing gap run is flushed at the end of each pair."""
    n = len(seqs)
    val = 0.0
    for i in range(n - 1):
        ai = seqs[i]
        for j in range(i + 1, n):
            bj = seqs[j]
            gl = 0
            for a, b in zip(ai, bj):
                ag, bg = (a == GAP), (b == GAP)
                if ag and bg:
                    continue
                if ag or bg:
                    gl += 1
                elif gl == 0:
                    val += sub(a, b)
                else:
                    val += sub(a, b) - (opn + ext * gl)
                    gl = 0
            if gl > 0:
                val -= opn + ext * gl
    return val


def gap_stats(seqs):
    width = len(seqs[0]) if seqs else 0
    cells = sum(row.count(GAP) for row in seqs)
    total = sum(len(row) for row in seqs) or 1
    return width, cells, 100.0 * cells / total


def validate(in_path, aln_path):
    """Return a list of problems (empty = genuine MSA). Alphabet-agnostic: the
    only reserved symbol is '-'. Checks the output has the same sequence set,
    equal row lengths, and de-gaps back to each exact input sequence."""
    in_ids, in_seqs = read_fasta(in_path)
    inp = dict(zip(in_ids, in_seqs))
    out_ids, out_seqs = read_fasta(aln_path)
    out = list(zip(out_ids, out_seqs))
    problems = []

    if set(inp) != {h for h, _ in out}:
        missing = sorted(set(inp) - {h for h, _ in out})
        extra = sorted({h for h, _ in out} - set(inp))
        if missing:
            problems.append(f"sequences missing from output: {missing}")
        if extra:
            problems.append(f"unexpected sequences in output: {extra}")

    lengths = {len(s) for _, s in out}
    if len(lengths) != 1:
        problems.append(f"aligned rows are not equal length: {sorted(lengths)}")
    elif next(iter(lengths)) == 0:
        problems.append("aligned rows have length 0")

    for hdr, row in out:
        if hdr in inp and row.replace(GAP, "") != inp[hdr]:
            problems.append(
                f"{hdr}: de-gapped row != input "
                f"(got {len(row.replace(GAP, ''))} residues, "
                f"expected {len(inp[hdr])})")
    return problems, out


def _cmd_score(args):
    sub, _ = load_matrix(args.matrix)
    _, seqs = read_fasta(args.aln)
    print(f"{get_value(seqs, sub, args.op, args.ep):.6f}")
    return 0


def _cmd_validate(args):
    problems, out = validate(args.input, args.aln)
    width = len(out[0][1]) if out else 0
    print(f"      sequences : {len(out)}    aln width : {width}")
    for hdr, row in out:
        print(f"        {hdr:<12} len={len(row):<4} "
              f"residues={len(row.replace(GAP, '')):<4} gaps={row.count(GAP)}")
    if problems:
        print("      VALIDATION PROBLEMS:")
        for p in problems:
            print(f"        - {p}")
        return 1
    print("      validation OK: genuine alignment of all input sequences.")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("score", help="print sum-of-pairs value of an alignment")
    s.add_argument("--matrix", required=True)
    s.add_argument("--op", type=float, required=True, help="gap-open penalty (opn_pen)")
    s.add_argument("--ep", type=float, required=True, help="gap-extend penalty (ext_pen)")
    s.add_argument("aln")
    s.set_defaults(func=_cmd_score)

    v = sub.add_parser("validate", help="check an alignment is a genuine MSA")
    v.add_argument("--input", required=True, help="original (unaligned) FASTA")
    v.add_argument("aln")
    v.set_defaults(func=_cmd_validate)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
