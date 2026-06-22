#!/usr/bin/env python3
"""Validate that an MSA output file is a genuine alignment of the given input.

Usage: validate_alignment.py <input.fa> <alignment.fa>

Checks: same set of sequences as input; all aligned rows equal length; each row
de-gaps ('-' removed) back to its exact input residues; and (for unequal-length
inputs) gaps are present. Exits 0 on success, 1 on any problem. Prints a report.
"""
import sys


def read_fasta(path):
    """List of (header, sequence), order preserved; sequence uppercased and
    whitespace-stripped. Header keeps everything after '>'."""
    recs, hdr, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(seq)))
                hdr = line[1:].strip()
                seq = []
            else:
                seq.append("".join(line.split()).upper())
        if hdr is not None:
            recs.append((hdr, "".join(seq)))
    return recs


def main():
    if len(sys.argv) != 3:
        print("usage: validate_alignment.py <input.fa> <alignment.fa>")
        return 2

    inp = dict(read_fasta(sys.argv[1]))
    out = read_fasta(sys.argv[2])
    problems = []

    # 1. same set of sequences
    in_names, out_names = set(inp), {h for h, _ in out}
    if in_names != out_names:
        missing, extra = in_names - out_names, out_names - in_names
        if missing:
            problems.append(f"sequences missing from output: {sorted(missing)}")
        if extra:
            problems.append(f"unexpected sequences in output: {sorted(extra)}")

    # 2. all aligned rows equal length (defining property of an MSA)
    lengths = {len(s) for _, s in out}
    if len(lengths) != 1:
        problems.append(f"aligned rows are not equal length: {sorted(lengths)}")
    aln_len = next(iter(lengths)) if len(lengths) == 1 else None
    if aln_len == 0:
        problems.append("aligned rows have length 0")

    # 3. de-gapping each row must recover its exact input residues
    ALLOWED = set("ACDEFGHIKLMNPQRSTVWYBZXU*")  # residues; '-' is the gap symbol
    for hdr, row in out:
        if hdr not in inp:
            continue
        degapped = row.replace("-", "")
        if degapped != inp[hdr]:
            problems.append(
                f"{hdr}: de-gapped row != input "
                f"(got {len(degapped)} residues, expected {len(inp[hdr])})")
        stray = set(row) - ALLOWED - {"-"}
        if stray:
            problems.append(f"{hdr}: unexpected characters in row: {sorted(stray)}")

    # 4. sanity: unequal-length inputs => a valid MSA must contain gaps
    total_gaps = sum(s.count("-") for _, s in out)
    max_in = max((len(s) for s in inp.values()), default=0)
    if aln_len is not None and aln_len < max_in:
        problems.append(f"alignment length {aln_len} < longest input {max_in}")
    if len(set(len(s) for s in inp.values())) > 1 and total_gaps == 0:
        problems.append("inputs differ in length but output has no gaps")

    # --- report ---
    print(f"      sequences   : {len(out)} (expected {len(inp)})")
    print(f"      aln length  : {aln_len}")
    print(f"      input range : {min((len(s) for s in inp.values()), default=0)}"
          f"-{max_in} residues")
    print(f"      total gaps  : {total_gaps}")
    print("      per-seq     :")
    for hdr, row in out:
        print(f"        {hdr:<8} len={len(row):<4} residues={len(row.replace('-',''))} "
              f"gaps={row.count('-')}")

    if problems:
        print("\n      VALIDATION PROBLEMS:")
        for p in problems:
            print(f"        - {p}")
        return 1
    print("\n      validation OK: genuine alignment of all input sequences.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
