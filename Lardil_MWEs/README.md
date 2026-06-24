# Lardil MWE — exact MSA of real linguistic data with the SeqAn warm start

A minimal working example: align a toy Lardil noun paradigm exactly with
`msa_aff`, using the **SeqAn** warm start (not MUSCLE, which rejects
non-amino-acid matrices). It is the first run of MSAMDD-Ling against real
linguistic data and doubles as the template for larger morphological paradigms.

## Run it

```sh
bash Lardil_MWEs/lardil_mwe.sh
```

from anywhere (paths resolve relative to the script). Pass a dataset as the first
argument (default the short toy set); the default scoring is **identity-anchored**
(the synchronic-morphology mode — see Scoring scheme). The full 33-form Lardil
paradigm lives at `Data/Lardil_plain.fasta`, but the exact solver **crashes above
12 forms** (see the Caveat), so cap it for an exact result:

```sh
MAXFORMS=12 bash Lardil_MWEs/lardil_mwe.sh Lardil_MWEs/Data/Lardil_plain.fasta
```

It builds the binaries if needed, then runs three steps. Output is namespaced per
input dataset under `Results/<dataset>/`, so runs on different inputs never
clobber each other (the short-set results sit in `Results/Lardil_short_plain/`).

In an interactive terminal you get live progress: a spinner during the build,
and during the exact solve a bar showing elapsed time against the `-time` budget
plus the solver's current optimality gap and incumbent (e.g.
`[solve] 0:42/10:00 [####------] 7%  gap 6.36%  inc 110.0`). The bar tracks time
used, not work done — an exact solver has no honest percent-complete; the gap
shrinking toward 0 is the real convergence signal. Piped/redirected runs print a
single static line instead.

| step | tool | output (under `Results/<dataset>/`) |
|------|------|--------|
| 1 | `lardil_prep.py` | `ascii.fa`, `matrix.mat`, `legend.tsv` |
| 2 | `seqan_warmstart` | `warmstart.fa` (heuristic seed handed to the exact solver) |
| 3 | `msa_aff` | `aligned.fa` (**the deliverable** — exact, provably-optimal when the solver closes) |

## Why the ASCII step

MSAMDD reads a sequence one **byte** at a time, so each phoneme must be a single
ASCII byte. The full Lardil set uses unicode symbols that are 2 bytes each and
would otherwise be read as two residues — retroflexes ʈ/ɻ, the velar nasal ŋ, a
laminal stop ṱ, and the stress-marked vowels á/é/ú. `lardil_prep.py` transliterates
losslessly to mnemonic single-byte ASCII (ʈ→`T`, ɻ→`R`, ŋ→`N`, ṱ→`D`, á→`A`, …;
plain-ASCII phonemes pass through), normalises headers, and writes the substitution
matrix over the resulting alphabet. `lardil_prep.py decode` turns an aligned ASCII
FASTA back into readable IPA, which the driver prints for you.

The stress-marked vowels (á/é/ú) are treated as **distinct phonemes** from their
plain counterparts — faithful to the 1:1 character↔phoneme encoding. Lardil stress
is positionally predictable (word-initial), so in this data they only ever occur on
the first vowel and never force a stressed/unstressed mismatch. To treat stress as
non-phonemic instead, fold them in `PHONEME_MAP` (á→`a`, é→`e`, ú→`u`).

## Scoring scheme

No natural-class structure: positive diagonal (identical phoneme), one uniform
negative off-diagonal (every substitution penalised equally), affine gaps set by
`OP`/`EP`. All values are integers, so they pass losslessly through the warm
start's integer score. The gap penalties are solver flags (`-op`/`-ep`), not
matrix entries.

The `|OFFDIAG|` vs `OP` relationship sets the behaviour, and the driver prints the
regime:

- **`|OFFDIAG| > OP` — identity-anchored (the default, `-2`/`1`).** A substitution
  costs more than a gap, so the optimiser only aligns *identical* segments and gaps
  the rest apart. This is the mode for **synchronic morphology** — line up the
  material shared across related inflected forms (confirmed with Erich as the
  intended default). On the toy set it anchors the shared `r`.
- **`|OFFDIAG| < OP` — positional correspondence (`OFFDIAG=-1 OP=2`).** A gap costs
  more than a substitution, so it lines up *slots* even where the segments differ.
  This is the mode for **diachronic / historical** alignment (cognate
  correspondence), a later phase of the project. On the toy set it gives the
  width-5 onset/vowel/`r`/suffix alignment.

## Caveat: the exact solve has a hard ceiling

The exact solve (step 3) is the point — it is what makes the alignment *provably
optimal* rather than heuristic. It has two limits:

- **Crashes above 12 forms.** At ≥13 forms `msa_aff` segfaults inside its Benders
  cut-separation — an upstream robustness bug, not a result. The script reports it
  and stops cleanly. So the full 33-form paradigm can't be solved as a single
  instance today; **12 is the current ceiling**.
- **Certification is slow even within range.** Runtime is wildly instance-
  sensitive; the 9-form toy found its (very likely) optimum but didn't *prove* it
  in 600 s. Exact MSA is NP-hard, so part of this ceiling is fundamental.

For an exact result, keep the instance small:

```sh
MAXFORMS=12 bash Lardil_MWEs/lardil_mwe.sh Lardil_MWEs/Data/Lardil_plain.fasta  # largest that runs
MAXFORMS=3  bash Lardil_MWEs/lardil_mwe.sh                                       # small -> closes fast
TIME_LIMIT=3600 MAXFORMS=12 bash Lardil_MWEs/lardil_mwe.sh Lardil_MWEs/Data/Lardil_plain.fasta
```

## Building msa_aff (CPLEX)

`seqan_warmstart` builds with just a C++14 compiler. `msa_aff` needs CPLEX, and
the Makefile finds it by globbing `BASISDIR` (default `/Applications` on macOS)
for `CPLEX_Studio*`. If that glob misses your install you'll see
`ilcplex/ilocplex.h file not found` and step 3 is skipped. Point it at your
install — easiest is to pin the studio dir:

```sh
ls -d /Applications/CPLEX_Studio* ~/CPLEX_Studio* 2>/dev/null   # find yours
BASISILOG=/Applications/CPLEX_Studio2211 bash Lardil_MWEs/lardil_mwe.sh
```

(macOS CPLEX ships x86-64, so `msa_aff` builds x86-64 and runs under Rosetta;
`seqan_warmstart` stays native arm64 — a Rosetta parent can exec a native child.)

## Knobs

`OP`, `EP`, `DIAG`, `OFFDIAG`, `TIME_LIMIT`, `MAXFORMS`, `REBUILD=1`
(force rebuild), `SKIP_BUILD=1` (use existing binaries), and CPLEX-location
passthrough `BASISILOG` / `BASISDIR` / `SYSTEM`. A first positional argument
overrides the input FASTA.
