# MSAMDD-Ling_MWE — Phase 4 validation

Minimal working examples that validate the SeqAn warm-start integration end to
end. Two questions:

1. **Did swapping MUSCLE for SeqAn change anything it must not?** The warm start
   only seeds the exact solver's first incumbent, so it may change *speed* but
   never the *optimum*. We confirm `msa_aff` reaches the same proven optimum from
   either warm start, and check that SeqAn's incumbents are in the same quality
   ballpark as MUSCLE's.
2. **Can the toolchain now do what MUSCLE could not?** We align a small
   morphological paradigm over a non-protein phoneme alphabet with a custom
   substitution matrix — input MUSCLE cannot consume — and solve it to a proven
   optimum.

Everything runs inside the **linux/amd64** Docker image (the bundled `muscle`
and CPLEX-linked `msa_aff` are x86-64 Linux binaries; they do not run natively
on Apple Silicon). You run it; the scripts are arranged so you only ever invoke
one wrapper.

## Prerequisites (one-time)

- Docker Desktop, with **Use Rosetta for x86_64/amd64 emulation** enabled.
- CPLEX Linux x86-64 installer at `MSAMDD_MWE/docker/cplex_linux.bin`
  (see `MSAMDD_MWE/docker/README.md` — same file the Phase 0 smoke test used).
- SeqAn submodule initialised: `git submodule update --init`.

The runner **reuses the `msamdd-ling:amd64` image from Phase 0** and bind-mounts
your live repo, so it does not rebuild the image (no apt/CPLEX, no network) on a
normal run. It builds the image only if it is missing or you pass `FORCE_BUILD=1`
— that one-time build is the only step needing the CPLEX `.bin` and a network
inside Docker.

## Run it

From the repository root (`MSAMDD-Ling/`):

```bash
# Everything (reuse the Phase-0 image, then all three validations):
bash MSAMDD-Ling_MWE/run_in_docker.sh

# Or one stage at a time:
bash MSAMDD-Ling_MWE/run_in_docker.sh incumbent     # warm-start quality, no solve (fast)
bash MSAMDD-Ling_MWE/run_in_docker.sh optimum       # same proven optimum, end-to-end
bash MSAMDD-Ling_MWE/run_in_docker.sh linguistic    # exact MSA over a phoneme alphabet
```

Tunables are environment variables, forwarded into the container:

```bash
# Prove optimality on a chosen set, with a higher per-arm time budget:
TIME_LIMIT=1800 ENTRIES="Data/bali1/1tgxA.fa Data/bali1/1aab.fa Data/bali1/1tvxA.fa" \
  bash MSAMDD-Ling_MWE/run_in_docker.sh optimum
```

Outputs (logs, per-arm alignments) land in `MSAMDD-Ling_MWE/results/` on your
Mac. That directory is git-ignored.

## Troubleshooting

**`ERROR: missing CPLEX installer ... cplex_linux.bin`** — only the one-time image
build needs it. If the image already exists, this should not appear; if it does,
the image was pruned. Put the Linux x86-64 CPLEX `.bin` at
`MSAMDD_MWE/docker/cplex_linux.bin` and re-run.

**Image build fails with `Could not resolve archive.ubuntu.com` / connection
reset** — transient Docker DNS/network flakiness under amd64 emulation, not a
repo problem. Check whether the image already exists (then no build is needed):

```bash
docker image ls msamdd-ling:amd64
```

If it is listed, just re-run the validation — the runner reuses it. If you must
build, the apt step now retries, and: re-run once (layers resume), restart Docker
Desktop, disable any VPN, or add a DNS server (e.g. `8.8.8.8`) in Docker Desktop
→ Settings → Docker Engine. Once built, later runs need no network.

## What each stage does

### `incumbent` — warm-start quality (no CPLEX solve, fast)
Runs `seqan_warmstart` and `muscle` standalone on a handful of BAliBASE entries
with the exact flags `msa_aff` passes, then scores both alignments with the
solver's own sum-of-pairs objective (`lib/score_sop.py`, a faithful port of
`get_value`). The printed score is exactly the incumbent `msa_aff` would load
from each tool, so "higher = better warm start". Width / gap% columns show
whether one tool over-inserts gaps.

### `optimum` — same proven optimum from either warm start (end-to-end)
The exact solver is one fixed binary; the only thing varied is which aligner
seeds it. `msa_aff` shells out to `./seqan_warmstart`, so the SeqAn arm uses the
real binary and the MUSCLE arm temporarily replaces it with a one-line shim that
`exec ./muscle "$@"` — i.e. the warm-start call runs MUSCLE with byte-for-byte
the flags the old `run_MUSCLE` used. Both arms are run to optimality; a fixed
MILP has one global optimum, so the objectives must match. They differ only in
the warm-start incumbent and the time to converge — which is the whole claim:
*warm start = speed, not optimum*.

The script reports, per entry, whether each arm actually **proved** optimality
(`Gap = 0`). Only proven-vs-proven pairs are a clean comparison; if an arm times
out the entry is flagged `INCONCLUSIVE` (raise `TIME_LIMIT` or pick smaller
`ENTRIES`). The exact affine solve can be slow and is **not** guaranteed to
finish inside the time budget on every instance — start with the small defaults
and scale up deliberately.

### `linguistic` — exact MSA over a phoneme alphabet (end-to-end)
`linguistic/paradigm.fa` is a six-form toy verb paradigm written in a 13-symbol
single-character phoneme alphabet; `linguistic/phon.mat` is a custom matrix
(positive diagonal, uniform negative off-diagonal — the synchronic-morphology
convention, no natural classes) regenerable with `make_matrix.py`. Gaps use
Erich's scheme (`-op 1 -ep 0`); all values are integers, so they round-trip
losslessly through `seqan_warmstart`'s integer score. The script:

1. runs MUSCLE on this input + matrix and **captures what it actually does** (it
   is built for the 20 amino acids and a protein matrix; a phoneme matrix with
   non-amino-acid labels is outside what it accepts);
2. runs `seqan_warmstart` on the same input — a valid warm-start MSA;
3. runs `msa_aff` end-to-end to a **proven-optimal** alignment over the phoneme
   alphabet with the custom matrix.

## Files

| Path | Purpose |
|------|---------|
| `run_in_docker.sh` | Host entry point: build/reuse the amd64 image, run a stage, mount results back. |
| `validate_all.sh` | In-container driver: build once, run the selected stage(s). |
| `lib/score_sop.py` | Authoritative `get_value` port: score / validate any alignment. |
| `balibase/incumbent_bench.py` | SeqAn-vs-MUSCLE warm-start incumbent comparison. |
| `balibase/optimum_compare.sh` | End-to-end same-optimum check (SeqAn arm vs MUSCLE-shim arm). |
| `linguistic/paradigm.fa` | Phoneme-coded morphological paradigm. |
| `linguistic/phon.mat` | Custom substitution matrix over the phoneme alphabet. |
| `linguistic/make_matrix.py` | Regenerates `phon.mat`. |
| `linguistic/linguistic_mwe.sh` | MUSCLE-refuses / SeqAn-aligns / msa_aff-optimal. |

## Note on the MUSCLE arm

The current solver short-circuits the intra-Benders refine call, so via the shim
MUSCLE supplies only the *initial* warm start (no `muscle -refine` mid-solve).
That cannot change the optimum — which is the property under test — but it means
the `optimum` stage is a same-optimum check, not a re-run of the old end-to-end
MUSCLE refine behaviour.
