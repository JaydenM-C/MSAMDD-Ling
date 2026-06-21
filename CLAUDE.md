# MSAMDD-Ling

Public fork of [aminhn/MSAMDD](https://github.com/aminhn/MSAMDD) adapted for 
multiple sequence alignment of linguistic data (initially morphological 
paradigms, later diachronic/phonological work). The main planned divergence 
from upstream is replacing the MUSCLE warm-start with SeqAn v2, to support 
arbitrary ASCII character sets and custom substitution matrices.

## Status

Early-stage. Build not yet verified locally; SeqAn v2 integration not yet 
started; no original code has landed. Treat the conventions below as 
provisional — update this file as real decisions get made.

## Relation to upstream

- The `upstream` remote points at `aminhn/MSAMDD` and is **read-only**. 
  Never push there.
- This fork's default branch (`master`, matching upstream) is protected. 
  Work on topic branches; open PRs into `master`; don't commit directly.
- Never fetch-merge from upstream without explicit confirmation in the 
  current turn. Upstream is a moving target and merges should be deliberate, 
  not ambient.

## SeqAn version

**v2, not v3.** The two libraries share a name and little else; APIs, 
templates, and idioms differ substantially, and v3 examples will not 
compile against v2. When generating SeqAn code, citing documentation, or 
recommending headers, explicitly target v2. If you find yourself reaching 
for `seqan3::` anything, stop.

## Code conventions

Match upstream MSAMDD style until something else is established. When 
in doubt:

- Look at adjacent code in the same file or module.
- Don't reformat code you're not otherwise touching.
- Keep modifications localised and reversible — this is a research fork 
  and individual experiments may be backed out.
- Prefer additive changes (new files, new functions) over invasive edits 
  to upstream code, where it doesn't cost clarity.

A real style guide, formatter config, and lint rules will be added here 
once the project has enough code to make those choices meaningful. Don't 
invent conventions pre-emptively.

## Build and test

TODO. Once a working build is verified locally, document here:

- Compiler and version
- SeqAn v2 install path or vendoring approach
- CMake/Make invocations and any required flags
- A minimal smoke test command and expected output
- Platform notes (macOS / Linux / HPC)

Until then, do not guess at build commands — ask, or read the upstream 
README and adapt.

## Git safety

- Never `git push` without explicit confirmation in the current turn. 
  Prior approval does not carry across turns or sessions.
- Never force-push. No `--force`, no `--force-with-lease`, regardless of 
  perceived safety.
- Commits are fine without asking. Write clear messages; reference the 
  specific upstream behaviour being changed where relevant.
- Never rewrite or squash existing history without confirmation.

## Scope reminder

Some project notes are kept outside this repository and should not be 
referenced from code, comments, or commit messages here. If a cross-
reference is needed, use this repo's commit hash as the anchor and point 
to it from the external notes — not the other way around.
