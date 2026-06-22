# Running the MSAMDD smoke test in Docker (Apple Silicon)

This builds `msa_aff` and runs the smoke test inside a **Linux x86-64** container.
That matters because the Makefile (`SYSTEM=x86-64_linux`), the CPLEX libraries, and
the bundled `muscle` warm-start binary are all x86-64. On your M1 the container runs
under emulation. Inside it, `muscle` runs for real — so this is also how you confirm
the MUSCLE warm-start actually works, not just that the solver tolerates its absence.

You provide one thing the repo can't: the proprietary **IBM CPLEX Linux installer**.
Everything else is automated.

---

## 0. One-time Docker Desktop setup (optional but recommended)

Open Docker Desktop → Settings → General → enable **"Use Rosetta for x86_64/amd64
emulation"** → Apply & Restart. This makes the emulated build noticeably faster.
Make sure Docker Desktop is **running** (whale icon in the menu bar) before continuing.

## 1. Download the CPLEX **Linux x86-64** installer

You have the macOS package; you need the Linux one from the same IBM account.

1. Go to the IBM Academic Initiative and sign in with your IBMid:
   <https://www.ibm.com/academic> → search "CPLEX Optimization Studio".
   (Direct download notes: <https://www.ibm.com/support/pages/downloading-ibm-ilog-cplex-optimization-studio-2211>)
2. Choose **"Get the no-cost academic edition"** → Software → **ILOG CPLEX
   Optimization Studio** → Download.
3. In the file list pick the **Linux x86-64** build — a file named like
   `cplex_studio2211.linux_x86_64.bin` (the version number may differ; match your
   22.x entitlement). **Do not** pick macOS, Windows, or Linux-on-ARM/Power.
4. Download it (~500 MB).

If the portal only offers what you already have, the macOS `.pkg` will **not** work
here — you specifically need the `*.linux_x86_64.bin`.

## 2. Drop the installer into this folder and rename it

Move the downloaded `.bin` into `MSAMDD_MWE/docker/` and rename it to exactly
`cplex_linux.bin`. From a Mac terminal, e.g.:

```bash
mv ~/Downloads/cplex_studio*.linux_x86_64.bin \
   "<path-to-repo>/MSAMDD-Ling/MSAMDD_MWE/docker/cplex_linux.bin"
```

This file is git-ignored and Docker keeps it out of the final image — it never ships.

## 3. Build and run

From the repo root (`MSAMDD-Ling/`), the one-shot wrapper does everything:

```bash
bash MSAMDD_MWE/docker/build_and_run.sh
```

That's equivalent to these two commands if you prefer to run them yourself:

```bash
# build (run from repo root)
docker build --platform linux/amd64 -f MSAMDD_MWE/docker/Dockerfile -t msamdd-ling:amd64 .

# run the smoke test, writing results back to the host
docker run --rm --platform linux/amd64 \
  -v "$PWD/MSAMDD_MWE/results:/work/MSAMDD_MWE/results" \
  msamdd-ling:amd64
```

The first build is slow (CPLEX install + compile under emulation — minutes, not
seconds). Rebuilds are cached. Results land in `MSAMDD_MWE/results/` on your Mac.

## 4. What success looks like

```
[2/4] probe MUSCLE warm-start binary (standalone)
      ./muscle       : ELF 64-bit LSB executable, x86-64 ...
      -> MUSCLE IS RUNNING: produced a valid warm-start alignment.
...
RESULT: PASS  -  toolchain works; valid alignment written to:
  MSAMDD_MWE/results/1tgxA_aligned.fa
MUSCLE warm-start status (standalone probe): working
```

`PASS` means the solver builds, runs, and emits a valid alignment. The MUSCLE line
tells you separately that the warm-start binary executed — which is the whole reason
for running under amd64.

## 5. Troubleshooting

- **`missing CPLEX installer ...`** — step 2 not done, or the file isn't named
  exactly `cplex_linux.bin`.
- **Build fails during CPLEX install** — your `.bin` may be a different
  architecture (re-download the **linux_x86_64** one) or a corrupt/partial
  download. The installer log prints to the build output.
- **`No valid Java virtual machine ... exit code 106`** — the CPLEX installer
  needs Java >=11. The Dockerfile already installs a JRE and points the installer
  at it (`LAX_VM /usr/bin/java`); if you still see this, the JRE apt-install above
  it likely failed (check the build output).
- **Build fails at `make msa_aff`** with CPLEX-not-found — confirm the install
  step printed a `CPLEX installed at: /opt/ibm/ILOG/CPLEX_Studio...` line. If your
  edition installed elsewhere, edit `USER_INSTALL_DIR` in `cplex_silent.properties`.
- **Runtime error about a missing `libcplex*.so`** — the build uses static CPLEX
  libs so this is unlikely, but if it happens, the CPLEX `bin/x86-64_linux`
  directory needs to be on `LD_LIBRARY_PATH`; tell me and I'll add it to the image.
- **Community Edition only** (no academic entitlement) — CPLEX CE caps model size
  (~1000 vars/constraints) and may reject even this small instance. Use the full
  academic edition.
- **Very slow** — enable Rosetta (step 0); emulated CPLEX is CPU-heavy.
