# syntax=docker/dockerfile:1
#
# Build + run the MSAMDD-Ling exact MSA toolchain in a Linux x86-64 environment.
#
# Everything in this stack is x86-64: the Makefile target (SYSTEM=x86-64_linux),
# the CPLEX libraries, and the bundled `muscle` ELF. So the image is pinned to
# linux/amd64. On Apple Silicon (M1/M2/M3) Docker Desktop runs it under emulation
# - enable "Use Rosetta for x86_64/amd64 emulation" in Docker Desktop settings for
# a big speedup. Inside this amd64 container the bundled muscle warm-start binary
# runs natively, so the MUSCLE step is actually exercised (not just tolerated).
#
# Build context MUST be the repo root (MSAMDD-Ling/). You provide the CPLEX Linux
# x86-64 installer at MSAMDD_MWE/docker/cplex_linux.bin (see docker/README.md).
# That .bin is NOT committed to the repo.

# ---- Stage 1: install CPLEX + build the solver -----------------------------
FROM --platform=linux/amd64 ubuntu:22.04 AS build

ENV DEBIAN_FRONTEND=noninteractive
# openjdk: the CPLEX installer (InstallAnywhere) needs a Java VM >=11. Its own
# bundled JRE fails to launch under emulation, so we provide one and point the
# installer at it explicitly with LAX_VM. Java is only needed to install CPLEX,
# not to run the solver, so it stays in this build stage (not the runtime image).
# Acquire::Retries hardens the fetch against transient mirror/DNS drops, which
# are common when this amd64 image is built under emulation on Apple Silicon.
RUN apt-get update -o Acquire::Retries=5 \
    && apt-get install -y --no-install-recommends -o Acquire::Retries=5 \
        build-essential g++ make file ca-certificates python3 \
        openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Silent CPLEX install. The installer reads the properties file via -f; LAX_VM
# forces it to use the JRE we just installed instead of its broken bundled one.
# CPLEX lands under /opt/ibm/ILOG/CPLEX_Studio, which the Makefile finds via
# `find $(BASISDIR) -maxdepth 1 -type d -name "CPLEX_Studio*"` (BASISDIR=/opt/ibm/ILOG).
COPY MSAMDD_MWE/docker/cplex_silent.properties /tmp/cplex_silent.properties
COPY MSAMDD_MWE/docker/cplex_linux.bin /tmp/cplex_linux.bin
RUN chmod +x /tmp/cplex_linux.bin \
    && /tmp/cplex_linux.bin LAX_VM /usr/bin/java -f /tmp/cplex_silent.properties \
    && rm -f /tmp/cplex_linux.bin /tmp/cplex_silent.properties \
    && echo "CPLEX installed at:" \
    && find /opt/ibm/ILOG -maxdepth 1 -type d -name 'CPLEX_Studio*' \
    && test -n "$(find /opt/ibm/ILOG -maxdepth 1 -type d -name 'CPLEX_Studio*')"

# Build the exact affine solver. Remove the installer copy that the bulk `COPY .`
# drags in, so it never reaches the final image.
WORKDIR /work
COPY . /work
RUN rm -f /work/MSAMDD_MWE/docker/cplex_linux.bin \
    && make msa_aff \
    && test -x /work/msa_aff \
    && ls -l /work/msa_aff

# ---- Stage 2: clean runtime image ------------------------------------------
# Carries the built binary, the repo, and CPLEX (harmless if msa_aff is statically
# linked; present in case any component is resolved at runtime). No installer.
FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update -o Acquire::Retries=5 \
    && apt-get install -y --no-install-recommends -o Acquire::Retries=5 \
        g++ make file python3 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/ibm/ILOG /opt/ibm/ILOG
COPY --from=build /work /work
WORKDIR /work

# Default action: run the smoke test (build is skipped since msa_aff is present).
CMD ["bash", "MSAMDD_MWE/smoke_test.sh"]
