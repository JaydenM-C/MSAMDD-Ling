# --- PLATFORM ---
#
# Detect the host OS to pick the CPLEX platform directory and toolchain flags,
# so one Makefile builds the CPLEX solver on Linux (incl. Eureka2 / Rocky 8) and
# on macOS. On macOS the build matches the INSTALLED CPLEX architecture: native
# arm64 against a 22.1.1+ Apple-Silicon CPLEX (lib/arm64_osx), or x86-64 under
# Rosetta 2 against an Intel CPLEX (lib/x86-64_osx). SYSTEM defaults to the host
# arch and ARCHFLAG follows it. The seqan_warmstart target is CPLEX-independent
# and always builds NATIVELY -- see below.
#
# Anything here can be overridden on the command line if your install differs:
#   make SYSTEM=x86-64_linux BASISDIR=$$HOME/cplex     # CPLEX in HPC home/scratch
#   make SYSTEM=x86-64_osx                             # force Rosetta build (Intel CPLEX on M1)
#   make BASISILOG=/Applications/CPLEX_Studio222       # pin a specific studio dir
#   make CXX=g++-13                                    # a non-default compiler

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  # Default SYSTEM to the host arch (arm64 CPLEX on Apple Silicon, x86-64 on Intel
  # / under Rosetta); ARCHFLAG follows SYSTEM so the two never disagree. Override
  # SYSTEM to force the other (e.g. an x86-64 CPLEX on an M1 -> SYSTEM=x86-64_osx).
  UNAME_M   := $(shell uname -m)
  ifeq ($(UNAME_M),arm64)
    SYSTEM  ?= arm64_osx
  else
    SYSTEM  ?= x86-64_osx
  endif
  ifeq ($(SYSTEM),arm64_osx)
    ARCHFLAG := -arch arm64
  else
    ARCHFLAG := -arch x86_64   # x86-64 macOS libs -> Rosetta 2 on Apple Silicon
  endif
  BASISDIR  ?= /Applications
  PIEFLAG   :=                 # -no-pie is GNU/Linux-only; omit on macOS
  LDLIBS_OS := -framework CoreFoundation -framework IOKit -framework Accelerate  # CPLEX on macOS; Accelerate supplies CBLAS (cblas_dgemm)
  SORTV     := sort            # BSD sort lacks -V; pin BASISILOG if multi-version
else
  SYSTEM    ?= x86-64_linux
  BASISDIR  ?= /opt/ibm/ILOG
  ARCHFLAG  := -m64
  PIEFLAG   := -no-pie
  LDLIBS_OS := -ldl
  SORTV     := sort -V
endif

LIBFORMAT  = static_pic

# --- DIRECTORIES ---

CXX ?= g++
CCC = $(CXX) -std=gnu++11 $(ARCHFLAG) $(PIEFLAG) -Iincludes
BASISILOG ?= $(shell find $(BASISDIR) -maxdepth 1 -type d -name "CPLEX_Studio*" | $(SORTV) | tail -1)
CONCERTDIR = $(BASISILOG)/concert
CPLEXDIR   = $(BASISILOG)/cplex

# --- FLAGS ---

# Arch flag (-m64 / -arch x86_64) now lives in CCC so it applies at link too.
CCOPT = -fPIC -fno-strict-aliasing -fexceptions -DIL_STD -Wno-deprecated-declarations -Wno-ignored-attributes
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

# --- OPTIMIZATION FLAGS ---

DEBUG_OPT = -DNDEBUG -O3
#DEBUG_OPT = -g3 -O0
#PROF = -pg
PROF =

CFLAGS += $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -I./include $(DEBUG_OPT) -c $(PROF)

LDFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread $(LDLIBS_OS)

# ---- COMPILE  ----
SRC_DIR_cnv   := src/src_cnv
OBJ_DIR_cnv   := obj/obj_cnv

SRC_DIRS_cnv  := $(shell find $(SRC_DIR_cnv) -type d)
OBJ_DIRS_cnv  := $(addprefix $(OBJ_DIR_cnv)/,$(SRC_DIRS_cnv))

SOURCES_cnv   := $(shell find $(SRC_DIR_cnv) -name '*.cpp')
OBJ_FILES_cnv := $(addprefix $(OBJ_DIR_cnv)/, $(SOURCES_cnv:.cpp=.o))

SRC_DIR_aff   := src/src_aff
OBJ_DIR_aff   := obj/obj_aff

SRC_DIRS_aff  := $(shell find $(SRC_DIR_aff) -type d)
OBJ_DIRS_aff  := $(addprefix $(OBJ_DIR_aff)/,$(SRC_DIRS_aff))

SOURCES_aff   := $(shell find $(SRC_DIR_aff) -name '*.cpp')
OBJ_FILES_aff := $(addprefix $(OBJ_DIR_aff)/, $(SOURCES_aff:.cpp=.o))

vpath %.cpp $(SRC_DIRS_cnv)
vpath %.cpp $(SRC_DIRS_aff)

# ---- TARGETS ----

EXECUTABLE1 = msa_cnv 
EXECUTABLE2 = msa_aff

EXECUTABLES = $(EXECUTABLE1) $(EXECUTABLE2)

# SeqAn v2 warm-start aligner (drop-in replacement for the external MUSCLE
# call). Built standalone at C++17 against the vendored SeqAn headers; it links
# no CPLEX/Concert and is independent of the gnu++11 msa_* targets above, so the
# two toolchains never mix in a single translation unit. Note: NO $(ARCHFLAG)
# here -- it builds for the host's native arch (arm64 on Apple Silicon), which is
# fine even when msa_* are x86-64: a Rosetta parent can exec a native child.
SEQAN_INC     = includes/seqan/include
WARMSTART     = seqan_warmstart
WARMSTART_SRC = src/seqan_warmstart/seqan_warmstart.cpp
WARMSTART_CXX = $(CXX) -std=gnu++17 -O3 -DNDEBUG -Wno-deprecated-declarations

all: $(EXECUTABLES) $(WARMSTART)

$(EXECUTABLE1): makedir $(SOURCES_cnv) $(OBJ_FILES_cnv) 
	$(CCC) $(OBJ_FILES_cnv) $(LDFLAGS) $(PROF) -o $@

$(EXECUTABLE2): makedir $(SOURCES_aff) $(OBJ_FILES_aff)
	$(CCC) $(OBJ_FILES_aff) $(LDFLAGS) $(PROF) -o $@

$(WARMSTART): $(WARMSTART_SRC)
	$(WARMSTART_CXX) -I $(SEQAN_INC) $< -o $@

$(OBJ_DIR_cnv)/%.o: %.cpp
	$(CCC) $(CFLAGS) $< -o $@
$(OBJ_DIR_aff)/%.o: %.cpp
	$(CCC) $(CFLAGS) $< -o $@

makedir: $(OBJ_DIRS_cnv)
makedir: $(OBJ_DIRS_aff)

$(OBJ_DIRS_cnv):
	@mkdir -p $@
$(OBJ_DIRS_aff):
	@mkdir -p $@

clean:
	@rm -rf $(OBJ_DIR_cnv)
	@rm -rf $(OBJ_DIR_aff)
	@rm -rf $(EXECUTABLE1)
	@rm -rf $(EXECUTABLE2)
	@rm -rf $(WARMSTART)


