## Notes

This directory contains the files needed to generate
a clang++ executable that is optimized for line-level profiling

Edit Makefile.conf so that the default build has other properties,
or use this to pass additional options to make:
```
make -j <n> EXTRA="-DFOO -DBAR"
```
Makefile.conf has details for build

we only need to tweak the flags so binaries keep line tables and frame pointers while staying optimized. That gives great source-level profiling in perf/Hotspot/pprof without ballooning debug size.

for Makefile.conf

added -fno-omit-frame-pointer; 
switched -g → -gline-tables-only for lean line info; 
kept other existing opts

Why these exact changes
	•	-fno-omit-frame-pointer: lets perf -g fp (or DWARF) recover stacks cheaply and reliably in large C++.
	•	-gline-tables-only (Clang): keeps line numbers + inlining info with a fraction of the size of full -g, perfect for source/line profiling.
	•	Keep -O3 for your “normal” binaries so hotspots are representative; use DEBUGOPT only if your top-level Make lets you switch modes.


workflow
cd cloudy/source/sys_llvm
make clean && make -j
# profile
perf record -F 400 --call-graph dwarf -- ./cloudy < path/to/input.in
perf report
# line-level for a hot symbol:
perf annotate --stdio --source --symbol 'SomeHotFunction'


final Makefile.conf below
-----------------------------
# source/sys_llvm/Makefile.conf

# avoid -menable-unsafe-fp-math as it breaks the vectorized math routines
# this implies that -ffast-math or -Ofast should not be used
CXX := clang++

ifeq "arm64-" "$(findstring arm64-, $(shell ${CXX} --version | grep Target))"
  ARCH = -mcpu=apple-m1
else
  ARCH = -march=native
endif

# Keep high optimization, but ensure reliable stacks for profilers:
#  -fno-omit-frame-pointer  → cheap & robust stack unwinding
#  -gline-tables-only       → line info without huge debug size (Clang)
#  -fno-limit-debug-info    → keep inline attribution sensible under -O3
OPT      = -O3 -fno-signed-zeros ${ARCH} -fasynchronous-unwind-tables \
           -fno-omit-frame-pointer -fno-limit-debug-info \
           -Wno-deprecated -Qunused-arguments

# For “debug” builds you still want optimization for realistic hotspots.
# Use -O2 and keep line tables; switch to full -g if you need variables.
DEBUGOPT = -O2 -fasynchronous-unwind-tables -fno-omit-frame-pointer

# Use line tables by default (great for perf/Hotspot/pprof annotate):
CXXFLAGS    = ${OPT} -Wall -W -std=c++17 -gline-tables-only
CXXFLAGSEXC = -fnon-call-exceptions

# Link with the same unwind/frame-pointer settings and line tables
LDFLAGS = ${OPT} -Wall -W -gline-tables-only
