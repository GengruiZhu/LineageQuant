#!/bin/bash
###############################################################################
# build.sh -- compile LineageQuant (seconds; touches no data)
#   bash build.sh                              # portable build (SSE4.2/POPCNT)
#   ARCHFLAGS=-march=native bash build.sh      # only if every compute node has the build machine's CPU
#   CXX=g++-12 bash build.sh
#   bash build.sh --static                     # fully static binary (no glibc / libstdc++ / zlib needed at run time):
#                                              # runs on any x86-64 Linux; this is what goes into a GitHub release
# Result: bin/lq, bin/lineagequant (the same, long name) and bin/lqcore (the engine name of 0.x / 1.0.0 scripts;
# `lqcore run` = one sample, `lq run` = the pipeline). libstdc++ is linked statically when possible.
# Needs GCC >= 7 (or clang >= 6) and zlib. Without them: conda install -c conda-forge gxx_linux-64 zlib make
#
# The new binaries are installed as bin/<name>.<time> and bin/<name> is switched to them atomically; older builds
# are kept (deleting or overwriting a binary that is running on another node of a network file system crashes it).
# Build on the cluster that runs it: a binary built elsewhere may need a newer glibc.
###############################################################################
cd "$(dirname "$0")" || exit 1
FULLSTATIC=0; [ "${1:-}" = --static ] && FULLSTATIC=1
ver () { "$1" -dumpfullversion -dumpversion 2>/dev/null | head -1 | cut -d. -f1; }
CAND=(${CXX:-} g++ x86_64-conda-linux-gnu-g++ x86_64-conda_cos7-linux-gnu-g++ x86_64-conda_cos6-linux-gnu-g++ clang++)
for c in "${CAND[@]}"; do
  [ -n "$c" ] && command -v "$c" >/dev/null || continue
  v=$(ver "$c"); [ "${v:-0}" -ge 7 ] && { CXXBIN=$c; break; }
  echo "  skip $c (version ${v:-?}, need >= 7)"
done
[ -n "${CXXBIN:-}" ] || { echo "x no GCC >= 7 found. Install one: conda install -c conda-forge gxx_linux-64 zlib make"; exit 1; }
echo "compiler: $CXXBIN $($CXXBIN -dumpfullversion -dumpversion)"
EXI=""; EXL=""
if [ -n "${CONDA_PREFIX:-}" ] && [ -f "$CONDA_PREFIX/include/zlib.h" ]; then
  EXI="-I$CONDA_PREFIX/include"; EXL="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib"
fi
AF=(); [ -n "${ARCHFLAGS:-}" ] && AF=(ARCHFLAGS="$ARCHFLAGS")
STAMP=$(date +%Y%m%d_%H%M%S)
BD=$(mktemp -d "$PWD/.build_XXXXXX") || { echo "x cannot create a build directory"; exit 1; }
cp -r src Makefile "$BD"/
# never reuse object files left in src/ by an in-place `make`: after the copy they look newer than the sources,
# so make would link stale objects into the new binary
rm -f "$BD"/src/*.o
cd "$BD" || exit 1
if [ $FULLSTATIC = 1 ]; then
  # needs the static C library (libc.a, from glibc-static / libc6-dev) and libz.a; libgomp's dlopen warning is harmless
  make -j4 CXX="$CXXBIN" "${AF[@]}" STATIC=-static EXTRA_CXXFLAGS="$EXI" EXTRA_LDFLAGS="$EXL" > build.log 2>&1 || {
    tail -20 build.log; cp build.log ../build.log
    echo "x fully static build failed: the static C library / libz.a is missing on this machine (e.g. dnf install glibc-static zlib-static,"
    echo "  apt install libc6-dev zlib1g-dev). A normal 'bash build.sh' builds the usual binary."; cd ..; rm -rf "$BD"; exit 1; }
elif ! make -j4 CXX="$CXXBIN" "${AF[@]}" EXTRA_CXXFLAGS="$EXI" EXTRA_LDFLAGS="$EXL" > build.log 2>&1; then
  echo "  static libstdc++ failed, linking dynamically"
  make clean >/dev/null 2>&1
  make -j4 CXX="$CXXBIN" "${AF[@]}" STATIC= EXTRA_CXXFLAGS="$EXI" EXTRA_LDFLAGS="$EXL" > build.log 2>&1 || {
    tail -30 build.log; cp build.log ../build.log; echo "x build failed (full log: build.log)"; cd ..; rm -rf "$BD"; exit 1; }
fi
grep -i "warning" build.log | head -5
cp build.log ../build.log
cd ..
mkdir -p bin
for n in lq lineagequant lqcore; do
  cp "$BD/bin/lineagequant" "bin/$n.$STAMP" && chmod 755 "bin/$n.$STAMP" || { echo "x cannot install bin/$n.$STAMP"; rm -rf "$BD"; exit 1; }
done
rm -rf "$BD"
"bin/lineagequant.$STAMP" version >/dev/null || { echo "x the new bin/lineagequant.$STAMP does not run"; exit 1; }
for n in lq lineagequant lqcore; do
  if [ -e "bin/$n" ] && [ ! -L "bin/$n" ]; then mv "bin/$n" "bin/$n.pre_$STAMP"; fi
  ln -sfn "$n.$STAMP" "bin/.$n.new" && mv -Tf "bin/.$n.new" "bin/$n"
  ls -1t bin/$n.2* 2>/dev/null | tail -n +21 | xargs -r rm -f   # keep the 20 most recent builds
done
./bin/lq version && ./bin/lqcore version && echo "ok $(pwd)/bin/lq -> $(readlink bin/lq)  (also bin/lineagequant, bin/lqcore)"
if [ $FULLSTATIC = 1 ]; then
  file -L bin/lq 2>/dev/null | grep -q "statically linked" && echo "   fully static: copy bin/lq (and the lineagequant / lqcore copies) to any x86-64 Linux"
fi
echo "   add to PATH:  export PATH=$(pwd)/bin:\$PATH"
