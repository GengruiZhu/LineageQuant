# LineageQuant 1.1.0 (C++17, OpenMP, zlib) -> bin/lq (= bin/lineagequant) + bin/lqcore    (bash build.sh does this safely)
#   make                      portable build (SSE4.2/POPCNT: runs on every node of the cluster)
#   make ARCHFLAGS=-march=native    tuned for the build machine (only if all compute nodes share its CPU)
#   make CXX=g++-12           choose a compiler (GCC >= 7 or clang >= 6 for C++17)
#   make STATIC=              link libstdc++ dynamically (if libstdc++.a is missing)
#   EXTRA_CXXFLAGS / EXTRA_LDFLAGS   e.g. -I$CONDA_PREFIX/include / -L$CONDA_PREFIX/lib
CXX       ?= g++
ARCHFLAGS ?= -msse4.2 -mpopcnt
STATIC    ?= -static-libstdc++ -static-libgcc
CXXFLAGS  := -O3 -std=c++17 -ffp-contract=off -Wall -Wextra -Wno-unused-parameter -fopenmp $(ARCHFLAGS) $(EXTRA_CXXFLAGS)
LDFLAGS   := -fopenmp $(STATIC) $(EXTRA_LDFLAGS)
LDLIBS    := -lz

SRC = src/main.cpp src/em.cpp src/run.cpp src/model.cpp src/lineage.cpp src/scan.cpp src/theta.cpp src/lqmap.cpp src/longread.cpp src/tools.cpp src/prep.cpp src/pipeline.cpp
HDR = src/binio.hpp src/prep.hpp src/proc.hpp src/lqmap.hpp src/common.hpp src/model.hpp src/sam.hpp src/lineage.hpp src/scan.hpp
OBJ = $(SRC:.cpp=.o)

all: bin/lineagequant bin/lq bin/lqcore

bin/lineagequant: $(OBJ)
	@mkdir -p bin
	$(CXX) $(OBJ) -o $@ $(LDFLAGS) $(LDLIBS)

# bin/lqcore: the engine name of 0.x / 1.0.0 scripts. A copy, not a link: the name decides what `run` means
# (lqcore run = one sample), and scripts often resolve links.
bin/lqcore: bin/lineagequant
	cp -f $< $@

# bin/lq: the short name (= lineagequant)
bin/lq: bin/lineagequant
	cp -f $< $@

src/%.o: src/%.cpp $(HDR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

PREFIX ?= /usr/local
install: all
	install -d $(PREFIX)/bin
	install -m 755 bin/lineagequant $(PREFIX)/bin/lineagequant
	install -m 755 bin/lq $(PREFIX)/bin/lq
	install -m 755 bin/lqcore $(PREFIX)/bin/lqcore

clean:
	rm -f src/*.o bin/lineagequant bin/lq bin/lqcore

.PHONY: all clean install
