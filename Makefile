# set of flags for compiling
CC = c++
CFLAGS = -O3 -Wall -shared -std=c++11 -fPIC -Wno-sign-compare -Wno-unknown-pragmas
SRC = PICAlgorithm.cpp interpolationMethods.cpp iterativePoissonSolvers.cpp
DEP = particle.h interpolationMethods.h iterativePoissonSolvers.h
DEBUGFLAG = -g -O0
NONUMAFLAG = -D PARTICLENONUMA

# sequential code target
main: $(SRC) $(DEP)
	$(CC) $(CFLAGS) $(shell python3 -m pybind11 --includes) $(SRC) -o PICAlgorithm$(shell python3-config --extension-suffix ) $(NONUMAFLAG)

# parallel code target with NUMA considerations
omp: $(SRC) $(DEP)
	$(CC) $(CFLAGS) -fopenmp $(shell python3 -m pybind11 --includes) $(SRC) -o PICAlgorithm$(shell python3-config --extension-suffix )

# parallel code target without NUMA considerations
ompNoNUMA: $(SRC) $(DEP)
	$(CC) $(CFLAGS) -fopenmp $(shell python3 -m pybind11 --includes) $(SRC) -o PICAlgorithm$(shell python3-config --extension-suffix ) $(NONUMAFLAG)

# debug code target
debug: $(SRC) $(DEP)
	$(CC) $(CFLAGS) $(shell python3 -m pybind11 --includes) $(SRC) -o PICAlgorithm$(shell python3-config --extension-suffix ) $(DEBUGFLAG)

# clean target 
clean:
	rm -f PICAlgorithm$(shell python3-config --extension-suffix )
