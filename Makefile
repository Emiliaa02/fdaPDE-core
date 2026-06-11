# LIBRARIES PATHS
FDAPDE_DIR := .

SYS_INCLUDE_DIR := /usr/include
ARCH_INCLUDE_DIR := /usr/include/x86_64-linux-gnu
NLOHMANN_DIR := /usr/include/nlohmann
EIGEN_DIR := /usr/include/eigen3

# FLAGS                  
CXXFLAGS := -std=c++20 -g -march=native -DFDAPDE_NO_DEBUG 
INCLUDES := -I$(EIGEN_DIR) -I$(FDAPDE_DIR) -I$(SYS_INCLUDE_DIR) -I$(ARCH_INCLUDE_DIR) -I$(NLOHMANN_DIR)
SHAPE ?= unit_square_16

# Run all the tests
.PHONY: all
all: generate_clipped_voronoi check_voronoi_definition computational_times

# Compare the total area of the Clipped Voronoi diagram and the area of the input mesh
.PHONY: generate_clipped_voronoi
generate_clipped_voronoi:
	g++ $(CXXFLAGS) $(INCLUDES) -O2 -o test/generate_clipped_voronoi \
	test/generate_clipped_voronoi.cpp
	./test/generate_clipped_voronoi $(SHAPE)
	python3 test/plots/plotting_clipped_voronoi.py --shape $(SHAPE)

# Check Voronoi definition 
.PHONY: check_voronoi_definition
check_voronoi_definition:
	g++ $(CXXFLAGS) $(INCLUDES) -O2 -o test/check_voronoi_definition \
	test/check_voronoi_definition.cpp 
	./test/check_voronoi_definition $(SHAPE)

# Computational Times 
.PHONY: computational_times
computational_times:
	g++ $(CXXFLAGS) $(INCLUDES) -O2 -o test/measure_computational_times \
	test/measure_computational_times.cpp
	./test/measure_computational_times
	python3 test/plots/plotting_computational_times.py

# Cleaning directories
.PHONY: clean
clean:
	@rm -f test/generate_clipped_voronoi  test/check_voronoi_definition test/measure_computational_times
	@rm -f test/plots/data/*.json
	@rm -f -r test/data/computational_times
	@rm -f test/plots/*.png