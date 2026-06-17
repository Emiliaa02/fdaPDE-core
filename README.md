<h1> A C++ Implementation of the Clipped Voronoi Diagram </h1>

This repository is a fork of the `stable` branch of the `fdaPDE` library and contains a `C++` implementation of Unrestricted and Clipped Voronoi diagrams, integrated within `fdaPDE`. 

## Dependencies
The requirements for running the code and reproducing the described pipelines are listed below:

- a C++20-compliant compiler;
- `make`;
- `Eigen` library (version 3.4.0);
- `nhlomann/json` header;
- `Python 3`, along with the `matplotlib` package.

## Installation
To install the project locally on your laptop, open the terminal and run:

```bash 
git clone https://github.com/Emiliaa02/fdaPDE-core.git 
cd fdaPDE-core
```

## Instructions to run tests

Before proceeding, it is necessary to edit the initial Makefile variables to reflect the correct local paths on your system.
As the construction of the Voronoi diagram builds upon an existing Delaunay mesh, it is necessary that the `test/data/mesh/` folder is filled with all data of the meshes to be used. In particular, for a specific mesh (e.g. `unit_square_16`), there must be a folder containing the files `boundary.csv`, `edges.csv`, `elements.csv` and `points.csv`. In case after cloning the repo these files are not present, manually download them from GitHub and put them in the right folder.

The following tests can be run once inside the `fdaPDE-core` folder. We will denote with *shape* the general shape one can choose to compute the Voronoi of. In the code, one can choose between `unit_square_16`, `unit_square_32`,`unit_square_64`, `unit_square_128`, `quasi_circle`, `c_shaped`, `north_italy` and `brain`. When it is not specified, the default *shape* is `unit_square_16`:

- To generate a clipped Voronoi diagram and compare its area with the starting Delaunay mesh area run:

```bash
make generate_clipped_voronoi SHAPE=shape
```

The plotted Voronoi diagram will be stored in `test/plots/`, whereas the data structures encoding the built Voronoi as DCEL object are stores as `.json` files inside `test/plots/data/`.

- To check that the Voronoi diagram of *shape* satisfies the definition of geodesic Voronoi run:

```bash
make check_voronoi_definition SHAPE=shape
```

This will check the Voronoi definition on 100,000 points lying inside the domain, and outputs the percentage of points satisfying it.

- To measure the computational times required by `Voronoi` and `VoronoiClipped` to build the diagram, run:

```bash
make computational_times
```

which will run the constructor both for clipped and unclipped Voronoi over all the four meshes of `unit_square`. The output plots representing the computational cost in terms of the number of seeds can be found in `test/plots/`.

- To run all these tests in once, do:

```bash
make all SHAPE=shape
```

- To clean the directories from `.json` and `.png` files, run:

```bash
make clean
```