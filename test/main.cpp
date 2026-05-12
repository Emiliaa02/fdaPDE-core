// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// #include <gtest/gtest.h>   // testing framework
// // include eigen now to avoid possible linking errors
// #include <Eigen/Dense>
// #include <Eigen/Sparse>

/*
// utils
#include "src/scalar_field_test.cpp"
#include "src/vector_field_test.cpp"
#include "src/matrix_field_test.cpp"
#include "src/type_erasure_test.cpp"
#include "src/binary_tree_test.cpp"
// geometry
#include "src/simplex_test.cpp"
// #include "src/triangulation_test.cpp"
#include "src/point_location_test.cpp"
#include "src/kd_tree_test.cpp"
#include "src/voronoi_test.cpp"
// linear_algebra
#include "src/kronecker_product_test.cpp"
#include "src/vector_space_test.cpp"
#include "src/binary_matrix_test.cpp"
*/

// #include "src/rand_linear_algebra_test.cpp"
/*
// finite_elements
#include "src/fem_operators_test.cpp"
#include "src/fem_pde_test.cpp"
#include "src/integration_test.cpp"
#include "src/lagrangian_basis_test.cpp"
// optimization
#include "src/optimization_test.cpp"
// splines
#include "src/spline_test.cpp"
// fspai
#include "src/fspai_test.cpp"
*/

// int main(/*int argc, char** argv*/) {
//     // // start testing
//     // testing::InitGoogleTest(&argc, argv);
//     // return RUN_ALL_TESTS();

//   return 0;
// }


#include <fdaPDE/geometry.h>
#include <fdaPDE/src/geometry/voronoi.h>
#include <fdaPDE/src/geometry/delaunay.h>
#include <fdaPDE/src/geometry/triangulation.h>
#include <json.hpp>
#include <limits>
#include <chrono>
#include <cmath>
using namespace fdapde;
// using json = nlohmann::json;



int main() {

// caricate mesh dall'esterno
Triangulation<2, 2> mesh("data/mesh/quasi_circle/points.csv", "data/mesh/quasi_circle/elements.csv", "data/mesh/quasi_circle/boundary.csv", true, true);

std::cout<<"Mesh measure:"<<mesh.measure()<<std::endl;

// auto start = std::chrono::high_resolution_clock::now();
Voronoi<2, 2> voronoi_obj(mesh);
// auto end = std::chrono::high_resolution_clock::now();
// auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
// std::cout << "Computational time: " << duration.count() << " us\n";

std::cout << "Asking for nodes" << std::endl;
int n_nodes = voronoi_obj.n_nodes();

std::cout << "Asking for edges" << std::endl;
int n_edges = voronoi_obj.n_edges();

std::cout << "Asking for cells" << std::endl;
int n_cells = voronoi_obj.n_cells();


std::cout << "Voronoi has " << n_nodes << " nodes, " << n_edges << " edges and " << n_cells << " cells." << std::endl;

// fate operazioni...
int i = 0;
float sum=0;
int unbdd=0;
for(auto it = voronoi_obj.cells_begin(); it != voronoi_obj.cells_end(); ++it) {
  // std::cout << "Cell: " << it->id() << std::endl;

	auto measure = it->measure();
  // std::cout << "Computed measure" << std::endl;

  if (measure != std::numeric_limits<double>::infinity()) {
	sum+=measure;  
  i += 1;
//   std::cout << "Measure of cell with ID " << it->id() << ": " << measure << std::endl;
}
  if(it->is_unbounded()){

    unbdd+=1;
  }
}

std::cout << "Printed area of " << i << " cells" << std::endl;
std::cout << "Area " << sum <<std::endl;
std::cout << "Unbounded " << unbdd << std::endl;

std::cout << "Exporting..." << std::endl;
voronoi_obj.export_to_json("plots/data/quasi_circle.json");
std::cout << "Finished exporting" << std::endl;

return 0;
}

