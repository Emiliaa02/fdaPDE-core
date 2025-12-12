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
using namespace fdapde;
// using json = nlohmann::json;



int main() {

// caricate mesh dall'esterno
// Triangulation<2, 2> mesh("data/mesh/c_shaped/points.csv", "data/mesh/c_shaped/elements.csv", "data/mesh/c_shaped/boundary.csv", true, true);

// fate operazioni...
// for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it) {
//   std::cout << it->measure() << std::endl;
// }

// Construct an object boundary for the Delaunay
std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> boundaries_entry;
Eigen::Matrix<double, Eigen::Dynamic, 2> first_side;
// Eigen::Matrix<double, Eigen::Dynamic, 2> second_side;
// Eigen::Matrix<double, Eigen::Dynamic, 2> third_side;
// Eigen::Matrix<double, Eigen::Dynamic, 2> fourth_side;
first_side.resize(4,2);
first_side(0, 0) = 0.0; first_side(0, 1) = 0.0;
first_side(1, 0) = 1.0; first_side(1, 1) = 0.0;
first_side(2, 0) = 1.0; first_side(2, 1) = 1.0;
first_side(3, 0) = 0.0; first_side(3, 1) = 1.0;

// first_side(6, 0) = 0.0; first_side(6, 1) = 1.0;
// first_side(7, 0) = 0.0; first_side(7, 1) = 0.0;
boundaries_entry.resize(1);
boundaries_entry[0] = first_side;
// boundaries_entry.push_back(second_side);
// boundaries_entry.push_back(third_side);
// boundaries_entry.push_back(fourth_side);

// DCEL<2,2> dcel_obj;
// std::vector<Eigen::Matrix<double, Eigen::Dynamic, 2>> holes;
//dcel_obj.from_triangulation(mesh, holes);

// BinaryVector<Dynamic> boundaries = mesh.boundary_nodes();
// std::cout<<"Boundary nodes"<< boundaries<<std::endl;
double min_angle = 30;
double max_area = 1;
Delaunay<2,2> delaunay_obj(boundaries_entry, min_angle, max_area, 10);
//delaunay_obj.flip();
DCEL<2,2> dcel_obj = delaunay_obj.dcel();
std::string filename = "delaunay_square";

dcel_obj.export_to_json(filename);

Triangulation<2, 2> mesh_obj = dcel_obj.to_triangulation<Triangulation<2, 2>>();
Voronoi<2,2> voronoi_obj(mesh_obj);

int n_nodes = voronoi_obj.n_nodes();
int n_edges = voronoi_obj.n_edges();
int n_cells = voronoi_obj.n_cells();


std::cout << "Voronoi has " << n_nodes << " nodes, " << n_edges << " edges and " << n_cells << " cells." << std::endl;

voronoi_obj.export_to_json("voronoi_circle.json");
 
return 0;
}

