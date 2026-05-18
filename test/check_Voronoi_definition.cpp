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

#include <fdaPDE/geometry.h>
#include <fdaPDE/src/geometry/voronoi.h>
#include <fdaPDE/src/geometry/voronoi_clipped.h>
#include <fdaPDE/src/geometry/delaunay.h>
#include <fdaPDE/src/geometry/triangulation.h>
#include <json.hpp>
#include <limits>
#include <chrono>
#include <cmath>
#include <exception>
using namespace fdapde;
using coords_t = Eigen::Matrix<double, 1, 2>;


int main() {

// Load a mesh
Triangulation<2, 2> mesh("data/mesh/north_italy/points.csv", "data/mesh/north_italy/elements.csv", "data/mesh/north_italy/boundary.csv", true, true);
// Construct the Voronoi object starting from the mesh
VoronoiClipped<2, 2> voronoi_obj(mesh);

//Ask to the user how many points he wants to sample
int n_points_grid;
std::cout<<"How many points do you want to sample inside the domain?"<<std::endl;
std::cin >> n_points_grid;

// Generate a set of random points inside the domain
Eigen::Matrix<double, Dynamic, Dynamic> sampled_points = mesh.sample(n_points_grid);
assert(sampled_points.cols() == 2);

int n_correct_ones = 0;
coords_t cur_point;
for (int i=0; i<sampled_points.rows(); ++i){
	cur_point = sampled_points.row(i);
	n_correct_ones += voronoi_obj.check_definition_on_single_point(cur_point);
}

std::cout<<"The total number of points in the grid is: "<< n_points_grid <<std::endl;
std::cout<<"The number of points that satisfy the Voronoi definition is: "<< n_correct_ones <<std::endl;

return 0;
}