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
#include <fdaPDE/src/geometry/dcel.h>
#include <fdaPDE/src/geometry/delaunay.h>
#include <fdaPDE/src/geometry/triangulation.h>
#include <json.hpp>
#include <limits>
#include <chrono>
#include <cmath>
#include <exception>
using namespace fdapde;


int main() {
    // Load a mesh
    Triangulation<2, 2> mesh("data/mesh/brain/points.csv", "data/mesh/brain/elements.csv", "data/mesh/brain/boundary.csv", true, true);
    // auto mesh = Triangulation<2, 2>::Square(0.0, 2.0, 256);

    // Print mesh measure
    std::cout<<"Mesh measure:"<<mesh.measure()<<std::endl;

    // auto start = std::chrono::high_resolution_clock::now();
    VoronoiClipped<2, 2> voronoi_obj(mesh);
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    // std::cout << "Computational time: " << duration.count() << " us\n";

    std::cout << "Voronoi has " << voronoi_obj.n_nodes() << " nodes, " << voronoi_obj.n_edges() << " edges and " << voronoi_obj.n_cells() << " cells." << std::endl;

    // Compute the measure for each Voronoi cell and the total Voronoi measure
    int i = 0;
    float sum=0;
    int unbdd=0;

    for (auto it = voronoi_obj.cells_begin(); it != voronoi_obj.cells_end(); ++it){
        if (it->is_unbounded()) {
            unbdd++;
            continue;
        }
        // compute the measure for each cell
        auto measure = it->measure();
        if (std::isfinite(measure)) {
            sum += measure;
            i++;
        }
    }

    // Print the total measure
    std::cout << "Printed area of " << i << " cells" << std::endl;
    std::cout << "Area: " << sum <<std::endl;
    std::cout << "Number of unbounded cell: " << unbdd << std::endl;

    // Export the Voronoi in a json file
    std::cout << "Exporting..." << std::endl;
    voronoi_obj.export_to_json("plots/data/brain.json");
    std::cout << "Finished exporting" << std::endl;

    return 0;
}