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


int main(int argc, char* argv[]) {

    // Read shape if present
    if (argc < 2){
        std::cerr << "Shape needed as input" << std::endl;
        return 1;
    }
    std::string shape = argv[1];

    // Load a mesh
    Triangulation<2, 2> mesh("test/data/mesh/" + shape + "/points.csv", 
                        "test/data/mesh/" + shape + "/elements.csv", 
                        "test/data/mesh/" + shape + "/boundary.csv", true, true);

    VoronoiClipped<2, 2> voronoi_obj(mesh);

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
    std::cout << "Mesh Area: " << mesh.measure() << "; Voronoi Area: " << sum << std::endl;
    if(mesh.measure()==sum){
        std::cout << "The two areas coincide!" << std::endl;
    }

    // Create the directory data
    std::filesystem::create_directories("test/plots/data");
    // Export the Voronoi in a json file
    std::cout << "\n\n\n" << std::endl;
    std::cout << "Exporting..." << std::endl;
    voronoi_obj.export_to_json("test/plots/data/" + shape + ".json");
    std::cout << "Finished exporting" << std::endl;

    return 0;
}