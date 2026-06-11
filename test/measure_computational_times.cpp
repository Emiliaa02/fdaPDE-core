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
#include <filesystem>
#include <exception>
using namespace fdapde;
namespace fs = std::filesystem;
using coords_t = Eigen::Matrix<double, 1, 2>;

int main(){
    std::cout << "\n\n\n" << std::endl;
    std::cout << "Measuring computational times..." << std::endl;
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

    // Create the directory "computational times" if it does not exist
    fs::path dir = "test/data/computational_times";

    if (!fs::exists(dir)) {
        fs::create_directories(dir);
    }

    // Files definition
    std::ofstream file_voronoi("test/data/computational_times/voronoi.txt");
    std::ofstream file_voronoi_clipped("test/data/computational_times/voronoi_clipped.txt");

    if (!file_voronoi || !file_voronoi_clipped) {
        std::cerr << "Failed to open output files.\n";
        return 1;
    }

    // VORONOI

    std::chrono::steady_clock::duration duration_1, duration_2, duration_3, duration_4;

    // Load first mesh
    {
        Triangulation<2, 2> mesh(
            "test/data/mesh/unit_square_16/points.csv",
            "test/data/mesh/unit_square_16/elements.csv",
            "test/data/mesh/unit_square_16/boundary.csv",
            true, true);

        start = std::chrono::steady_clock::now();
        Voronoi<2, 2> voronoi_obj(mesh);
        end = std::chrono::steady_clock::now();
        duration_1 = end - start;
    }

    // Load second mesh
    {
        Triangulation<2, 2> mesh(
            "test/data/mesh/unit_square_32/points.csv",
            "test/data/mesh/unit_square_32/elements.csv",
            "test/data/mesh/unit_square_32/boundary.csv",
            true, true);

        start = std::chrono::steady_clock::now();
        Voronoi<2, 2> voronoi_obj(mesh);
        end = std::chrono::steady_clock::now();
        duration_2 = end - start;
    }

    // Load third mesh
    {
        Triangulation<2, 2> mesh(
            "test/data/mesh/unit_square_64/points.csv",
            "test/data/mesh/unit_square_64/elements.csv",
            "test/data/mesh/unit_square_64/boundary.csv",
            true, true);

        start = std::chrono::steady_clock::now();
        Voronoi<2, 2> voronoi_obj(mesh);
        end = std::chrono::steady_clock::now();
        duration_3 = end - start;
    }

    // Load fourth mesh
    {
        Triangulation<2, 2> mesh(
            "test/data/mesh/unit_square_128/points.csv",
            "test/data/mesh/unit_square_128/elements.csv",
            "test/data/mesh/unit_square_128/boundary.csv",
            true, true);

        start = std::chrono::steady_clock::now();
        Voronoi<2, 2> voronoi_obj(mesh);
        end = std::chrono::steady_clock::now();
        duration_4 = end - start;
    }

    file_voronoi
        << std::chrono::duration_cast<std::chrono::microseconds>(duration_1).count() << '\n'
        << std::chrono::duration_cast<std::chrono::microseconds>(duration_2).count() << '\n'
        << std::chrono::duration_cast<std::chrono::microseconds>(duration_3).count() << '\n'
        << std::chrono::duration_cast<std::chrono::microseconds>(duration_4).count() << '\n';

    // VORONOI CLIPPED

    // Load first mesh
    {
        Triangulation<2, 2> mesh(
            "test/data/mesh/unit_square_16/points.csv",
            "test/data/mesh/unit_square_16/elements.csv",
            "test/data/mesh/unit_square_16/boundary.csv",
            true, true);

        start = std::chrono::steady_clock::now();
        VoronoiClipped<2, 2> voronoi_obj(mesh);
        end = std::chrono::steady_clock::now();
        duration_1 = end - start;
    }

    // Load second mesh
    {
        Triangulation<2, 2> mesh(
            "test/data/mesh/unit_square_32/points.csv",
            "test/data/mesh/unit_square_32/elements.csv",
            "test/data/mesh/unit_square_32/boundary.csv",
            true, true);

        start = std::chrono::steady_clock::now();
        VoronoiClipped<2, 2> voronoi_obj(mesh);
        end = std::chrono::steady_clock::now();
        duration_2 = end - start;
    }

    // Load third mesh
    {
        Triangulation<2, 2> mesh(
            "test/data/mesh/unit_square_64/points.csv",
            "test/data/mesh/unit_square_64/elements.csv",
            "test/data/mesh/unit_square_64/boundary.csv",
            true, true);

        start = std::chrono::steady_clock::now();
        VoronoiClipped<2, 2> voronoi_obj(mesh);
        end = std::chrono::steady_clock::now();
        duration_3 = end - start;
    }

    // Load fourth mesh
    {
        Triangulation<2, 2> mesh(
            "test/data/mesh/unit_square_128/points.csv",
            "test/data/mesh/unit_square_128/elements.csv",
            "test/data/mesh/unit_square_128/boundary.csv",
            true, true);

        start = std::chrono::steady_clock::now();
        VoronoiClipped<2, 2> voronoi_obj(mesh);
        end = std::chrono::steady_clock::now();
        duration_4 = end - start;
    }

    file_voronoi_clipped
        << std::chrono::duration_cast<std::chrono::microseconds>(duration_1).count() << '\n'
        << std::chrono::duration_cast<std::chrono::microseconds>(duration_2).count() << '\n'
        << std::chrono::duration_cast<std::chrono::microseconds>(duration_3).count() << '\n'
        << std::chrono::duration_cast<std::chrono::microseconds>(duration_4).count() << '\n';

    std::cout << "Finished measuring computational times..." << std::endl;

    return 0;
}