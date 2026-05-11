#ifndef __FDAPDE_VORONOI_H__
#define __FDAPDE_VORONOI_H__


namespace fdapde{

template <int LocalDim, int EmbedDim>
class Voronoi {

    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;

    using dcel_t = DCEL<local_dim, embed_dim>;
    using simplex_t = Simplex<local_dim, embed_dim>;
    using triangulation_t = Triangulation<local_dim, embed_dim>;
    using triangle_t = Triangle<triangulation_t>;
    using cell_t = dcel_t::cell_t;
    using coords_t = Eigen::Matrix<double, 1, embed_dim>;
    
    public:

    Voronoi(const Triangulation<local_dim, embed_dim>& mesh) {

    // Number of faces in the mesh 
    int n_mesh_faces = mesh.n_cells();
    // Number of vertices in the mesh
    int n_mesh_vertices = mesh.n_nodes();
    // Centroids lookup
    std::vector<typename simplex_t::NodeType> centroid_lookup(n_mesh_faces+1);
    // Visited or not
    std::vector<bool> visited_centroids(n_mesh_faces+1, false);
    // Old node ID : new node ID
    std::vector<int> old2new(n_mesh_faces);
    // counter to set the IDs
    int counter = 0;
    // Same centroids lookup
    std::vector<int> first_id_lookup(n_mesh_faces);
    // Lookup for midpoints coordinates
    std::map<int, typename simplex_t::NodeType> midpoints_lookup;
    std::map<int, typename simplex_t::NodeType> midpoints_lookup_raw;
    // Map that associates a pair of centroids to the halfedge that connects them
    std::multimap<std::pair<int, int>, typename dcel_t::halfedge_t*> vertexes2halfedge;
    int cur_id = 0;
    // Structure to store cells without finite halfedges
    // std::list<int> not_created_cells;
    std::set<int> not_created_cells;
    // Set to store d_centroids out of the boundary
    std::set<int> out_of_boundary_d_centroids;
    // Map to store infinity halfedge -> midpoint ID (for understanding infinity edges intersections)
    std::map<int, int> infty_halfedges_midpoints;

    std::vector<std::vector<int>> node_neighbors_lookup = v_vertex_computation_(mesh, 
                                                                                dcel_, 
                                                                                old2new, 
                                                                                centroid_lookup,
                                                                                midpoints_lookup,
                                                                                midpoints_lookup_raw,
                                                                                cur_id,
                                                                                out_of_boundary_d_centroids);


    // Define infinity v_vertex and add it to DCEL
    int infty_id = cur_id;
    typename simplex_t::NodeType v_vertex_infty = simplex_t::NodeType::Constant(std::numeric_limits<double>::infinity());
    centroid_lookup[infty_id] = v_vertex_infty;
    typename dcel_t::node_t* infty_node = dcel_.insert_node(typename dcel_t::node_t(infty_id, false, v_vertex_infty));

    auto start = std::chrono::high_resolution_clock::now();
    // Imagine to have the correspondence cell_id: centroid coordinates
    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it) {  // O(N)

        // Retrieve the old ID of the current cell
        int old_cell_id = it->id();
        // Retrieve new ID of the current cell
        int cell_id = old2new.at(old_cell_id);  // O(logN)
        // Order the neighbouring cells in a counter-clockwise way
        std::vector<int> ordered_neigh_ids = order_neighbours_cclw_(cell_id, node_neighbors_lookup, old2new, 
                                                                    centroid_lookup, midpoints_lookup);  // O(1)
        
        // Create and connect the halfedges
        if(!visited_centroids[cell_id]){
            create_and_connect_halfedges_(dcel_, cell_id, infty_id, infty_node, centroid_lookup, 
                                        midpoints_lookup_raw, midpoints_lookup, ordered_neigh_ids, vertexes2halfedge, infty_halfedges_midpoints);  // O(logN)
            // Set the centroid as visited
            visited_centroids[cell_id] = true;
        }

        auto neighbor_simplexes = it->neighbors();
        // Create cells
        create_cells_(dcel_, mesh, cell_id, old_cell_id, old2new,
                    centroid_lookup, neighbor_simplexes,vertexes2halfedge, not_created_cells);  // O(N)
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout<<"computational times: "<< duration.count() << std::endl;
    // Add to the DCEL structure the not created cells
    add_not_created_cells_(infty_id, not_created_cells, mesh, old2new, vertexes2halfedge);

    // Connect the infinity halfedges
    connect_infty_halfedges_(infty_id);

    // Update halfedges-cells structure in DCEL
    this->dcel_.update_halfedges_with_cells();

    std::map<int, std::vector<typename dcel_t::node_t*>> supplementary_points;
    std::set<int> on_boundary;

    // PER PROVARE SE RUNNA
    clip_voronoi(mesh, infty_id, out_of_boundary_d_centroids, old2new, vertexes2halfedge, midpoints_lookup, infty_halfedges_midpoints);

    std::cout << "\nFinished constructor" << std::endl;

    }


    // iterators
    using cell_iterator = std::list<cell_t>::iterator;
    using const_cell_iterator = std::list<cell_t>::const_iterator;

    cell_iterator cells_begin() { return this->dcel_.cells_begin(); }
    cell_iterator cells_end() { return this->dcel_.cells_end(); }

    const_cell_iterator cells_cbegin() { return this->dcel_.cells_cbegin(); }
    const_cell_iterator cells_cend() { return this->dcel_.cells_cend(); }

    // observers
    // matrix of coordinates (n_nodes x EmbedDim)
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes() const {
        // Simply calls nodes() from dcel_t
        return dcel_.nodes();
    }

    int n_nodes() { return dcel_.n_nodes(); }

    // matrix of int (n_edge x 2) where i-th row: [id_node_1 i-th edge, id_node_2 i-th edge]
    Eigen::Matrix<int, Eigen::Dynamic, 2> edges() const {
        // calls edges() from dcel_t
        return dcel_.edges();
    }

    int n_edges() { return dcel_.n_edges(); }

    int n_cells() { return dcel_.n_cells(); }

    void export_to_json(const std::string& filename){
        this->dcel_.export_to_json(filename);
    }

    
    bool areAligned(const coords_t& A,
                const coords_t& B,
                const coords_t& C)
            {
                float x1 = A(0), y1 = A(1);
                float x2 = B(0), y2 = B(1);
                float x3 = C(0), y3 = C(1);

                const float EPS = 1e-6f;

                // 1. Check collinearity (cross product)
                float cross = (x2 - x1) * (y3 - y1) -
                            (y2 - y1) * (x3 - x1);

                if (std::fabs(cross) > EPS)
                    return false;

                // 2. Check if B is between A and C (dot product)
                float dot = (x2 - x1) * (x2 - x3) +
                            (y2 - y1) * (y2 - y3);

                return (dot <= EPS);
            }


    // std::map<int, int> n_intersections_halfedges_boundary(std::set<int>& out_of_boundary_d_centroids, int infty_id){

    //     std::map<int, int> n_intersections_halfedges_boundary;
    //     for(auto it = this->dcel_.halfedges_begin(); it != this->dcel_.halfedges_end(); ++it){
    //         int count = 0;
    //         if(it->node()->id() == infty_id || out_of_boundary_d_centroids.find(it->node()->id()) != out_of_boundary_d_centroids.end()){
    //             count += 1;
    //         }
    //         if(it->twin()->node()->id() == infty_id || out_of_boundary_d_centroids.find(it->twin()->node()->id()) != out_of_boundary_d_centroids.end()){
    //             count += 1;
    //         }
    //         n_intersections_halfedges_boundary[it->id()] = count;
    //     }
    //     return n_intersections_halfedges_boundary;
    // }

    void clip_voronoi(const Triangulation<local_dim, embed_dim>& mesh, int infty_id,
                      std::set<int>& out_of_boundary_d_centroids, std::vector<int>& d_centroid2v_vertex,
                      std::multimap<std::pair<int, int>, typename dcel_t::halfedge_t*>& vertexes2halfedge,
                      std::map<int, typename simplex_t::NodeType>&  midpoints_lookup,
                      std::map<int, int>& infty_halfedges_midpoints){

        std::map<int, std::vector<typename dcel_t::node_t*>> supplementary_points;
        std::map<std::pair<double, double>, std::vector<coords_t>> adjacent_points_map = adjacent_points(mesh);
        std::map<int, typename dcel_t::node_t*> cell_centroid_coords;
        std::set<int> on_boundary;
        std::vector<int> node_ids_to_remove;

        boundary_cells_detection(mesh, supplementary_points, on_boundary,
                                out_of_boundary_d_centroids, d_centroid2v_vertex,
                                infty_id, node_ids_to_remove, cell_centroid_coords, infty_halfedges_midpoints, midpoints_lookup);

        // Loop over cells of the DCEL structure
        for (auto cur_cell_it = this->dcel_.cells_begin(); cur_cell_it != this->dcel_.cells_end(); cur_cell_it++){
            if(on_boundary.find((*cur_cell_it).id()) == on_boundary.end()){ continue; };
            typename dcel_t::halfedge_t* start = nullptr;
            typename dcel_t::halfedge_t* cur_halfedge = nullptr;
            typename dcel_t::halfedge_t* prev_halfedge = nullptr;

            // Create a list of its points
            std::list<typename dcel_t::node_t*> pts;
            // Get nodes
            std::vector<typename dcel_t::node_t*> nodes = cur_cell_it -> cell_nodes();
            for (auto cell_node = nodes.cbegin(); cell_node != nodes.cend(); cell_node++){
                auto to_exclude_it = std::find(node_ids_to_remove.begin(), node_ids_to_remove.end(), (*cell_node)->id());
                if(to_exclude_it == node_ids_to_remove.end()){
                    pts.push_back(*cell_node);
                }
            }

            // Get supplementary points
            auto check_supp_pt = supplementary_points.find(cur_cell_it->id());
            if (check_supp_pt != supplementary_points.end()){
                // Insert supplementary points in DCEL structure, if they are not yet included 
                for(typename dcel_t::node_t* supp_pt : check_supp_pt->second){
                    auto node_ptr = this->dcel_.find_node_by_coords(supp_pt->coords());
                    if(node_ptr == nullptr){
                        typename dcel_t::node_t* inserted_node = dcel_.insert_node(*supp_pt);
                        pts.push_back(inserted_node);
                    }
                    else{pts.push_back(node_ptr);};     
                }
            }

            // Compute the center of mass to perform the ordering
            typename simplex_t::NodeType center_of_mass = compute_center_of_mass(pts);

            // Put in counterclockwise order the points
            pts.sort([&](dcel_t::node_t* pa, dcel_t::node_t* pb)
                {
                    double ang_a = std::atan2(pa->coords()(1) - center_of_mass(1),
                                            pa->coords()(0) - center_of_mass(0));

                    double ang_b = std::atan2(pb->coords()(1) - center_of_mass(1),
                                            pb->coords()(0) - center_of_mass(0));

                    return ang_a < ang_b;
                }
            );

            // Check if the centroid has to be included into the polygon
            auto centroid_it = cell_centroid_coords.find((*cur_cell_it).id());
            if(centroid_it != cell_centroid_coords.end()){
                // std::cout << "Found centroid" << std::endl;
                double centroid_x = centroid_it-> second -> coords()(0);
                double centroid_y = centroid_it-> second -> coords()(1);
                std::pair<double, double> centroid_xy = std::make_pair(centroid_x, centroid_y);
                if(adjacent_points_map.find(centroid_xy) == adjacent_points_map.end()){
                    std::cout<<"row 263"<<std::endl;
                }
                std::vector<coords_t> adjacent_nodes_on_boundary = adjacent_points_map.at(centroid_xy);
                std::vector<typename dcel_t::node_t*> on_boundary_pts;

                for (auto pt_: pts){
                    for (coords_t ad_node : adjacent_nodes_on_boundary){
                        // std::cout << "Printing points: " << std::endl;
                        // std::cout << ad_node << std::endl;
                        // std::cout << pt_->coords() << std::endl;
                        // std::cout << centroid_it-> second -> coords() << std::endl;
                        // std::cout << "Condition on dot: " << areAligned(ad_node, pt_->coords(), centroid_it-> second -> coords()) << std::endl;
                        // std::cout << "Condition 2a: " << (ad_node(0) != centroid_x) << std::endl;
                        // std::cout << "Condition 2b: " << (ad_node(0) != centroid_x) << std::endl;
                        if((ad_node(0) != centroid_x | ad_node(1) != centroid_y) & areAligned(ad_node, pt_->coords(), centroid_it-> second -> coords())){
                            on_boundary_pts.push_back(pt_);
                            break;
                        }
                    }
                }


                auto first_position  = std::find(pts.begin(), pts.end(), on_boundary_pts[0]);
                auto second_position = std::find(pts.begin(), pts.end(), on_boundary_pts[1]);

                if (std::distance(pts.begin(), first_position) >
                    std::distance(pts.begin(), second_position))
                {
                    std::swap(first_position, second_position);
                }

                if (first_position != pts.end() && second_position != pts.end()){
                    if(first_position == pts.begin() && second_position == std::prev(pts.end())){
                        if(cell_centroid_coords.find((*cur_cell_it).id()) == cell_centroid_coords.end()){
                            std::cout<<"row 297"<<std::endl;
                        }
                        pts.push_back(cell_centroid_coords.at((*cur_cell_it).id()));
                        typename dcel_t::node_t* inserted_node = dcel_.insert_node(*cell_centroid_coords.at((*cur_cell_it).id()));
                    }
                    else{
                        if(cell_centroid_coords.find((*cur_cell_it).id()) == cell_centroid_coords.end()){
                            std::cout<<"row 304"<<std::endl;
                        }
                        pts.insert(std::next(first_position), cell_centroid_coords.at((*cur_cell_it).id()));
                        typename dcel_t::node_t* inserted_node = dcel_.insert_node(*cell_centroid_coords.at((*cur_cell_it).id()));
                    }
                }
            }
            
            if(pts.size() < 2){
                std::cout << "LESS THAN 2 POINTS" << std::endl;
                continue;
            }

            // Loop points one after the other
            for(auto pt = pts.begin(); pt != pts.end(); ++pt){
                // Next point
                auto pt_next = std::next(pt);
                if(pt_next == pts.end()) {
                    pt_next = pts.begin();
                }
                // Look for cur->next halfedge
                auto search_it = vertexes2halfedge.find({(*pt)->id(), (*pt_next)->id()});
                // If cur -> next halfedge already exists, retrieve it...
                if (search_it != vertexes2halfedge.end()){
                    cur_halfedge = search_it -> second;
                    cur_cell_it -> set_halfedge(cur_halfedge);
                }
                // ... otherwise create it
                else{
                    // Create new halfedge
                    typename dcel_t::halfedge_t* new_halfedge = dcel_.emplace_halfedge(*pt);
                    // Set it as halfedge for pt
                    (*pt) -> set_halfedge(new_halfedge);
                    // Add halfedge to the node_t structure
                    (*pt) -> add_halfedge((*pt_next)->id(), new_halfedge);

                    // If its twin already exists, connect them
                    auto twin_it = vertexes2halfedge.find({(*pt_next)->id(), (*pt)->id()});
                    if(twin_it != vertexes2halfedge.end()){
                        typename dcel_t::halfedge_t* twin_halfedge = twin_it->second;
                        new_halfedge -> set_twin(twin_halfedge);
                        twin_halfedge -> set_twin(new_halfedge);
                    }
                    vertexes2halfedge.insert({std::make_pair((*pt)->id(), (*pt_next)->id()), new_halfedge});

                    // Now it is the current halfedge
                    cur_halfedge = new_halfedge;
                    cur_cell_it -> set_halfedge(cur_halfedge);
                }

                // If this is the first loop, initialize start
                if(pt == pts.begin()){
                    start = cur_halfedge;
                }
                // If previous halfedge is not null, connect them
                if (prev_halfedge){
                    prev_halfedge -> set_next(cur_halfedge);
                    cur_halfedge -> set_prev(prev_halfedge);
                }
                // Previous halfedge becomes current halfedge
                prev_halfedge = cur_halfedge;
        }
        if (start){
            cur_halfedge -> set_next(start);
            start -> set_prev(cur_halfedge);
        }

        // Declare cell as clipped
        cur_cell_it->clipped();
    }
    // Remove the nodes outside the domain and the corrensponding halfedges
    this->dcel_.remove_nodes(node_ids_to_remove);
}   


    // What the following function does is understand which v_cells are on the boundary
    // and compile supplementary_points[v_cell_id]: for each v_cell, it assigns, if it is on boundary,
    // the supplementary points we will have to consider in order to have it clipped.
    // In the clipping phase, all we have to do is:
    // - Order supplementary points + already existing points in counterclockwise order;
    // - Following the order, if halfedge n1->n2 already exists, OK
    // - Else, create n1->n2, removing the appropriate halfedges in the DCEL
    // void boundary_cells_detection(const Triangulation<local_dim, embed_dim>& mesh, 
    //                             std::map<int, std::vector<typename dcel_t::node_t*>>& supplementary_points,
    //                             std::set<int>& on_boundary,
    //                             std::set<int>& out_of_boundary_d_centroids,
    //                             std::vector<int>& d_centroid2v_vertex,
    //                             int infty_id, std::vector<int>& node_ids_to_remove,
    //                             std::map<int, typename dcel_t::node_t*>& cell_centroid_coords){

    //     std::set<int> entered_first_time;
    //     int counter = infty_id + 1;
    //     node_ids_to_remove.push_back(infty_id);
    //     for(int d_centroid_id : out_of_boundary_d_centroids){
    //         int v_vertex_id = d_centroid2v_vertex.at(d_centroid_id);
    //         node_ids_to_remove.push_back(v_vertex_id);
    //     }

    //     // Compute the v_cell2halfedges structure
    //     std::map<int, std::vector<int>> v_cell_halfedges = compute_v_cell2halfedges_(infty_id);
    //     // Loop over mesh boundary edges of the Delaunay
    //     for (auto d_boundary_edge = mesh.boundary_edges_begin(); d_boundary_edge != mesh.boundary_edges_end(); ++d_boundary_edge){
    //         // Find closest v_centroid to d_edge_midpoint
    //         auto current_midpoint = d_boundary_edge->compute_midpoint();
    //         auto closest_v_centroid = find_closest_v_centroid(mesh, current_midpoint);

    //         // Initialize list of v_cells to check
    //         std::list<typename dcel_t::cell_t*> to_check_v_cells{closest_v_centroid};
    //         // Initialize map of checked halfedges
    //         std::map<int, bool> already_checked_halfedges;

    //         // While list to check is not empty
    //         while (to_check_v_cells.size() > 0){
    //             typename dcel_t::cell_t* cur_v_centroid = to_check_v_cells.front(); 
    //             to_check_v_cells.pop_front();
    //             std::vector<int> cur_v_cell_halfedges = v_cell_halfedges.at(cur_v_centroid->id());

    //             // Loop over cur_v_cell halfedges
    //             for(int cur_halfedge_id: cur_v_cell_halfedges){                  
    //                 typename dcel_t::halfedge_t* cur_halfedge;
    //                 bool found_halfedge = false;
                
    //                 for(auto it = this->dcel_.halfedges_begin(); it != this->dcel_.halfedges_end(); ++it){
    //                     if(it->id() == cur_halfedge_id){
    //                         cur_halfedge = &(*it);
    //                         found_halfedge = true;
    //                         break;
    //                     }
    //                 }
    //                 // Check if the current halfedge exists
    //                 if (!found_halfedge) {
    //                     throw std::runtime_error("Halfedge not found");
    //                 }

    //                 // Check if the twin halfedge exists
    //                 auto twin = cur_halfedge->twin();
    //                 if (twin == nullptr) {
    //                     throw std::runtime_error("Twin is null");
    //                 }
    //                 if (twin->cell() == nullptr) {
    //                     throw std::runtime_error("Twin cell is null");
    //                 }

    //                 // If we already know halfedge intersects boundary edge (necessarily on the midpoint)...
    //                 auto check_it = already_checked_halfedges.find(cur_halfedge_id);
    //                 if((check_it == already_checked_halfedges.end() || !check_it->second) && intersection_d_edges_.find(cur_halfedge_id) != intersection_d_edges_.end()
    //                    && intersection_d_edges_.at(cur_halfedge_id) == d_boundary_edge->id()){
                        
    //                     // Cur_v_cell is a boundary cell
    //                     on_boundary.insert(cur_v_centroid->id());

    //                     // Store the site of the Voronoi cell only one time
    //                     if(entered_first_time.find(cur_v_centroid->id()) == entered_first_time.end()){
    //                         typename simplex_t::NodeType centroid_coords;
    //                         centroid_coords(0) = mesh.node(cur_v_centroid->id())(0);
    //                         centroid_coords(1) = mesh.node(cur_v_centroid->id())(1);
    //                         // Check if the current site is already a DCEL node
    //                         auto existing_node = dcel_.find_node_by_coords(centroid_coords);
    //                         if (existing_node != nullptr) {
    //                             cell_centroid_coords[cur_v_centroid->id()] = existing_node;
    //                             // supplementary_points[cur_v_centroid.id()].push_back(existing_node);
    //                         } else {
    //                             typename dcel_t::node_t* v_centroid_coords = new typename dcel_t::node_t();
    //                             v_centroid_coords->set_coords(centroid_coords);
    //                             v_centroid_coords->set_id(counter++);
    //                             cell_centroid_coords[cur_v_centroid->id()] = v_centroid_coords;
    //                             // supplementary_points[cur_v_centroid->id()].push_back(v_centroid_coords);
    //                         }
    //                     }  
    //                     entered_first_time.insert(cur_v_centroid->id());

    //                     // Add to cur_v_cell the midpoint
    //                     if(supplementary_points.find(cur_v_centroid->id()) == supplementary_points.end()){
    //                         supplementary_points[cur_v_centroid->id()];
    //                     }
    //                     typename dcel_t::node_t* v_current_midpoint;
    //                     auto node_ptr = dcel_.find_node_by_coords(current_midpoint);
    //                     if(node_ptr == nullptr){
    //                         v_current_midpoint = new typename dcel_t::node_t();
    //                         v_current_midpoint->set_coords(current_midpoint);
    //                         v_current_midpoint->set_id(counter++);
    //                     }
    //                     else{
    //                         v_current_midpoint = node_ptr;
    //                     }
    //                     v_current_midpoint -> set_boundary(true);
    //                     supplementary_points[cur_v_centroid->id()].push_back(v_current_midpoint);

    //                     // Add to list neighbouring v_cell
    //                     to_check_v_cells.push_back(cur_halfedge->twin()->cell());

    //                     already_checked_halfedges[cur_halfedge_id] = true;
                        
    //                 }
    //                 // ... otherwise
    //                 else if (check_it == already_checked_halfedges.end() || !check_it->second){
    //                     // Check whether halfedge and boundary edge intersect in a point
    //                     bool check_intersection= false;
    //                     typename simplex_t::NodeType intersection_point;
    //                     find_intersection(mesh, *cur_halfedge, *d_boundary_edge, check_intersection, intersection_point);

    //                     // If you found it...
    //                     if(check_intersection){ 
    //                         typename dcel_t::node_t* v_intersection_point;
    //                         auto intersection_ptr =  dcel_.find_node_by_coords(intersection_point);
    //                         if(intersection_ptr == nullptr){
    //                             v_intersection_point = new typename dcel_t::node_t();
    //                             v_intersection_point->set_coords(intersection_point);
    //                             v_intersection_point->set_id(counter++);
    //                         }
    //                         else{
    //                             v_intersection_point = intersection_ptr;
    //                         }
    //                         // cur_v_cell is boundary cell
    //                         on_boundary.insert(cur_v_centroid->id());

    //                         // Add to cur_v_cell the point
    //                         if(supplementary_points.find(cur_v_centroid->id()) == supplementary_points.end()){
    //                             supplementary_points[cur_v_centroid->id()];
    //                         }
    //                         supplementary_points[cur_v_centroid->id()].push_back(v_intersection_point);
    //                         v_intersection_point->set_boundary(true);

    //                         // Add to list neighbouring v_cell
    //                         to_check_v_cells.push_back(cur_halfedge->twin()->cell());
    //                         already_checked_halfedges[cur_halfedge_id] = true;
    //                     }
    //                 }
                
    //             }
                
    //         }
    //     }
    //     std::cout << "Finished" << std::endl;

    // }










//////////////////// NUOVA IMPLEMENTAZIONE //////////////////////////////

    void boundary_cells_detection(const Triangulation<local_dim, embed_dim>& mesh, 
                                std::map<int, std::vector<typename dcel_t::node_t*>>& supplementary_points,
                                std::set<int>& on_boundary,
                                std::set<int>& out_of_boundary_d_centroids,
                                std::vector<int>& d_centroid2v_vertex,
                                int infty_id, std::vector<int>& node_ids_to_remove,
                                std::map<int, typename dcel_t::node_t*>& cell_centroid_coords,
                                std::map<int, int>& infty_halfedges_midpoints,
                                std::map<int, typename simplex_t::NodeType>& midpoints_lookup){

        node_ids_to_remove.push_back(infty_id);
        for(int d_centroid_id : out_of_boundary_d_centroids){
            int v_vertex_id = d_centroid2v_vertex.at(d_centroid_id);
            node_ids_to_remove.push_back(v_vertex_id);
        }

        // Take one boundary d_edge 
        auto d_boundary_edge = mesh.boundary_edges_begin();
        if (mesh.boundary_edges_begin() == mesh.boundary_edges_end()){
            std::cerr << "Mesh does not have boundary edges, maybe there is no need for clipping..." << std::endl;
        }

        // Consider one of its two endpoints, and the corresponding v_cell
        Eigen::Matrix<int, Dynamic, 1> node_ids = d_boundary_edge->node_ids();
        int d_vertex_id = node_ids(0);
        Eigen::Matrix<double, 2, 1> d_vertex = mesh.node(d_vertex_id);
        if(cells_.find(d_vertex_id) == cells_.end()){
            std::cout<<"row 582"<<std::endl;
        }
        typename dcel_t::cell_t v_cell = cells_.at(d_vertex_id);

        // Set counter for IDs
        int counter = infty_id + 1;

        // Cells queue
        std::list<int> to_check_v_cells;

        // Already checked data structure
        std::set<std::tuple<int, int, int>> already_checked;

        // Initialize to_check_v_cells
        to_check_v_cells.push_back(d_vertex_id);


        // While queue is not empty
        while (to_check_v_cells.size() > 0){

            // Get v_cell and d_vertex ID
            d_vertex_id = to_check_v_cells.front();
            d_vertex = mesh.node(d_vertex_id);
            if(cells_.find(d_vertex_id) == cells_.end()){
                std::cout<<"row 606"<<std::endl;
            }
            v_cell = cells_.at(d_vertex_id); 
            to_check_v_cells.pop_front();

            // Neighbouring d_vertexes
            std::vector<int> neigh_d_vertexes_ids =  mesh.node_one_ring(d_vertex_id);

            // Loop over d_edges connected to that d_vertex
            for (int neigh_d_vertex_id : neigh_d_vertexes_ids){

                coords_t neigh_coords = mesh.node(neigh_d_vertex_id);
                std::vector<typename dcel_t::halfedge_t*> v_cell_halfedges =  v_cell.cell_edges();
                bool intersected_once = false;
                // Loop over v_halfedges connected to that v_centroid
                for (typename dcel_t::halfedge_t* v_cell_halfedge : v_cell_halfedges){

                    // If the triple (v_centroid, neigh_d_vertex, v_halfedge) is already present, skip
                    if(already_checked.find(std::make_tuple(d_vertex_id, neigh_d_vertex_id, v_cell_halfedge->id())) != already_checked.end()){
                        continue;
                    }

                    // Find halfedge points
                    coords_t he_p1, he_p2;
                    if (v_cell_halfedge->node()->id()==infty_id){
                        if(infty_halfedges_midpoints.find(v_cell_halfedge->id())== infty_halfedges_midpoints.end()){
                            std::cout<<"row 632"<<std::endl;
                        }
                        if(midpoints_lookup.find(infty_halfedges_midpoints.at(v_cell_halfedge->id()))== midpoints_lookup.end()){
                            std::cout<<"row 635"<<std::endl;
                        }
                        he_p1 = midpoints_lookup.at(infty_halfedges_midpoints.at(v_cell_halfedge->id()));
                        he_p2 = v_cell_halfedge->twin()->node()->coords();
                    }
                    else if (v_cell_halfedge->twin()->node()->id()==infty_id){
                        if(infty_halfedges_midpoints.find(v_cell_halfedge->id())== infty_halfedges_midpoints.end()){
                            std::cout<<"row 2"<<std::endl;
                        }
                        he_p1 = v_cell_halfedge->node()->coords();
                        he_p2 = midpoints_lookup.at(infty_halfedges_midpoints.at(v_cell_halfedge->id()));
                    }
                    else{
                        he_p1 = v_cell_halfedge->node()->coords();
                        he_p2 = v_cell_halfedge->twin()->node()->coords();
                    }

                    // Useful structures
                    bool check_intersection = false;
                    typename simplex_t::NodeType intersection_point;

                    // Check intersection
                    find_intersection_by_points(
                        d_vertex,
                        neigh_coords,
                        he_p1,
                        he_p2,
                        check_intersection,
                        intersection_point
                    );

                    // If you find hafedge-edge intersection, add it to v_cell supplementary_points and twin_v_cell supplementary points
                    if (check_intersection){
                        on_boundary.insert(d_vertex_id);
                        typename dcel_t::node_t* v_intersection_point;
                        v_intersection_point = new typename dcel_t::node_t();
                        v_intersection_point->set_coords(intersection_point);
                        v_intersection_point->set_id(counter++);
                        if(supplementary_points.find(d_vertex_id) == supplementary_points.end()){
                            supplementary_points[d_vertex_id];
                        }
                        if(supplementary_points.find(neigh_d_vertex_id) == supplementary_points.end()){
                            supplementary_points[neigh_d_vertex_id];
                        }
                        supplementary_points[d_vertex_id].push_back(v_intersection_point);
                        supplementary_points[neigh_d_vertex_id].push_back(v_intersection_point);
                        v_intersection_point->set_boundary(true);
                        intersected_once = true;
                    }

                    // Add to already_checked the triple (v_centroid, neigh_d_vertex, v_halfedge) and corresponding triple for the neighbouring v_cell
                    already_checked.insert({d_vertex_id, neigh_d_vertex_id, v_cell_halfedge->id()});
                    already_checked.insert({neigh_d_vertex_id, d_vertex_id, v_cell_halfedge->twin()->id()});
                }
                if(intersected_once){
                    // Add to queue twin_v_cell
                    to_check_v_cells.push_back(neigh_d_vertex_id);
                }
            }
        }

        // Here you might do a loop over all cells: if their centroid is at the boundary and is not already present as a node in the DCEL, add it and set that it is a boundary node
        for(auto cell_it = dcel_.cells_begin(); cell_it != dcel_.cells_end(); ++cell_it){
            int cur_cell_id = cell_it -> id();
            auto cur_centroid_coords = mesh.node(cur_cell_id);
            typename dcel_t::node_t* new_centroid_node = dcel_.find_node_by_coords(cur_centroid_coords);
            if(new_centroid_node == nullptr){
                new_centroid_node = new typename dcel_t::node_t();
                new_centroid_node->set_coords(cur_centroid_coords);
                new_centroid_node->set_id(counter++);
            }
            if(mesh.is_node_on_boundary(cur_cell_id)){
                new_centroid_node->set_boundary(true);
                cell_centroid_coords[cur_cell_id] = new_centroid_node;
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////


    typename dcel_t::cell_t* find_closest_v_centroid(const Triangulation<local_dim, embed_dim>& mesh, 
                                                    Eigen::Matrix<double, 2, 1> midpoint){
            
        // Phisical coordinates of mesh vertices
        Eigen::Matrix<double, Dynamic, Dynamic> nodes = mesh.nodes();

        int best_v_centroid_id = 0;
        double best_distance = euclidean_distance(midpoint(0), midpoint(1), nodes(0,0), nodes(0,1));

        for(int i=1; i<nodes.rows(); ++i){
            double current_distance = euclidean_distance(midpoint(0), midpoint(1), nodes(i,0), nodes(i,1));
            if(current_distance < best_distance){
                best_distance = current_distance;
                best_v_centroid_id = i;
            }
        }

        typename dcel_t::cell_t* best_cell;

        for(auto it = this->dcel_.cells_begin(); it != this->dcel_.cells_end(); ++it){
            if(it->id() == best_v_centroid_id){
                best_cell = (*it);
                return best_cell;
            }
        }
        throw std::runtime_error("Best cell ID is not included among Voronoi cell IDs");
    }

    typename simplex_t::NodeType compute_center_of_mass(std::list<typename dcel_t::node_t*> points){
        typename simplex_t::NodeType center_of_mass;
		center_of_mass(0) = 0;
        center_of_mass(1) = 0;
        for(auto pt : points){
            auto pt_coords = pt->coords();
            center_of_mass(0) += pt_coords(0);
            center_of_mass(1) += pt_coords(1);
        }
        center_of_mass(0) = center_of_mass(0) / points.size();
        center_of_mass(1) = center_of_mass(1) / points.size();
        return center_of_mass;
    }


    double euclidean_distance(double x1, double y1, double x2, double y2) {
        return std::sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
    }

    void find_intersection_by_points(
        Eigen::Matrix<double, 2, 1> p1,
        Eigen::Matrix<double, 2, 1> p2,
        Eigen::Matrix<double, 2, 1> q1,
        Eigen::Matrix<double, 2, 1> q2,
        bool& check_intersection,
        typename simplex_t::NodeType& intersection_point
    ){

        // Compute differences
        double dx1 = p2(0) - p1(0);
        double dy1 = p2(1) - p1(1);
        double dx2 = q2(0) - q1(0);
        double dy2 = q2(1) - q1(1);

        // Compute determinant
        double det = dx1 * dy2 - dy1 * dx2;
        if (det != 0.0) {
            // Compute parameters t and s
            double t = ((q1(0) - p1(0)) * dy2 - (q1(1) - p1(1)) * dx2) / det;
            double s = ((q1(0) - p1(0)) * dy1 - (q1(1) - p1(1)) * dx1) / det;

            // Check if intersection is within both segments
            if (t >= 0 && t <= 1 && s >= 0 && s <= 1) {
                intersection_point(0) = p1(0) + t * dx1;
                intersection_point(1) = p1(1) + t * dy1;
                check_intersection = true;
            }
        }

    }

    void find_intersection(const Triangulation<local_dim, embed_dim>& mesh,
                            typename dcel_t::halfedge_t cur_halfedge, 
                            typename triangle_t::EdgeType d_boundary_edge,
                            bool& check_intersection,
                            typename simplex_t::NodeType& intersection_point){
        
        Eigen::Matrix<double, 2, 1> p1 = cur_halfedge.node()->coords();
        Eigen::Matrix<double, 2, 1>  p2 = cur_halfedge.twin()->node()->coords();
        Eigen::Matrix<int, Dynamic, 1> node_ids = d_boundary_edge.node_ids();
        Eigen::Matrix<double, 2, 1> q1 = mesh.node(node_ids[0]);
        Eigen::Matrix<double, 2, 1> q2 = mesh.node(node_ids[1]);

        // Compute differences
        double dx1 = p2(0) - p1(0);
        double dy1 = p2(1) - p1(1);
        double dx2 = q2(0) - q1(0);
        double dy2 = q2(1) - q1(1);

        // Compute determinant
        double det = dx1 * dy2 - dy1 * dx2;
        if (det != 0.0) {
            // Compute parameters t and s
            double t = ((q1(0) - p1(0)) * dy2 - (q1(1) - p1(1)) * dx2) / det;
            double s = ((q1(0) - p1(0)) * dy1 - (q1(1) - p1(1)) * dx1) / det;

            // Check if intersection is within both segments
            if (t >= 0 && t <= 1 && s >= 0 && s <= 1) {
                intersection_point(0) = p1(0) + t * dx1;
                intersection_point(1) = p1(1) + t * dy1;
                check_intersection = true;
            }
        }
    }

    // std::map<std::pair<double, double>, std::vector<coords_t>> adjacent_points(const Triangulation<local_dim, embed_dim>& mesh){
    //     std::map<std::pair<double, double>, std::vector<coords_t>> adjacent_points_map;
    //     const Eigen::Matrix<double, Dynamic, Dynamic> mesh_nodes = mesh.nodes();
    //     for(int id = 0; id < mesh.n_nodes(); ++id){
    //         if(mesh.is_node_on_boundary(id)){
    //             double centroid_x = mesh.node(id)(0);
    //             double centroid_y = mesh.node(id)(1);
    //             std::pair<double, double> centroid_xy = std::make_pair(centroid_x, centroid_y);
    //             adjacent_points_map[centroid_xy];
    //             std::vector<int> adjacent_nodes_ids = mesh.node_one_ring(id);

    //             for(int node_id: adjacent_nodes_ids){
    //                 if(mesh.is_node_on_boundary(node_id)){   // TODO: here you should check that the entire edge is actually at the boundary, otherwise you get issues (see unit square)
    //                     coords_t node_coords = mesh.node(node_id);
    //                     adjacent_points_map[centroid_xy].push_back(node_coords);
    //                 }
    //             }
    //         }
    //     }
    //     return adjacent_points_map;
    // }

    std::map<std::pair<double, double>, std::vector<coords_t>> adjacent_points(const Triangulation<local_dim, embed_dim>& mesh){

        std::map<std::pair<double, double>, std::vector<coords_t>> adjacent_points_map;
        
        // Loop over boundary edges
        for (auto d_boundary_edge = mesh.boundary_edges_begin(); d_boundary_edge != mesh.boundary_edges_end(); ++d_boundary_edge){

            // Define first and second endpoints
            Eigen::Matrix<int, Dynamic, 1> endpoints_matrix = d_boundary_edge -> node_ids();
            int first_endpoint_id = endpoints_matrix(0);
            int second_endpoint_id = endpoints_matrix(1);
            coords_t first_endpoint = mesh.node(first_endpoint_id);
            coords_t second_endpoint = mesh.node(second_endpoint_id);
            std::pair<double, double> first_endpoint_pair = std::make_pair(first_endpoint(0), first_endpoint(1));
            std::pair<double, double> second_endpoint_pair = std::make_pair(second_endpoint(0), second_endpoint(1));

            // If first endpoint already in map, add second endpoint
            if(adjacent_points_map.find(first_endpoint_pair) != adjacent_points_map.end()){
                adjacent_points_map[first_endpoint_pair].push_back(second_endpoint);
            }
            // Else, create its entry and add second endpoint
            else{
                adjacent_points_map[first_endpoint_pair];
                adjacent_points_map[first_endpoint_pair].push_back(second_endpoint);
            }

            // Repeat for second endpoint
            if(adjacent_points_map.find(second_endpoint_pair) != adjacent_points_map.end()){
                adjacent_points_map[second_endpoint_pair].push_back(first_endpoint);
            }
            else{
                adjacent_points_map[second_endpoint_pair];
                adjacent_points_map[second_endpoint_pair].push_back(first_endpoint);
            }
        }

        return adjacent_points_map;
    }

    protected:
    dcel_t dcel_;

    // Internal function for computing v_vertex neighbours
    std::vector<int> compute_neighbours_(const Triangulation<local_dim, embed_dim>& mesh,
        const triangle_t& d_simplex, 
        std::map<int, typename simplex_t::NodeType>& midpoints_lookup,
        std::map<int, typename simplex_t::NodeType>& midpoints_lookup_raw,
        int& cur_midpoint_id,
        bool print,
        std::set<int>& out_of_boundary_d_centroids){
                                        
        // Get neighbours
        Eigen::Matrix<int, Eigen::Dynamic, 1> neigh = d_simplex.neighbors();

        // Store them inside a vector
        std::vector<int> neigh_vec(neigh.data(), neigh.data() + neigh.size());

        // Midpoint IDs
        std::vector<int> mp_ids;

        // Loop over d_simplex edges
        for (auto d_simplex_edge = d_simplex.edges_begin(); d_simplex_edge != d_simplex.edges_end(); d_simplex_edge++){

                // If edge on boundary ...
                if (d_simplex_edge -> on_boundary()){

                    // ... compute midpoint ...
                    auto mp = d_simplex_edge -> compute_midpoint();

                    // ... and create a new midpoint index
                    midpoints_lookup_raw[cur_midpoint_id] = mp;

                    // update midpoint id to edge map id
                    midpoint_to_edge_[cur_midpoint_id] = d_simplex_edge->id();

                    // Retrieve centroid
                    auto centroid_ = d_simplex.circumcenter();

                    // Get d_simplex vertexes ids
                    Eigen::Matrix<int, Eigen::Dynamic, 1> d_simplex_vertexes_ids = d_simplex.node_ids();
                    Eigen::Matrix<int, Dynamic, 1> d_simplex_edge_ids = d_simplex_edge->node_ids();

                    typename simplex_t::NodeType in_vertex_coords;
                    std::vector<typename simplex_t::NodeType> boundary_vertex_coords;

                    // Loop over d_vertexes of this d_simplex
                    for(auto d_simplex_vertex_id : d_simplex_vertexes_ids){

                        // Vertex coords
                        auto d_simplex_vertex_coords = mesh.node(d_simplex_vertex_id);

                        if(d_simplex_vertex_id == d_simplex_edge_ids(0) || d_simplex_vertex_id == d_simplex_edge_ids(1)){
                            boundary_vertex_coords.push_back(d_simplex_vertex_coords);
                        }
                        else{
                            in_vertex_coords = d_simplex_vertex_coords;
                        }

                    }

                    // Create vectors
                    typename simplex_t::NodeType v0, v1;
                    v0 = boundary_vertex_coords[0] - in_vertex_coords;
                    v1 = boundary_vertex_coords[1] - in_vertex_coords;

                    // Compute scalar product
                    float sp = v0(0)*v1(0) + v0(1)*v1(1);

                    // If scalar product negative, angle is > 90 -> change midpoint position
                    if(sp < 0){
                        mp(0) = 2*centroid_(0) - mp(0);
                        mp(1) = 2*centroid_(1) - mp(1);
                        out_of_boundary_d_centroids.insert(d_simplex.id());
                    }

                    midpoints_lookup[cur_midpoint_id] = mp;

                    // Add in the added midpoint ids
                    mp_ids.push_back(cur_midpoint_id);

                    cur_midpoint_id -= 1;
                }
            }

        // Replace inside neigh_vec -1s with actual midpoint IDs 
        std::size_t ii=0;
        std::size_t jj=0;
        for (int& neigh_tmp_id : neigh_vec){
            if (neigh_tmp_id == -1){
                neigh_vec[jj] = mp_ids[ii];
                ii += 1;
            }
            jj += 1;
        }

        return neigh_vec;

    }

    // Retrieve v_vertex corresponding to d_centroid, if present
    int v_vertex_from_d_centroid_(
        const std::vector<typename simplex_t::NodeType>& v_vertex_lookup, 
        const typename simplex_t::NodeType& d_centroid){

            for (auto ii = 0; ii < v_vertex_lookup.size(); ++ii){
                float thresh = 1e-9;
                if (std::abs(v_vertex_lookup[ii](0)-d_centroid(0))<thresh and std::abs(v_vertex_lookup[ii](1)-d_centroid(1))<thresh){
                    return ii;
                }
            }

            return -1;
    }
    
    // Compute v_vertexes
    std::vector<std::vector<int>> v_vertex_computation_(const Triangulation<local_dim, embed_dim>& mesh, 
                            dcel_t& dcel_, 
                            std::vector<int>& d_centroid2v_vertex, 
                            std::vector<typename simplex_t::NodeType>& v_vertex_lookup,
                            std::map<int, typename simplex_t::NodeType>& midpoints_lookup,
                            std::map<int, typename simplex_t::NodeType>& midpoints_lookup_raw,
                            int& cur_id,
                            std::set<int>& out_of_boundary_d_centroids){

        // Number of faces in the mesh
        int n_d_simplexes = mesh.n_cells();

        // Decremental number for midpoint IDs
        int cur_midpoint_id = -1;

        // Neighbours structure
        std::vector<std::vector<int>> v_vertex_neighbours_lookup(n_d_simplexes);

        // Helper function associating to each v_vertex_id the first d_centroid_id corresponding to it
        std::vector<int> first_d_centroid(n_d_simplexes);

        // Loop over all d_simplexes
        for (auto d_simplex_it = mesh.cells_begin(); d_simplex_it != mesh.cells_end(); ++d_simplex_it){

            // Retrieve d_centroid ID
            int d_centroid_id = d_simplex_it -> id();

            // Compute d_centroid
            typename simplex_t::NodeType d_centroid = d_simplex_it->circumcenter();

            // Compute neighbours
            std::vector<int> neigh_vec = compute_neighbours_(mesh, *(d_simplex_it), midpoints_lookup, midpoints_lookup_raw, cur_midpoint_id, false, out_of_boundary_d_centroids);

            // Check whether d_centroid is already associated to a v_vertex
            int v_vertex_id = v_vertex_from_d_centroid_(v_vertex_lookup, d_centroid);

            // If d_centroid already associated to a v_vertex, just update neighbours
            if (v_vertex_id != -1){

                // Insert current d_centroid neighbours to corresponding v_vertex neighbours
                v_vertex_neighbours_lookup[v_vertex_id].insert(v_vertex_neighbours_lookup[v_vertex_id].end(),
                                                                neigh_vec.begin(), neigh_vec.end());

                // Remove current d_centroid from corresponding v_vertex neighbours
                v_vertex_neighbours_lookup[v_vertex_id].erase(
                    std::remove(v_vertex_neighbours_lookup[v_vertex_id].begin(), v_vertex_neighbours_lookup[v_vertex_id].end(), d_centroid_id),
                    v_vertex_neighbours_lookup[v_vertex_id].end());

                // Remove from corresponding v_vertex neighbours the first found d_centroid corresponding to it
                v_vertex_neighbours_lookup[v_vertex_id].erase(
                    std::remove(v_vertex_neighbours_lookup[v_vertex_id].begin(), v_vertex_neighbours_lookup[v_vertex_id].end(), first_d_centroid[v_vertex_id]),
                    v_vertex_neighbours_lookup[v_vertex_id].end());
            }

            // Else if d_centroid is not already associated to a v_vertex, create its v_vertex and neighbours
            else{

                // Create v_vertex ID
                v_vertex_id = cur_id;

                // Create a node for v_vertex
                typename dcel_t::node_t v_vertex(v_vertex_id, false, d_centroid);

                // Insert in DCEL
                dcel_.insert_node(v_vertex);

                // Update all structures
                v_vertex_neighbours_lookup[cur_id] = neigh_vec;
                first_d_centroid[cur_id] = d_centroid_id;
                v_vertex_lookup[cur_id] = d_centroid;

                // Increment current ID
                cur_id += 1;
            }

            // Update d_centroid2v_vertex
            d_centroid2v_vertex[d_centroid_id] = v_vertex_id;

        }

        return v_vertex_neighbours_lookup;

    }

    int retrieve_v_vertex_id_(
                            int d_centroid_id,
                            std::vector<int>& d_centroid2v_vertex){

        // If positive, than it is a physical d_centroid
        if (d_centroid_id >= 0) return d_centroid2v_vertex[d_centroid_id];

        // Else it is a midpoint
        return d_centroid_id;
    }

    typename simplex_t::NodeType id2point_(
                                        int id,
                                        std::vector<typename simplex_t::NodeType>& v_vertex_lookup,
                                        std::map<int, typename simplex_t::NodeType>& midpoints_lookup
    ){
        // If positive, then it is a physical d_centroid
        if (id >= 0) return v_vertex_lookup.at(id);

        // Else it is a midpoint
        return midpoints_lookup.at(id);
    }

    std::vector<int> order_neighbours_cclw_(
                                        int v_vertex_id,
                                        std::vector<std::vector<int>>& v_vertex_neighbours_lookup,
                                        std::vector<int>& d_centroid2v_vertex,
                                        std::vector<typename simplex_t::NodeType>& v_vertex_lookup,
                                        std::map<int, typename simplex_t::NodeType>& midpoints_lookup
    ){

        // Retrieve v_vertex_coords
        typename simplex_t::NodeType v_vertex = id2point_(v_vertex_id, v_vertex_lookup, midpoints_lookup);

        // Create ordered neighbours structure
        std::vector<int> ordered_neigh_ids;
        ordered_neigh_ids.reserve(v_vertex_neighbours_lookup.at(v_vertex_id).size());

        // Loop over neighbours
        for (auto d_centroid_neigh_id : v_vertex_neighbours_lookup.at(v_vertex_id)) {  // Access O(1)

            int neigh_id = retrieve_v_vertex_id_(d_centroid_neigh_id, d_centroid2v_vertex);  // O(1)
            
            ordered_neigh_ids.push_back(neigh_id);
        }

        std::sort(ordered_neigh_ids.begin(), ordered_neigh_ids.end(),
            [&](int id_a, int id_b)
        {
            typename simplex_t::NodeType pa = id2point_(id_a, v_vertex_lookup, midpoints_lookup);
            typename simplex_t::NodeType pb = id2point_(id_b, v_vertex_lookup, midpoints_lookup);

            double ang_a = std::atan2(pa.y() - v_vertex.y(),
                                    pa.x() - v_vertex.x());
            double ang_b = std::atan2(pb.y() - v_vertex.y(),
                                    pb.x() - v_vertex.x());

            return ang_a > ang_b;
        });

        return ordered_neigh_ids;
    }

    void create_and_connect_halfedges_(dcel_t& dcel_, 
                            int v_vertex_id, int infty_id,
                            typename dcel_t::node_t* infty_node,
                            std::vector<typename simplex_t::NodeType>& v_vertex_lookup,
                            std::map<int, typename simplex_t::NodeType>& midpoints_lookup_raw,
                            std::map<int, typename simplex_t::NodeType>& midpoints_lookup,
                            const std::vector<int>& ordered_neigh_ids,
                            std::multimap<std::pair<int, int>, typename dcel_t::halfedge_t*>& vertexes2halfedge,
                            std::map<int, int>& infty_halfedges_midpoints) {  // --> O(logN)
                 
        int count_existing_neigh = 0;
        // Keep the first halfedge
        typename dcel_t::halfedge_t* first_halfedge=nullptr;
        // Pointers to e1 and e2 initialized as null_pointer
        typename dcel_t::halfedge_t* e1 = nullptr;
        typename dcel_t::halfedge_t* e2 = nullptr;
        // Loop over neighbouring cells
        for (auto neigh_id : ordered_neigh_ids) {  // On average should be O(1)
            // flag for infinity node
            bool is_infinity = false;
            // id for computing the midpoint id
            int mid_id = 0;
            // Create cell node
            typename dcel_t::node_t* v_vertex;
            // Find nodes in the DCEL structure
            v_vertex = dcel_.find_node(dcel_.nodes().row(v_vertex_id));  // O(1)
            // Create twin node 
            typename dcel_t::node_t* twin_node;
            // Create halfedge centroid -> neigh
            typename dcel_t::halfedge_t* cell_halfedge = new dcel_t::halfedge_t();
            // Create twin halfedge (neigh -> centroid)
            typename dcel_t::halfedge_t* twin_halfedge = new dcel_t::halfedge_t();

            if(neigh_id >= 0){
                // Retrieve new ID for neighbour
                if (neigh_id == v_vertex_id){
                    continue;
                }
                // Retrieve centroid of the neighbouring cell
                typename simplex_t::NodeType neigh_v_vertex = v_vertex_lookup.at(neigh_id);
                twin_node = dcel_.find_node(dcel_.nodes().row(neigh_id));  // O(1)
            }
            else{
                // The centroid of the neighbouring cell is the infinity node
                mid_id = neigh_id;
                neigh_id = infty_id;
                twin_node = infty_node;
                is_infinity = true;
            }

            // If neighbouring cell was not visited, its node needs an halfedge
            if ((twin_node->neighID2halfedge(v_vertex_id)==nullptr) or (neigh_id == infty_id)){
                // Set IDs based on counter
                cell_halfedge = dcel_.emplace_halfedge(v_vertex);  // O(1)
                twin_halfedge = dcel_.emplace_halfedge(twin_node);  // O(1)

                // Set halfedge to the cell centroid
                v_vertex->set_halfedge(cell_halfedge);
                // Set halfedge to the neighbouring centroid
                twin_node->set_halfedge(twin_halfedge);

                cell_halfedge->set_twin(twin_halfedge);
                twin_halfedge->set_twin(cell_halfedge);

                // Update lookup
                v_vertex->add_halfedge(neigh_id, cell_halfedge);
                twin_node->add_halfedge(v_vertex_id, twin_halfedge);
                vertexes2halfedge.insert({std::make_pair(v_vertex->id(), twin_node->id()), cell_halfedge});
                vertexes2halfedge.insert({std::make_pair(twin_node->id(), v_vertex->id()), twin_halfedge});

                // Update intersection d_edges map
                if((is_infinity) and midpoints_lookup_raw.at(mid_id)==midpoints_lookup.at(mid_id)){ // O(logN), anche se secondo me sarà tipo in media O(log(sqrt(N))), essendo i punti medi solo al bordo
                    intersection_d_edges_[cell_halfedge->id()]= midpoint_to_edge_.at(mid_id);
                    intersection_d_edges_[twin_halfedge->id()] = midpoint_to_edge_.at(mid_id);
                    infty_halfedges_midpoints[cell_halfedge->id()] = mid_id;
                    infty_halfedges_midpoints[twin_halfedge->id()] = mid_id;
                }
            }
            // If neighbouring cell was visited, only retrieve the halfedges
            else{
                delete cell_halfedge;
                delete twin_halfedge;
                twin_halfedge = (twin_node->neighID2halfedge(v_vertex_id));
                cell_halfedge = (twin_halfedge->twin());
            }

            // Set previous and next
            e2 = cell_halfedge;  

            if (count_existing_neigh != 0){
                e1->twin()->set_next(e2);
                e2->set_prev(e1->twin());
            }
            else{
                first_halfedge = cell_halfedge;
            }

            e1 = cell_halfedge;
            count_existing_neigh += 1;
    }

    // Complete the association of previous and next
    if (count_existing_neigh != 1)
    {
    e2->twin()->set_next(first_halfedge);
    first_halfedge->set_prev(e2->twin());
    }
        
    }

    void create_cells_(dcel_t& dcel_, 
                        const Triangulation<local_dim, embed_dim>& mesh, 
                        int v_vertex_id,
                        int d_centroid_id,
                        std::vector<int>& d_centroid2v_vertex,
                        std::vector<typename simplex_t::NodeType>& v_vertex_lookup,
                        const Eigen::Matrix<int, Eigen::Dynamic, 1>& neighbor_simplexes,
                        std::multimap<std::pair<int, int>, typename dcel_t::halfedge_t*>& vertexes2halfedge,
                        std::set<int>& not_created_cells) {  // --> O(N)
        
        int count_neigh_tria = 0;
        int diff_id = 0;
        // Retrieve the IDs of the vertices of the current cell
        triangle_t curr_tria(d_centroid_id, &mesh); 
        Eigen::Matrix<int, Dynamic, 1> d_vertex_ids = curr_tria.node_ids();
        // Retrieve the centroid of the current cell
        typename simplex_t::NodeType v_vertex = v_vertex_lookup.at(v_vertex_id);  // O(logN)
        // Helper structure for cells creation
        
        for (int neigh_d_centroid_id : neighbor_simplexes){  // O(1) on average
            if(neigh_d_centroid_id != -1){
                count_neigh_tria += 1;
                // Retrieve the triangle in the mesh using the old id
                triangle_t tria(neigh_d_centroid_id, &mesh);  // O(?)
                auto neigh_d_vertex_ids = tria.node_ids();

                // Find the common and different IDs between these vertices and the cell ones
                std::vector<int> vec_1(d_vertex_ids.data(), d_vertex_ids.data()+d_vertex_ids.size());
                std::vector<int> vec_2(neigh_d_vertex_ids.data(), neigh_d_vertex_ids.data()+neigh_d_vertex_ids.size());
                // Sort the vectors
                std::sort(vec_1.begin(), vec_1.end());
                std::sort(vec_2.begin(), vec_2.end());
                std::vector<int> common;

                // Take the common elements between the two
                std::set_intersection(vec_1.begin(), vec_1.end(),
                                    vec_2.begin(), vec_2.end(),
                                    std::back_inserter(common));

                // Assert that the common elements are two since they are adjacent cells
                assert(common.size()==2);

                // Create the cells and associate the correct halfedges to them
                for(int i=0; i<common.size(); ++i){  // O(1)

                    int idx = common[i];

                    // Retrieve the corrensponding coordinates
                    typename simplex_t::NodeType d_vertex = mesh.node(idx); // O(1)

                    int neigh_v_vertex_id = d_centroid2v_vertex.at(neigh_d_centroid_id);
                    typename simplex_t::NodeType neigh_v_vertex = v_vertex_lookup.at(neigh_v_vertex_id);  // O(logN)
                    // Create the new cell
                    cell_t curr_cell(idx);

                    // create segments around the vertex
                    typename simplex_t::NodeType u = v_vertex - d_vertex;
                    typename simplex_t::NodeType v = neigh_v_vertex - d_vertex;
                    // Compute the cross product to assign the halfedges properly 
                    double cross = u(0)*v(1)-u(1)*v(0);
                    typename dcel_t::halfedge_t* current_cell_halfedge = vertexes2halfedge.find(std::make_pair(v_vertex_id, neigh_v_vertex_id))->second;  // O(logN)
                    typename dcel_t::halfedge_t* current_twin_halfedge = vertexes2halfedge.find(std::make_pair(neigh_v_vertex_id, v_vertex_id))->second;

                    if (mesh.is_node_on_boundary(idx)){
                        curr_cell.set_unbounded();
                    }

                    // Add the new cell to the list if it is not yet added 
                    // if(std::find(cells_.begin(), cells_.end(), curr_cell) == cells_.end())
                    if (cells_.find(curr_cell.id()) == cells_.end()){  // O(logN)
                        if(cross < -1e-9){
                            // cross product < tolerance means that neigh_centroid is clock-wise with respect to centroid, so we want neigh -> centroid
                            curr_cell.set_halfedge(current_twin_halfedge);
                            // cells_.push_back(curr_cell);
                            cells_[curr_cell.id()] = curr_cell;
                            dcel_.insert_cell(curr_cell);
                            // not_created_cells.remove(idx);
                            not_created_cells.erase(idx);
                        }
                        else if(cross > 1e-9){
                            // cross product > tolerance means that neigh_centroid is CCW with respect to centroid, so we want centroid -> neigh
                            curr_cell.set_halfedge(current_cell_halfedge);
                            // cells_.push_back(curr_cell);
                            cells_[curr_cell.id()] = curr_cell;
                            dcel_.insert_cell(curr_cell);
                            // not_created_cells.remove(idx);
                            not_created_cells.erase(idx);
                        }
                        else{
                            // if(std::find(not_created_cells.begin(), not_created_cells.end(), idx) == not_created_cells.end())
                            if (not_created_cells.find(idx) == not_created_cells.end()){ // O(N)
                            not_created_cells.insert(idx);}   
                        }
                        }

                }

                // Take the different id between the two
                std::vector<int> tmp_diff;
                // Element that is in vec1 but not in vec2
                std::set_difference(
                    vec_1.begin(), vec_1.end(),
                    vec_2.begin(), vec_2.end(),
                    std::back_inserter(tmp_diff)
                );

                // Assert that the different element is one since they are adjacent cells
                assert(tmp_diff.size()==1);
                diff_id = tmp_diff[0];
            }
    }

    if(count_neigh_tria == 1){
        cell_t curr_cell_diff(diff_id);
        // if (std::find(cells_.begin(), cells_.end(), curr_cell_diff) == cells_.end())
        if (cells_.find(curr_cell_diff.id()) == cells_.end()){  // O(N)
            not_created_cells.insert(diff_id);
        }
    }
    }

    std::map<int, std::vector<int>> compute_v_cell2halfedges_(int infty_id){
        std::map<int, std::vector<int>> v_cell2halfedges;
        for(auto cell_it = dcel_.cells_cbegin(); cell_it != dcel_.cells_cend(); ++cell_it){
            std::vector<int> halfedges_vector;
            if(!cell_it->is_unbounded()){
                for(auto v: cell_it->cell_edges()){
                    halfedges_vector.push_back(v->id());
                }
            }
            else{
                for(auto v: cell_it->cell_edges_with_infty(infty_id, false)){
                    halfedges_vector.push_back(v->id());
                }
            }

            v_cell2halfedges[cell_it->id()] = halfedges_vector;
        }

        return v_cell2halfedges;
    }


    void add_not_created_cells_(int infty_id, std::set<int>& not_created_cells,
                                const Triangulation<local_dim, embed_dim>& mesh,
                                std::vector<int>& d_centroid2v_vertex,
                                std::multimap<std::pair<int, int>, typename dcel_t::halfedge_t*>& vertexes2halfedge){

        std::map<int, std::vector<int>> v_cell2halfedges = compute_v_cell2halfedges_(infty_id);

        if (not_created_cells.size()>0){

            for(int internal_id : not_created_cells){
                cell_t curr_cell(internal_id);
                // if (std::find(cells_.begin(), cells_.end(), curr_cell) != cells_.end())
                if (cells_.find(curr_cell.id()) != cells_.end()){
                    continue;
                }
                if (mesh.is_node_on_boundary(internal_id)){
                    curr_cell.set_unbounded();
                }
                
                // Retreive one neighbour
                int neighbour_v_cell = 0;
                if(internal_id != mesh.node_one_ring(internal_id)[0]){
                    neighbour_v_cell = mesh.node_one_ring(internal_id)[0];
                }
                else{
                    neighbour_v_cell = mesh.node_one_ring(internal_id)[1];
                }
                // Vector with cell halfedges
                std::vector<int> v_cell_halfedges = v_cell2halfedges[neighbour_v_cell];
                
                std::vector<int> v_vertexes = mesh.node_patch(internal_id);
                int v_vertex_id = d_centroid2v_vertex[v_vertexes[0]];
                
                bool breaking_check = false;

                for(int halfedge_id : v_cell_halfedges){
                    auto to_infty_list = vertexes2halfedge.equal_range(std::make_pair(v_vertex_id, infty_id));
                    auto from_infty_list = vertexes2halfedge.equal_range(std::make_pair(infty_id, v_vertex_id));
                    for(auto to_infty = to_infty_list.first; to_infty != to_infty_list.second; ++to_infty){
                        if(to_infty->second and (halfedge_id == to_infty->second->id())){
                            curr_cell.set_halfedge(to_infty->second->twin());
                            breaking_check = true;
                            break;
                        }
                    }

                    if(breaking_check) break;

                    for(auto from_infty = from_infty_list.first; from_infty != from_infty_list.second; ++from_infty){
                        if(from_infty->second and (halfedge_id == from_infty->second->id())){
                            curr_cell.set_halfedge(from_infty->second->twin());
                            breaking_check = true;
                            break;
                        }
                    }

                    if(breaking_check) break;
                }

                // cells_.push_back(curr_cell);
                cells_[curr_cell.id()] = curr_cell;
                dcel_.insert_cell(curr_cell);
            }
        }
    }


    void connect_infty_halfedges_(int infty_id){
        std::map<int, std::vector<int>> v_cell2halfedges = compute_v_cell2halfedges_(infty_id);

        std::map<int, typename dcel_t::halfedge_t*> halfedges_map;
        for(auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it){
            halfedges_map[it->id()] = &(*it);
        }

        for (const auto& [v_cell_id, halfedges_ids] : v_cell2halfedges) {
            typename dcel_t::halfedge_t* from_infty;
            typename dcel_t::halfedge_t* to_infty;
            for (int halfedge_id : halfedges_ids) {
                typename dcel_t::halfedge_t* cur_halfedge = halfedges_map[halfedge_id];
                if (cur_halfedge->node()->id()==infty_id){
                    from_infty = cur_halfedge;
                }
                else if (cur_halfedge->twin()->node()->id()==infty_id){
                    to_infty = cur_halfedge;
                }
            }

            to_infty->set_next(from_infty);
            from_infty->set_prev(to_infty);
        }
    }   

    // List of cells
    // std::list<typename dcel_t::cell_t> cells_;
    std::map<int, typename dcel_t::cell_t> cells_;
    // Midpoint id to egde id
    std::map<int, int> midpoint_to_edge_;
    // Intersection d_edges
    std::map<int, int> intersection_d_edges_;
};

}


#endif