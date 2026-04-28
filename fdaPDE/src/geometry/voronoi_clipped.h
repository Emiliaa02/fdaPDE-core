#ifndef __FDAPDE_VORONOI_CLIPPED_H__
#define __FDAPDE_VORONOI_CLIPPED_H__


namespace fdapde{

template <int LocalDim, int EmbedDim>
class Voronoi_Clipped : public Voronoi<LocalDim, EmbedDim>{

    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;

    using dcel_t = DCEL<local_dim, embed_dim>;
    using simplex_t = Simplex<local_dim, embed_dim>;
    using triangulation_t = Triangulation<local_dim, embed_dim>;
    using triangle_t = Triangle<triangulation_t>;
    using cell_t = dcel_t::cell_t;
    
    public:

    Voronoi_Clipped(const Triangulation<local_dim, embed_dim>& mesh) : Voronoi<LocalDim, EmbedDim>(mesh){
        clip_voronoi(mesh, infty_id, out_of_boundary_d_centroids, old2new, vertexes2halfedge);
        std::cout << "\nFinished constructor" << std::endl;
    }


    private:

    void clip_voronoi(const Triangulation<local_dim, embed_dim>& mesh, int infty_id,
                      std::set<int>& out_of_boundary_d_centroids, std::vector<int>& d_centroid2v_vertex,
                      std::multimap<std::pair<int, int>, typename dcel_t::halfedge_t*>& vertexes2halfedge){

        std::map<int, std::vector<typename dcel_t::node_t*>> supplementary_points;
        std::set<int> on_boundary;
        std::vector<int> node_ids_to_remove;

        boundary_cells_detection(mesh, supplementary_points, on_boundary,
                                out_of_boundary_d_centroids, d_centroid2v_vertex,
                                infty_id, node_ids_to_remove);

        // Loop over cells of the DCEL structure
        for (auto cur_cell_it = this->dcel_.cells_begin(); cur_cell_it != this->dcel_.cells_end(); cur_cell_it++){
            if(on_boundary.find((*cur_cell_it).id()) == on_boundary.end()){ continue; };
            typename dcel_t::halfedge_t* start = nullptr;
            typename dcel_t::halfedge_t* cur_halfedge = nullptr;
            typename dcel_t::halfedge_t* prev_halfedge = nullptr;

            // Create a list of its points
            std::vector<typename dcel_t::node_t*> pts;
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
            std::sort(pts.begin(), pts.end(),
                [&](dcel_t::node_t* pa, dcel_t::node_t* pb)
            {

                double ang_a = std::atan2(pa->coords()(1) - center_of_mass(1),
                                        pa->coords()(0) - center_of_mass(0));
                double ang_b = std::atan2(pb->coords()(1) - center_of_mass(1),
                                        pb->coords()(0) - center_of_mass(0));

                return ang_a < ang_b;
            });
            
            if(pts.size() < 2){
                std::cout << "LESS THAN 2 POINTS" << std::endl;
                continue;
            }
            // Loop points one after the other
            for(auto pt = pts.begin(); pt != pts.end(); ++pt){
                // Next point
                auto pt_next = std::next(pt);
                if (pt_next == pts.end()) {
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
    void boundary_cells_detection(const Triangulation<local_dim, embed_dim>& mesh, 
                                std::map<int, std::vector<typename dcel_t::node_t*>>& supplementary_points,
                                std::set<int>& on_boundary,
                                std::set<int>& out_of_boundary_d_centroids,
                                std::vector<int>& d_centroid2v_vertex,
                                int infty_id, std::vector<int>& node_ids_to_remove){

        std::set<int> entered_first_time;
        int counter = infty_id + 1;
        node_ids_to_remove.push_back(infty_id);
        for(int d_centroid_id : out_of_boundary_d_centroids){
            int v_vertex_id = d_centroid2v_vertex.at(d_centroid_id);
            node_ids_to_remove.push_back(v_vertex_id);
        }
        // Loop over mesh boundary edges of the Delaunay
        for (auto d_boundary_edge = mesh.boundary_edges_begin(); d_boundary_edge != mesh.boundary_edges_end(); ++d_boundary_edge){
            // Find closest v_centroid to d_edge_midpoint
            auto current_midpoint = d_boundary_edge->compute_midpoint();
            auto closest_v_centroid = find_closest_v_centroid(mesh, current_midpoint);

            // Initialize list of v_cells to check
            std::list<typename dcel_t::cell_t> to_check_v_cells{closest_v_centroid};
            // Initialize map of checked halfedges
            std::map<int, bool> already_checked_halfedges;
            // Compute the v_cell2halfedges structure
            std::map<int, std::vector<int>> v_cell_halfedges = compute_v_cell2halfedges_(infty_id);

            // While list to check is not empty
            while (to_check_v_cells.size() > 0){
                typename dcel_t::cell_t cur_v_centroid = to_check_v_cells.front(); 
                to_check_v_cells.pop_front();
                std::vector<int> cur_v_cell_halfedges = v_cell_halfedges.at(cur_v_centroid.id());

                // Loop over cur_v_cell halfedges
                for(int cur_halfedge_id: cur_v_cell_halfedges){                  
                    typename dcel_t::halfedge_t cur_halfedge;
                    bool found_halfedge = false;
                
                    for(auto it = this->dcel_.halfedges_begin(); it != this->dcel_.halfedges_end(); ++it){
                        if(it->id() == cur_halfedge_id){
                            cur_halfedge = *it;
                            found_halfedge = true;
                            break;
                        }
                    }
                    // Check if the current halfedge exists
                    if (!found_halfedge) {
                        throw std::runtime_error("Halfedge not found");
                    }

                    // Check if the twin halfedge exists
                    auto twin = cur_halfedge.twin();
                    if (twin == nullptr) {
                        throw std::runtime_error("Twin is null");
                    }
                    if (twin->cell() == nullptr) {
                        throw std::runtime_error("Twin cell is null");
                    }

                    // If we already know halfedge intersects boundary edge (necessarily on the midpoint)...
                    auto check_it = already_checked_halfedges.find(cur_halfedge_id);
                    if((check_it == already_checked_halfedges.end() || !check_it->second) && intersection_d_edges_.find(cur_halfedge_id) != intersection_d_edges_.end()
                       && intersection_d_edges_.at(cur_halfedge_id) == d_boundary_edge->id()){
                        
                        // Cur_v_cell is a boundary cell
                        on_boundary.insert(cur_v_centroid.id());

                        // Store the site of the Voronoi cell only one time
                        if(entered_first_time.find(cur_v_centroid.id()) == entered_first_time.end()){
                            typename simplex_t::NodeType centroid_coords;
                            centroid_coords(0) = mesh.node(cur_v_centroid.id())(0);
                            centroid_coords(1) = mesh.node(cur_v_centroid.id())(1);
                            // Check if the current site is already a DCEL node
                            auto existing_node = dcel_.find_node_by_coords(centroid_coords);
                            if (existing_node != nullptr) {
                                supplementary_points[cur_v_centroid.id()].push_back(existing_node);
                            } else {
                                typename dcel_t::node_t* v_centroid_coords = new typename dcel_t::node_t();
                                v_centroid_coords->set_coords(centroid_coords);
                                v_centroid_coords->set_id(counter++);
                                supplementary_points[cur_v_centroid.id()].push_back(v_centroid_coords);
                            }
                        }  
                        entered_first_time.insert(cur_v_centroid.id());

                        // Add to cur_v_cell the midpoint
                        if(supplementary_points.find(cur_v_centroid.id()) == supplementary_points.end()){
                            supplementary_points[cur_v_centroid.id()];
                        }
                        typename dcel_t::node_t* v_current_midpoint;
                        auto node_ptr = dcel_.find_node_by_coords(current_midpoint);
                        if(node_ptr == nullptr){
                            v_current_midpoint = new typename dcel_t::node_t();
                            v_current_midpoint->set_coords(current_midpoint);
                            v_current_midpoint->set_id(counter++);
                        }
                        else{
                            v_current_midpoint = node_ptr;
                        }
                        supplementary_points[cur_v_centroid.id()].push_back(v_current_midpoint);

                        // Add to list neighbouring v_cell
                        to_check_v_cells.push_back(*(cur_halfedge.twin()->cell()));

                        already_checked_halfedges[cur_halfedge_id] = true;
                        
                    }
                    // ... otherwise
                    else if (check_it == already_checked_halfedges.end() || !check_it->second){
                        // Check whether halfedge and boundary edge intersect in a point
                        bool check_intersection= false;
                        typename simplex_t::NodeType intersection_point;
                        find_intersection(mesh, cur_halfedge, *d_boundary_edge, check_intersection, intersection_point);

                        // If you found it...
                        if(check_intersection){ 
                            typename dcel_t::node_t* v_intersection_point;
                            auto intersection_ptr =  dcel_.find_node_by_coords(intersection_point);
                            if(intersection_ptr == nullptr){
                                v_intersection_point = new typename dcel_t::node_t();
                                v_intersection_point->set_coords(intersection_point);
                                v_intersection_point->set_id(counter++);
                            }
                            else{
                                v_intersection_point = intersection_ptr;
                            }
                            // cur_v_cell is boundary cell
                            on_boundary.insert(cur_v_centroid.id());

                            // Add to cur_v_cell the point
                            if(supplementary_points.find(cur_v_centroid.id()) == supplementary_points.end()){
                                supplementary_points[cur_v_centroid.id()];
                            }
                            supplementary_points[cur_v_centroid.id()].push_back(v_intersection_point);

                            // Add to list neighbouring v_cell
                            to_check_v_cells.push_back(*(cur_halfedge.twin()->cell()));
                            already_checked_halfedges[cur_halfedge_id] = true;
                        }
                    }
                
                }
                
            }
        }

    }


    typename dcel_t::cell_t find_closest_v_centroid(const Triangulation<local_dim, embed_dim>& mesh, 
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

        typename dcel_t::cell_t best_cell;

        for(auto it = this->dcel_.cells_begin(); it != this->dcel_.cells_end(); ++it){
            if(it->id() == best_v_centroid_id){
                best_cell = (*it);
                return best_cell;
            }
        }
        throw std::runtime_error("Best cell ID is not included among Voronoi cell IDs");
    }


	typename simplex_t::NodeType compute_center_of_mass(std::vector<typename dcel_t::node_t*> points){
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


    }
}
#endif