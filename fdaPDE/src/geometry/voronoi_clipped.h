#ifndef __FDAPDE_VORONOI_CLIPPED_H__
#define __FDAPDE_VORONOI_CLIPPED_H__


namespace fdapde{

template <int LocalDim, int EmbedDim>
class VoronoiClipped : public Voronoi<LocalDim, EmbedDim>{

    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;

    using dcel_t = DCEL<local_dim, embed_dim>;
    using simplex_t = Simplex<local_dim, embed_dim>;
    using triangulation_t = Triangulation<local_dim, embed_dim>;
    using triangle_t = Triangle<triangulation_t>;
    using cell_t = dcel_t::cell_t;
    using voronoi_t = Voronoi<local_dim, embed_dim>;
    using coords_t = Eigen::Matrix<double, 1, embed_dim>;
    using node_t = dcel_t::node_t;
    
    public:

    // Constructor for the Clipped Voronoi structure
    VoronoiClipped(const Triangulation<local_dim, embed_dim>& mesh) : Voronoi<LocalDim, EmbedDim>(mesh){  // --> O(N^2)
        // Voronoi constructor O(N^2)
        clip_voronoi(mesh, this->infty_id, this->out_of_boundary_d_centroids, this->d_centroid2v_vertex, this->vertexes2halfedge,
                    this->midpoints_lookup, this->infty_halfedges_midpoints);  // O(N^2)
        std::cout << "\nFinished Voronoi clipped constructor" << std::endl;
    }

    // Check if Voronoi definition holds for single point
    const bool check_definition_on_single_point(const coords_t& point){

        assert(point.size() == embed_dim);

        // Find its closest v_centroid
        int closest_v_centroid = find_closest_v_centroid(point);
        if (closest_v_centroid < 0 || closest_v_centroid >= this->v_cell2v_centroid.size()) {
            std::cerr << "Invalid closest_v_centroid: " << closest_v_centroid
                    << ", v_cell2v_centroid size: " << this->v_cell2v_centroid.size()
                    << ", point: " << point << std::endl;
            return false;
        }

        // Check that point is within associated v_cell
        bool def_holds = is_point_in_cell(point, closest_v_centroid);

        // If definition does not hold, display node
        if (!def_holds){
            std::cerr << "Point does not satisfy definition of Voronoi: " << point  << ". Cell ID: " << 
            closest_v_centroid << " corresponding centroid: " << this->v_cell2v_centroid.at(closest_v_centroid) << std::endl;
        }

        return def_holds;
    }

    private:

    std::vector<std::set<int>> partition_boundary_edges(const Triangulation<local_dim, embed_dim>& mesh){ // -> O(N*logN)
       
        // Already visited cells
        std::set<int> visited;
        // Connected boundary-edge components (outer boundary and holes)
        std::vector<std::set<int>> boundaries;
        // Adjacent points map
        std::map<int, std::vector<int>> adj = adjacent_points(mesh); // RIVEDIIIIIIIIIIIII

        for (auto& kv : adj) // O(N)
        {
            int start = kv.first;

            if (visited.count(start))
                continue;

            // Component of the boundary
            std::set<int> component;
            std::queue<int> q;
            q.push(start); // O(1)

            visited.insert(start); // O(logN)
            while (!q.empty()) 
            {
                int v = q.front();
                q.pop();
                component.insert(v); 

                for (int nb : adj[v])
                {
                    if (!visited.count(nb))
                    {
                        visited.insert(nb);
                        q.push(nb);
                    }
                }
            }
            boundaries.push_back(component);
        }
        return boundaries;
    }


    void clip_voronoi(const Triangulation<local_dim, embed_dim>& mesh, int infty_id,
                      std::set<int>& out_of_boundary_d_centroids, std::vector<int>& d_centroid2v_vertex,
                      std::multimap<std::pair<int, int>, typename dcel_t::halfedge_t*>& vertexes2halfedge,
                      std::map<int, typename simplex_t::NodeType>&  midpoints_lookup,
                      std::map<int, int>& infty_halfedges_midpoints){

        std::map<int, std::vector<typename dcel_t::node_t*>> supplementary_points;
        std::map<int, std::vector<int>> adjacent_points_map = adjacent_points(mesh); // O(N*logN) (A dire il vero O(sqrt(N) * log(sqrt(N))))
        std::map<int, typename dcel_t::node_t*> cell_centroid_coords;
        std::set<int> on_boundary;
        std::set<int> node_ids_to_remove;
        std::map<int, std::pair<int, int>> edge_vertexes2intersection;
        
        boundary_cells_detection(mesh, supplementary_points, on_boundary,
                                out_of_boundary_d_centroids, d_centroid2v_vertex,
                                infty_id, node_ids_to_remove, cell_centroid_coords, infty_halfedges_midpoints,
                                midpoints_lookup, edge_vertexes2intersection); // O(N^2)


        // Loop over cells of the DCEL structure
        for (auto cur_cell_it = this->dcel_.cells_begin(); cur_cell_it != this->dcel_.cells_end(); cur_cell_it++){ // O(N)
            if(on_boundary.find((*cur_cell_it).id()) == on_boundary.end()){ continue; }; // O(logN)
            typename dcel_t::halfedge_t* start = nullptr;
            typename dcel_t::halfedge_t* cur_halfedge = nullptr;
            typename dcel_t::halfedge_t* prev_halfedge = nullptr;

            // Create a list of its points
            std::list<typename dcel_t::node_t*> pts;
            // Get nodes
            std::vector<typename dcel_t::node_t*> nodes = cur_cell_it -> cell_nodes(); // O(1)
            for (auto cell_node = nodes.cbegin(); cell_node != nodes.cend(); cell_node++){ // O(1)
                auto to_exclude_it = node_ids_to_remove.find((*cell_node)->id());  // O(log N)
                if(to_exclude_it == node_ids_to_remove.end()){
                    pts.push_back(*cell_node); // O(1)
                }
            }


            // Get supplementary points
            auto check_supp_pt = supplementary_points.find(cur_cell_it->id()); // O(logN)
            if (check_supp_pt != supplementary_points.end()){
                // Insert supplementary points in DCEL structure, if they are not yet included 
                for(typename dcel_t::node_t* supp_pt : check_supp_pt->second){ // O(1)                    
                    auto node_ptr = this->nodes_map.find(supp_pt->id()); // O(logN)
                    if(node_ptr == this->nodes_map.end()){
                        typename dcel_t::node_t* inserted_node = this->dcel_.insert_node(*supp_pt); // O(1)
                        pts.push_back(inserted_node); // O(1)
                    }
                    else{pts.push_back(node_ptr->second);}; // O(1)
                }
            }

            // Compute the center of mass to perform the ordering
            typename simplex_t::NodeType center_of_mass = compute_center_of_mass(pts); // O(1)

            // If pts has only size 2, insert centroid directly
            bool only_two_pts = false;
            if (pts.size()==2){
                pts.push_back(cell_centroid_coords.at((*cur_cell_it).id())); // O(logN)
                typename dcel_t::node_t* inserted_node = this->dcel_.insert_node(*cell_centroid_coords.at((*cur_cell_it).id())); // O(logN)
                only_two_pts = true;
            }

            // Put in counterclockwise order the points
            pts.sort([&](dcel_t::node_t* pa, dcel_t::node_t* pb)
                {
                    double ang_a = std::atan2(pa->coords()(1) - center_of_mass(1),
                                            pa->coords()(0) - center_of_mass(0));

                    double ang_b = std::atan2(pb->coords()(1) - center_of_mass(1),
                                            pb->coords()(0) - center_of_mass(0));

                    return ang_a < ang_b;
                }
            ); // O(1)

            // Check if the centroid has to be included into the polygon
            if (!only_two_pts){
                auto centroid_it = cell_centroid_coords.find((*cur_cell_it).id()); // O(logN)
                if(centroid_it != cell_centroid_coords.end()){

                    std::vector<int> adjacent_nodes_on_boundary = adjacent_points_map.at((*cur_cell_it).id()); // O(logN)
                    std::vector<typename dcel_t::node_t*> on_boundary_pts;

                    for (auto pt_: pts){
                        for (int ad_node : adjacent_nodes_on_boundary){
                            bool areAligned = false;
                            auto areAligned_it = edge_vertexes2intersection.find(pt_->id()); // O(logN)
                            if(areAligned_it != edge_vertexes2intersection.end()){

                                areAligned = (edge_vertexes2intersection.at(pt_->id()) == std::make_pair((*cur_cell_it).id(), ad_node)) 
                                            || (edge_vertexes2intersection.at(pt_->id()) == std::make_pair(ad_node, (*cur_cell_it).id())); // O(logN)
                            }
                            if((ad_node != (*cur_cell_it).id()) && areAligned){
                                on_boundary_pts.push_back(pt_); // O(1)
                                break;
                            }
                        }
                    }

                    if (on_boundary_pts.size() < 2){
                        std::cout<<"Cell ID: "<<(*cur_cell_it).id()<<std::endl;
                        throw std::runtime_error("on_boundary_pts has size < 2");
                    }
                    
                    auto first_position  = std::find(pts.begin(), pts.end(), on_boundary_pts[0]); // O(1)
                    auto second_position = std::find(pts.begin(), pts.end(), on_boundary_pts[1]); // O(1)

                    if (std::distance(pts.begin(), first_position) >
                        std::distance(pts.begin(), second_position))
                    {
                        std::swap(first_position, second_position);
                    }

                    if (first_position != pts.end() && second_position != pts.end()){
                        if(first_position == pts.begin() && second_position == std::prev(pts.end())){
                            pts.push_back(cell_centroid_coords.at((*cur_cell_it).id()));
                            typename dcel_t::node_t* inserted_node = this->dcel_.insert_node(*cell_centroid_coords.at((*cur_cell_it).id()));
                        }
                        else{
                            pts.insert(std::next(first_position), cell_centroid_coords.at((*cur_cell_it).id()));
                            typename dcel_t::node_t* inserted_node = this->dcel_.insert_node(*cell_centroid_coords.at((*cur_cell_it).id()));
                        }
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
                auto search_it = vertexes2halfedge.find({(*pt)->id(), (*pt_next)->id()}); // O(logN)
                // If cur -> next halfedge already exists, retrieve it...
                if (search_it != vertexes2halfedge.end()){
                    cur_halfedge = search_it -> second;
                    cur_cell_it -> set_halfedge(cur_halfedge);
                }
                // ... otherwise create it
                else{
                    // Create new halfedge
                    typename dcel_t::halfedge_t* new_halfedge = this->dcel_.emplace_halfedge(*pt);
                    // Set it as halfedge for pt
                    (*pt) -> set_halfedge(new_halfedge);
                    // Add halfedge to the node_t structure
                    (*pt) -> add_halfedge((*pt_next)->id(), new_halfedge);

                    // If its twin already exists, connect them
                    auto twin_it = vertexes2halfedge.find({(*pt_next)->id(), (*pt)->id()}); // O(logN)
                    if(twin_it != vertexes2halfedge.end()){
                        typename dcel_t::halfedge_t* twin_halfedge = twin_it->second;
                        new_halfedge -> set_twin(twin_halfedge);
                        twin_halfedge -> set_twin(new_halfedge);
                    }
                    vertexes2halfedge.insert({std::make_pair((*pt)->id(), (*pt_next)->id()), new_halfedge}); // O(logN)

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
        this->dcel_.remove_nodes(node_ids_to_remove); // O(N)
    }   


    void boundary_cells_detection(const Triangulation<local_dim, embed_dim>& mesh, 
                                std::map<int, std::vector<typename dcel_t::node_t*>>& supplementary_points,
                                std::set<int>& on_boundary,
                                std::set<int>& out_of_boundary_d_centroids,
                                std::vector<int>& d_centroid2v_vertex,
                                int infty_id, std::set<int>& node_ids_to_remove,
                                std::map<int, typename dcel_t::node_t*>& cell_centroid_coords,
                                std::map<int, int>& infty_halfedges_midpoints,
                                std::map<int, typename simplex_t::NodeType>& midpoints_lookup,
                                std::map<int, std::pair<int, int>>& edge_vertexes2intersection){ // --> O(N^2)

        // Fill the structure containing the node ids to remove during the clipping (i.e. the infinity node and the ones out of the boundary)
        node_ids_to_remove.insert(infty_id); // O(log N)
        for(int d_centroid_id : out_of_boundary_d_centroids){ // O(1)
            int v_vertex_id = d_centroid2v_vertex.at(d_centroid_id); // O(1)
            node_ids_to_remove.insert(v_vertex_id); // O(logN)
        }

        if (mesh.boundary_edges_begin() == mesh.boundary_edges_end()){
            std::cerr << "Mesh does not have boundary edges, maybe there is no need for clipping..." << std::endl;
        }

        // Already checked data structure
        std::set<std::tuple<int, int, int>> already_checked;
        // Vector containing the differents components of the boundary
        std::vector<std::set<int>> boundaries = partition_boundary_edges(mesh); // O(N*logN)
        // Set counter for IDs
        int counter = infty_id + 1;

        for(std::set<int> boundary: boundaries){ // O(BUCHI)
            // Find the boundary_edge_iterator corrensponding to the first ID of the current component
            typename triangulation_t::boundary_edge_iterator d_boundary_edge = find_boundary_edge(mesh, boundary); // O(N* logN)
            // Consider one of its two endpoints, and the corresponding v_cell
            Eigen::Matrix<int, Dynamic, 1> node_ids = d_boundary_edge->node_ids();
            int d_vertex_id = node_ids(0);
            Eigen::Matrix<double, 2, 1> d_vertex = mesh.node(d_vertex_id);
            typename dcel_t::cell_t* v_cell = this->cells_map[d_vertex_id]; // O(1)

            // Cells queue
            std::list<int> to_check_v_cells;
            // Initialize to_check_v_cells
            to_check_v_cells.push_back(d_vertex_id); // O(1)

            // While queue is not empty
            while (to_check_v_cells.size() > 0){
                // Propagate towards the boundary to fill the supplementary points structure
                propagate_boundary_intersections(mesh, to_check_v_cells, already_checked,  infty_halfedges_midpoints, supplementary_points, 
                                                midpoints_lookup, on_boundary, edge_vertexes2intersection, infty_id, counter); // O(logN)

            }
        }
        // If the cell centroid is at the boundary and is not already present as a node in the DCEL, add it and set that it as a boundary node
        add_boundary_voronoi_site(mesh, cell_centroid_coords, counter); // O(N^2)
    }


    void add_boundary_voronoi_site(const Triangulation<local_dim, embed_dim>& mesh, 
                                std::map<int, typename dcel_t::node_t*>& cell_centroid_coords, int& counter){ // -> O(N^2)
        // If the cell centroid is at the boundary and is not already present as a node in the DCEL, add it and set that it as a boundary node
        for(auto cell_it = this->dcel_.cells_begin(); cell_it != this->dcel_.cells_end(); ++cell_it){ // O(N)
            int cur_cell_id = cell_it -> id();
            auto cur_centroid_coords = mesh.node(cur_cell_id);
            typename dcel_t::node_t* new_centroid_node;
            typename dcel_t::node_t* centroid_node = this->dcel_.find_node_by_coords(cur_centroid_coords); // O(N)  - BOTTLENECK
            if(centroid_node == nullptr){
                new_centroid_node = new typename dcel_t::node_t();
                new_centroid_node->set_coords(cur_centroid_coords);
                new_centroid_node->set_id(counter++);
            }
            else{
                new_centroid_node = centroid_node;
            }
            if(mesh.is_node_on_boundary(cur_cell_id)){
                new_centroid_node->set_boundary(true);
                cell_centroid_coords[cur_cell_id] = new_centroid_node; // O(logN)
            }
        }
    }


    triangulation_t::boundary_edge_iterator find_boundary_edge(const Triangulation<local_dim, embed_dim>& mesh,
                                            const std::set<int>& boundary){

        bool found_boundary_edge = false;
        auto d_boundary_edge = mesh.boundary_edges_begin(); 
        for(auto delaunay_edge = mesh.boundary_edges_begin(); delaunay_edge != mesh.boundary_edges_end(); ++delaunay_edge){ // O(N)
            int first_node_id = delaunay_edge->node_ids()(0); 
            // int second_node_id = delaunay_edge->node_ids()(1);
            if(boundary.find(first_node_id) != boundary.end()){ //O(log N)
                d_boundary_edge = delaunay_edge;
                found_boundary_edge = true;
                break;
            }
        }
        if(!found_boundary_edge){
            throw std::runtime_error("Mesh does not contain this boundary...");
        }

        return d_boundary_edge;
    }


    void propagate_boundary_intersections(const Triangulation<local_dim, embed_dim>& mesh, std::list<int>& to_check_v_cells,  
                                        std::set<std::tuple<int, int, int>>& already_checked,  std::map<int, int>& infty_halfedges_midpoints,
                                        std::map<int, std::vector<typename dcel_t::node_t*>>& supplementary_points, 
                                        std::map<int, typename simplex_t::NodeType>& midpoints_lookup, std::set<int>& on_boundary, 
                                        std::map<int, std::pair<int, int>>& edge_vertexes2intersection, int infty_id, int& counter){  // --> O(logN)
        // Get v_cell and d_vertex ID
        int d_vertex_id = to_check_v_cells.front();
        Eigen::Matrix<double, 2, 1> d_vertex = mesh.node(d_vertex_id); // O(1)
        typename dcel_t::cell_t* v_cell = this->cells_map[d_vertex_id];  // O(logN)
        to_check_v_cells.pop_front();

        // Neighbouring d_vertexes
        std::vector<int> neigh_d_vertexes_ids =  mesh.node_one_ring(d_vertex_id);  // O(1)
        // Loop over d_edges connected to that d_vertex
        for (int neigh_d_vertex_id : neigh_d_vertexes_ids){  // O(1)
            bool skippa = true;
            for (auto delaunay_edge = mesh.boundary_edges_begin(); delaunay_edge != mesh.boundary_edges_end(); ++delaunay_edge){  // O(BoundaryEdges) (sqrt(N) average)
                bool condition_to_not_skip =  (delaunay_edge->node_ids()(0) == d_vertex_id && delaunay_edge->node_ids()(1) == neigh_d_vertex_id) ||
                    (delaunay_edge->node_ids()(1) == d_vertex_id && delaunay_edge->node_ids()(0) == neigh_d_vertex_id);
                if(condition_to_not_skip){
                    skippa = false;
                }
            }
            if (skippa){
                continue;
            }
            coords_t neigh_coords = mesh.node(neigh_d_vertex_id);
            std::vector<typename dcel_t::halfedge_t*> v_cell_halfedges =  v_cell->cell_edges();  // O(1)
            bool intersected_once = false;
            // Loop over v_halfedges connected to that v_centroid
            for (typename dcel_t::halfedge_t* v_cell_halfedge : v_cell_halfedges){  // O(1)

                // If the triple (v_centroid, neigh_d_vertex, v_halfedge) is already present, skip
                if(already_checked.find(std::make_tuple(d_vertex_id, neigh_d_vertex_id, v_cell_halfedge->id())) != already_checked.end()){ // O(logN)
                    continue;
                }

                // Find halfedge points
                coords_t he_p1, he_p2;
                if (v_cell_halfedge->node()->id()==infty_id){
                    he_p1 = midpoints_lookup.at(infty_halfedges_midpoints.at(v_cell_halfedge->id())); // O(logN)
                    he_p2 = v_cell_halfedge->twin()->node()->coords();
                }
                else if (v_cell_halfedge->twin()->node()->id()==infty_id){
                    he_p1 = v_cell_halfedge->node()->coords();
                    he_p2 = midpoints_lookup.at(infty_halfedges_midpoints.at(v_cell_halfedge->id()));  // O(logN)
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
                );  // O(1)

                // Find ID of twin cell (most of the times it will coincide with neigh_d_vertex_id, but not always)
                int twin_cell_id = v_cell_halfedge->twin()->cell()->id();

                // If you find hafedge-edge intersection, add it to v_cell supplementary_points and twin_v_cell supplementary points
                if (check_intersection){

                    on_boundary.insert(d_vertex_id); // O(logN)
                    on_boundary.insert(neigh_d_vertex_id); // O(logN)
                    on_boundary.insert(twin_cell_id); // O(logN)
                    typename dcel_t::node_t* v_intersection_point;
                    auto intersection_ptr =  this->dcel_.find_node_by_coords(intersection_point);  // O(N)
                    if(intersection_ptr==nullptr){
                        v_intersection_point = new typename dcel_t::node_t();
                        v_intersection_point->set_coords(intersection_point);
                        v_intersection_point->set_id(counter++);
                    }
                    else{
                        v_intersection_point = intersection_ptr;
                    }
                    if(supplementary_points.find(d_vertex_id) == supplementary_points.end()){ // O(logN)
                        supplementary_points[d_vertex_id];
                    }
                    if(supplementary_points.find(twin_cell_id) == supplementary_points.end()){  // O(logN)
                        supplementary_points[twin_cell_id];
                    }
                    supplementary_points[d_vertex_id].push_back(v_intersection_point);  // O(logN)
                    supplementary_points[twin_cell_id].push_back(v_intersection_point);  // O(logN)
                    this->nodes_map[v_intersection_point->id()] = v_intersection_point;  // O(logN)
                    this->dcel_.insert_node(*v_intersection_point);  // O(1)
                    v_intersection_point->set_boundary(true);
                    intersected_once = true;

                    to_check_v_cells.push_back(twin_cell_id); // O(1)
                    to_check_v_cells.push_back(neigh_d_vertex_id);  // O(1)
                    edge_vertexes2intersection[v_intersection_point->id()] = std::make_pair(d_vertex_id, neigh_d_vertex_id);  // O(logN)
                }

                // Add to already_checked the triple (v_centroid, neigh_d_vertex, v_halfedge) and corresponding triple for the neighbouring v_cell
                already_checked.insert({d_vertex_id, neigh_d_vertex_id, v_cell_halfedge->id()});  // O(logN)
                already_checked.insert({neigh_d_vertex_id, d_vertex_id, v_cell_halfedge->twin()->id()});  // O(logN)
            }
        }
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


    void find_intersection_by_points(const Eigen::Matrix<double, 2, 1>& p1,
                                    const Eigen::Matrix<double, 2, 1>& p2,
                                    const Eigen::Matrix<double, 2, 1>& q1,
                                    const Eigen::Matrix<double, 2, 1>& q2,
                                    bool& check_intersection,
                                    typename simplex_t::NodeType& intersection_point){  // --> O(1)
        check_intersection = false;

        const double dx1 = p2(0) - p1(0);
        const double dy1 = p2(1) - p1(1);

        const double dx2 = q2(0) - q1(0);
        const double dy2 = q2(1) - q1(1);

        const double det = dx1 * dy2 - dy1 * dx2;

        // geometric scale of the problem
        const double len1 = std::sqrt(dx1*dx1 + dy1*dy1);
        const double len2 = std::sqrt(dx2*dx2 + dy2*dy2);

        const double EPS = 1e-10;
        const double tol = EPS * len1 * len2;

        // almost collinear segments
        if (std::abs(det) <= tol) {
            return;
        }

        const double rx = q1(0) - p1(0);
        const double ry = q1(1) - p1(1);

        const double t = (rx * dy2 - ry * dx2) / det;
        const double s = (rx * dy1 - ry * dx1) / det;

        const double param_tol = 1e-9;

        if (t >= -param_tol && t <= 1.0 + param_tol &&
            s >= -param_tol && s <= 1.0 + param_tol)
        {
            double tt = std::clamp(t, 0.0, 1.0);

            intersection_point(0) = p1(0) + tt * dx1;
            intersection_point(1) = p1(1) + tt * dy1;

            check_intersection = true;
        }
    }


    std::map<int, std::vector<int>> adjacent_points(const Triangulation<local_dim, embed_dim>& mesh){  // --> O(N*logN) (A dire il vero O(sqrt(N) * log(sqrt(N))))

        std::map<int, std::vector<int>> adjacent_points_map;
        
        // Loop over boundary edges
        for (auto d_boundary_edge = mesh.boundary_edges_begin(); d_boundary_edge != mesh.boundary_edges_end(); ++d_boundary_edge){  // O(sqrt(N))

            // Define first and second endpoints
            Eigen::Matrix<int, Dynamic, 1> endpoints_matrix = d_boundary_edge -> node_ids();
            int first_endpoint_id = endpoints_matrix(0);
            int second_endpoint_id = endpoints_matrix(1);

            // If first endpoint already in map, add second endpoint
            if(adjacent_points_map.find(first_endpoint_id) != adjacent_points_map.end()){  // O(logN)
                adjacent_points_map[first_endpoint_id].push_back(second_endpoint_id);  // O(logN)
            }
            // Else, create its entry and add second endpoint
            else{
                adjacent_points_map[first_endpoint_id];  // O(logN)
                adjacent_points_map[first_endpoint_id].push_back(second_endpoint_id);  // O(logN)
            }

            // Repeat for second endpoint
            if(adjacent_points_map.find(second_endpoint_id) != adjacent_points_map.end()){  // O(logN)
                adjacent_points_map[second_endpoint_id].push_back(first_endpoint_id);  // O(logN)
            }
            else{
                adjacent_points_map[second_endpoint_id];  // O(logN)
                adjacent_points_map[second_endpoint_id].push_back(first_endpoint_id);  // O(logN)
            }
        }

        return adjacent_points_map;
    }


    // Function to check if a point lies in a cell
    const bool is_point_in_cell(const coords_t& point, int cell_id){
        
        // Retrieve cell
        auto cell = this->cells_map[cell_id];

        // Compute its points
        std::vector<node_t*> cell_nodes_vector = cell->cell_nodes();

        // Collect them inside a proper matrix
        Eigen::Matrix<double, Dynamic, embed_dim> cell_nodes(cell_nodes_vector.size(), 2);
        for (int i=0; i < cell_nodes_vector.size(); ++i){
            cell_nodes.row(i) = cell_nodes_vector[i]->coords();
        }

        // Check condition
        return internals::point_in_2d_polygon(cell_nodes, point);
    }


    // Find Voronoi centroid closest to a given point
    const int find_closest_v_centroid(const coords_t& point){
        // Closest v_centroid
        int min_id = 0;

        // Point coordinates
        double x_point = point(0);
        double y_point = point(1);

        // Prohibited centroids
        std::set<int> prohibited_v_centroids;

        while(true){
            double min_distance = std::numeric_limits<double>::max();
            for (auto v_cell_id = 0; v_cell_id < this->v_cell2v_centroid.size(); ++v_cell_id){
                if(prohibited_v_centroids.find(v_cell_id) != prohibited_v_centroids.end()){
                    continue;
                }
                // Store data
                typename simplex_t::NodeType v_centroid = this->v_cell2v_centroid[v_cell_id];

                // Centroid coordinates
                double x_centroid = v_centroid(0);
                double y_centroid = v_centroid(1);

                // Compute differences
                double dx = x_point - x_centroid;
                double dy = y_point - y_centroid;

                // Compute distance
                double distance = std::sqrt(dx * dx + dy * dy);

                // Check if it "beats the one before"
                if (distance < min_distance){
                    min_distance = distance;
                    min_id = v_cell_id;
                }
            }

            // Loop over halfedges
            bool found_one = false;
            for (auto halfedge = this->dcel_.halfedges_begin(); halfedge != this->dcel_.halfedges_end(); ++halfedge){
                if (!(halfedge->on_boundary())) continue;
                // Verify if there is intersection
                Eigen::Matrix<double, 2, 1> p1 = this->v_cell2v_centroid[min_id];
                Eigen::Matrix<double, 2, 1> p2 = point;
                Eigen::Matrix<double, 2, 1> q1 = halfedge->node()->coords();
                Eigen::Matrix<double, 2, 1> q2 = halfedge->next()->node()->coords();
                bool check_intersection = false;
                typename simplex_t::NodeType intersection_point;
                find_intersection_by_points(p1, p2, q1, q2, check_intersection, intersection_point);
                if (check_intersection && 
                    (std::abs(intersection_point(0)-p1(0))>1e-6 && std::abs(intersection_point(1)-p1(1))>1e-6)){
                    // std::cout << "Maleeeee" << std::endl;
                    prohibited_v_centroids.insert(min_id);
                    found_one = true;
                    break;
                }
            }
            if (!found_one){
                // std::cout << "Beneeeee" << std::endl;
                return min_id;
            }
        }
        std::cerr << "Did not find any valid v_centroid" << std::endl;
    }


    };
}
#endif