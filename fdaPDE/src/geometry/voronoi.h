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
    // counter to set the IDs
    int counter = 0;
    // Lookup for midpoints coordinates
    std::map<int, typename simplex_t::NodeType> midpoints_lookup_raw;
    // Current ID
    int cur_id = 0;
    // Structure to store cells without finite halfedges
    std::set<int> not_created_cells;
    // Resize d_centroid2v_vertex
    d_centroid2v_vertex.resize(n_mesh_faces);
    // Resize v_cell2v_centroid
    v_cell2v_centroid.resize(n_mesh_vertices);
    
    std::vector<std::vector<int>> node_neighbors_lookup = v_vertex_computation_(mesh, 
                                                                                dcel_, 
                                                                                centroid_lookup,
                                                                                midpoints_lookup_raw,
                                                                                cur_id); 

    auto start = std::chrono::high_resolution_clock::now();

    // Define infinity v_vertex and add it to DCEL
    this->infty_id = cur_id;
    typename simplex_t::NodeType v_vertex_infty = simplex_t::NodeType::Constant(std::numeric_limits<double>::infinity());
    centroid_lookup[this->infty_id] = v_vertex_infty; 
    typename dcel_t::node_t* infty_node = dcel_.insert_node(typename dcel_t::node_t(this->infty_id, false, v_vertex_infty)); 
    nodes_map[this->infty_id] = infty_node; 

    // Imagine to have the correspondence cell_id: centroid coordinates
    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it) {  

        // Retrieve the old ID of the current cell
        int old_cell_id = it->id();
        // Retrieve new ID of the current cell
        int cell_id = this->d_centroid2v_vertex.at(old_cell_id);  
        // Order the neighbouring cells in a counter-clockwise way
        std::vector<int> ordered_neigh_ids = order_neighbours_cclw_(cell_id, node_neighbors_lookup, 
                                                                    centroid_lookup); 
        
        // Create and connect the halfedges
        if(!visited_centroids[cell_id]){
            create_and_connect_halfedges_(dcel_, cell_id, infty_node, centroid_lookup, 
                                        midpoints_lookup_raw, ordered_neigh_ids);  
            // Set the centroid as visited
            visited_centroids[cell_id] = true;
        }

        auto neighbor_simplexes = it->neighbors(); 
        // Create cells
        create_cells_(dcel_, mesh, cell_id, old_cell_id, 
                    centroid_lookup, neighbor_simplexes,  not_created_cells);  
    }
    // Add to the DCEL structure the not created cells
    add_not_created_cells_(not_created_cells, mesh); 

    // Connect the infinity halfedges
    connect_infty_halfedges_(); 

    // Update halfedges-cells structure in DCEL
    this->dcel_.update_halfedges_with_cells(); 

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout<<"computational times: "<< duration.count() << std::endl;

    std::cout << "\nFinished Voronoi constructor" << std::endl;
    }


    // Cell iterators
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

    // number of nodes
    int n_nodes() { return dcel_.n_nodes(); }

    // matrix of int (n_edge x 2) where i-th row: [id_node_1 i-th edge, id_node_2 i-th edge]
    Eigen::Matrix<int, Eigen::Dynamic, 2> edges() const {
        // calls edges() from dcel_t
        return dcel_.edges();
    }

    // number of edges
    int n_edges() { return dcel_.n_edges(); }

    // number of cells
    int n_cells() { return dcel_.n_cells(); }

    void export_to_json(const std::string& filename){
        this->dcel_.export_to_json(filename);
    }

    private:
    // Internal function for computing v_vertex neighbours
    std::vector<int> compute_neighbours_(const Triangulation<local_dim, embed_dim>& mesh,
        const triangle_t& d_simplex, std::map<int, typename simplex_t::NodeType>& midpoints_lookup_raw,
        int& cur_midpoint_id, bool print){
                                        
        // Get neighbours
        Eigen::Matrix<int, Eigen::Dynamic, 1> neigh = d_simplex.neighbors();
        // Store them inside a vector
        std::vector<int> neigh_vec(neigh.data(), neigh.data() + neigh.size()); 
        // Midpoint IDs
        std::map<int, int> mp_ids;
        // Incremental counter for mp_ids
        int counter = 0;

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
                        this->out_of_boundary_d_centroids.insert(d_simplex.id());
                    }

                    this->midpoints_lookup[cur_midpoint_id] = mp; 

                    // Add in the added midpoint ids
                    mp_ids[counter] = cur_midpoint_id;
                    counter += 1;
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
    int v_vertex_from_d_centroid_(const std::vector<typename simplex_t::NodeType>& v_vertex_lookup, 
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
                            std::vector<typename simplex_t::NodeType>& v_vertex_lookup,
                            std::map<int, typename simplex_t::NodeType>& midpoints_lookup_raw,
                            int& cur_id){  

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
            std::vector<int> neigh_vec = compute_neighbours_(mesh, *(d_simplex_it), midpoints_lookup_raw, cur_midpoint_id, false); 
            
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
                typename dcel_t::node_t* inserted_node = dcel_.insert_node(v_vertex);  
                nodes_map[v_vertex_id] = inserted_node; 

                // Update all structures
                v_vertex_neighbours_lookup[cur_id] = neigh_vec; 
                first_d_centroid[cur_id] = d_centroid_id; 
                v_vertex_lookup[cur_id] = d_centroid; 

                // Increment current ID
                cur_id += 1;
            }

            // Update d_centroid2v_vertex
            this->d_centroid2v_vertex[d_centroid_id] = v_vertex_id; 

        }

        return v_vertex_neighbours_lookup;
    }

    int retrieve_v_vertex_id_(int d_centroid_id){

        // If positive, than it is a physical d_centroid
        if (d_centroid_id >= 0) return (this->d_centroid2v_vertex[d_centroid_id]);
        // Else it is a midpoint
        return d_centroid_id;
    }

    typename simplex_t::NodeType id2point_(
                                        int id,
                                        std::vector<typename simplex_t::NodeType>& v_vertex_lookup
    ){ 
        // If positive, then it is a physical d_centroid
        if (id >= 0) return v_vertex_lookup.at(id);

        // Else it is a midpoint
        return this->midpoints_lookup.at(id); 
    } 

    std::vector<int> order_neighbours_cclw_(
                                        int v_vertex_id,
                                        std::vector<std::vector<int>>& v_vertex_neighbours_lookup,
                                        std::vector<typename simplex_t::NodeType>& v_vertex_lookup
    ){ 

        // Retrieve v_vertex_coords
        typename simplex_t::NodeType v_vertex = id2point_(v_vertex_id, v_vertex_lookup); 

        // Create ordered neighbours structure
        std::vector<int> ordered_neigh_ids;
        ordered_neigh_ids.reserve(v_vertex_neighbours_lookup.at(v_vertex_id).size()); 

        // Loop over neighbours
        for (auto d_centroid_neigh_id : v_vertex_neighbours_lookup.at(v_vertex_id)) {  

            int neigh_id = retrieve_v_vertex_id_(d_centroid_neigh_id);  
            
            ordered_neigh_ids.push_back(neigh_id); 
        }

        std::sort(ordered_neigh_ids.begin(), ordered_neigh_ids.end(),
            [&](int id_a, int id_b)
        {
            typename simplex_t::NodeType pa = id2point_(id_a, v_vertex_lookup);
            typename simplex_t::NodeType pb = id2point_(id_b, v_vertex_lookup);

            double ang_a = std::atan2(pa.y() - v_vertex.y(),
                                    pa.x() - v_vertex.x());
            double ang_b = std::atan2(pb.y() - v_vertex.y(),
                                    pb.x() - v_vertex.x());

            return ang_a > ang_b;
        }); 

        return ordered_neigh_ids;
    }

    void create_and_connect_halfedges_(dcel_t& dcel_, 
                            int v_vertex_id,
                            typename dcel_t::node_t* infty_node,
                            std::vector<typename simplex_t::NodeType>& v_vertex_lookup,
                            std::map<int, typename simplex_t::NodeType>& midpoints_lookup_raw,
                            const std::vector<int>& ordered_neigh_ids) { 
                 
        int count_existing_neigh = 0;
        // Keep the first halfedge
        typename dcel_t::halfedge_t* first_halfedge=nullptr;
        // Pointers to e1 and e2 initialized as null_pointer
        typename dcel_t::halfedge_t* e1 = nullptr;
        typename dcel_t::halfedge_t* e2 = nullptr;
        // Loop over neighbouring cells
        for (auto neigh_id : ordered_neigh_ids) {  
            // flag for infinity node
            bool is_infinity = false;
            // id for computing the midpoint id
            int mid_id = 0;
            // Create cell node
            typename dcel_t::node_t* v_vertex;
            // Find nodes in the DCEL structure
            v_vertex = nodes_map.at(v_vertex_id); 
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
                twin_node = nodes_map.at(neigh_id); 
            }
            else{
                // The centroid of the neighbouring cell is the infinity node
                mid_id = neigh_id;
                neigh_id = this->infty_id;
                twin_node = infty_node;
                is_infinity = true;
            }

            // If neighbouring cell was not visited, its node needs an halfedge
            if ((twin_node->neighID2halfedge(v_vertex_id)==nullptr) or (neigh_id == this->infty_id)){
                // Set IDs based on counter
                cell_halfedge = dcel_.emplace_halfedge(v_vertex);  
                twin_halfedge = dcel_.emplace_halfedge(twin_node);  

                // Set halfedge to the cell centroid
                v_vertex->set_halfedge(cell_halfedge);
                // Set halfedge to the neighbouring centroid
                twin_node->set_halfedge(twin_halfedge);

                cell_halfedge->set_twin(twin_halfedge);
                twin_halfedge->set_twin(cell_halfedge);

                // Update lookup
                v_vertex->add_halfedge(neigh_id, cell_halfedge);
                twin_node->add_halfedge(v_vertex_id, twin_halfedge);
                this->vertexes2halfedge.insert({std::make_pair(v_vertex->id(), twin_node->id()), cell_halfedge}); 
                this->vertexes2halfedge.insert({std::make_pair(twin_node->id(), v_vertex->id()), twin_halfedge}); 

                // Update intersection d_edges map
                if((is_infinity) and midpoints_lookup_raw.at(mid_id)==this->midpoints_lookup.at(mid_id)){ 
                    intersection_d_edges_[cell_halfedge->id()]= midpoint_to_edge_.at(mid_id); 
                    intersection_d_edges_[twin_halfedge->id()] = midpoint_to_edge_.at(mid_id); 
                }
                if (is_infinity){
                    this->infty_halfedges_midpoints[cell_halfedge->id()] = mid_id; 
                    this->infty_halfedges_midpoints[twin_halfedge->id()] = mid_id; 
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
                        std::vector<typename simplex_t::NodeType>& v_vertex_lookup,
                        const Eigen::Matrix<int, Eigen::Dynamic, 1>& neighbor_simplexes,
                        std::set<int>& not_created_cells) {  
        
        int count_neigh_tria = 0;
        int diff_id = 0;
        // Retrieve the IDs of the vertices of the current cell
        triangle_t curr_tria(d_centroid_id, &mesh); 
        Eigen::Matrix<int, Dynamic, 1> d_vertex_ids = curr_tria.node_ids();
        // Retrieve the centroid of the current cell
        typename simplex_t::NodeType v_vertex = v_vertex_lookup.at(v_vertex_id);  
        
        for (int neigh_d_centroid_id : neighbor_simplexes){  
            if(neigh_d_centroid_id != -1){
                count_neigh_tria += 1;
                // Retrieve the triangle in the mesh using the old id
                triangle_t tria(neigh_d_centroid_id, &mesh); 
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
                for(int i=0; i<common.size(); ++i){  

                    int idx = common[i];

                    // Retrieve the corrensponding coordinates
                    typename simplex_t::NodeType d_vertex = mesh.node(idx); 

                    int neigh_v_vertex_id = this->d_centroid2v_vertex.at(neigh_d_centroid_id); 
                    typename simplex_t::NodeType neigh_v_vertex = v_vertex_lookup.at(neigh_v_vertex_id);  
                    // Create the new cell
                    cell_t curr_cell(idx);

                    // Store v_cell - v_centroid information
                    v_cell2v_centroid[idx] = d_vertex;  

                    // create segments around the vertex
                    typename simplex_t::NodeType u = v_vertex - d_vertex;
                    typename simplex_t::NodeType v = neigh_v_vertex - d_vertex;
                    // Compute the cross product to assign the halfedges properly 
                    double cross = u(0)*v(1)-u(1)*v(0);
                    typename dcel_t::halfedge_t* current_cell_halfedge = this->vertexes2halfedge.find(std::make_pair(v_vertex_id, neigh_v_vertex_id))->second;  
                    typename dcel_t::halfedge_t* current_twin_halfedge = this->vertexes2halfedge.find(std::make_pair(neigh_v_vertex_id, v_vertex_id))->second;  

                    if (mesh.is_node_on_boundary(idx)){ 
                        curr_cell.set_unbounded();
                    }

                    // Add the new cell to the list if it is not yet added 
                    if(cells_map.find(curr_cell.id()) == cells_map.end()){   
                        if(cross < -1e-9){
                            // cross product < tolerance means that neigh_centroid is clock-wise with respect to centroid, so we want neigh -> centroid
                            curr_cell.set_halfedge(current_twin_halfedge);
                            typename dcel_t::cell_t* inserted_cell = dcel_.insert_cell(curr_cell); 
                            cells_map[curr_cell.id()] = inserted_cell;
                            not_created_cells.erase(idx); 
                        }
                        else if(cross > 1e-9){
                            // cross product > tolerance means that neigh_centroid is CCW with respect to centroid, so we want centroid -> neigh
                            curr_cell.set_halfedge(current_cell_halfedge);
                            typename dcel_t::cell_t* inserted_cell = dcel_.insert_cell(curr_cell);
                            cells_map[curr_cell.id()] = inserted_cell;
                            not_created_cells.erase(idx);
                        }
                        else{
                            if (not_created_cells.find(idx) == not_created_cells.end()){ 
                                not_created_cells.insert(idx); 
                            }    
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
        if(cells_map.find(diff_id) == cells_map.end()){ 
            not_created_cells.insert(diff_id); 
        }
    }
    }

    std::map<int, std::vector<int>> compute_v_cell2halfedges_(){
        std::map<int, std::vector<int>> v_cell2halfedges;
        for(auto cell_it = dcel_.cells_cbegin(); cell_it != dcel_.cells_cend(); ++cell_it){ 
            std::vector<int> halfedges_vector;
            if(!cell_it->is_unbounded()){
                for(auto v: cell_it->cell_edges()){
                    halfedges_vector.push_back(v->id()); 
                }
            }
            else{
                for(auto v: cell_it->cell_edges_with_infty(this->infty_id, false)){ 
                    halfedges_vector.push_back(v->id()); 
                }
            }

            v_cell2halfedges[cell_it->id()] = halfedges_vector; 
        }

        return v_cell2halfedges;
    }


    void add_not_created_cells_(std::set<int>& not_created_cells,
                                const Triangulation<local_dim, embed_dim>& mesh){ 

        std::map<int, std::vector<int>> v_cell2halfedges = compute_v_cell2halfedges_(); 

        if (not_created_cells.size()>0){

            for(int internal_id : not_created_cells){ 
                cell_t curr_cell(internal_id);
                if(cells_map.find(internal_id) != cells_map.end()){ 
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
                int v_vertex_id = this->d_centroid2v_vertex[v_vertexes[0]]; 
                
                bool breaking_check = false;

                for(int halfedge_id : v_cell_halfedges){ 
                    auto to_infty_list = this->vertexes2halfedge.equal_range(std::make_pair(v_vertex_id, this->infty_id)); 
                    auto from_infty_list = this->vertexes2halfedge.equal_range(std::make_pair(this->infty_id, v_vertex_id)); 
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

                typename dcel_t::cell_t* inserted_cell = dcel_.insert_cell(curr_cell); 
                cells_map[internal_id] = inserted_cell;  

                // Store v_cell - v_centroid information
                v_cell2v_centroid[curr_cell.id()] = mesh.node(curr_cell.id()); 
            }
        }
    }


    void connect_infty_halfedges_(){ 
        std::map<int, std::vector<int>> v_cell2halfedges = compute_v_cell2halfedges_(); 

        std::map<int, typename dcel_t::halfedge_t*> halfedges_map;
        for(auto it = dcel_.halfedges_begin(); it != dcel_.halfedges_end(); ++it){ 
            halfedges_map[it->id()] = &(*it); 
        }

        for (const auto& [v_cell_id, halfedges_ids] : v_cell2halfedges) { 
            typename dcel_t::halfedge_t* from_infty;
            typename dcel_t::halfedge_t* to_infty;
            for (int halfedge_id : halfedges_ids) {
                typename dcel_t::halfedge_t* cur_halfedge = halfedges_map[halfedge_id]; 
                if (cur_halfedge->node()->id()==this->infty_id){
                    from_infty = cur_halfedge;
                }
                else if (cur_halfedge->twin()->node()->id()==this->infty_id){
                    to_infty = cur_halfedge;
                }
            }

            to_infty->set_next(from_infty);
            from_infty->set_prev(to_infty);
        }
    }   


    protected:
    // DCEL data structure
    dcel_t dcel_;
    // Map of cells
    std::map<int, typename dcel_t::cell_t*> cells_map;
    // Midpoint id to egde id
    std::map<int, int> midpoint_to_edge_;
    // Intersection d_edges
    std::map<int, int> intersection_d_edges_;
    // Infty id
    int infty_id;
    // Set to store d_centroids out of the boundary
    std::set<int> out_of_boundary_d_centroids;
    // Conversion from Delaunay centroid ids to Voronoi vertex ids
    std::vector<int> d_centroid2v_vertex;
    // Multimap that associated an halfedge to the corrensponding vertexes
    std::multimap<std::pair<int, int>, typename dcel_t::halfedge_t*> vertexes2halfedge;
    std::map<int, typename simplex_t::NodeType>  midpoints_lookup;
    // Map to store infinity halfedge -> midpoint ID (for understanding infinity edges intersections)
    std::map<int, int> infty_halfedges_midpoints;
    // Map associating to each v_cell the corresponding v_centroid
    std::vector<typename simplex_t::NodeType> v_cell2v_centroid;
    // map containing the nodes
    std::map<int, typename dcel_t::node_t*> nodes_map;

};

}


#endif