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
    // using delaunay_t = Delaunay<local_dim, embed_dim>;
    using cell_t = dcel_t::cell_t;
    
    public:

    // Voronoi(Matrix seed)
    // Ho scelto questi ingressi basandomi sul constructor in delaunay.h che usa random generated points e no refinement
    // fare il delaunay dei seed O(n log(n)) --> I think that we can rely on the constructor in delaunay.h with random generated points, without no refinement
    // chiamo quello sotto (faccio il duale)
    // Voronoi(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& boundaries, int N=0, const std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>>& holes = {{}}) 
    // : Voronoi(typename myDelaunay(boundaries, N, holes))
    // {};


    Voronoi(const Triangulation<local_dim, embed_dim>& mesh) {

    // Number of faces in the mesh 
    int n_mesh_faces = mesh.n_cells();
    // Number of vertices in the mesh
    int n_mesh_vertices = mesh.n_nodes();
    // Centroids lookup
    std::vector<typename simplex_t::NodeType> centroid_lookup(n_mesh_faces+1);
    // Visited or not
    std::vector<bool> visited_centroids(n_mesh_faces+1, false);
    // Boundary cell or not
    std::vector<bool> on_boundary(n_mesh_faces+1);
    // Old node ID : new node ID
    std::vector<int> old2new(n_mesh_faces);
    // counter to set the IDs
    int counter = 0;
    // Same centroids lookup
    std::vector<int> first_id_lookup(n_mesh_faces);
    // Lookup for midpoints coordinates
    std::map<int, typename simplex_t::NodeType> midpoints_lookup;
    // Map that associates a pair of centroids to the halfedge that connects them
    std::map<std::pair<int, int>, typename dcel_t::halfedge_t*> vertexes2halfedge;

    int cur_id = 0;

    std::vector<std::vector<int>> node_neighbors_lookup = v_vertex_computation_(mesh, 
                                                                                dcel_, 
                                                                                old2new, 
                                                                                centroid_lookup,
                                                                                midpoints_lookup,
                                                                                cur_id);

    // Define infinity v_vertex and add it to DCEL
    int infty_id = cur_id;
    typename simplex_t::NodeType v_vertex_infty = simplex_t::NodeType::Constant(std::numeric_limits<double>::infinity());
    centroid_lookup[infty_id] = v_vertex_infty;
    typename dcel_t::node_t* infty_node =
    dcel_.insert_node(
        typename dcel_t::node_t(infty_id, false, v_vertex_infty)
    );

// ==========================================================================================================================================
    std::unordered_set<int> not_created_cells;
    // Imagine to have the correspondence cell_id: centroid coordinates
    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it) {

        // Retrieve the IDs of the vertices of the current cell
        Eigen::Matrix<int, Dynamic, 1> ids_node = it->node_ids();
        // Retrieve the old ID of the current cell
        int old_cell_id = it->id();
        // Retrieve new ID of the current cell
        int cell_id = old2new.at(old_cell_id);

        // Check if the cell is on the boundary
        if (it->on_boundary()){
            on_boundary[cell_id] = true;
        }

        // Retrieve circumcenter
        typename simplex_t::NodeType centroid = centroid_lookup.at(cell_id);

        // Keep the first halfedge
        typename dcel_t::halfedge_t* first_halfedge=nullptr;

        // Pointers to e1 and e2 initialized as null_pointer
        typename dcel_t::halfedge_t* e1 = nullptr;
        typename dcel_t::halfedge_t* e2 = nullptr;

        int count_existing_neigh = 0;
        int diff_id = 0;
        // Order the neighbouring cells in a counter-clockwise way
        std::vector<int> ordered_neigh_ids = order_neighbours_cclw_(cell_id, node_neighbors_lookup,
        old2new, centroid_lookup, midpoints_lookup);

        // std::cout << "Neighbours of 2 " << std::endl;
        // for(auto elimina : node_neighbors_lookup.at(2)){
        //     std::cout << "Old " << elimina << std::endl;
        //     if (elimina >= 0){
        //         std::cout << "New " << old2new.at(elimina) << std::endl;
        //     }
        //     else{
        //         std::cout << "New " << midpoints_lookup.at(elimina) << std::endl;
        //     }
        // }


        // ordered_neigh_ids.reserve(node_neighbors_lookup.at(cell_id).size());

        // for (auto old_neigh_cell_id : node_neighbors_lookup.at(cell_id)) {

        //     int new_neigh_id = retrieve_v_vertex_id_(old_neigh_cell_id, old2new);
            
        //     ordered_neigh_ids.push_back(new_neigh_id);
        // }

        // std::sort(ordered_neigh_ids.begin(), ordered_neigh_ids.end(),
        //     [&](int id_a, int id_b)
        // {
        //     typename simplex_t::NodeType pa = id2point_(id_a, midpoints_lookup, centroid_lookup);
        //     typename simplex_t::NodeType pb = id2point_(id_b, midpoints_lookup, centroid_lookup);

        //     double ang_a = std::atan2(pa.y() - centroid.y(),
        //                             pa.x() - centroid.x());
        //     double ang_b = std::atan2(pb.y() - centroid.y(),
        //                             pb.x() - centroid.x());

        //     return ang_a < ang_b;
        // });

        // Loop over neighbouring cells
        for (auto neigh_cell_id : ordered_neigh_ids) {
            // Define the variable for the neigh cell ID
            // int neigh_cell_id;
            // Create cell node
            typename dcel_t::node_t* cell_node;
            // Find nodes in the DCEL structure
            cell_node = dcel_.find_node(dcel_.nodes().row(cell_id));
            // Create twin node 
            typename dcel_t::node_t* twin_node;
            // Create halfedge centroid -> neigh
            typename dcel_t::halfedge_t* cell_halfedge = new dcel_t::halfedge_t();
            // Create twin halfedge (neigh -> centroid)
            typename dcel_t::halfedge_t* twin_halfedge = new dcel_t::halfedge_t();

            if(neigh_cell_id >= 0){
                // Retrieve new ID for neighbour
                // neigh_cell_id = old2new.at(old_neigh_cell_id);
                if (neigh_cell_id == cell_id){
                    continue;
                }
                // Retrieve centroid of the neighbouring cell
                typename simplex_t::NodeType neigh_centroid = centroid_lookup.at(neigh_cell_id);
                twin_node = dcel_.find_node(dcel_.nodes().row(neigh_cell_id));
            }
            else{
                // The centroid of the neighbouring cell is the infinity node
                neigh_cell_id = infty_id;
                twin_node = infty_node;
            }

            // If neighbouring cell was not visited, its node needs an halfedge
            // if (!visited_centroids.at(neigh_cell_id))
            if (twin_node->neighID2halfedge(cell_id)==nullptr){
                // Set IDs based on counter
                cell_halfedge = dcel_.emplace_halfedge(cell_node);
                twin_halfedge = dcel_.emplace_halfedge(twin_node);

                // Set halfedge to the cell centroid
                cell_node->set_halfedge(cell_halfedge);
                // Set halfedge to the neighbouring centroid
                twin_node->set_halfedge(twin_halfedge);

                cell_halfedge->set_twin(twin_halfedge);
                twin_halfedge->set_twin(cell_halfedge);

                // Update lookup
                cell_node->add_halfedge(neigh_cell_id, cell_halfedge);
                twin_node->add_halfedge(cell_id, twin_halfedge);

                vertexes2halfedge[std::make_pair<int, int>(cell_node->id(), twin_node->id())] = cell_halfedge;
                vertexes2halfedge[std::make_pair<int, int>(twin_node->id(), cell_node->id())] = twin_halfedge;
            }
            // If neighbouring cell was not visited, only retrieve the halfedges
            else{
                delete cell_halfedge;
                delete twin_halfedge;
                twin_halfedge = (twin_node->neighID2halfedge(cell_id));
                cell_halfedge = (twin_halfedge->twin());
            }
        
            // Set the current cell as visited
            visited_centroids[cell_id] = true;
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
                
    // Complete the association of previous and next
    if (count_existing_neigh != 1)
    {
    e2->twin()->set_next(first_halfedge);
    first_halfedge->set_prev(e2->twin());
    }
    }























    auto neighbor_simplexes = it->neighbors();

    // NELL'ALTRO CONSIDERIAMO TUTTI I NEIGHS
    int count_neigh_tria = 0;
    for (int old_neigh_cell_id : neighbor_simplexes){
        if(old_neigh_cell_id != -1){
                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                count_neigh_tria += 1;
                // Retrieve the triangle in the mesh using the old id
                triangle_t tria(old_neigh_cell_id, &mesh);  
                auto ids_neigh_node = tria.node_ids();

                // Find the common and different IDs between these vertices and the cell ones
                std::vector<int> vec_1(ids_node.data(), ids_node.data()+ids_node.size());
                std::vector<int> vec_2(ids_neigh_node.data(), ids_neigh_node.data()+ids_neigh_node.size());

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

                // Associate the halfedges to the cells. IDEA: consideriamo le coordinare dei due centroidi e dei due vertici
                // e assegniamo gli halfedges in senso antiorario
                for(int i=0; i<common.size(); ++i){
                    // find the element into the vector 
                    auto found_element = std::find(vec_2.begin(), vec_2.end(), common[i]);
                    // check 
                    assert(found_element != vec_2.end());
                    // find the corrensponding index
                    //int local_idx = std::distance(vec_2.begin(), found_element);
                    //typename simplex_t::NodeType ver_coords = tria.node(local_idx);
                    //SECONDO ME ERA SBAGLIATO QUELLO
                    typename simplex_t::NodeType ver_coords = mesh.node(common[i]);

                    // orientation check
                    int new_neigh_id = old2new.at(old_neigh_cell_id);
                    typename simplex_t::NodeType neigh_centroid_coords = centroid_lookup.at(new_neigh_id);
                    // Create the new cell
                    cell_t curr_cell(common[i]);

                    // create segments around the vertex
                    typename simplex_t::NodeType u = centroid - ver_coords;
                    typename simplex_t::NodeType v = neigh_centroid_coords - ver_coords;
                    // cross product
                    double cross = u(0)*v(1)-u(1)*v(0);

                    // if(cell_id == new_neigh_id){
                    //     not_created_cells.insert(common[i]);
                    //     continue;
                    // }

                    typename dcel_t::halfedge_t* current_cell_halfedge = vertexes2halfedge[std::make_pair(cell_id, new_neigh_id)];
                    typename dcel_t::halfedge_t* current_twin_halfedge = vertexes2halfedge[std::make_pair(new_neigh_id, cell_id)];

                    if (mesh.is_node_on_boundary(common[i])){
                        curr_cell.set_unbounded();
                    }
                    // if(current_cell_halfedge == nullptr & current_twin_halfedge == nullptr){
                    //     std::cout<<"Cell id: "<< cell_id << "  Neigh cell id: "<<new_neigh_id<<std::endl;
                    //     std::cout<<"Cell coords: "<< centroid << "  Neigh cell coords: "<<neigh_centroid_coords<<std::endl;
                    // }
                    // Add the new cell to the list if it is not yet added 

                    if (std::find(cells_.begin(), cells_.end(), curr_cell) == cells_.end()){
                        if(cross < -1e-9){
                            // cross product < 0 means that neigh_centroid is clock-wise with respect to centroid, so we want neigh -> centroid
                            curr_cell.set_halfedge(current_twin_halfedge);
                            cells_.push_back(curr_cell);
                            not_created_cells.erase(common[i]);
                        }
                        else if(cross > 1e-9){
                            // cross product > 0 means that neigh_centroid is CCW with respect to centroid, so we want centroid -> neigh
                            curr_cell.set_halfedge(current_cell_halfedge);
                            cells_.push_back(curr_cell);
                            not_created_cells.erase(common[i]);
                        }
                        else{
                            not_created_cells.insert(common[i]);
                        }
                        }
                }

                // Take the different id between the two
                std::vector<int> tmp_diff;

                // Element that is in vec_1 but not in vec_2
                std::set_difference(
                    vec_1.begin(), vec_1.end(),
                    vec_2.begin(), vec_2.end(),
                    std::back_inserter(tmp_diff)
                );

                // Assert that the different element is one since they are adjacent cells
                assert(tmp_diff.size()==1);
                diff_id = tmp_diff[0];
                

                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
    }

    if(count_neigh_tria == 1){
        cell_t curr_cell_diff(diff_id);
        curr_cell_diff.set_unbounded();
        if (std::find(cells_.begin(), cells_.end(), curr_cell_diff) == cells_.end()){
            cells_.push_back(curr_cell_diff);
        }
    }

    }

    for(auto it = not_created_cells.begin(); it != not_created_cells.end(); ++it){
        cell_t curr_cell(*it);
        if (std::find(cells_.begin(), cells_.end(), curr_cell) != cells_.end()){
            continue;
        }
        if (mesh.is_node_on_boundary(*it)){
            curr_cell.set_unbounded();
        }
        cells_.push_back(curr_cell);
    }

// ==========================================================================================================================================
    int n_cells = dcel_.n_cells();
    std::cout << "\nFinished constructor" << std::endl;

    }














































    // iteratori sulle celle
    // in futuro vorremmo poter fare
    // for(auto it = voronoi.cells_begin(); it != voronoi.cells_end(); ++it) {
    //       it->measure(); // misura della cella
    // }

    

    // iterators
    using cell_iterator = std::list<cell_t>::iterator;
    using const_cell_iterator = std::list<cell_t>::const_iterator;

    cell_iterator cells_begin() { return cells_.begin(); }
    cell_iterator cells_end() { return cells_.end(); }

    const_cell_iterator cells_cbegin() { return cells_.cbegin(); }
    const_cell_iterator cells_cend() { return cells_.cend(); }

    // DA IMPLEMENTARE SUBITO: (per plotting)

    // matrice dei nodi
    //   una matrice di coordinate (n_nodi x EmbedDim)

    // observers
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes() const {
        // Simply calls nodes() from dcel_t
        return dcel_.nodes();
    }

    int n_nodes() { return dcel_.n_nodes(); }

    // matrici degli edge
    //   una matrice di interi (n_edge x 2) dove riga i-esima: 
    //          [id_nodo_1 i-esimo edge, id_nodo_2 i-esimo edge]

    Eigen::Matrix<int, Eigen::Dynamic, 2> edges() const {
        // calls edges() from dcel_t
        return dcel_.edges();
    }

    int n_edges() { return dcel_.n_edges(); }

    int n_cells() { return cells_.size(); }

    void export_to_json(const std::string& filename){
        this->dcel_.export_to_json(filename);
    }



    private:
    dcel_t dcel_;

    // Internal function for computing v_vertex neighbours
    std::vector<int> compute_neighbours_(const triangle_t& d_simplex, 
        std::map<int, typename simplex_t::NodeType>& midpoints_lookup,
        int& cur_midpoint_id,
        bool print){
                                        
        // Get neighbours
        Eigen::Matrix<int, Eigen::Dynamic, 1> neigh = d_simplex.neighbors();

        // Store them inside a vector
        std::vector<int> neigh_vec(neigh.data(), neigh.data() + neigh.size());

        if(print){
            for (int elim : neigh_vec){
                std::cout << elim << std::endl;
            }
        }

        // Midpoint IDs
        std::vector<int> mp_ids;

        // Loop over d_simplex edges
        for (auto d_simplex_edge = d_simplex.edges_begin(); d_simplex_edge != d_simplex.edges_end(); d_simplex_edge++){

                // If edge on boundary ...
                if (d_simplex_edge -> on_boundary()){

                    // ... compute midpoint ...
                    auto mp = d_simplex_edge -> compute_midpoint();

                    // ... and create a new midpoint index
                    midpoints_lookup[cur_midpoint_id] = mp;

                    // Add in the added midpoint ids
                    mp_ids.push_back(cur_midpoint_id);

                    if(print){std::cout << "Found one on boundary to which we assigned ID: " << cur_midpoint_id << std::endl;}

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

        if (print){
            std::cout << "After changes" << std::endl;
        for (int id_elimm : neigh_vec){
            std::cout << id_elimm << std::endl;
        }
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
            std::vector<int> neigh_vec = compute_neighbours_(*(d_simplex_it), midpoints_lookup, cur_midpoint_id, false);

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

                if (v_vertex_id==2){
                    std::cout << "\nIf the first time:" << std::endl;
                    for(auto neigh_elim : neigh_vec){
                        std::cout << neigh_elim << std::endl;
                    }
                }
                first_d_centroid[cur_id] = d_centroid_id;
                v_vertex_lookup[cur_id] = d_centroid;

                // Increment current ID
                cur_id += 1;
            }

            // if (v_vertex_id==2){
            //     std::cout << "\nIt's 2 " << std::endl;
            //     std::vector<int>elim = compute_neighbours_(*(d_simplex_it), midpoints_lookup, cur_midpoint_id, true);
            //     std::cout << "\nThis is neigh: " << std::endl;
            //     for(auto neigh_elim : v_vertex_neighbours_lookup[2]){
            //         std::cout << neigh_elim << std::endl;
            //     }
            // }

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
        std::vector<int> ordered_neigh_ids(v_vertex_neighbours_lookup.at(v_vertex_id).size());

        // Loop over neighbours
        for (auto d_centroid_neigh_id : v_vertex_neighbours_lookup.at(v_vertex_id)) {

            int neigh_id = retrieve_v_vertex_id_(d_centroid_neigh_id, d_centroid2v_vertex);
            
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

            return ang_a < ang_b;
        });

        return ordered_neigh_ids;
    }
   

    // List of cells
    std::list<typename dcel_t::cell_t> cells_;
};
}


#endif