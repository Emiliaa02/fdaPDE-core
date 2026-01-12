#ifndef __FDAPDE_VORONOI_H__
#define __FDAPDE_VORONOI_H__


// PROMEMORIA: CI MANCA DA AGGIUNGERE GLI HALFEDGE ALLE CELLE --> VA CAPITO COME PERCHE' AD UNA VA ASSOCIATO L'HALFEDGE E AD UNA IL TWIN MA NON SO COME


namespace fdapde{

template <int LocalDim, int EmbedDim>
class Voronoi {

    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;

    using dcel_t = DCEL<local_dim, embed_dim>;
    using simplex_t = Simplex<local_dim, embed_dim>;
    using triangulation_t = Triangulation<local_dim, embed_dim>;
    using triangle_t = Triangle<triangulation_t>;
    // using myDelaunay = Delaunay<local_dim, embed_dim>;
    
    public:

    // Voronoi(Matrix seed)
    // Ho scelto questi ingressi basandomi sul constructor in delaunay.h che usa random generated points e no refinement
    // fare il delaunay dei seed O(n log(n)) --> I think that we can rely on the constructor in delaunay.h with random generated points, without no refinement
    // chiamo quello sotto (faccio il duale)
    // Voronoi(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& boundaries, int N=0, const std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>>& holes = {{}}) 
    // : Voronoi(typename myDelaunay(boundaries, N, holes))
    // {};

    Voronoi(const Triangulation<local_dim, embed_dim>& mesh) {

    // controlla in triangulation se esiste la funzione che ti da tutti i vertici 
    // della mesh (che diventano centroidi in dcel_t)
    // vedi riga 95 di trinagulation.h, in teoria semplicemente chiamando mesh->nodes() 
    // dovresti avere una matrice di coordinate di nodi

    // calcola il duale (controllare gli articoli se quello sotto è il modo migliore, se esistono 
    // algoritmi più efficienti e noti in lettaratura, implementate quelli!!)


    //   - calcolare i centroidi
    // Per farlo, all'interno del loop deve:
    //      - Trovare il centroide
    //      - Inserirlo nel dcel_t
    // Problema: il tipo di nodo che ti restituisce cell.circumenter() è Eigen::Matrix<double, embed_dim, 1>,
    // mentre quello che vuole il dcel_t è un oggetto del suo tipo node_t, il quale richiede di conoscere anche l'halfedge:
    // node_t(int id, halfedge_t* halfedge, bool boundary, const CoordsType& coords)
    // IDEA:
    // - Loopa su tutte le celle e trova il centroide, e salva questa informazione, associando cell_id: centroide
    // - Secondo loop sulle celle, per ciascuna cella: retrieva il centroide, loopa sulle neighbours, collega i due centroidi
    // (considera che il minimal constructor del halfedge_t è halfedge_t(int id, node_t* node, bool sub= false))
    // NB secondo me conviene che l'ID del centroide coincida con quello della cella di Delaunay, e l'ID del
    // site del Voronoi con quello del nodo di Delaunay che gli coincide.

    // Number of faces in triangulations
    int n_mesh_faces = mesh.n_cells();
    // Number of vertices in the mesh
    int n_mesh_vertices = mesh.n_nodes();

    // Lookup
    std::vector<typename simplex_t::NodeType> centroid_lookup(n_mesh_faces+1);
    // Visited or not
    std::vector<bool> visited_centroids(n_mesh_faces+1);
    // boundary cell or not
    std::vector<bool> on_boundary(n_mesh_faces+1);
    // Old node ID : new node ID
    std::vector<int> old2new(n_mesh_faces);
    // cell with halfedge
    std::vector<int> cell2half(n_mesh_vertices);  // 0 for cell node e 1 for twin node
    // counter to set the IDs
    int counter = 0;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Loop to compute centroids
    int cur_id = 0;
    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it) {

        // Retrieve ID 
        int cell_id = it->id();

        // for this, take a look at simplex.h (simplex_t::NodeType is a Eigen::Matrix<double, embed_dim, 1>)
        typename simplex_t::NodeType centroid = it->circumcenter();

        // Check that the centroid is not already present
        bool already_present = false;
        int same_centroid_idx;
        for (auto ii = 0; ii < centroid_lookup.size(); ++ii){
            float thresh = 0.000000001;
            if (std::abs(centroid_lookup[ii](0)-centroid(0))<thresh and std::abs(centroid_lookup[ii](1)-centroid(1))<thresh){
                already_present = true;
                same_centroid_idx = ii;
                break;
            }
        }

        if (!already_present){

            // New ID
            old2new[cell_id] = cur_id;

            // Insert inside the DCEL the node, without the halfedge
            typename dcel_t::node_t cell_node(cur_id, false, centroid);
            dcel_.insert_node(cell_node);

            // Add to lookup
            centroid_lookup[cur_id] = centroid;
            visited_centroids[cur_id] = false; 
            on_boundary[cur_id] = false;

            cur_id += 1;
        }
        else{

            // New ID
            old2new[cell_id] = same_centroid_idx;
        }

    }

    int infty_id = cur_id;

    typename simplex_t::NodeType centroid_infty = simplex_t::NodeType::Constant(std::numeric_limits<double>::infinity());

    centroid_lookup[infty_id] = centroid_infty;
    visited_centroids[infty_id] = false;
    on_boundary[infty_id] = false;

    typename dcel_t::node_t infty_node(infty_id, false, centroid_infty);
    dcel_.insert_node(infty_node);

    // Shrink the vectors capacity if less elements than estimated were used
    centroid_lookup.shrink_to_fit();
    visited_centroids.shrink_to_fit();
    on_boundary.shrink_to_fit();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ==========================================================================================================================================

    // Imagine to have the correspondence cell_id: centroid coordinates
    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it) {

        // Retrieve the IDs of the vertices of the current cell
        Eigen::Matrix<int, Dynamic, 1> ids_node = it->node_ids();

        // Retrieve old ID 
        int old_cell_id = it->id();

        // Retrieve new ID
        int cell_id = old2new.at(old_cell_id);
        
        // check if the cell is on the boundary
        if (it->on_boundary()){
            on_boundary[cell_id] = true;
        }

        // Retrieve circumcenter
        typename simplex_t::NodeType centroid = centroid_lookup.at(cell_id);

        // Get adjajent cells
        Eigen::Matrix<int, Eigen::Dynamic, 1> cell_neighbours = it->neighbors();

        // Keep the previous and the first halfedge
        // typename dcel_t::halfedge_t* passed_halfedge=nullptr;
        typename dcel_t::halfedge_t* first_halfedge=nullptr;

        // Pointers to e1 and e2 initialized as null_pointer
        typename dcel_t::halfedge_t* e1 = nullptr;
        typename dcel_t::halfedge_t* e2 = nullptr;

        // Make the orientation consistent
        // double determinant = 1;
        // std::vector<typename simplex_t::NodeType> nodes_cell;    
        // bool compute_det = true;
        // for (auto old_neigh_cell_id : cell_neighbours){
            
        //     if(old_neigh_cell_id == -1){
        //         compute_det = false;
        //         break;
        //     }

        //     int new_neigh_id = old2new.at(old_neigh_cell_id);

        //     typename simplex_t::NodeType neigh_centroid_coords = centroid_lookup.at(new_neigh_id);
        //     typename simplex_t::NodeType centered_neigh_coords;
        //     for(int i = 0; i < embed_dim; i++){
        //         centered_neigh_coords(i) = neigh_centroid_coords(i) - centroid(i);
        //     }
        //     nodes_cell.push_back(centered_neigh_coords);
        // }
        // if(compute_det){
        //     determinant = (nodes_cell[1](0) - nodes_cell[0](0))*(nodes_cell[2](1) - nodes_cell[0](1)) - (nodes_cell[1](1) - nodes_cell[0](1))*(nodes_cell[2](0) - nodes_cell[0](0));
        //     if(determinant < 0){
        //         int tmp = cell_neighbours(0);
        //         cell_neighbours(0) = cell_neighbours(1);
        //         cell_neighbours(1) = tmp;
        //     }
        // }

        bool prev_was_minus_1 = false;
        int count_existing_neigh = 0;
        int diff_id = 0;
        
        // Loop over neighbouring cells (ASSUMING THEY ARE LOOPED IN CLOCKWISE ORDER)
        for (auto old_neigh_cell_id : cell_neighbours) {

            int neigh_cell_id;

            // Create cell node
            typename dcel_t::node_t* cell_node;

            // Find nodes
            cell_node = dcel_.find_node(dcel_.nodes().row(cell_id));

            // Create twin node 
            typename dcel_t::node_t* twin_node;

            // Create halfedge centroid -> neigh
            typename dcel_t::halfedge_t* cell_halfedge = new dcel_t::halfedge_t();

            // Create twin halfedge (neigh -> centroid)
            typename dcel_t::halfedge_t* twin_halfedge = new dcel_t::halfedge_t();

            if(old_neigh_cell_id != -1){
                // Retrieve new ID for neighbour
                neigh_cell_id = old2new.at(old_neigh_cell_id);
                if (neigh_cell_id == cell_id){
                    continue;
                }

                // Retrieve centroid of the neighbouring cell
                typename simplex_t::NodeType neigh_centroid = centroid_lookup.at(neigh_cell_id);
                
                twin_node = dcel_.find_node(dcel_.nodes().row(neigh_cell_id));
            }
            else{

                // prev_was_minus_1 = true;

                // It is the infinity node
                neigh_cell_id = infty_id;

                twin_node = &infty_node;
            }


            // If neighbouring cell was not visited, its node needs an halfedge
            if (!visited_centroids.at(neigh_cell_id)){
                // Set IDs based on counter
                cell_halfedge->set_id(counter++);
                twin_halfedge->set_id(counter++);

                // Set halfedge to the cell centroid
                cell_node->set_halfedge(cell_halfedge);
                // Set halfedge to the neighbouring centroid
                twin_node->set_halfedge(twin_halfedge);

                // Add node and twin to the halfedges
                cell_halfedge->set_node(cell_node);
                twin_halfedge->set_node(twin_node);
                cell_halfedge->set_twin(twin_halfedge);
                twin_halfedge->set_twin(cell_halfedge);

                // Add the halfedges (true because, being twins, cells for the two halfedges are different)
                dcel_.insert_edge(cell_halfedge, twin_halfedge, true);

                // Update lookup
                cell_node->add_halfedge(neigh_cell_id, cell_halfedge);
                twin_node->add_halfedge(cell_id, twin_halfedge);
            }

            // Else, only retrieve them
            else{
                // Retrieve halfedge
                delete cell_halfedge;
                delete twin_halfedge;
                twin_halfedge = (twin_node->neighID2halfedge(cell_id));
                cell_halfedge = (twin_halfedge->twin());

            }
        
            // set the cell as visited
            visited_centroids[cell_id] = true;

            e2 = cell_halfedge;  
            if (count_existing_neigh != 0){
                e1->twin()->set_next(e2);
                e2->set_prev(e1->twin());
            }
            else{
                first_halfedge = cell_halfedge;
            }

            // prev_was_minus_1 = false;
            e1 = cell_halfedge;
            count_existing_neigh += 1;
                        
            if(old_neigh_cell_id != -1){

                //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                // Retrieve the triangle in the mesh using the old id
                triangle_t tria(old_neigh_cell_id, &mesh);  
                auto ids_neigh_node = tria.node_ids();

                // Find the different id between this vertices and the cell ones
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
                    auto it = std::find(vec_2.begin(), vec_2.end(), common[i]);
                    // check 
                    assert(it != vec_2.end());
                    // find the corrensponding index
                    int local_idx = std::distance(vec_2.begin(), it);
                    typename simplex_t::NodeType ver_coords = tria.node(local_idx);

                    // orientation check
                    int new_neigh_id = old2new.at(old_neigh_cell_id);
                    typename simplex_t::NodeType neigh_centroid_coords = centroid_lookup.at(new_neigh_id);
                    // create segments around the vertex
                    typename simplex_t::NodeType u = centroid - ver_coords;
                    typename simplex_t::NodeType v = neigh_centroid_coords - ver_coords;
                    // cross product
                    double cross = u(0)*v(1)-u(1)*v(0);
                    int new_ids = old2new.at(common[i]);
                    cell_t curr_cell(new_ids);
                    // Add the new cell to the list if it is not yet added 
                    if (std::find(cells_.begin(), cells_.end(), curr_cell) == cells_.end()){
                        if(cross < 0){
                            cell2half[new_ids] = 0;
                            // cross product < 0 means that neigh_centroid is clock-wise with respect to centroid, so we want neigh -> centroid
                            curr_cell.set_halfedge(twin_halfedge);
                            cells_.push_back(curr_cell);
                        }
                        else if(cross > 0){
                            cell2half[new_ids] = 1;
                            // cross product > 0 means that neigh_centroid is CCW with respect to centroid, so we want centroid -> neigh
                            curr_cell.set_halfedge(cell_halfedge);
                            cells_.push_back(curr_cell);
                        }
                        else{
                            std::cout<<"Undetermined: cross product is zero"<<std::endl;
                            std::cout << "centroid coordinates" << centroid << std::endl;
                            std::cout << "neigh_centroid_coords coordinates" << neigh_centroid_coords << std::endl;
                            std::cout << "Vertex coordinates" << ver_coords << std::endl; 
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

    
    if (count_existing_neigh != 1)
    {
    e2->twin()->set_next(first_halfedge);
    first_halfedge->set_prev(e2->twin());
    }


    if(count_existing_neigh == 1){
        int new_diff_id = old2new.at(diff_id);
        cell_t curr_cell_diff(new_diff_id);
        if (std::find(cells_.begin(), cells_.end(), curr_cell_diff) == cells_.end()){
            cells_.push_back(curr_cell_diff);
        }
    }

    }


// ==========================================================================================================================================



    //   - "collegare" i centroidi
    //    -- costruire dcel_t in modo tale che codifichi il voronoi
    //    -- identificare le celle unbounded

    int n_cells = dcel_.n_cells();
    //std::list<cell_t> cells_(n_cells);

    // for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it){
    //     cell_t whatever;
    //     whatever.set_cell(&(*it));
    //     cells_.push_back(whatever);
        
    //     // it->set_unbounded(on_boundary.at(it->id()));
    // }

    std::cout << "\nFinished constructor" << std::endl;

    }















































    // struttura per codificare una cella di voronoi
    struct cell_t : public dcel_t::cell_t{

    // accedere ai vertici della cella
    std::vector<typename dcel_t::node_t*> cell_nodes() const {
        std::vector<typename dcel_t::node_t*> cell_nodes;
        std::vector<typename dcel_t::halfedge_t*> edges = cell_edges();
        for(auto it = edges.cbegin(); it != edges.cend(); ++it) {
            cell_nodes.push_back((*it)->node());
        }

        return cell_nodes;
    }   
    // accedere agli edge della cella
    std::vector<typename dcel_t::halfedge_t*> cell_edges() const {
        std::vector<typename dcel_t::halfedge_t*> cell_edges;
        typename dcel_t::halfedge_t* start = this->halfedge();
        // for (typename dcel_t::halfedge_t::circulator it(start); it; ++it) {
        //     cell_edges.push_back(&(*it));
        // }
        typename dcel_t::halfedge_t* he = start;

        if (!start) return cell_edges;

        do {
            cell_edges.push_back(he);
            he = he->next();
        } while(he != start);
        
        return cell_edges;
    }

    void set_unbounded() {unbounded_ = true; }

    bool is_unbounded() const { return unbounded_; }

    void set_cell(typename dcel_t::cell_t* cell_ptr){cell_ = cell_ptr;}

    typename dcel_t::cell_t* cell_;
    bool unbounded_; // cella unbounded o no?

    double measure(){
        if(is_unbounded()){
            return std::numeric_limits<double>::infinity();
        }
        std::vector<typename dcel_t::node_t*> points_list = cell_nodes();
        Eigen::Matrix<double, Eigen::Dynamic, embed_dim> points_coords(points_list.size(), embed_dim);
        for(int i=0; i<points_list.size(); ++i){
            points_coords.row(i) = points_list[i] -> coords();
        }
        return internals::signed_measure_2d_polygon(points_coords);
    }
    };
















































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
    
   

    // List of cells
    std::list<cell_t> cells_;
};
}


#endif