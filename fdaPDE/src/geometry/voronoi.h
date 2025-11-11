#ifndef __FDAPDE_VORONOI_H__
#define __FDAPDE_VORONOI_H__


// PROMEMORIA: CI MANCA DA AGGIUNGERE GLI HALFEDGE ALLE CELLE --> VA CAPITO COME PERCHE' AD UNA VA ASSOCIATO L'HALFEDGE E AD UNA IL TWIN MA NON SO COME


namespace fdapde{

template <int LocalDim, int EmbedDim>
class Voronoi {
    public:

    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;

    using myDCEL = DCEL<local_dim, embed_dim>;
    using mySimplex = Simplex<local_dim, embed_dim>;
    using myTriangulation = Triangulation<local_dim, embed_dim>;
    // using myDelaunay = Delaunay<local_dim, embed_dim>;

    // Voronoi(Matrix seed)
    // Ho scelto questi ingressi basandomi sul constructor in delaunay.h che usa random generated points e no refinement
    // fare il delaunay dei seed O(n log(n)) --> I think that we can rely on the constructor in delaunay.h with random generated points, without no refinement
    // chiamo quello sotto (faccio il duale)
    // Voronoi(const std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>& boundaries, int N=0, const std::vector<std::vector<Eigen::Matrix<double, Eigen::Dynamic, embed_dim>>>& holes = {{}}) 
    // : Voronoi(typename myDelaunay(boundaries, N, holes))
    // {};

    Voronoi(const Triangulation<local_dim, embed_dim>& mesh) {

    // controlla in triangulation se esiste la funzione che ti da tutti i vertici 
    // della mesh (che diventano centroidi in myDCEL)
    // vedi riga 95 di trinagulation.h, in teoria semplicemente chiamando mesh->nodes() 
    // dovresti avere una matrice di coordinate di nodi

    // calcola il duale (controllare gli articoli se quello sotto è il modo migliore, se esistono 
    // algoritmi più efficienti e noti in lettaratura, implementate quelli!!)


    //   - calcolare i centroidi
    // Per farlo, all'interno del loop deve:
    //      - Trovare il centroide
    //      - Inserirlo nel myDCEL
    // Problema: il tipo di nodo che ti restituisce cell.circumenter() è Eigen::Matrix<double, embed_dim, 1>,
    // mentre quello che vuole il myDCEL è un oggetto del suo tipo node_t, il quale richiede di conoscere anche l'halfedge:
    // node_t(int id, halfedge_t* halfedge, bool boundary, const CoordsType& coords)
    // IDEA:
    // - Loopa su tutte le celle e trova il centroide, e salva questa informazione, associando cell_id: centroide
    // - Secondo loop sulle celle, per ciascuna cella: retrieva il centroide, loopa sulle neighbours, collega i due centroidi
    // (considera che il minimal constructor del halfedge_t è halfedge_t(int id, node_t* node, bool sub= false))
    // NB secondo me conviene che l'ID del centroide coincida con quello della cella di Delaunay, e l'ID del
    // site del Voronoi con quello del nodo di Delaunay che gli coincide.

    // Lookup
    std::map<int, typename mySimplex::NodeType> centroid_lookup;
    // Visited or not
    std::map<int, bool> visited_centroids;
    // boundary cell or not
    std::map<int, bool> on_boundary;
    // Old node ID : new node ID
    std::map<int, int> old2new;

    // Loop to compute centroids
    int cur_id = 0;
    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it) {

        // Retrieve ID 
        int cell_id = it->id();

        // New ID
        old2new[cell_id] = cur_id;
        
        // for this, take a look at simplex.h (mySimplex::NodeType is a Eigen::Matrix<double, embed_dim, 1>)
        typename mySimplex::NodeType centroid = it->circumcenter();

        // Insert inside the DCEL the node, without the halfedge
        typename myDCEL::node_t cell_node(cur_id, false, centroid);
        dcel_.insert_node(cell_node);

        // Add to lookup
        centroid_lookup[cur_id] = centroid;
        visited_centroids[cur_id] = false; 
        on_boundary[cur_id] = false;
        cur_id += 1;
    }

// ==========================================================================================================================================

    // Imagine to have the correspondence cell_id: centroid coordinates
    for(auto it = mesh.cells_begin(); it != mesh.cells_end(); ++it) {
        std::cout << "New mesh cell\n" << std::endl;

        // Retrieve the IDs of the vertices of the current cell
        Eigen::Matrix<int, Dynamic, 1> ids_node = it->node_ids();

        // Retrieve old ID 
        int old_cell_id = it->id();

        // Retrieve new ID
        int cell_id = old2new.at(old_cell_id);

        // std::cout << "\nCell with ID: " << cell_id << std::endl;
        
        // check if the cell is on the boundary
        if (it->on_boundary()){
            on_boundary[cell_id] = true;
        }

        // Retrieve circumcenter
        typename mySimplex::NodeType centroid = centroid_lookup.at(cell_id);

        // Get adjajent cells
        Eigen::Matrix<int, Eigen::Dynamic, 1> cell_neighbours = it->neighbors();
        
        int counter = 0;

        // Keep the previous and the first halfedge
        // typename myDCEL::halfedge_t* passed_halfedge=nullptr;
        typename myDCEL::halfedge_t* first_halfedge=nullptr;

        // Pointers to e1 and e2 initialized as null_pointer
        typename myDCEL::halfedge_t* e1 = nullptr;
        typename myDCEL::halfedge_t* e2 = nullptr;

        // Loop over neighbouring cells (ASSUMING THEY ARE LOOPED IN CLOCKWISE ORDER)
        for (auto old_neigh_cell_id : cell_neighbours) {
            

        //     // Retrieve the triangle in the mesh using the old id
        //     Triangle tria(old_neigh_cell_id, mesh);  // questo secondo me è sbagliato ma non so come si scrive
        //     Eigen::Matrix<int, Dynamic, 1> ids_neigh_node = tria.node_ids();

        //     // Find the different id between this vertices and the cell ones
        //     std::vector<int> vec_1(ids_node.data(), ids_node.data()+ids_node.size());
        //     std::vector<int> vec_2(ids_neigh_node.data(), ids_neigh_node.data()+ids_neigh_node.size());

        //     // Sort the vectors
        //     std::sort(vec_1.begin(), vec_1.end());
        //     std::sort(vec_2.begin(), vec_2.end());

        //     std::vector<int> common;

        //     // Take the common elements between the two
        //     std::set_intersection(vec_1.begin(), vec_1.end(),
        //                         vec_2.begin(), vec_2.end(),
        //                         std::back_inserter(common));

        //     // Assert that the common elements are two since they are adjacent cells
        //     assert(common.size()==2);

        //     // Add the new cell to the list if it is not yet added 
        //     for(int ids : common){
        //         cell_t curr_cell(ids);
        //         if (std::find(cells_.begin(), cells_.end(), curr_cell) == cells_.end()){
        //             cells_.push_back(curr_cell);
        //         }
        //     }

            if (old_neigh_cell_id==-1){
                // std::cout << "\nAbout to break..." << std::endl;
                break;}

            if (old_neigh_cell_id != -1){
                // Retrieve new ID for neighbour
                int neigh_cell_id = old2new.at(old_neigh_cell_id);
                std::cout << "My ID" << cell_id << " neighbour ID " << neigh_cell_id << std::endl;
            
                // Retrieve centroid of the neighbouring cell
                // std::cout << "\nRetrieving centroid...";
                typename mySimplex::NodeType neigh_centroid = centroid_lookup.at(neigh_cell_id); 
                
                // Create cell node
                typename myDCEL::node_t* cell_node;

                // Create twin node 
                typename myDCEL::node_t* twin_node;

                // Create halfedge centroid -> neigh
                typename myDCEL::halfedge_t* cell_halfedge = new myDCEL::halfedge_t();

                // Create twin halfedge (neigh -> centroid)
                typename myDCEL::halfedge_t* twin_halfedge = new myDCEL::halfedge_t();

                // // Set IDs based on counter
                // cell_halfedge.set_id(counter++);
                // twin_halfedge.set_id(counter++);

                // Find nodes
                cell_node = dcel_.find_node(dcel_.nodes().row(cell_id));
                twin_node = dcel_.find_node(dcel_.nodes().row(neigh_cell_id));


                // If neighbouring cell was not visited, its node needs an halfedge
                if (! visited_centroids.at(neigh_cell_id)){
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

                // Else, only the cell centroid needs an halfedge
                else{
                    // Retrieve halfedge
                    std::cout << "I am asking for ID " << cell_id << std::endl;
                    delete cell_halfedge;
                    delete twin_halfedge;
                    twin_halfedge = (twin_node->neighID2halfedge(cell_id));
                    cell_halfedge = (twin_halfedge->twin());

                }
              
                // set the cell as visited
                visited_centroids[cell_id] = true;
                // visited_centroids[neigh_cell_id] = true;

                // Provo a fare qui
                e2 = cell_halfedge;  
                if (e1 != nullptr){
                    e1->twin()->set_next(e2);
                    e2->set_prev(e1->twin());
                }
                else{
                    first_halfedge = cell_halfedge;
                }
                e1 = e2;
        }
    }

    first_halfedge->twin()->set_next(e2);
    e2->set_prev(first_halfedge->twin());

    }




// ==========================================================================================================================================



    //   - "collegare" i centroidi
    //    -- costruire myDCEL in modo tale che codifichi il voronoi
    //    -- identificare le celle unbounded

    int n_cells = dcel_.n_cells();
    std::list<cell_t> cells_(n_cells);

    // std::cout << "\nPopulating list of cells...";
    // for (auto it = dcel_.cells_begin(); it != dcel_.cells_end(); ++it){
    //     std::cout << "\nInside loop";
    //     cell_t whatever;
    //     whatever.set_cell(&(*it));
    //     cells_.push_back(whatever);
        
    //     // it->set_unbounded(on_boundary.at(it->id()));
    // }

    std::cout << "\nFinished constructor" << std::endl;

    }

    // struttura per codificare una cella di voronoi
    struct cell_t : public myDCEL::cell_t{

    // double measure() const { return internals::signed_measure_2d_polygon(); };
    // accedere ai vertici della cella
    std::list<typename myDCEL::node_t*> cell_nodes() const {
        std::list<typename myDCEL::node_t*> cell_nodes;

        std::list<typename myDCEL::halfedge_t*> edges = cell_edges();
        for(const auto it = edges.cbegin(); it != edges.cend(); ++it) {
            cell_nodes.push_back((*it)->node());
        }

        return cell_nodes;
    }   
    // accedere agli edge della cella
    std::list<typename myDCEL::halfedge_t*> cell_edges() const {
        std::list<typename myDCEL::halfedge_t*> cell_edges;

        typename myDCEL::halfedge_t* start = cell_->halfedge();

        for (typename myDCEL::halfedge_t::circulator it(start); it; ++it) {
            cell_edges.push_back(&(*it));
        }

        return cell_edges;

        // myDCEL::halfedge_t* he = start;

        // do {
        //     cell_edges.push_back(he);
        //     he = he->next();
        // } while(he != start)
    }

    void set_unbounded() {unbounded_ = true; }

    bool is_unbounded() const { return unbounded_; }

    void set_cell(typename myDCEL::cell_t* cell_ptr){cell_ = cell_ptr;}

    typename myDCEL::cell_t* cell_;
    bool unbounded_; // cella unbounded o no?
    };

    // iteratori sulle celle
    // in futuro vorremmo poter fare
    // for(auto it = voronoi.cells_begin(); it != voronoi.cells_end(); ++it) {
    //       it->measure(); // misura della cella
    // }

    // iterators
    using cell_iterator = std::list<cell_t>::iterator;
    using const_cell_iterator = std::list<cell_t>::const_iterator;

    cell_iterator cells_begin(void) { return cells_.begin(); }
    cell_iterator cells_end(void) { return cells_.end(); }

    const_cell_iterator cells_cbegin(void) { return cells_.cbegin(); }
    const_cell_iterator cells_cend(void) { return cells_.cend(); }

    // DA IMPLEMENTARE SUBITO: (per plotting)

    // matrice dei nodi
    //   una matrice di coordinate (n_nodi x EmbedDim)

    // observers
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodes() const {
        // Simply calls nodes() from myDCEL
        return dcel_.nodes();
    }

    int n_nodes(void) { return dcel_.n_nodes(); }

    // matrici degli edge
    //   una matrice di interi (n_edge x 2) dove riga i-esima: 
    //          [id_nodo_1 i-esimo edge, id_nodo_2 i-esimo edge]

    Eigen::Matrix<int, Eigen::Dynamic, 2> edges() const {
        // calls edges() from myDCEL
        return dcel_.edges();
    }

    int n_edges(void) { return dcel_.n_edges(); }

    int n_cells(void) { return cells_.size(); }

    private:
    myDCEL dcel_;

    // List of cells
    std::list<cell_t> cells_;
};
}


#endif