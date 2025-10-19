#ifndef __FDAPDE_VORONOI_H__
#define __FDAPDE_VORONOI_H__




namespace fdapde{

template <int LocalDim, int EmbedDim>
class Voronoi {

    static constexpr int local_dim = LocalDim;
    static constexpr int embed_dim = EmbedDim;

    using myDCEL = DCEL<local_dim, embed_dim>;
    using mySimplex = Simplex<local_dim, embed_dim>;

    // Voronoi(Matrix seed) {

    // // fare il delaunay dei seed O(n log(n))
    // // chiamo quello sotto (faccio il duale)

    // }

    Voronoi(const Triangulation<local_dim, embed_dim>& mesh) {

    // controlla in triangulation se esiste la funzione che ti da tutti i vertici 
    // della mesh (che diventano centroidi in myDCEL)
    // vedi riga 95 di trinagulation.h, in teoria semplicemente chiamando mesh->nodes() 
    // dovresti avere una matrice di coordinate di nodi

    // calcola il duale (controllare gli articoli se quello sotto è il modo migliore, se esistono 
    // algoritmi più efficienti e noti in lettaratura, implementate quelli!!)

    // Inizializza il myDCEL
    myDCEL voronoi_dcel;

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

    // Loop to compute centroids
    for(const auto& cell : mesh) {

        // Retrieve ID 
        int cell_id = mesh.id();
        
        // for this, take a look at simplex.h (mySimplex::NodeType is a Eigen::Matrix<double, embed_dim, 1>)
        typename mySimplex::NodeType centroid = cell.circumcenter();

        // Add to lookup
        centroid_lookup[cell_id] = centroid;
        visited_centroids[cell_id] = false; 
        on_boundary[cell_id] = false;
    }

    // Imagine to have the correspondence cell_id: centroid coordinates
    for(const auto& cell : mesh) {

        // Retrieve ID 
        int cell_id = mesh.id();
        
        // check if the cell is on the boundary
        if (cell.on_boundary()){
            on_boundary[cell_id] = true;
        }

        // Retrieve circumcenter
        typename mySimplex::NodeType centroid = centroid_lookup.at(cell_id);

        // Get adjajent cells
        Eigen::Matrix<int, Dynamic, 1> cell_neighbours = cell.neighbors();
        
        int counter = 0;

        // Loop over neighbouring cells
        for (int neigh_cell_id : cell_neighbours) {

            // Retrieve centroid of the neighbouring cell
            typename mySimplex::NodeType neigh_centroid = centroid_lookup.at(neigh_cell_id); 

            // Create halfededge
            typename myDCEL::halfedge_t cell_halfedge;
            cell_halfedge.set_id(counter++);
            // Create node or retrieve it, if it is already in the myDCEL (the part of checking is not implemented)
            // Check if already 
            typename myDCEL::node_t cell_node;
            if (! visited_centroids.at(cell_id)){
                cell_node = myDCEL::node_t(cell_id, cell_halfedge, false, centroid);
            }
            else{
                cell_node = *voronoi_dcel.find_node(voronoi_dcel.nodes()[cell_id]);
            }

            // Create twin halfedge
            typename myDCEL::halfedge_t twin_halfedge;
            twin_halfedge.set_id(counter++);
            // Create twin node or retrieve it, if it is already in the myDCEL (the part of checking is not implemented)
            typename myDCEL::node_t twin_node;
            if (! visited_centroids.at(neigh_cell_id)){
                twin_node = typename myDCEL::node_t(neigh_cell_id, twin_halfedge, false, neigh_centroid);
            }
            else{
                twin_node = *voronoi_dcel.find_node(voronoi_dcel.nodes().row(neigh_cell_id));
            }

            // Add node and twin to the halfedges
            cell_halfedge.set_node(&cell_node);
            twin_halfedge.set_node(&twin_node);
            cell_halfedge.set_twin(&twin_halfedge);
            twin_halfedge.set_twin(&cell_halfedge);

            // Add the nodes IF THEY DO NOT ALREADY EXIST
            voronoi_dcel.insert_node(cell_node);
            voronoi_dcel.insert_node(twin_node);

            // Add the halfedges
            voronoi_dcel.insert_edge(&cell_halfedge, &twin_halfedge);

        }

        // set the cell as visited
        visited_centroids[cell_id] = true;

    }

    //   - "collegare" i centroidi
    //    -- costruire myDCEL in modo tale che codifichi il voronoi
    //    -- identificare le celle unbounded

    int n_cells = voronoi_dcel.n_cells();
    std::list<cell_t> cells_(n_cells);

    for (auto it = voronoi_dcel.cells_cbegin(); it != voronoi_dcel.cells_cend(); ++it){
        cell_t whatever;
        whatever.set_cell(&(*it));
        cells_.push_back(whatever);
        
        // it->set_unbounded(on_boundary.at(it->id()));
    }

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

    void set_cell(const typename myDCEL::cell_t* cell_ptr){cell_ = cell_ptr;}

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

    private:
    myDCEL dcel_;

    // List of cells
    std::list<cell_t> cells_;
};
}


#endif