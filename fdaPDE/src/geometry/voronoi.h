template <int LocalDim, int EmbedDim>
class Voronoi {

    Voronoi(Matrix seed) {

    // fare il delaunay dei seed O(n log(n))
    // chiamo quello sotto (faccio il duale)

    }

    Voronoi(const Triangulation& mesh) {

    // controlla in triangulation se esiste la funzione che ti da tutti i vertici 
    // della mesh (che diventano centroidi in DCEL)
    // vedi riga 95 di trinagulation.h, in teoria semplicemente chiamando mesh->nodes() 
    // dovresti avere una matrice di coordinate di nodi

    // calcola il duale (controllare gli articoli se quello sotto è il modo migliore, se esistono 
    // algoritmi più efficienti e noti in lettaratura, implementate quelli!!)

    // Inizializza il DCEL
    DCEL voronoi_dcel();

    //   - calcolare i centroidi
    // Per farlo, all'interno del loop deve:
    //      - Trovare il centroide
    //      - Inserirlo nel DCEL
    // Problema: il tipo di nodo che ti restituisce cell.circumenter() è Eigen::Matrix<double, embed_dim, 1>,
    // mentre quello che vuole il DCEL è un oggetto del suo tipo node_t, il quale richiede di conoscere anche l'halfedge:
    // node_t(int id, halfedge_t* halfedge, bool boundary, const CoordsType& coords)
    // IDEA:
    // - Loopa su tutte le celle e trova il centroide, e salva questa informazione, associando cell_id: centroide
    // - Secondo loop sulle celle, per ciascuna cella: retrieva il centroide, loopa sulle neighbours, collega i due centroidi
    // (considera che il minimal constructor del halfedge_t è halfedge_t(int id, node_t* node, bool sub= false))
    // NB secondo me conviene che l'ID del centroide coincida con quello della cella di Delaunay, e l'ID del
    // site del Voronoi con quello del nodo di Delaunay che gli coincide.

    // Lookup
    std::map<int, Simplex::NodeType> centroid_lookup();

    // Loop to compute centroids
    for(const auto& cell : mesh) {

        // Retrieve ID 
        int cell_id = mesh.id();
        
        // for this, take a look at simplex.h (Simplex::NodeType is a Eigen::Matrix<double, embed_dim, 1>)
        Simplex::NodeType centroid = cell.circumcenter();

        // Add to lookup
        centroid_lookup[cell_id] = centroid;
    }

    // Imagine to have the correspondence cell_id: centroid coordinates
    for(const auto& cell : mesh) {

        // Retrieve ID 
        int cell_id = mesh.id();
        
        // Retrieve circumcenter
        Simplex::NodeType centroid = centroid_lookup.at(cell_id);

        // Get adjajent cells
        Eigen::Matrix<int, Dynamic, 1> cell_neighbours = cell.neighbors();

        // Loop over neighbouring cells
        for (int neigh_cell_id : cell_neighbours) {

            // Retrieve centroid of the neighbouring cell
            Simplex::NodeType neigh_centroid = ; // Non ho trovato il metodo per farlo in trianulation.h

            // Create halfededge
            DCEL::halfedge_t cell_halfedge();
            // Create node or retrieve it, if it is already in the DCEL (the part of checking is not implemented)
            DCEL::node_t cell_node(cell_id, cell_halfedge, false, centroid);

            // Create twin halfedge
            DCEL::helfedge_t twin_halfedge()
            // Create twin node or retrieve it, if it is already in the DCEL (the part of checking is not implemented)
            DCEL::node_t twin_node(neigh_cell_id, twin_halfedge, false, neigh_centroid);

            // Add node and twin to the halfedges
            cell_halfedge.set_node(cell_node);
            twin_halfedge.set_node(twin_node);
            cell_halfedge.set_twin(twin_halfedge);
            twin_halfedge.set_twin(cell_halfedge);

            // Add the nodes IF THEY DO NOT ALREADY EXIST
            voronoi_dcel.insert_node(cell_node)
            voronoi_dcel.insert_node(twin_node)

            // Add the halfedges
            voronoi_dcel.insert_edge(cell_halfedge, twin_halfedge)

        }

    }

    //   - "collegare" i centroidi
    //    -- costruire DCEL in modo tale che codifichi il voronoi
    //    -- identificare le celle unbounded

    }

    // struttura per codificare una cella di voronoi
    struct cell_t : public DCEL::cell_t{

    double measure() const { internals::measure_2d_polygon(); };
    // accedere ai vertici della cella
    std::list<DCEL::nodes_t*> cell_nodes() const {
        std::list<DCEL::nodes_t*> cell_nodes;

        std::list<DCEL::halfedge_t*> edges = cell_edges();
        for(const auto it = edges.cbegin(); it != edges.cend(); ++it) {
            cell_nodes.push_back((*it)->node());
        }

        return cell_nodes;
    }   
    // accedere agli edge della cella
    std::list<DCEL::halfedge_t*> cell_edges() const {
        std::list<DCEL::halfedge_t*> cell_edges;

        DCEL::halfedge_t* start = cell_->halfedge();
        if (!start) return;  // safety check

        for (DCEL::halfedge_t::circulator it(start); it; ++it) {
            cell_edges.push_back(&(*it));
        }

        return cell_edges;

        // DCEL::halfedge_t* he = start;

        // do {
        //     cell_edges.push_back(he);
        //     he = he->next();
        // } while(he != start)
    }

    bool is_unbounded() const { return unbounded_; }

    DCEL::cell_t* cell_;
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
        // Simply calls nodes() from DCEL
        return dcel_.nodes();
    }

    int n_nodes(void) { return dcel_.n_nodes(); }

    // matrici degli edge
    //   una matrice di interi (n_edge x 2) dove riga i-esima: 
    //          [id_nodo_1 i-esimo edge, id_nodo_2 i-esimo edge]

    Eigen::Matrix<int, Eigen::Dynamic, 2> edges() const {
        // calls edges() from DCEL
        return dcel_.edges();
    }

    int n_edges(void) { return dcel_.n_edges(); }

    private:
    DCEL<LocalDim, EmbedDim> dcel_;

    // List of cells
    std::list<cell_t> cells_;
};
