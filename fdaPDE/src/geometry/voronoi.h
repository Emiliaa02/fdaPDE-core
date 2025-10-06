template <int LocalDim, int EmbedDim>
class Voronoi {

    Voronoi(Matrix seed) {

    // fare il delaunay dei seed O(n log(n))
    // chiamo quello sotto (faccio il duale)

    }

    Voronoi(const Triangulation& mesh) {

    // calcola il duale (controllare gli articoli se quello sotto è il modo migliore, se esistono 
    // algoritmi più efficienti e noti in lettaratura, implementate quelli!!)
    //   - calcolare i centroidi
    for(const auto& cell : mesh) {
    cell.centroid();
    }

    //   - "collegare" i centroidi
    //    -- costruire DCEL in modo tale che codifichi il voronoi
    //    -- identificare le celle unbounded

    }

    // struttura per codificare una cella di voronoi
    struct cell_t {

    double measure() const { internals::measure_2d_polygon(); };
    // accedere ai vertici della cella
    // accedere agli edge della cella
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

    private:
    DCEL<LocalDim, EmbedDim> dcel_;

    // List of cells
    std::list<cell_t> cells_;
};
