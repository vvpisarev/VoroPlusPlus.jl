This page summarizes the low-level API of the package (direct exports from Voro++).

# Export Status of `container` Class Methods

This table shows which public methods from `container` and its base classes (`container_base`, `wall_list`, `voro_base`) are exported to Julia via CxxWrap.jl.

| Method | Signature | Exported As | Status |
|--------|-----------|-------------|--------|
| **Constructor** |
| `container` | `container(double ax_, double bx_, ...)` | `RawContainer(...)` | ✅ Exported |
| **Inherited from base classes** |
| `point_inside` | `bool point_inside(double x, double y, double z)` | `__cxxwrap_isinside(con, x, y, z)` | ✅ Exported |
| `draw_domain_gnuplot` | `void draw_domain_gnuplot(FILE *fp = stdout)` | — | ❌ Not exported |
| `draw_domain_gnuplot` | `void draw_domain_gnuplot(const char *filename)` | `draw_domain_gnuplot(con, filename)` | ✅ Exported |
| `draw_domain_pov` | `void draw_domain_pov(FILE *fp = stdout)` | — | ❌ Not exported |
| `draw_domain_pov` | `void draw_domain_pov(const char *filename)` | `draw_domain_pov(con, filename)` | ✅ Exported |
| `total_particles` | `int total_particles()` | `total_particles(con)` | ✅ Exported |
| `add_wall` | `void add_wall(wall *w)` | — | ❌ Not exported (generic) |
| `add_wall` | `void add_wall(wall &w)` | — | ❌ Not exported (generic) |
| `add_wall` | `void add_wall(wall_list &wl)` | — | ❌ Not exported |
| `point_inside_walls` | `bool point_inside_walls(double x, double y, double z)` | `__cxxwrap_point_inside_walls(con, x, y, z)` | ✅ Exported |
| `apply_walls` | `template<class c_class> bool apply_walls(...)` | — | ❌ Not exported |
| **Particle Management** |
| `clear` | `void clear()` | `__cxxwrap_clear!(con)` | ✅ Exported |
| `put` | `void put(int n, double x, double y, double z)` | `__cxxwrap_put!(con, id, x, y, z)` | ✅ Exported |
| `put` | `void put(particle_order &vo, int n, double x, double y, double z)` | `__cxxwrap_put!(con, ord, id, x, y, z)` | ✅ Exported |
| `import` | `void import(FILE *fp = stdin)` | `import!(con, fptr)` | ✅ Exported |
| `import` | `void import(particle_order &vo, FILE *fp = stdin)` | — | ❌ Not exported |
| `import` | `void import(const char* filename)` | `import!(con, path)` | ✅ Exported |
| `import` | `void import(particle_order &vo, const char *filename)` | — | ❌ Not exported |
| `compute_all_cells` | `void compute_all_cells()` | `compute_all_cells(con)` | ✅ Exported |
| `sum_cell_volumes` | `double sum_cell_volumes()` | `sum_cell_volumes(con)` | ✅ Exported |
| **Cell Computation** |
| `find_voronoi_cell` | `bool find_voronoi_cell(double x, double y, double z, ...)` | `__cxxwrap_find_cell(con, x, y, z)` | ✅ Exported (lambda) |
| `compute_cell` | `template<class v_cell, class c_loop> bool compute_cell(...)` | `__cxxwrap_compute_cell!(vc, con, cl)` | ✅ Exported (lambda, specific types) |
| **Particle Output (text)** |
| `draw_particles` | `template<class c_loop> void draw_particles(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_particles` | `void draw_particles(FILE *fp = stdout)` | `draw_particles(con, fptr)` | ✅ Exported |
| `draw_particles` | `void draw_particles(const char *filename)` | `draw_particles(con, path)` | ✅ Exported |
| **Particle Output (POV-Ray)** |
| `draw_particles_pov` | `template<class c_loop> void draw_particles_pov(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_particles_pov` | `void draw_particles_pov(FILE *fp = stdout)` | `draw_particles_pov(con, fptr)` | ✅ Exported |
| `draw_particles_pov` | `void draw_particles_pov(const char *filename)` | `draw_particles_pov(con, path)` | ✅ Exported |
| **Cell Output (Gnuplot)** |
| `draw_cells_gnuplot` | `template<class c_loop> void draw_cells_gnuplot(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_cells_gnuplot` | `void draw_cells_gnuplot(FILE *fp = stdout)` | `draw_cells_gnuplot(con, fptr)` | ✅ Exported |
| `draw_cells_gnuplot` | `void draw_cells_gnuplot(const char *filename)` | `draw_cells_gnuplot(con, path)` | ✅ Exported |
| **Cell Output (POV-Ray)** |
| `draw_cells_pov` | `template<class c_loop> void draw_cells_pov(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_cells_pov` | `void draw_cells_pov(FILE *fp = stdout)` | `draw_cells_pov(con, fptr)` | ✅ Exported |
| `draw_cells_pov` | `void draw_cells_pov(const char *filename)` | `draw_cells_pov(con, path)` | ✅ Exported |

# Export Status of `container_poly` Class Methods

This table shows which public methods from `container_poly` and its base classes (`container_base`, `wall_list`, `voro_base`, `radius_poly`) are exported to Julia via CxxWrap.jl.

| Method | Signature | Exported As | Status |
|--------|-----------|-------------|--------|
| **Constructor** |
| `container_poly` | `container_poly(double ax_, double bx_, double ay_, double by_, double az_, double bz_, int nx_, int ny_, int nz_, bool xperiodic_, bool yperiodic_, bool zperiodic_, int init_mem)` | `RawContainerPoly(...)` | ✅ Exported |
| `point_inside` | `bool point_inside(double x, double y, double z)` | `__cxxwrap_isinside(con, x, y, z)` | ✅ Exported |
| `draw_domain_gnuplot` | `void draw_domain_gnuplot(FILE *fp = stdout)` | — | ❌ Not exported |
| `draw_domain_gnuplot` | `void draw_domain_gnuplot(const char *filename)` | — | ❌ Not exported |
| `draw_domain_pov` | `void draw_domain_pov(FILE *fp = stdout)` | — | ❌ Not exported |
| `draw_domain_pov` | `void draw_domain_pov(const char *filename)` | — | ❌ Not exported |
| `total_particles` | `int total_particles()` | — | ❌ Not exported |
| `add_wall` | `void add_wall(wall *w)` | — | ❌ Not exported (generic) |
| `add_wall` | `void add_wall(wall &w)` | — | ❌ Not exported (generic) |
| `add_wall` | `void add_wall(wall_list &wl)` | — | ❌ Not exported |
| `point_inside_walls` | `bool point_inside_walls(double x, double y, double z)` | `__cxxwrap_point_inside_walls(con, x, y, z)` | ✅ Exported (lambda) |
| `apply_walls` | `template<class c_class> bool apply_walls(...)` | — | ❌ Not exported |
| **Particle Management** |
| `clear` | `void clear()` | `__cxxwrap_clear!(con)` | ✅ Exported |
| `put` | `void put(int n, double x, double y, double z, double r)` | `__cxxwrap_put!(con, id, x, y, z, r)` | ✅ Exported |
| `put` | `void put(particle_order &vo, int n, double x, double y, double z, double r)` | `__cxxwrap_put!(con, ord, id, x, y, z, r)` | ✅ Exported |
| `import` | `void import(FILE *fp = stdin)` | `import!(con, fptr)` | ✅ Exported |
| `import` | `void import(particle_order &vo, FILE *fp = stdin)` | — | ❌ Not exported |
| `import` | `void import(const char* filename)` | `import!(con, path)` | ✅ Exported |
| `import` | `void import(particle_order &vo, const char *filename)` | — | ❌ Not exported |
| `compute_all_cells` | `void compute_all_cells()` | `compute_all_cells(con)` | ✅ Exported |
| `sum_cell_volumes` | `double sum_cell_volumes()` | `sum_cell_volumes(con)` | ✅ Exported |
| **Cell Computation** |
| `find_voronoi_cell` | `bool find_voronoi_cell(double x, double y, double z, ...)` | `__cxxwrap_find_cell(con, x, y, z)` | ✅ Exported (lambda) |
| `compute_cell` | `template<class v_cell, class c_loop> bool compute_cell(...)` | — | ❌ Not exported |
| `compute_cell` | `template<class v_cell> bool compute_cell(v_cell &c, int ijk, int q)` | — | ❌ Not exported |
| **Particle Output (text)** |
| `draw_particles` | `template<class c_loop> void draw_particles(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_particles` | `void draw_particles(FILE *fp = stdout)` | `draw_particles(con, fptr)` | ✅ Exported |
| `draw_particles` | `void draw_particles(const char *filename)` | `draw_particles(con, path)` | ✅ Exported |
| **Particle Output (POV-Ray)** |
| `draw_particles_pov` | `template<class c_loop> void draw_particles_pov(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_particles_pov` | `void draw_particles_pov(FILE *fp = stdout)` | `draw_particles_pov(con, fptr)` | ✅ Exported |
| `draw_particles_pov` | `void draw_particles_pov(const char *filename)` | `draw_particles_pov(con, path)` | ✅ Exported |
| **Cell Output (Gnuplot)** |
| `draw_cells_gnuplot` | `template<class c_loop> void draw_cells_gnuplot(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_cells_gnuplot` | `void draw_cells_gnuplot(FILE *fp = stdout)` | `draw_cells_gnuplot(con, fptr)` | ✅ Exported |
| `draw_cells_gnuplot` | `void draw_cells_gnuplot(const char *filename)` | `draw_cells_gnuplot(con, path)` | ✅ Exported |
| **Cell Output (POV-Ray)** |
| `draw_cells_pov` | `template<class c_loop> void draw_cells_pov(c_loop &vl, FILE *fp)` | — | ❌ Not exported |
| `draw_cells_pov` | `void draw_cells_pov(FILE *fp = stdout)` | `draw_cells_pov(con, fptr)` | ✅ Exported |
| `draw_cells_pov` | `void draw_cells_pov(const char *filename)` | `draw_cells_pov(con, path)` | ✅ Exported |

## Public Data Members

| Member | Type | Exported As | Status |
|--------|------|-------------|--------|
| `max_radius` | `double` | — | ❌ Not exported |

# Export Status of `voronoicell_neighbor` Class Methods

This table shows which public methods from `voronoicell_neighbor` and its base class `voronoicell_base` are exported to Julia via CxxWrap.jl.

| Method | Signature | Exported As | Status |
|--------|-----------|-------------|--------|
| **Constructors/Destructor** |
| `voronoicell_neighbor` | `voronoicell_neighbor()` | `VoronoiCell()` | ✅ Exported |
| `voronoicell_neighbor` | `voronoicell_neighbor(double max_len_sq_)` | `VoronoiCell(max_len_sq)` | ✅ Exported |
| `voronoicell_neighbor` | `template<class c_class> voronoicell_neighbor(c_class &con)` | `VoronoiCell(con)` | ✅ Exported |
| **Assignment** |
| `operator=` | `void operator=(voronoicell_neighbor &c)` | `__cxxwrap_copy!(dest, src)` | ✅ Exported (lambda) |
| **Plane Cutting (with neighbor tracking)** |
| `nplane` | `bool nplane(double x, double y, double z, double rsq, int p_id)` | `__cxxwrap_nplane!(vc, x, y, z, rsq, p_id)` | ✅ Exported |
| `nplane` | `bool nplane(double x, double y, double z, int p_id)` | `__cxxwrap_nplane!(vc, x, y, z, p_id)` | ✅ Exported |
| `plane` | `bool plane(double x, double y, double z, double rsq)` | — | ❌ Not exported |
| `plane` | `bool plane(double x, double y, double z)` | — | ❌ Not exported |
| **Initialization** |
| `init` | `void init(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)` | `__cxxwrap_init!(vc, xmin, xmax, ymin, ymax, zmin, zmax)` | ✅ Exported |
| `init_octahedron` | `void init_octahedron(double l)` | `__cxxwrap_init_octahedron!(vc, l)` | ✅ Exported |
| `init_tetrahedron` | `void init_tetrahedron(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3)` | `__cxxwrap_init_tetrahedron!(vc, ...)` | ✅ Exported |
| **Translation** |
| `translate` | `void translate(double x, double y, double z)` | `__cxxwrap_translate!(vc, x, y, z)` | ✅ Exported |
| **Output (POV-Ray)** |
| `draw_pov` | `void draw_pov(double x, double y, double z, FILE *fp = stdout)` | `__cxxwrap_draw_pov(vc, x, y, z, fptr)` | ✅ Exported |
| `draw_pov` | `void draw_pov(double x, double y, double z, const char *filename)` | `__cxxwrap_draw_pov(vc, x, y, z, filename)` | ✅ Exported |
| `draw_pov_mesh` | `void draw_pov_mesh(double x, double y, double z, FILE *fp = stdout)` | `__cxxwrap_draw_pov_mesh(vc, x, y, z, fptr)` | ✅ Exported |
| `draw_pov_mesh` | `void draw_pov_mesh(double x, double y, double z, const char *filename)` | `__cxxwrap_draw_pov_mesh(vc, x, y, z, filename)` | ✅ Exported |
| **Output (Gnuplot)** |
| `draw_gnuplot` | `void draw_gnuplot(double x, double y, double z, FILE *fp = stdout)` | `__cxxwrap_draw_gnuplot(vc, x, y, z, fptr)` | ✅ Exported |
| `draw_gnuplot` | `void draw_gnuplot(double x, double y, double z, const char *filename)` | `__cxxwrap_draw_gnuplot(vc, x, y, z, filename)` | ✅ Exported |
| **Geometric Properties** |
| `volume` | `double volume()` | `volume(vc)` | ✅ Exported |
| `max_radius_squared` | `double max_radius_squared()` | `max_radius_squared(vc)` | ✅ Exported |
| `total_edge_distance` | `double total_edge_distance()` | `total_edge_distance(vc)` | ✅ Exported |
| `surface_area` | `double surface_area()` | `surface_area(vc)` | ✅ Exported |
| `centroid` | `void centroid(double &cx, double &cy, double &cz)` | `__cxxwrap_centroid(vc)` → `(x, y, z)` | ✅ Exported (lambda) |
| `number_of_faces` | `int number_of_faces()` | `number_of_faces(vc)` | ✅ Exported |
| `number_of_edges` | `int number_of_edges()` | `number_of_edges(vc)` | ✅ Exported |
| **Vertex Information** |
| `vertex_orders` | `void vertex_orders(std::vector<int> &v)` | `__cxxwrap_get_vertex_orders!(v, vc)` | ✅ Exported (lambda) |
| `output_vertex_orders` | `void output_vertex_orders(FILE *fp = stdout)` | — | ❌ Not exported |
| `vertices` | `void vertices(std::vector<double> &v)` | `__cxxwrap_vertices!(v, vc)` | ✅ Exported (lambda) |
| `vertices` | `void vertices(double x, double y, double z, std::vector<double> &v)` | `__cxxwrap_vertices!(v, vc, x, y, z)` | ✅ Exported (lambda) |
| `output_vertices` | `void output_vertices(FILE *fp = stdout)` | — | ❌ Not exported |
| `output_vertices` | `void output_vertices(double x, double y, double z, FILE *fp = stdout)` | — | ❌ Not exported |
| **Face Information** |
| `face_areas` | `void face_areas(std::vector<double> &v)` | `__cxxwrap_face_areas!(v, vc)` | ✅ Exported (lambda) |
| `output_face_areas` | `void output_face_areas(FILE *fp = stdout)` | — | ❌ Not exported |
| `face_orders` | `void face_orders(std::vector<int> &v)` | `__cxxwrap_face_orders!(v, vc)` | ✅ Exported (lambda) |
| `output_face_orders` | `void output_face_orders(FILE *fp = stdout)` | — | ❌ Not exported |
| `face_freq_table` | `void face_freq_table(std::vector<int> &v)` | `__cxxwrap_face_freq_table!(v, vc)` | ✅ Exported (lambda) |
| `output_face_freq_table` | `void output_face_freq_table(FILE *fp = stdout)` | — | ❌ Not exported |
| `face_vertices` | `void face_vertices(std::vector<int> &v)` | `__cxxwrap_face_vertices!(v, vc)` | ✅ Exported (lambda) |
| `output_face_vertices` | `void output_face_vertices(FILE *fp = stdout)` | — | ❌ Not exported |
| `face_perimeters` | `void face_perimeters(std::vector<double> &v)` | `__cxxwrap_face_perimeters!(v, vc)` | ✅ Exported (lambda) |
| `output_face_perimeters` | `void output_face_perimeters(FILE *fp = stdout)` | — | ❌ Not exported |
| **Normal Vectors** |
| `normals` | `void normals(std::vector<double> &v)` | `__cxxwrap_normals!(v, vc)` | ✅ Exported (lambda) |
| `output_normals` | `void output_normals(FILE *fp = stdout)` | — | ❌ Not exported |
| **Custom Output** |
| `output_custom` | `void output_custom(const char *format, FILE *fp = stdout)` | — | ❌ Not exported |
| `output_custom` | `void output_custom(const char *format, int i, double x, double y, double z, double r, FILE *fp)` | — | ❌ Not exported |
| **Plane Cutting** |
| `plane_intersects` | `bool plane_intersects(double x, double y, double z, double rsq)` | `__cxxwrap_plane_intersects(vc, x, y, z, rsq)` | ✅ Exported |
| `plane_intersects_guess` | `bool plane_intersects_guess(double x, double y, double z, double rsq)` | — | ❌ Not exported |
| **Neighbor Information (Virtual)** |
| `neighbors` | `virtual void neighbors(std::vector<int> &v)` | `__cxxwrap_get_neighbors!(v, vc)` | ✅ Exported (lambda) |
| `output_neighbors` | `virtual void output_neighbors(FILE *fp = stdout)` | — | ❌ Not exported |
| **Utilities** |
| `cycle_up` | `inline int cycle_up(int a, int p)` | `__cycle_up(vc, a, p)` | ✅ Exported |
| `cycle_down` | `inline int cycle_down(int a, int p)` | `__cycle_down(vc, a, p)` | ✅ Exported |
| `print_edges` | `void print_edges()` | `print_edges(vc)` | ✅ Exported |

## Public Data Members

| Member | Type | Exported As | Status |
|--------|------|-------------|--------|
| `current_vertices` | `int` | `__get_current_vertices(vc)` | ✅ Exported |
| `p` | `int` | `__get_p(vc)` | ✅ Exported |
| `up` | `int` | `__get_up(vc)` | ✅ Exported |
| `ed` | `int **` | `__cxxwrap_get_ed(vc)` | ✅ Exported |
| `nu` | `int *` | `__cxxwrap_get_nu(vc)` | ✅ Exported |
| `pts` | `double *` | `__cxxwrap_get_pts(vc)` | ✅ Exported |
| `mne` | `int **` | `__cxxwrap_get_mne(vc)` | ✅ Exported |
| `ne` | `int **` | `__cxxwrap_get_ne(vc)` | ✅ Exported |
| `ed[i][j]` | `int` | `__get_ed_ij(vc, i, j)` / `__set_ed_ij!(vc, i, j, k)` | ✅ Exported |
| `tol` | `double` | `__get_tol(vc)` | ✅ Exported |
| `tol_cu` | `double` | — | ❌ Not exported |
| `big_tol` | `double` | — | ❌ Not exported |
