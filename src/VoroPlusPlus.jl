module VoroPlusPlus
    using CxxWrap
    using Preferences
    using Printf: Format, format

    # Functions for Wall Class
    export Wall
    export point_inside

    # Functios for Wall_List Class
    export Wall_List
    export add_wall

    # Functions for Container Class
    export Container
    export add_point!, import!
    export draw_particles, draw_cells_gnuplot
    export draw_particles_pov
    export print_custom!, draw_cells_pov!
    export find_voronoi_cell
    export clear, compute_all_cells
    export sum_cell_volumes, point_inside!
    export region_count
    export initialize_search, frac_pos
    export region_index, draw_domain_gnuplot
    export draw_domain_pov, total_particles
    export add_wall!
    export point_inside_walls
    #Defined outside type Container, as an anonymus function
    export compute_cell!, compute_ghost_cell!
    export apply_walls!

    # Functions for Container Iterator(c_loop_all) Class
    export Container_Iterator
    export start!, next!, pos
    export ci_x, ci_y, ci_z, ci_pid

    # Functions for Container Iterator_Subset(c_loop_subset) Class
    export Container_Iterator_Subset
    export cis_start, cis_next, cis_pos
    export ci_x, ci_y, ci_z, ci_pid
    export cis_setup_sphere, cis_setup_box, cis_setup_intbox
    export cis_x, cis_y, cis_z, cis_pid

    # Functions for Particle_Order Class
    export Particle_Order
    export po_add

    # Functions for Container Iterator_Order(c_loop_order) Class
    export Container_Iterator_Order
    export cio_start, cio_next, cio_pos
    export cio_x, cio_y, cio_z, cio_pid
    
    # Functions for Container_Poly Class
    export Container_Poly
    export conp_add_point!, conp_import!
    export conp_draw_particles, conp_draw_cells_gnuplot
    export conp_draw_particles_pov
    export conp_print_custom!, conp_draw_cells_pov!

    # Functions for VoronoiCell Class
    export VoronoiCell
    export init!, init_l_shape!, add_plane!
    export volume, check_relations, check_duplicates
    export max_radius_squared, number_of_edges
    export total_edge_distance
    export number_of_faces, surface_area
    export draw_gnuplot!, draw_pov, draw_pov_mesh
    export num_vertices, root_vertex
    export init_octahedron, plane_intersects
    export centroid!, nplane, init_tetrahedron
    #Inherited from voronoicell_base
    export init_base
    export init_octahedron_base, init_tetrahedron_base
    export translate, plane_intersects_guess
    export construct_relations
    export print_edges, cycle_up, cycle_down
    #Implemented with get and set functions
    export vertex_positions, vertex_positions!
    # Implemented with extern C
    export draw_gnuplot, output_vertices, output_vertex_orders
    export output_face_perimeters, output_face_freq_table
    export output_face_orders, output_face_areas
    export output_normals, output_face_vertices

    # Functions for VoronoiCell_Neighbor Class
    export VoronoiCell_Neighbor
    export init, init_octahedron, init_tetrahedron
    export nplane_rsq, nplane, plane_rsq, plane
    export check_facets, print_edges_neighbors

    # Functions for Containter Periodic Poly (conprdply) Class
    export Container_Periodic_Poly
    export conprdply_add_point!
    export conprdply_compute_ghost_cell
    export conprdply_draw_particles
    export conprdply_draw_cells_gnuplot
    export conprdply_draw_domain_gnuplot

    # Functions for Wall_Sphere Class
    export Wall_Sphere
    export point_inside_wph
    export cut_cell_vc, cut_cell_vcn

   
    # Refactoring for Ref substitution
    export get_centroid
    export get_pos
    export find_voro_cell


    ####################################################3

    function set_wrapper_path(path::AbstractString="/usr/lib")
        # Set it in our runtime values, as well as saving it to disk
        @set_preferences!("VORO_JL_WRAPPER_PATH" => path)
        @info("Voro++ wrapper path set; restart your Julia session for this change to take effect!")
    end

    const VORO_JL_WRAPPER_PATH = @load_preference("VORO_JL_WRAPPER_PATH")

    @static if VORO_JL_WRAPPER_PATH !== nothing
        @wrapmodule(() -> joinpath(VORO_JL_WRAPPER_PATH, "libvoro++wrap"))

        include("config.jl")
        include("container.jl")
        include("cell.jl")
        include("cell_iter.jl")
        include("container_prd.jl")

        function __init__()
            @initcxx
        end
    end

end # module VoroPlusPlus
