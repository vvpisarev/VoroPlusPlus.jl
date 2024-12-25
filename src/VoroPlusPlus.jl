module VoroPlusPlus
    using CxxWrap
    using Preferences
    using Printf: Format, format

    # Functions for Container Class
    export Container
    export add_point!, import!
    export draw_particles, draw_cells_gnuplot
    export draw_particles_pov
    export print_custom!, draw_cells_pov!
    export find_voronoi_cell

    #Functions for Container Iterator(c_loop_all) Class
    export Container_Iterator
    export start!, next!, pos
    export compute_cell!


    # Functions for Container_Poly Class
    export Container_Poly
    export conp_add_point!, conp_import!
    export conp_draw_particles, conp_draw_cells_gnuplot
    export conp_draw_particles_pov
    export conp_print_custom!, conp_draw_cells_pov!


    # Functions for VoronoiCell Class
    export VoronoiCell
    export init!, add_plane!
    export num_vertices, root_vertex, volume
    export vertex_positions, vertex_positions!
    export init_l_shape!, check_relations, check_duplicates
    export draw_gnuplot, draw_pov, draw_pov_mesh
    export output_vertices, output_vertex_orders
    export output_face_perimeters
    export max_radius_squared, number_of_edges, total_edge_distance
    export number_of_faces, surface_area
    export output_face_freq_table, output_face_orders
    export output_face_areas, output_normals
    export output_face_vertices
    #export centroid
    export init_octahedron, plane_intersects

    #Functions for Containter Periodic Poly (conprdply) Class
    export Container_Periodic_Poly
    export conprdply_add_point!
    export conprdply_compute_ghost_cell
    export conprdply_draw_particles
    export conprdply_draw_cells_gnuplot
    export conprdply_draw_domain_gnuplot


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
