module VoroPlusPlus
    using Base: OneTo, @propagate_inbounds
    using CxxWrap
    using LinearAlgebra: Symmetric, dot, eigen
    using Preferences
    using Printf: Format, format
    using StaticArrays

    # Functions for Wall Class
    export Wall
    export point_inside

    # Functios for Wall_List Class
    export Wall_List
    export add_wall

    # ordering

    # Functions for Container Class
    export AbstractContainer
    export Container
    export container, polydisperse_container

    export bounding_box
    export periodicity
    export add_point!
    export import_particles, import_polydisperse_particles
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
    export compute_cell!
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
    export voronoicell_box, voronoicell_tetrahedron, voronoicell_octahedron
    export reset_to_box!, reset_to_tetrahedron!, reset_to_octahedron!
    export cut_by_particle_position!
    export translate!
    export volume, centroid
    export number_of_vertices, number_of_edges, number_of_faces
    export vertex_positions, vertex_positions!
    export get_neighbors!, append_neighbors!, neighbors
    export normals, get_normals!
    export get_face_perimeters!, face_perimeters
    export get_face_vertices!, face_vertices
    export get_face_orders!, face_orders
    export total_edge_distance
    export surface_area
    export draw_pov, draw_pov_mesh
    export root_vertex
    export check_relations, check_duplicates
    export plane_intersects
    # Implemented with extern C
    export draw_gnuplot

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
    export get_pos
    export find_voro_cell


    ####################################################3

    function set_wrapper_path(path::AbstractString="/usr/lib")
        # Set it in our runtime values, as well as saving it to disk
        if isfile(joinpath(path, "libvoro++wrap.so"))
            @set_preferences!("VORO_JL_WRAPPER_PATH" => joinpath(path, "libvoro++wrap.so"))
            @info("Voro++ wrapper path set; restart your Julia session for this change to take effect!")
        elseif isfile(joinpath(path, "lib", "libvoro++wrap.so"))
            @set_preferences!("VORO_JL_WRAPPER_PATH" => joinpath(path, "lib", "libvoro++wrap.so"))
            @info("Voro++ wrapper path set; restart your Julia session for this change to take effect!")
        else
            mesg = "Could not find valid wrapper library at " *
                    "$(joinpath(path, "libvoro++wrap.so")) or " *
                    "$(joinpath(path, "lib", "libvoro++wrap.so")). " *
                    "The setting has no effect."
            @warn mesg
        end
    end

    const VORO_JL_WRAPPER_PATH = @load_preference("VORO_JL_WRAPPER_PATH")

    @static if VORO_JL_WRAPPER_PATH !== nothing
        include("typedefs.jl")
        @wrapmodule(() -> VORO_JL_WRAPPER_PATH)

        include("config.jl")
        include("nested_ptr.jl")
        include("container.jl")
        include("particle_info.jl")
        include("file_import.jl")
        include("cell.jl")
        include("iteration.jl")
        #include("container_prd.jl")

        function __init__()
            @initcxx
        end
    end

end # module VoroPlusPlus
