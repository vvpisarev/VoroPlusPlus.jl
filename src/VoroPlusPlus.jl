module VoroPlusPlus
    using CxxWrap
    using Preferences
    using Printf: Format, format

    export Container
    export add_point!, import!
    export draw_particles, draw_cells_gnuplot
    export VoronoiCell
    export init!, add_plane!
    export num_vertices, root_vertex, volume
    export vertex_positions, vertex_positions!
    export init_l_shape!, check_relations, check_duplicates
    export draw_gnuplot

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

        function __init__()
            @initcxx
        end
    end

end # module VoroPlusPlus
