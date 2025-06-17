@testset "Single Cell" begin
    rng = copy(Random.default_rng())
    Random.seed!(rng, 99887766)
    # Initialize the Voronoi cell to be a cube of side length 2, centered
    # on the origin
    vc = voronoicell_box((-1, -1, -1), (1, 1, 1))
    @info "Vertex positions"
    foreach(println, vertex_positions(vc))

    @info "" volume(vc)

    # Cut the cell by 250 random planes which are all a distance 1 away
    # from the origin, to make an approximation to a sphere
    for _ in 1:250
        x = 2 * rand(rng) - 1
        y = 2 * rand(rng) - 1
        z = 2 * rand(rng) - 1
        rsq = x * x + y * y + z * z
        if rsq > 0.01 && rsq < 1
            r = inv(sqrt(rsq))
            x *= r
            y *= r
            z *= r
            cut_by_particle_position!(vc, (x, y, z, 1.0))
        end
    end

    @info "" volume(vc) 4/3 * pi / 8
    @test 0 < volume(vc) - 4/3 * pi / 8 < 0.05
    # Output the Voronoi cell to a file, in the gnuplot format
    draw_gnuplot("single_cell.gnu", vc)
    buf = IOBuffer()
    draw_gnuplot(buf, vc)
    cell_str = String(take!(buf))
    # open("single_cell_jl.gnu", "w") do io; write(io, cell_str); end
    @test read("single_cell.gnu", String) == cell_str

    vertices_matrix = vertex_positions(Matrix, vc)
    vertices_vec = vertex_positions(vc)
    vertices_stdvec = vertex_positions!(StdVector{Float64}(), vc)
    @test eachcol(vertices_matrix) == vertices_vec
    @test reinterpret(Float64, vertices_vec) == vertices_stdvec
    @test reinterpret(Float64, vertices_vec) == vertex_positions!(Float64[], vc)
    @test vertices_vec == vertex_positions!(SVector{3,Float64}[], vc, (0, 0, 0))
    @test vertices_vec == vertex_positions(vc, (0, 0, 0))
end

@testset "Platonic Solids" begin

    Phi = 0.5 * (1 + sqrt(5.0))
    phi = 0.5 * (1 - sqrt(5.0))

    vc = voronoicell_box((-2, -2, -2), (2, 2, 2))

    # Create a tetrahedron
    cut_by_particle_position!(vc, (1, 1, 1))
    cut_by_particle_position!(vc, (1, -1, -1))
    cut_by_particle_position!(vc, (-1, 1, -1))
    cut_by_particle_position!(vc, (-1, -1, 1))
    @info "Tetrahedron volume" volume(vc)
    @test isapprox(volume(vc), 9.0; atol=1e-8)
    tet = mktemp() do _, io
        draw_gnuplot(io, vc)
        seek(io, 0)
        read(io, String)
    end

    @info "" tet

    # Create a cube. Since this is the default shape
    # we don't need to do any plane cutting.
    reset_to_box!(vc, (-1, -1, -1), (1, 1, 1))
    @info "Cube volume" volume(vc)
    @test isapprox(volume(vc), 8.0; atol=1e-8)

    # Create an octahedron
    reset_to_box!(vc, (-2, -2, -2), (2, 2, 2))
    cut_by_particle_position!(vc, (1, 1, 1))
    cut_by_particle_position!(vc, (-1, 1, 1))
    cut_by_particle_position!(vc, (1, -1, 1))
    cut_by_particle_position!(vc, (-1, -1, 1))
    cut_by_particle_position!(vc, (1, 1, -1))
    cut_by_particle_position!(vc, (-1, 1, -1))
    cut_by_particle_position!(vc, (1, -1, -1))
    cut_by_particle_position!(vc, (-1, -1, -1))
    @info "Octahedron volume" volume(vc)
    @test isapprox(volume(vc), 4.5; atol=1e-8)

    # Create a dodecahedron
    reset_to_box!(vc, (-2, -2, -2), (2, 2, 2))
    cut_by_particle_position!(vc, (0, Phi, 1))
    cut_by_particle_position!(vc, (0, -Phi, 1))
    cut_by_particle_position!(vc, (0, Phi, -1))
    cut_by_particle_position!(vc, (0, -Phi, -1))
    cut_by_particle_position!(vc, (1, 0, Phi))
    cut_by_particle_position!(vc, (-1, 0, Phi))
    cut_by_particle_position!(vc, (1, 0, -Phi))
    cut_by_particle_position!(vc, (-1, 0, -Phi))
    cut_by_particle_position!(vc, (Phi, 1, 0))
    cut_by_particle_position!(vc, (-Phi, 1, 0))
    cut_by_particle_position!(vc, (Phi, -1, 0))
    cut_by_particle_position!(vc, (-Phi, -1, 0))

    r_in = sqrt(1 + Phi^2) / 2
    a = 2 * sqrt(3 - Phi) / Phi^2 * r_in
    vol_ref = 5 * Phi^3 / (6 - 2 * Phi) * a^3
    @info "Dodecahedron volume" volume(vc)
    @test isapprox(volume(vc), vol_ref; atol=1e-8)

    # Create an icosahedron
    reset_to_box!(vc, (-2, -2, -2), (2, 2, 2))
    cut_by_particle_position!(vc, (1, 1, 1))
    cut_by_particle_position!(vc, (-1, 1, 1))
    cut_by_particle_position!(vc, (1, -1, 1))
    cut_by_particle_position!(vc, (-1, -1, 1))
    cut_by_particle_position!(vc, (1, 1, -1))
    cut_by_particle_position!(vc, (-1, 1, -1))
    cut_by_particle_position!(vc, (1, -1, -1))
    cut_by_particle_position!(vc, (-1, -1, -1))
    cut_by_particle_position!(vc, (0, phi, Phi))
    cut_by_particle_position!(vc, (0, phi, -Phi))
    cut_by_particle_position!(vc, (0, -phi, Phi))
    cut_by_particle_position!(vc, (0, -phi, -Phi))
    cut_by_particle_position!(vc, (Phi, 0, phi))
    cut_by_particle_position!(vc, (Phi, 0, -phi))
    cut_by_particle_position!(vc, (-Phi, 0, phi))
    cut_by_particle_position!(vc, (-Phi, 0, -phi))
    cut_by_particle_position!(vc, (phi, Phi, 0))
    cut_by_particle_position!(vc, (phi, -Phi, 0))
    cut_by_particle_position!(vc, (-phi, Phi, 0))
    cut_by_particle_position!(vc, (-phi, -Phi, 0))

    r_in = sqrt(3) / 2
    a = 2 * sqrt(3) / Phi^2 * r_in
    vol_ref = 5 * Phi^2 / 6 * a^3
    @info "Icosahedron volume" volume(vc)
    @test isapprox(volume(vc), vol_ref; atol=1e-8)
end

@testset "Random points" begin
    rng = copy(Random.default_rng())
    Random.seed!(rng, 9876)

    (x_min, x_max) = (y_min, y_max) = (z_min, z_max) = -1, 1
    cvol = (x_max - x_min) * (y_max - y_min) * (z_max - z_min)

    n_x = n_y = n_z = 6

    nparticles = 20

    # Create a container with the geometry given above, and make it
    # non-periodic in each of the three coordinates. Allocate space for
    # eight particles within each computational block
    con = VoroPlusPlus.container(
        ;
        bounds = ((x_min, y_min, z_min), (x_max, y_max, z_max)),
        nblocks = (n_x, n_y, n_z),
        periodic = (false, false, false),
        particles_per_block = 8,
    )

    # Randomly add particles into the container
    for i in 1:nparticles
        x = x_min + rand(rng) * (x_max - x_min)
        y = y_min + rand(rng) * (y_max - y_min)
        z = z_min + rand(rng) * (z_max - z_min)
        add_point!(con, i - 1, (x, y, z))
    end

    # Sum up the volumes, and check that this matches the container volume
    vvol = sum(con) do (part, cell)
        volume(cell)
    end

    vvol_unsafe = sum(VoroPlusPlus.Unsafe(con)) do (part, cell)
        volume(cell)
    end
    @info "Container volume" cvol
    @info "Voronoi volume" vvol
    @info "Difference" cvol - vvol
    @test isapprox(vvol, cvol; atol=1e-8)
    @test isapprox(vvol_unsafe, cvol; atol=1e-8)
    neigh_vec = Int32[]
    neigh_std_vec = StdVector{Int32}()
    for (part, cell) in con
        @test get_neighbors!(neigh_vec, cell) == get_neighbors!(neigh_std_vec, cell)
    end

    # Output the particle positions in gnuplot format
    #@test draw_particles(con, "random_points_p.gnu") === nothing

    # Output the Voronoi cells in gnuplot format
    #@test draw_cells_gnuplot(con, "random_points_v.gnu") === nothing
end

@testset "File import" begin
    # Set up constants for the container geometry
    x_min, x_max = -5.0, 5.0
    y_min, y_max = -5.0, 5.0
    z_min, z_max = 0.0, 10.0

    # Set up the number of blocks that the container is divided into
    n_x = n_y = n_z = 6

    con = VoroPlusPlus.container(
        ;
        bounds=((x_min, y_min, z_min), (x_max, y_max, z_max)),
        nblocks=(n_x, n_y, n_z),
        periodic=(false, false, false),
        particles_per_block=8,
    )

    @test bounding_box(con) == ((x_min, y_min, z_min), (x_max, y_max, z_max))
    @test periodicity(con) == (false, false, false)

    #@test VoroPlusPlus.import!(con, "./data/pack_ten_cube") === nothing
    #@test VoroPlusPlus.draw_cells_gnuplot(con, "pack_ten_cube.gnu") === nothing
end

# @testset "Convex Test" begin
#     v = VoronoiCell()
#     @test v isa VoronoiCell

#     @test init_l_shape!(v) === nothing

#     @test draw_gnuplot(v, 0, 0, 0, "single_cell.gnu") === nothing
#     lp = ls = -1
#     l = u = 1e-20

#     @test cut_by_particle_position!(v, -1, 3, 0, 0.5)
#     @test draw_gnuplot(v, 0, 0, 0, "single_cell2.gnu") === nothing
#     @test cut_by_particle_position!(v, -1, 3, 0.4, 0.53);
#     @test cut_by_particle_position!(v, -1, 3, -0.4, 0.54);
#     print("cr")
#     check_relations(v)
#     check_duplicates(v)
#     print("fi")

#     suc = true

#     fmt = Format("%s lp=%d ls=%d l=%g u=%g up=%d\n")
#     format(stdout, fmt, suc ? "True" : "False", lp, ls, l, u, root_vertex(v))

#     @test draw_gnuplot(v, 0, 0, 0, "single_cell3.gnu") === nothing
# end
