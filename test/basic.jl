@testset "Single Cell" begin
    v = VoronoiCell()

    # Initialize the Voronoi cell to be a cube of side length 2, centered
    # on the origin
    init!(v, -1, 1, -1, 1, -1, 1)

    # Cut the cell by 250 random planes which are all a distance 1 away
    # from the origin, to make an approximation to a sphere
    for _ in 1:250
        x = 2 * rand() - 1
        y = 2 * rand() - 1
        z = 2 * rand() - 1
        rsq = x * x + y * y + z * z
        if rsq > 0.01 && rsq < 1
            r = 1 / sqrt(rsq)
            x *= r
            y *= r
            z *= r
            add_plane!(v, x, y, z, 1)
        end
    end

    @info "" volume(v) 4/3 * pi / 8
    @test 0 < volume(v) - 4/3 * pi / 8 < 0.05
    # Output the Voronoi cell to a file, in the gnuplot format
    @test draw_gnuplot(v, 0, 0, 0, "single_cell.gnu") === nothing
end

@testset "Platonic Solids" begin

    Phi = 0.5 * (1 + sqrt(5.0))
    phi = 0.5 * (1 - sqrt(5.0))

    v = VoronoiCell(-2, 2, -2, 2, -2, 2)

    # Create a tetrahedron
    add_plane!(v, 1, 1, 1)
    add_plane!(v, 1, -1, -1)
    add_plane!(v, -1, 1, -1)
    add_plane!(v, -1, -1, 1)
    @info "Tetrahedron volume" volume(v)
    @test isapprox(volume(v), 9.0; atol=1e-8)

    # Create a cube. Since this is the default shape
    # we don't need to do any plane cutting.
    init!(v, -1, 1, -1, 1, -1, 1)
    @info "Cube volume" volume(v)
    @test isapprox(volume(v), 8.0; atol=1e-8)

    # Create an octahedron
    init!(v, -2, 2, -2, 2, -2, 2)
    add_plane!(v, 1, 1, 1)
    add_plane!(v, -1, 1, 1)
    add_plane!(v, 1, -1, 1)
    add_plane!(v, -1, -1, 1)
    add_plane!(v, 1, 1, -1)
    add_plane!(v, -1, 1, -1)
    add_plane!(v, 1, -1, -1)
    add_plane!(v, -1, -1, -1)
    @info "Octahedron volume" volume(v)
    @test isapprox(volume(v), 4.5; atol=1e-8)

    # Create a dodecahedron
    init!(v, -2, 2, -2, 2, -2, 2)
    add_plane!(v, 0, Phi, 1)
    add_plane!(v, 0, -Phi, 1)
    add_plane!(v, 0, Phi, -1)
    add_plane!(v, 0, -Phi, -1)
    add_plane!(v, 1, 0, Phi)
    add_plane!(v, -1, 0, Phi)
    add_plane!(v, 1, 0, -Phi)
    add_plane!(v, -1, 0, -Phi)
    add_plane!(v, Phi, 1, 0)
    add_plane!(v, -Phi, 1, 0)
    add_plane!(v, Phi, -1, 0)
    add_plane!(v, -Phi, -1, 0)

    r_in = sqrt(1 + Phi^2) / 2
    a = 2 * sqrt(3 - Phi) / Phi^2 * r_in
    vol_ref = 5 * Phi^3 / (6 - 2 * Phi) * a^3
    @info "Dodecahedron volume" volume(v)
    @test isapprox(volume(v), vol_ref; atol=1e-8)

    # Create an icosahedron
    init!(v, -2, 2, -2, 2, -2, 2)
    add_plane!(v, 1, 1, 1)
    add_plane!(v, -1, 1, 1)
    add_plane!(v, 1, -1, 1)
    add_plane!(v, -1, -1, 1)
    add_plane!(v, 1, 1, -1)
    add_plane!(v, -1, 1, -1)
    add_plane!(v, 1, -1, -1)
    add_plane!(v, -1, -1, -1)
    add_plane!(v, 0, phi, Phi)
    add_plane!(v, 0, phi, -Phi)
    add_plane!(v, 0, -phi, Phi)
    add_plane!(v, 0, -phi, -Phi)
    add_plane!(v, Phi, 0, phi)
    add_plane!(v, Phi, 0, -phi)
    add_plane!(v, -Phi, 0, phi)
    add_plane!(v, -Phi, 0, -phi)
    add_plane!(v, phi, Phi, 0)
    add_plane!(v, phi, -Phi, 0)
    add_plane!(v, -phi, Phi, 0)
    add_plane!(v, -phi, -Phi, 0)

    r_in = sqrt(3) / 2
    a = 2 * sqrt(3) / Phi^2 * r_in
    vol_ref = 5 * Phi^2 / 6 * a^3
    @info "Icosahedron volume" volume(v)
    @test isapprox(volume(v), vol_ref; atol=1e-8)
end

@testset "Ramdom points" begin
    (x_min, x_max) = (y_min, y_max) = (z_min, z_max) = -1, 1
    cvol = (x_max - x_min) * (y_max - y_min) * (z_max - z_min)

    n_x = n_y = n_z = 6

    nparticles = 20

    # Create a container with the geometry given above, and make it
    # non-periodic in each of the three coordinates. Allocate space for
    # eight particles within each computational block
    con = Container(
        ;
        bounds = (x_min, x_max, y_min, y_max, z_min, z_max),
        nblocks = (n_x, n_y, n_z),
        periodic = (false, false, false),
        particles_per_block = 8,
    )

    # Randomly add particles into the container
    for i in 1:nparticles
        x = x_min + rand() * (x_max - x_min)
        y = y_min + rand() * (y_max - y_min)
        z = z_min + rand() * (z_max - z_min)
        add_point!(con, i - 1, x, y, z)
    end

    # Sum up the volumes, and check that this matches the container volume
    vvol = sum(volume, con)
    @info "Container volume" cvol
    @info "Voronoi volume" vvol
    @info "Difference" cvol - vvol
    @test isapprox(vvol, cvol; atol=1e-8)
    @test isapprox(sum(volume, VoroPlusPlus.Unsafe(con)), cvol; atol=1e-8)

    # Output the particle positions in gnuplot format
    @test draw_particles(con, "random_points_p.gnu") === nothing

    # Output the Voronoi cells in gnuplot format
    @test draw_cells_gnuplot(con, "random_points_v.gnu") === nothing
end

@testset "File import" begin
    # Set up constants for the container geometry
    x_min, x_max = -5.0, 5.0
    y_min, y_max = -5.0, 5.0
    z_min, z_max = 0.0, 10.0

    # Set up the number of blocks that the container is divided into
    n_x = n_y = n_z = 6

    con = VoroPlusPlus.Container(
        x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, false, false, false, 8
    )

    @test con isa Container

    @test VoroPlusPlus.import!(con, "./data/pack_ten_cube") === nothing
    @test VoroPlusPlus.draw_cells_gnuplot(con, "pack_ten_cube.gnu") === nothing
end

# @testset "Convex Test" begin
#     v = VoronoiCell()
#     @test v isa VoronoiCell

#     @test init_l_shape!(v) === nothing

#     @test draw_gnuplot(v, 0, 0, 0, "single_cell.gnu") === nothing
#     lp = ls = -1
#     l = u = 1e-20

#     @test add_plane!(v, -1, 3, 0, 0.5)
#     @test draw_gnuplot(v, 0, 0, 0, "single_cell2.gnu") === nothing
#     @test add_plane!(v, -1, 3, 0.4, 0.53);
#     @test add_plane!(v, -1, 3, -0.4, 0.54);
#     print("cr")
#     check_relations(v)
#     check_duplicates(v)
#     print("fi")

#     suc = true

#     fmt = Format("%s lp=%d ls=%d l=%g u=%g up=%d\n")
#     format(stdout, fmt, suc ? "True" : "False", lp, ls, l, u, get_up(v))

#     @test draw_gnuplot(v, 0, 0, 0, "single_cell3.gnu") === nothing
# end
