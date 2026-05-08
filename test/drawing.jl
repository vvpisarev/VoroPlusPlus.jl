@testset "Gnuplot output" begin
    rng = copy(Random.default_rng())
    Random.seed!(rng, 99887766)

    con = VoroPlusPlus.read_particles(
        "./data/pack_six_cube",
        ;
        bounds=((-3, -3, 0), (3, 3, 6)),
    )

    vc, = VoroPlusPlus.eachcell(con)

    # Output the Voronoi cell to a file, in the gnuplot format
    VoroPlusPlus.__cxxwrap_draw_gnuplot(vc[], 0.0, 0.0, 0.0, "./image_data/single_cell.gnu")
    cell_str = let buf = IOBuffer()
        draw_gnuplot(buf, vc)
        String(take!(buf))
    end
    @test read("./image_data/single_cell.gnu", String) == cell_str

    dx, dy, dz = rand(rng, Float64, 3)

    VoroPlusPlus.__cxxwrap_draw_gnuplot(vc[], dx, dy, dz, "./image_data/single_cell_displaced.gnu")
    cell_str = let buf = IOBuffer()
        draw_gnuplot(buf, vc, (dx, dy, dz))
        String(take!(buf))
    end
    @test read("./image_data/single_cell_displaced.gnu", String) == cell_str

    VoroPlusPlus.draw_particles(VoroPlusPlus.__raw(con), "./image_data/particles.gnu")
    VoroPlusPlus.draw_cells_gnuplot(VoroPlusPlus.__raw(con), "./image_data/cells.gnu")

    VoroPlusPlus.draw_gnuplot("./image_data/julia.*.gnu", con; domain=true)

    VoroPlusPlus.draw_gnuplot("./image_data/only", con; domain=false, particles=false)
    VoroPlusPlus.draw_gnuplot("./image_data/only", con; domain=false, cells=false)

    @test count('\n', read("./image_data/julia.domain.gnu", String)) == 20
    @test read("./image_data/particles.gnu", String) == read("./image_data/julia.pts.gnu", String)
    @test read("./image_data/cells.gnu", String) == read("./image_data/julia.cells.gnu", String)

    @test read("./image_data/particles.gnu", String) == read("./image_data/only_pts.gnu", String)
    @test read("./image_data/cells.gnu", String) == read("./image_data/only_cells.gnu", String)
end

@testset "POV-Ray output" begin
    rng = copy(Random.default_rng())
    Random.seed!(rng, 99887766)

    con = VoroPlusPlus.read_particles(
        "./data/pack_six_cube",
        ;
        bounds=((-3, -3, 0), (3, 3, 6)),
    )

    vc, = VoroPlusPlus.eachcell(con)

    # Output the Voronoi cell to a file, in the POV-Ray format
    VoroPlusPlus.__cxxwrap_draw_pov(vc[], 0.0, 0.0, 0.0, "./image_data/single_cell.pov")
    cell_str = let buf = IOBuffer()
        draw_pov(buf, vc)
        String(take!(buf))
    end
    @test read("./image_data/single_cell.pov", String) == cell_str

    # Output the Voronoi cell to a file, in the POV-Ray Mesh format
    VoroPlusPlus.__cxxwrap_draw_pov_mesh(vc[], 0.0, 0.0, 0.0, "./image_data/single_cell_mesh.pov")
    cell_str = let buf = IOBuffer()
        draw_pov_mesh(buf, vc)
        String(take!(buf))
    end
    @test read("./image_data/single_cell_mesh.pov", String) == cell_str

    dx, dy, dz = rand(rng, Float64, 3)

    VoroPlusPlus.__cxxwrap_draw_pov(vc[], dx, dy, dz, "./image_data/single_cell_displaced.pov")
    cell_str = let buf = IOBuffer()
        draw_pov(buf, vc, (dx, dy, dz))
        String(take!(buf))
    end
    @test read("./image_data/single_cell_displaced.pov", String) == cell_str

    VoroPlusPlus.__cxxwrap_draw_pov_mesh(vc[], dx, dy, dz, "./image_data/single_cell_displaced_mesh.pov")
    cell_str = let buf = IOBuffer()
        draw_pov_mesh(buf, vc, (dx, dy, dz))
        String(take!(buf))
    end
    @test read("./image_data/single_cell_displaced_mesh.pov", String) == cell_str

    VoroPlusPlus.draw_pov("./image_data/pack_six_cube.*.pov", con; domain=true, particles=true, cells=true)

    @test isfile("./image_data/pack_six_cube.domain.pov")
    @test isfile("./image_data/pack_six_cube.pts.pov")
    @test isfile("./image_data/pack_six_cube.cells.pov")
end
