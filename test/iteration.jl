@testset "Iteration over containers" begin
    x_min, x_max = -5.0, 5.0
    y_min, y_max = -5.0, 5.0
    z_min, z_max = 0.0, 10.0

    con = VoroPlusPlus.read_particles(
        "data/pack_six_cube"
        ;
        bounds=((x_min, y_min, z_min), (x_max, y_max, z_max)),
        periodic=(false, false, false),
    )

    part_pos = SVector{3,Float64}[]
    cell_vols = Float64[]
    for (p, cell) in con
        push!(part_pos, p.pos)
        push!(cell_vols, volume(cell))
    end

    @test all(zip(part_pos, eachparticle(con))) do (pos, part)
        pos == part.pos
    end
    @test all(zip(cell_vols, eachcell(con))) do (vol, cell)
        vol == volume(cell)
    end

    @test all(zip(part_pos, cell_vols, VoroPlusPlus.Unsafe(con))) do (p1, v1, (part, cell))
        p1 == part.pos && v1 == volume(cell)
    end

    @test all(zip(part_pos, eachparticle(VoroPlusPlus.Unsafe(con)))) do (pos, part)
        pos == part.pos
    end
    @test all(zip(cell_vols, eachcell(VoroPlusPlus.Unsafe(con)))) do (vol, cell)
        vol == volume(cell)
    end

    cells_u = [cell for (_, cell) in VoroPlusPlus.Unsafe(con)]
    @test all(cells_u) do cell; cell === cells_u[1]; end

    cells_u = collect(eachcell(VoroPlusPlus.Unsafe(con)))
    @test all(cells_u) do cell; cell === cells_u[1]; end
end
