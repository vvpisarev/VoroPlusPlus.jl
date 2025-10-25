@testset "Radical" begin

    x_min = -3
    x_max = 3
    y_min = -3
    y_max = 3
    z_min = 0
    z_max = 6

    n_x = 3
    n_y = 3
    n_z = 3

    # Container Test
    con = Container(
            ;
            bounds = (x_min, x_max, y_min, y_max, z_min, z_max),
            nblocks = (n_x, n_y, n_z),
            periodic = (false, false, false),
            particles_per_block = 8,
        )

    @test con isa Container
	@test VoroPlusPlus.import!(con, "./data/pack_six_cube") === nothing
	@test draw_cells_gnuplot(con, "pack_six_cube.gnu") === nothing
	@test draw_cells_pov!(con, "pack_six_cube_v.pov") === nothing
	@test draw_particles_pov(con, "pack_six_cube_p.pov") === nothing


    # Container_Poly Test
    conp = Container_Poly(
            ;
            bounds = (x_min, x_max, y_min, y_max, z_min, z_max),
            nblocks = (n_x, n_y, n_z),
            periodic = (false, false, false),
            particles_per_block = 8,
        )

    @test conp isa Container_Poly
	@test conp_import!(conp, "./data/pack_six_cube_poly") === nothing
	@test conp_draw_cells_gnuplot(conp, "pack_six_cube_poly.gnu") === nothing
	@test conp_draw_cells_pov!(conp, "pack_six_cube_poly_v.pov") === nothing
	@test conp_draw_particles_pov(conp, "pack_six_cube_poly_p.pov") === nothing



end
