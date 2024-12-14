
@testset "Ghost Test" begin
    
    x_min = -1
    x_max = 1
    y_min = -1
    y_max = 1
    z_min = -1
    z_max = 1
    cvol = (x_max - x_min)*(y_max - y_min)*(x_max - x_min)

    bx = 2.0
    bxy = 0.5
    by = 2.0
    bxz = 0.0
    byz = 0.0
    bz = 2.0
    
    n_x = 6
    n_y = 6
    n_z = 6
    m = 8
    
    particles = 20
	
    #conprdply = Container_Periodic_Poly( bx, bxy, by, bxz, byz, bz, n_x, n_y, n_z, m )
	#conprdply = Container_Periodic_Poly( (bx, bxy, by, bxz, byz, bz), (n_x, n_y, n_z), m )
    #conprdply = Container_Periodic_Poly(
    #    ;
    #    vectors = (2.0, 0.5, 2.0, 0.0, 0.0, 2.0),
    #    nblocks = (n_x, n_y, n_z),
    #    initial_memory = 8,
    #)

    conprdply = Container_Periodic_Poly(
            ;
            v_data = (bx, bxy, by, bxz, byz, bz),
            nblocks = (n_x, n_y, n_z),
            initial_memory = m,
        )
    
    @test conprdply isa Container_Periodic_Poly

    @test conprdply_add_point!(conprdply, 0, 1.18038, 1.78877, 0.0, 1.0) === nothing
    @test conprdply_add_point!(conprdply, 1, 0.59688, 0.823295, 0.0, 1.0) === nothing
    @test conprdply_add_point!(conprdply, 2, 1.67045, 0.536459, 0.0, 1.0) === nothing
    @test conprdply_add_point!(conprdply, 3, 0.60794, 1.95479, 0.0, 1.0) === nothing

	@test conprdply_draw_particles(conprdply, "ghost_test_p.gnu") === nothing
	@test conprdply_draw_cells_gnuplot(conprdply, "ghost_test_v.gnu") === nothing

    fp = "ghost_test_c.gnu"
	c = VoronoiCell()

	if ( conprdply_compute_ghost_cell( conprdply, c, 1.56, 0.67, 0.0, 1.0 ) )
        @test draw_gnuplot(fp, c, (1.56, 0.67, 0)) === nothing
    end

	@test conprdply_draw_domain_gnuplot(conprdply, "ghost_test_d.gnu") === nothing

end