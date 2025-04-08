@testset "Cell Statistics" begin

    v = VoronoiCell()

    # Initialize the Voronoi cell to be a cube of side length 2, centered
	# on the origin
	init!(v,-1,1,-1,1,-1,1)

    # Remove one edge of the cell with a single plane cut
	add_plane!(v,1,1,0,2)

    ## Output the Voronoi cell to a file in gnuplot format
	@test draw_gnuplot("simple_cell.gnu", v) === nothing

    # Output vertex-based statistics
	@info "Total vertices      : " num_vertices(v)
	@info "Vertex positions    : --> file: cell_statistics.txt"
    @test output_vertices("cell_statistics.txt", v, 0.0, 0.0, 0.0) === nothing
    #Base.Libc.flush_cstdio()
	@info "Vertex orders       : --> file: cell_statistics.txt"
    @test output_vertex_orders("cell_statistics.txt", v) === nothing
    #Base.Libc.flush_cstdio()
	@info "Max rad. sq. vertex : " 0.25*max_radius_squared(v)

	# Output edge-based statistics
	@info "Total edges         : " number_of_edges(v)
	@info "Total edge distance : " total_edge_distance(v)
	@info "Face perimeters     : --> file: cell_statistics.txt" 
    @test output_face_perimeters("cell_statistics.txt", v) === nothing

	# Output face-based statistics
	@info "Total faces         : " number_of_faces(v)
	@info "Surface area        : " surface_area(v)
	@info "Face freq. table    : --> file: cell_statistics.txt"
    @test output_face_freq_table("cell_statistics.txt", v) === nothing
	@info "Face orders         : --> file: cell_statistics.txt"
    @test output_face_orders("cell_statistics.txt", v) === nothing
	@info "Face areas          : --> file: cell_statistics.txt"
    @test output_face_areas("cell_statistics.txt", v) === nothing
	@info "Face normals        : --> file: cell_statistics.txt"
    @test output_normals("cell_statistics.txt", v) === nothing
	@info "Face vertices       : --> file: cell_statistics.txt"
    @test output_face_vertices("cell_statistics.txt", v) === nothing

    #Output volume and centroid
    @info "Volume              : " volume(v)
    x = Ref(0.0)
    y = Ref(0.0)
    z = Ref(0.0)
    centroid!(v, x, y, z)
	@info "Centroid vector     : " (x[], y[], z[])

end


@testset "Custom Output" begin
    
    x_min = -3
    x_max = 3
    y_min = -3
    y_max = 3
    z_min = 0
    z_max = 6

    n_x = 3
    n_y = 3
    n_z = 3

    con = Container(
            ;
            bounds = (x_min, x_max, y_min, y_max, z_min, z_max),
            nblocks = (n_x, n_y, n_z),
            periodic = (false, false, false),
            particles_per_block = 8,
        )

    @test con isa Container
	@test VoroPlusPlus.import!(con, "./data/pack_six_cube") === nothing
	@test print_custom!(con,
		"ID=%i, pos=(%x,%y,%z), vertices=%w, edges=%g, faces=%s",
		"packing.custom1") === nothing
	@test print_custom!(con, "%i %q %s %F %a %f %t %l %n","packing.custom2") === nothing
	@test print_custom!(con, "%i %q %v %c","packing.custom3") === nothing
	@test draw_cells_pov!(con, "pack_six_cube_v.pov") === nothing

end


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