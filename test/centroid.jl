@testset "Cell Statistics" begin

    v = VoronoiCell()

    # Initialize the Voronoi cell to be a cube of side length 2, centered
	# on the origin
	init!(v,-1,1,-1,1,-1,1)

    # Remove one edge of the cell with a single plane cut
	add_plane!(v,1,1,0,2)

    x = Ref(0.0)
    y = Ref(0.0)
    z = Ref(0.0)
    centroid!(v, x, y, z)
	@info "Centroid vector with Ref     : " (x[], y[], z[])

    c = Centroid(0.0, 0.0, 0.0)
    @test __centroid(v) === (x[], y[], z[])
    @info "Centroid vector with extern C : " (c.x, c.y, c.z)
    
end