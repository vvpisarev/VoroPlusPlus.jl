@testset "Cell Statistics" begin
    # Initialize the Voronoi cell to be a cube of side length 2, centered
    # on the origin
    vc = voronoicell_box((-1, -1, -1), (1, 1, 1))

    # Remove one edge of the cell with a single plane cut
    cut_by_particle_position!(vc, (1.0, 1.0, 0.0, 2.0))

    @test number_of_vertices(vc) == 10

    @test vertex_positions!(StdVector{Float64}(), vc) == vertex_positions!(Float64[], vc)
    @test reinterpret(Float64, vertex_positions(vc)) == vertex_positions!(Float64[], vc)
    for vert in vertex_positions(vc)
        @test hypot(vert...) <= hypot((vert .- (1.0, 1.0, 0.0))...)
    end

    @test reinterpret(Float64, normals(vc)) == get_normals!(Float64[], vc)
end
