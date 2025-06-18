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

    vertices_matrix = vertex_positions(Matrix, vc)
    vertices_vec = vertex_positions(vc)
    vertices_stdvec = vertex_positions!(StdVector{Float64}(), vc)
    @test eachcol(vertices_matrix) == vertices_vec
    @test reinterpret(Float64, vertices_vec) == vertices_stdvec
    @test reinterpret(Float64, vertices_vec) == vertex_positions!(Float64[], vc)
    @test vertices_vec == vertex_positions!(SVector{3,Float64}[], vc, (1, 1, 1)) .- Ref(SVector(1, 1, 1))
    @test vertices_vec == vertex_positions(vc, (1, 1, 1)) .- Ref(SVector(1, 1, 1))
    @test vertices_stdvec == vertex_positions!(StdVector{Float64}(), vc, (0, 0, 0))
    @test vertices_matrix == vertex_positions!(zero(vertices_matrix), vc, (1, 1, 1)) .- 1

    @test get_face_perimeters!(Float64[], vc) == get_face_perimeters!(StdVector{Float64}(), vc)
end
