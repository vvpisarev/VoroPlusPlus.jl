# Tests for region-restricted iteration (`VoroPlusPlus.Subset`, .known_issues #18)

struct NotASubsetContainer <: VoroPlusPlus.AbstractContainer end

function __ids(s)
    return sort!([Int(p.id) for p in eachparticle(s)])
end

"""Brute-force sphere query with the minimum-image convention on periodic axes only."""
function __bf_sphere(tess, center, r)
    per = periodicity(tess)
    (xl, yl, zl), (xh, yh, zh) = bounding_box(tess)
    lens = (xh - xl, yh - yl, zh - zl)
    ids = Int[]
    for p in eachparticle(tess)
        d = p.pos .- SVector{3,Float64}(center)
        ok = true
        @inbounds for k in 1:3
            if per[k]
                d = setindex(d, d[k] - lens[k] * round(d[k] / lens[k]), k)
            end
        end
        sum(abs2, d) <= r^2 && push!(ids, Int(p.id))
    end
    return sort!(ids)
end

"""Brute-force box query: periodic axes wrap, non-periodic axes clip against the box."""
function __bf_box(tess, lower, upper)
    per = periodicity(tess)
    (xl, yl, zl), (xh, yh, zh) = bounding_box(tess)
    lens = (xh - xl, yh - yl, zh - zl)
    l = SVector{3,Float64}(lower)
    u = SVector{3,Float64}(upper)
    ids = Int[]
    for p in eachparticle(tess)
        ok = true
        @inbounds for k in 1:3
            if per[k]
                ok &= mod(p.pos[k] - l[k], lens[k]) <= u[k] - l[k]
            else
                ok &= l[k] <= p.pos[k] <= u[k]
            end
        end
        ok && push!(ids, Int(p.id))
    end
    return sort!(ids)
end

@testset "Subset iteration: non-periodic containers" begin
    Random.seed!(20261125)
    bounds = ((-2.0, -2.0, 0.0), (2.0, 2.0, 1.0))
    pos = [
        SVector(-2.0 + 4.0 * rand(), -2.0 + 4.0 * rand(), 0.0 + 1.0 * rand())
        for _ in 1:1500
    ]
    tess = voronoi_tessellation(pos; bounds=bounds, periodic=(false, false, false))

    # interior sphere
    s = sphere_region(tess, (0.3, -0.2, 0.5), 0.3)
    @test __ids(s) == __bf_sphere(tess, (0.3, -0.2, 0.5), 0.3)
    @test !isempty(__ids(s))

    # sphere crossing the z = 1 face: result is the clipped set
    s = sphere_region(tess, (0.0, 0.0, 0.98), 0.25)
    @test __ids(s) == __bf_sphere(tess, (0.0, 0.0, 0.98), 0.25)

    # box region, also crossing faces
    s = box_region(tess, (-0.5, -0.5, 0.4), (0.5, 0.5, 1.4))
    @test __ids(s) == __bf_box(tess, (-0.5, -0.5, 0.4), (0.5, 0.5, 1.4))

    # empty region
    s = sphere_region(tess, (0.0, 0.0, -3.0), 0.1)
    @test isempty(collect(s))
    @test isempty(collect(eachcell(s)))
    @test Base.iterate(s) === nothing

    # bounds_test=false is a superset of bounds_test=true
    s_true = box_region(tess, (-0.25, -0.25, 0.25), (0.25, 0.25, 0.75); bounds_test=true)
    s_false = box_region(tess, (-0.25, -0.25, 0.25), (0.25, 0.25, 0.75); bounds_test=false)
    ids_true = __ids(s_true)
    ids_false = __ids(s_false)
    @test ids_true ⊆ ids_false
    @test length(ids_false) >= length(ids_true)
end

@testset "Subset iteration: fully periodic unit box" begin
    Random.seed!(20261126)
    pos = [SVector(rand(), rand(), rand()) for _ in 1:2000]
    tess = voronoi_tessellation(
        pos;
        bounds=((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
        periodic=(true, true, true),
    )

    queries = [
        ((0.98, 0.5, 0.5), 0.15),   # crosses the max face of x
        ((0.02, 0.5, 0.5), 0.15),   # crosses the min face of x (the upstream bug case)
        ((0.5, 0.03, 0.5), 0.15),   # crosses the min face of y
        ((0.02, 0.02, 0.98), 0.15), # two min faces (used to segfault)
        ((0.98, 0.97, 0.99), 0.15), # three max faces
        ((0.02, 0.5, 0.98), 0.15),  # min face of x, max face of z
        ((0.01, 0.01, 0.01), 0.15), # corner: all three min faces
        ((0.5, 0.5, 0.5), 0.2),     # interior
    ]
    for (center, r) in queries
        s = sphere_region(tess, center, r)
        ids = __ids(s)
        @test ids == __bf_sphere(tess, center, r)
        @test length(ids) == length(unique(ids)) # no duplicates
    end

    # box regions crossing min and max faces
    for (lower, upper) in [
        ((-0.1, 0.4, 0.4), (0.05, 0.6, 0.6)),  # min face of x
        ((0.9, 0.4, 0.4), (1.05, 0.6, 0.6)),   # max face of x
        ((-0.05, -0.05, 0.4), (0.05, 0.05, 0.6)), # two min faces
        ((0.02, 0.4, 0.95), (0.1, 0.6, 1.06)), # min x, max z
    ]
        s = box_region(tess, lower, upper)
        @test __ids(s) == __bf_box(tess, lower, upper)
    end

    # far query centers behave as their periodic image
    s = sphere_region(tess, (1000.02, 0.5, 0.5), 0.15)
    @test __ids(s) == __bf_sphere(tess, (0.02, 0.5, 0.5), 0.15)
    s = box_region(tess, (-99.1, 0.4, 0.4), (-98.95, 0.6, 0.6))
    @test __ids(s) == __bf_box(tess, (0.9, 0.4, 0.4), (1.05, 0.6, 0.6))
end

@testset "Subset iteration: partially periodic container" begin
    Random.seed!(20261127)
    bounds = ((-2.0, -2.0, 0.0), (2.0, 2.0, 1.0)) # periodic in x only, non-cubic
    pos = [SVector(-2.0 + 4.0 * rand(), -2.0 + 4.0 * rand(), rand()) for _ in 1:1500]
    tess = voronoi_tessellation(pos; bounds=bounds, periodic=(true, false, false))

    for (center, r) in [
        ((-1.98, 0.1, 0.5), 0.2),  # wraps across the x min face
        ((1.95, -0.3, 0.5), 0.15), # wraps across the x max face
        ((0.0, 0.0, 0.97), 0.15),  # crosses non-periodic z max: clipped
        ((0.0, 0.0, -0.05), 0.1),  # crosses non-periodic z min: clipped
    ]
        s = sphere_region(tess, center, r)
        @test __ids(s) == __bf_sphere(tess, center, r)
    end
    s = box_region(tess, (-2.1, -0.3, -0.1), (-1.9, 0.3, 1.1))
    @test __ids(s) == __bf_box(tess, (-2.1, -0.3, -0.1), (-1.9, 0.3, 1.1))
end

@testset "Subset iteration: polydisperse periodic container" begin
    Random.seed!(20261128)
    tess = polydisperse_container(
        ;
        bounds=((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
        nblocks=(8, 8, 8),
        periodic=(true, true, true),
    )
    for i in 1:1200
        add_point!(tess, i, (rand(), rand(), rand()), 0.02 + 0.03 * rand())
    end
    for (center, r) in [
        ((0.03, 0.5, 0.5), 0.15),  # min face
        ((0.01, 0.02, 0.99), 0.12), # corner
        ((0.5, 0.5, 0.5), 0.25),   # interior
    ]
        s = sphere_region(tess, center, r)
        ids = __ids(s)
        @test ids == __bf_sphere(tess, center, r)
        # particles carry radius information
        @test all(p -> p.radius !== nothing && p.radius > 0, collect(eachparticle(s)))
    end
end

@testset "Subset iteration: cells" begin
    Random.seed!(20261129)
    pos = [SVector(rand(), rand(), rand()) for _ in 1:1500]
    tess = voronoi_tessellation(
        pos;
        bounds=((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
        periodic=(true, true, true),
    )
    full = Dict{Int,Float64}()
    for (p, c) in tess
        full[Int(p.id)] = volume(c)
    end

    for (center, r) in [((0.02, 0.5, 0.5), 0.15), ((0.98, 0.97, 0.5), 0.15)]
        s = sphere_region(tess, center, r)
        ids = Int[]
        volsum = 0.0
        for (p, c) in s
            @test isvalid(c)
            v = volume(c)
            push!(ids, Int(p.id))
            volsum += v
            @test v ≈ full[Int(p.id)] rtol=1e-12
        end
        @test sort(ids) == __bf_sphere(tess, center, r)
        volsum_bf = sum(full[id] for id in ids)
        @test volsum ≈ volsum_bf rtol=1e-12
    end
end

@testset "Subset iteration: reuse and views" begin
    Random.seed!(20261130)
    pos = [SVector(rand(), rand(), rand()) for _ in 1:1200]
    tess = voronoi_tessellation(
        pos;
        bounds=((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
        periodic=(true, true, true),
    )
    s = sphere_region(tess, (0.02, 0.5, 0.5), 0.15)
    ids_1 = [Int(p.id) for (p, c) in s]
    vols_1 = [volume(c) for (p, c) in s]
    ids_2 = [Int(p.id) for (p, c) in s]
    vols_2 = [volume(c) for (p, c) in s]
    @test ids_1 == ids_2
    @test vols_1 == vols_2

    # eachparticle/eachcell agree with the pair iteration
    @test [Int(p.id) for p in eachparticle(s)] == ids_1
    @test [volume(c) for c in eachcell(s)] == vols_1

    # changing the region between passes is picked up
    VoroPlusPlus.reset_to_sphere!(s, (0.5, 0.5, 0.5), 0.1)
    ids_3 = [Int(p.id) for (p, c) in s]
    @test sort(ids_3) == __bf_sphere(tess, (0.5, 0.5, 0.5), 0.1)
    @test ids_3 != ids_1
    VoroPlusPlus.reset_to_box!(s, (0.4, 0.4, 0.4), (0.6, 0.6, 0.6))
    ids_box = __ids(s)
    @test ids_box == __bf_box(tess, (0.4, 0.4, 0.4), (0.6, 0.6, 0.6))
    @test ids_box != ids_3

    # two independent Subset objects can be iterated independently (interleaved)
    s1 = sphere_region(tess, (0.02, 0.5, 0.5), 0.15)
    s2 = sphere_region(tess, (0.98, 0.5, 0.5), 0.15)
    ids1 = [Int(p.id) for (p, c) in s1]
    ids2 = [Int(p.id) for (p, c) in s2]
    i1 = iterate(s1)
    i2 = iterate(s2)
    @test i1 !== nothing && i2 !== nothing
    seen1 = Int[]
    seen2 = Int[]
    while i1 !== nothing
        push!(seen1, Int(i1[1].first.id))
        if i2 !== nothing
            push!(seen2, Int(i2[1].first.id))
            i2 = iterate(s2, i2[2])
        end
        i1 = iterate(s1, i1[2])
    end
    while i2 !== nothing
        push!(seen2, Int(i2[1].first.id))
        i2 = iterate(s2, i2[2])
    end
    @test seen1 == ids1
    @test seen2 == ids2
    @test __ids(s1) == __bf_sphere(tess, (0.02, 0.5, 0.5), 0.15)
    @test __ids(s2) == __bf_sphere(tess, (0.98, 0.5, 0.5), 0.15)

    @test eltype(s) == Pair{eltype(eachparticle(s)), VoroPlusPlus.CheckedVoronoiCell}
    @test Base.IteratorSize(s) == Base.SizeUnknown()
end

@testset "Subset iteration: input validation" begin
    Random.seed!(20261131)
    pos = [SVector(rand(), rand(), rand()) for _ in 1:500]
    tess = voronoi_tessellation(
        pos;
        bounds=((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
        periodic=(true, true, true),
    )
    s = sphere_region(tess, (0.5, 0.5, 0.5), 0.1)

    for bad_center in [
        (NaN, 0.5, 0.5), (0.5, Inf, 0.5), (1e300, 0.5, 0.5), (0.5, -Inf, NaN),
    ]
        @test_throws ArgumentError VoroPlusPlus.reset_to_sphere!(s, bad_center, 0.1)
    end
    for bad_radius in [0.0, -0.1, NaN, Inf, 1e300]
        @test_throws ArgumentError VoroPlusPlus.reset_to_sphere!(s, (0.5, 0.5, 0.5), bad_radius)
    end
    @test_throws ArgumentError VoroPlusPlus.reset_to_sphere!(s, (0.5, 0.5), 0.1)
    @test_throws ArgumentError VoroPlusPlus.reset_to_sphere!(s, (0.5, 0.5, 0.5, 0.0), 0.1)

    # region spanning a whole period on a periodic axis
    @test_throws ArgumentError VoroPlusPlus.reset_to_sphere!(s, (0.5, 0.5, 0.5), 0.5)  # 2r == L
    @test_throws ArgumentError VoroPlusPlus.reset_to_sphere!(s, (0.5, 0.5, 0.5), 0.6)  # 2r > L
    @test_throws ArgumentError box_region(tess, (0.0, 0.4, 0.4), (1.0, 0.6, 0.6))     # w == L
    @test_throws ArgumentError box_region(tess, (0.0, 0.4, 0.4), (1.3, 0.6, 0.6))     # w > L
    @test_throws ArgumentError box_region(tess, (0.4, 0.4, 0.4), (0.4, 0.6, 0.6))     # empty
    @test_throws ArgumentError box_region(tess, (0.6, 0.4, 0.4), (0.4, 0.6, 0.6))     # inverted
    @test_throws ArgumentError box_region(tess, (0.0, NaN, 0.4), (0.2, 0.6, 0.6))

    # still alive and iterating after all those rejected queries
    @test __ids(s) == __bf_sphere(tess, (0.5, 0.5, 0.5), 0.1)

    # container classes without a `c_loop_subset` constructor are rejected early
    # (NOTE: `ContainerTriclinic` itself must not be constructed here: its upstream
    # `unitcell` constructor calls exit() for small boxes in this wrapper build)
    @test_throws ArgumentError VoroPlusPlus.Subset(NotASubsetContainer())
    @test_throws ArgumentError sphere_region(NotASubsetContainer(), (0.5, 0.5, 0.5), 0.1)

    # a fresh Subset cannot be iterated before a region is set
    fresh = VoroPlusPlus.Subset(tess)
    @test_throws ArgumentError collect(eachparticle(fresh))
end

@testset "Subset iteration: periodic min-face regression" begin
    # These are the exact queries from .known_issues #18: the raw `ContainerSubsetIterator`
    # returns 0 or 1 particles for them and can segfault; the public API must return the
    # full minimum-image set.
    Random.seed!(20261201)
    pos = [SVector(rand(), rand(), rand()) for _ in 1:2000]
    tess = voronoi_tessellation(
        pos;
        bounds=((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
        periodic=(true, true, true),
    )
    for center in [(0.02, 0.5, 0.5), (0.5, 0.03, 0.5), (0.02, 0.02, 0.98)]
        s = sphere_region(tess, center, 0.15)
        ids = __ids(s)
        @test !isempty(ids)
        @test ids == __bf_sphere(tess, center, 0.15)
    end
    s = box_region(tess, (-0.1, 0.4, 0.4), (0.05, 0.6, 0.6))
    @test !isempty(__ids(s))
    @test __ids(s) == __bf_box(tess, (-0.1, 0.4, 0.4), (0.05, 0.6, 0.6))

    # positions are reported in the primary domain; minimum-image helpers work
    s = sphere_region(tess, (0.02, 0.5, 0.5), 0.15)
    (xl, yl, zl), (xh, yh, zh) = bounding_box(tess)
    for p in eachparticle(s)
        @test all(xl .<= p.pos .<= [xh, yh, zh])
        @test VoroPlusPlus.minimum_image_distance(tess, p.pos, (0.02, 0.5, 0.5)) <= 0.15 + 1e-12
    end
end