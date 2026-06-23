struct ContainerDim
    # minimum and maximum dimensions for each coordinate
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    z_min::Float64
    z_max::Float64
    # number of blocks for each dimension
    n_x::Int
    n_y::Int
    n_z::Int
end

struct ExtendedContainer{T<:Union{Container,ContainerPoly}}
    con::T # outer container
    # inner bounds
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    z_min::Float64
    z_max::Float64
end

# Voronoi Tesselation
struct ParallelTessellation{T<:ExtendedContainer}
    # minimum and maximum dimensions for each coordinate
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    z_min::Float64
    z_max::Float64
    domain::Array{T, 3} #grid of subcontainers
    skin_distance::Float64
    owner::Vector{Int32} # vector owner of each particle ID
    periodic::NTuple{3,Bool}
end

const SMALL_PRIMES = @SVector Int32[
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97
]

function factorize(n::Integer)
    primes = SMALL_PRIMES
    @assert 1 <= n <= last(primes)^2
    factors = Int32[]
    for p in primes
        p * p > n && break
        while true
            f, m = fldmod(n, p)
            if m == 0
                push!(factors, p)
                n = f
            else
                break
            end
        end
    end
    if n > 1
        push!(factors, n)
    end
    return factors
end

function generate_containers(
    x_min::Float64,  # global bounds
    x_max::Float64,
    y_min::Float64,
    y_max::Float64,
    z_min::Float64,
    z_max::Float64,
    periodic::Ntuple{3,Bool},
    d_skin::Float64, # d_skin size
    nparticles::Integer, # total no. of particles
    ntasks::Int32 = Int32(1), # default number of tasks
)
    lx, ly, lz = (x_max, y_max, z_max) .- (x_min, y_min, z_min)
    subdomain_size = [lx, ly, lz]
    splits = Int32[1, 1, 1]
    factors = factorize(ntasks)
    while !isempty(factors)
        p = pop!(factors)
        split_dim = argmax(subdomain_size)
        subdomain_size[split_dim] /= p
        splits[split_dim] *= p
    end
    nnx, nny, nnz = splits
    range_x = range(x_min; stop=x_max, length=nnx+1)
    range_y = range(y_min; stop=y_max, length=nny+1)
    range_z = range(z_min; stop=z_max, length=nnz+1)
    ilscale = cbrt(nparticles / (OPT_PART_PER_BLOCK * lx * ly * lz))
    ox_min, oy_min, oz_min = (x_min, y_min, z_min) .- periodic .* d_skin
    ox_max, oy_max, oz_max = (x_max, y_max, z_max) .+ periodic .* d_skin
    containers = map(CartesianIndices((nnx, nny, nnz))) do ind
        ix, iy, iz = Tuple(ind)
        x_lo = max(ox_min, range_x[ix] - d_skin)
        x_hi = min(ox_max, range_x[ix+1] + d_skin)
        x_in_lo, x_in_hi = range_x[ix], range_x[ix+1]

        y_lo = max(oy_min, range_y[iy] - d_skin)
        y_hi = min(oy_max, range_y[iy+1] + d_skin)
        y_in_lo, y_in_hi = range_y[iy], range_y[iy+1]

        z_lo = max(oz_min, range_z[iz] - d_skin)
        z_hi = min(oz_max, range_z[iz+1] + d_skin)
        z_in_lo, z_in_hi = range_z[iz], range_z[iz+1]
        nx, ny, nz = floor.(
            Int32,
            (x_hi - x_lo, y_hi - y_lo, z_hi - z_lo) .* (ilscale + 1),
        )
        con = container(
            ;
            bounds = ((x_lo, y_lo, z_lo), (x_hi, y_hi, z_hi)),
            nblocks = (nx, ny, nz),
            periodic = (false, false, false),
            particles_per_block = 8,
        )
        return ExtendedContainer(
            con,
            x_in_lo,
            x_in_hi,
            y_in_lo,
            y_in_hi,
            z_in_lo,
            z_in_hi,
        )
    end

    return containers
end

"""
    parallel_voronoi_tessellation(pos; [id,] bounds, periodic=(false, false, false), ntasks)

Allocate space for a container of Voronoi cells and add coordinates stored in `pos`.
# Arguments
* `pos`: vector of positions
# Keywords
* `id`: use those IDs instead of indices of `pos`
* `bounds`: limits of the bounding box `((xmin, ymin, zmin), (xmax, ymax, zmax))`
* `periodic::NTuple{3,Bool}`: periodicity in each axis. Default: `(false, false, false)`
* `ordering`: `UnspecifiedOrder()` or `InsertionOrder()`. Default: `UnspecifiedOrder()`.
"""
function parallel_voronoi_tessellation(
    pos::AbstractVector
    ;
    bounds,
    id::AbstractVector{<:Integer}=OneTo(length(pos)),
    periodic=(false, false, false),
    ntasks::Integer=2,
)

    eachindex(bounds) == OneTo(2) &&
    eachindex(bounds[1]) == eachindex(bounds[2]) == OneTo(3) ||
    throw(DimensionMismatch("Bounds must be a 2-tuple of length-3 arrays or tuples"))

    eachindex(periodic) == OneTo(3) ||
    throw(DimensionMismatch("Periodic flags must be a length-3 boolean array or tuple"))

    eachindex(pos) == eachindex(id) ||
    throw(DimensionMismatch("Position and ID vectors must have the same dimensions"))

    ((ax, ay, az), (bx, by, bz)) = bounds
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.((ax, bx, ay, by, az, bz))
    px, py, pz = (Bool(p) for p in periodic)
    ppb = Int32(8)

    # container length in earch coordinate
    lx, ly, lz = Float64.((x_max - x_min, y_max - y_min, z_max - z_min))

    # container volume
    vol = lx * ly * lz

    # number of paticles
    npts = length(coords)

    # mean interparticle distance
    r_ave = cbrt(vol / npts)
    d_skin = 3 * r_ave

    containers = generate_containers(
        x_min, x_max, y_min, y_max, z_min, z_max, d_skin, npts, Int32(ntasks)
    )
    owner = zeros(Int32, npts) # task that owns each particle

    tessellation = ParallelTessellation(
        x_min, x_max, y_min, y_max, z_min, z_max,
        containers,
        d_skin,
        owner,
        periodic
    )

    nx, ny, nz = size(containers)
    dx, dy, dz = (lx, ly, lz) ./ (nx, ny, nz)

    # particle
    p = 0

    # sort particles into bins

    #vector of vectors(number of task)
    workload = [Int32[] for _ in CartesianIndices(containers)]
    ghosts = SVector{3,Float64}[]
    for ind in eachindex(coords)
        p += 1
        x, y, z = coords[ind]
        (ix::Int, xoff), (iy::Int, yoff), (iz::Int, zoff) =
            fldmod.((x, y, z) .- (x_min, y_min, z_min), (dx, dy, dz))
        ix += one(ix)
        iy += one(iy)
        iz += one(iz)
        i_con = LinearIndices(workload)[ix, iy, iz]
        push!(workload[i_con], p)
        owner[p] = i_con
        for dix in -1:1
            jx = ix + dix
            jx in axes(containers, 1) || continue
            xoffj = xoff - dix * dx
            -d_skin < xoffj <= dx + d_skin || continue
            for diy in -1:1
                jy = iy + diy
                jy in axes(containers, 2) || continue
                yoffj = yoff - diy * dy
                -d_skin < yoffj <= dy + d_skin || continue
                for diz in -1:1
                    (dix, diy, diz) == (0, 0, 0) && continue
                    jz = iz + diz
                    jz in axes(containers, 3) || continue
                    zoffj = zoff - diz * dz
                    -d_skin < zoffj <= dz + d_skin || continue
                    push!(workload[jx, jy, jz], p)
                end
            end
        end
    end


    @sync for i_con in 1:ntasks
        @spawn let
            work = workload[i_con]
            con = tessellation.domain[i_con].con
            for i_part in work
                ind = (i_part - 1) * 3
                x, y, z = coords[ind+1], coords[ind+2], coords[ind+3]
                add_point!(con, i_part, x, y, z)
            end
        end

        #=cell = VoronoiCell()
        for i_part in ghosts[i_con]
            ind = (i_part - 1) * 3
            x, y, z = coords[ind+1], coords[ind+2], coords[ind+3]
            compute_ghost_cell!(cell, con, x, y, z)
        end=#

        # After particle assigment, each thread calculates tesselation and generates files for particles and cells
        # draw_particles(tessellation.domain[i_con], "parallel/particles_$(i_con).gnu")
        # draw_cells_gnuplot(tessellation.domain[i_con], "parallel/voro_cells_$(i_con).gnu")

    end

    return tessellation
end

############################################################
function Base.foreach(fn, vt::ParallelTessellation)
    @sync for con_id in eachindex(vt.domain)
        @spawn let
            pid = Ref(Int32(0))
            x = Ref(0.0)
            y = Ref(0.0)
            z = Ref(0.0)
            r = Ref(0.0)
            c = VoronoiCell()
            con = vt.domain[con_id].con
            cla = ContainerIterator(con)
            if ( start!(cla) )
                while true
                    pos(cla, pid, x, y, z, r)
                    vt.owner[pid[]] == con_id && compute_cell!(c, con, cla) && fn(c)
                    next!(cla) || break
                end
            end
        end
    end
    return nothing
end

"""
    mapreduce(fn, op, vt::ParallelTessellation; init)

Apply function `fn` to each Voronoi cell in `vt`, and then reduce the result using the
    binary function `op`. `op` is assumed to be commutative and transitive, and the order of
    application of `fn` must not have effect on the result. `init` must be neutral element
    for `op` and has to be provided as it defines the type of output.
"""
function Base.mapreduce(fn, op, vt::ParallelTessellation; init)
    locals = Channel{typeof(init)}(length(vt.domain))
    @sync for con_id in LinearIndices(vt.domain)
        @spawn let
            val_local = init
            pid = Ref(Int32(0))
            x = Ref(0.0)
            y = Ref(0.0)
            z = Ref(0.0)
            r = Ref(0.0)
            c = VoronoiCell()
            con = vt.domain[con_id].con
            cla = Container_Iterator(con)
            if start!(cla)
                while true
                    pos(cla, pid, x, y, z, r)
                    if vt.owner[pid[]] == con_id && compute_cell!(c, con, cla)
                        val_local = op(val_local, fn(c))
                    end
                    next!(cla) || break
                end
                put!(locals, val_local)
            end
        end
    end
    val = init
    while !isempty(locals)
        val = op(val, take!(locals))
    end
    return val
end

"""
    sum(fn, vt::ParallelTessellation; init)

Sum the results of calling function `fn` on each cell in `vt`. The order of application of
    `fn` must not have effect on the result. `init` must be zero and has to be provided as
    it defines the type of output.
"""
function Base.sum(fn, vt::ParallelTessellation; init)
    locals = Channel{typeof(init)}(length(vt.domain))
    @sync for con_id in LinearIndices(vt.domain)
        @spawn let
            val_local = init
            pid = Ref(Int32(0))
            x = Ref(0.0)
            y = Ref(0.0)
            z = Ref(0.0)
            r = Ref(0.0)
            c = VoronoiCell()
            con = vt.domain[con_id].con
            cla = Container_Iterator(con)
            if start!(cla)
                while true
                    pos(cla, pid, x, y, z, r)
                    if vt.owner[pid[]] == con_id && compute_cell!(c, con, cla)
                        val_local += fn(c)
                    end
                    next!(cla) || break
                end
                put!(locals, val_local)
            end
        end
    end
    val = init
    while !isempty(locals)
        val += take!(locals)
    end
    return val
end

function GetVoroVolume(vt::VoronoiTessellation)
    return sum(volume, vt; init=0.0)

    # domain_vol = [0.0 for _ in vt.domain]

    # Threads.@threads for con_id in eachindex(vt.domain)
    #     pid = Ref(Int32(0))
    #     x = Ref(0.0)
    #     y = Ref(0.0)
    #     z = Ref(0.0)
    #     r = Ref(0.0)
    #     c = VoronoiCell()
    #     con = vt.domain[con_id].con
    #     cla = Container_Iterator(con)
    #     pvol = 0.0
    #     if ( start!(cla) )
    #         pos(cla, pid, x, y, z, r)
    #         if vt.owner[pid[]] == con_id
    #             if ( compute_cell!(c, con, cla) )
    #                 pvol += volume(c)
    #             end
    #         end
    #         while ( next!(cla) )
    #             pos(cla, pid, x, y, z, r)
    #             if vt.owner[pid[]] == con_id
    #                 if ( compute_cell!(c, con, cla) )
    #                     pvol += volume(c)
    #                 end
    #             end
    #         end
    #     end
    #     domain_vol[con_id] = pvol
    # end

    # return sum(domain_vol)

end
