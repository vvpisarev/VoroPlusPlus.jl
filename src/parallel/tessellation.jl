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
struct ParallelTessellation{T<:ExtendedContainer} <: AbstractTessellation{T}
    # minimum and maximum dimensions for each coordinate
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    z_min::Float64
    z_max::Float64
    domain::Array{T, 3} #grid of subcontainers
    skin_distance::Float64
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
    periodic::NTuple{3,Bool},
    d_skin::Float64, # d_skin size
    nparticles::Integer, # total no. of particles
    ntasks::Int32 = Int32(1), # default number of tasks
)
    lx, ly, lz = (x_max, y_max, z_max) .- (x_min, y_min, z_min)
    subdomain_size = @MVector [lx, ly, lz]
    splits = @MVector Int32[1, 1, 1]
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

function __prepare_coords(pos, ((xmin, ymin, zmin), (xmax, ymax, zmax)), periodic, rneigh)
    function __wrap(x, xmin, xmax, lx)
        xmin < x <= xmax ? x : xmin + mod(x - xmin, lx)
    end

    trunc_pos = SVector{3,Float64}[]
    id = Int32[]
    px, py, pz = periodic
    lx, ly, lz = (xmax, ymax, zmax) .- (xmin, ymin, zmin)

    for (ind, (x, y, z)) in pairs(pos)
        xx = px ? __wrap(x, xmin, xmax, lx) : x
        yy = py ? __wrap(y, ymin, ymax, ly) : y
        zz = pz ? __wrap(z, zmin, zmax, lz) : z

        if (xmin < xx <= xmax) & (ymin < yy <= ymax) & (zmin < zz <= zmax)
            push!(trunc_pos, (xx, yy, zz))
            push!(id, ind)
        end
    end

    for dim in 1:3
        __extend_ghosts!(trunc_pos, id, dim, periodic, (xmin, ymin, zmin), rneigh, (lx, ly, lz))
    end
    return trunc_pos, id
end

function __extend_ghosts!(coord, id, dim, pbc, origin, rneigh, boxdim)
    !pbc[dim] && return coord, id
    ncoord = length(coord)
    dr = boxdim .* ((1, 2, 3) .== dim)
    @inbounds for k in 1:ncoord
        r = coord[k]
        xoff = r[dim] - origin[dim]
        if xoff < rneigh
            push!(coord, r .+ dr)
            push!(id, id[k])
        end
        if boxdim[dim] - xoff < rneigh
            push!(coord, r .- dr)
            push!(id, id[k])
        end
    end
    return coord, id
end

"""
    parallel_voronoi_tessellation(pos; bounds, periodic=(false, false, false), ntasks)

Allocate space for a container of Voronoi cells and add coordinates stored in `pos`.
# Arguments
* `pos`: vector of positions
# Keywords
* `bounds`: limits of the bounding box `((xmin, ymin, zmin), (xmax, ymax, zmax))`
* `periodic::NTuple{3,Bool}`: periodicity in each axis. Default: `(false, false, false)`
* `ntasks`: number of tasks for parallel computation.
"""
function parallel_voronoi_tessellation(
    pos::AbstractVector
    ;
    bounds,
    periodic=(false, false, false),
    ntasks::Integer=2,
)

    eachindex(bounds) == OneTo(2) &&
    eachindex(bounds[1]) == eachindex(bounds[2]) == OneTo(3) ||
    throw(DimensionMismatch("Bounds must be a 2-tuple of length-3 arrays or tuples"))

    eachindex(periodic) == OneTo(3) ||
    throw(DimensionMismatch("Periodic flags must be a length-3 boolean array or tuple"))

    ((ax, ay, az), (bx, by, bz)) = bounds
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.((ax, bx, ay, by, az, bz))
    px, py, pz = periodic

    # container length in earch coordinate
    lx, ly, lz = Float64.((x_max - x_min, y_max - y_min, z_max - z_min))

    # container volume
    vol = lx * ly * lz

    # number of paticles
    npts = length(pos)

    # mean interparticle distance
    r_ave = cbrt(vol / npts)
    d_skin = 3 * r_ave

    containers = generate_containers(
        x_min, x_max, y_min, y_max, z_min, z_max, periodic, d_skin, npts, Int32(ntasks)
    )

    tessellation = ParallelTessellation(
        x_min, x_max, y_min, y_max, z_min, z_max,
        containers,
        d_skin,
        (px, py, pz)
    )

    nx, ny, nz = size(containers)
    dx, dy, dz = (lx, ly, lz) ./ (nx, ny, nz)

    # sort particles into bins

    #vector of vectors(number of task)
    workload = [Int32[] for _ in CartesianIndices(containers)]
    trunc_pos, id = __prepare_coords(pos, bounds, periodic, d_skin)
    for ind in eachindex(trunc_pos)
        x, y, z = trunc_pos[ind]
        (ix::Int, xoff), (iy::Int, yoff), (iz::Int, zoff) =
            fldmod.((x, y, z) .- (x_min, y_min, z_min), (dx, dy, dz))
        ix += one(ix)
        iy += one(iy)
        iz += one(iz)
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
                    jz = iz + diz
                    jz in axes(containers, 3) || continue
                    zoffj = zoff - diz * dz
                    -d_skin < zoffj <= dz + d_skin || continue
                    push!(workload[jx, jy, jz], ind)
                end
            end
        end
    end


    @sync for i_con in 1:ntasks
        @spawn let
            work = workload[i_con]
            con = tessellation.domain[i_con].con
            for i_part in work
                __cxxwrap_put!(con, id[i_part], trunc_pos[i_part]...)
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
            cell = VoronoiCell()
            (; con, x_min, x_max, y_min, y_max, z_min, z_max) = vt.domain[con_id]
            itor = ContainerSubsetIterator(con)
            __cxxwrap_setup_box!(itor, nextfloat(x_min), x_max, nextfloat(y_min), y_max, nextfloat(z_min), z_max, true)
            if ( __cxxwrap_start!(itor) )
                while true
                    pinfo = __cxxwrap_particle_info(itor)
                    particle = Particle(con, pinfo)
                    __cxxwrap_compute_cell!(cell, con, itor) && fn((particle, cell))
                    __cxxwrap_inc!(itor) || break
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
            cell = VoronoiCell()
            (; con, x_min, x_max, y_min, y_max, z_min, z_max) = vt.domain[con_id]
            itor = ContainerSubsetIterator(con)
            __cxxwrap_setup_box!(itor, nextfloat(x_min), x_max, nextfloat(y_min), y_max, nextfloat(z_min), z_max, true)
            if __cxxwrap_start!(itor)
                while true
                    pinfo = __cxxwrap_particle_info(itor)
                    particle = Particle(con, pinfo)
                    if __cxxwrap_compute_cell!(cell, con, itor)
                        val_local = op(val_local, fn((particle, cell)))
                    end
                    __cxxwrap_inc!(itor) || break
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

function Base.map!(fn, dest, vt::ParallelTessellation)
    @sync for con_id in LinearIndices(vt.domain)
        @spawn let
            cell = VoronoiCell()
            (; con, x_min, x_max, y_min, y_max, z_min, z_max) = vt.domain[con_id]
            itor = ContainerSubsetIterator(con)
            __cxxwrap_setup_box!(itor, nextfloat(x_min), x_max, nextfloat(y_min), y_max, nextfloat(z_min), z_max, true)
            if __cxxwrap_start!(itor)
                while true
                    pinfo = __cxxwrap_particle_info(itor)
                    particle = Particle(con, pinfo)
                    id = particle.id
                    if __cxxwrap_compute_cell!(cell, con, itor)
                        dest[id] = fn((particle, cell))
                    end
                    __cxxwrap_inc!(itor) || break
                end
            end
        end
    end
    return dest
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
            cell = VoronoiCell()
            (; con, x_min, x_max, y_min, y_max, z_min, z_max) = vt.domain[con_id]
            itor = ContainerSubsetIterator(con)
            __cxxwrap_setup_box!(itor, nextfloat(x_min), x_max, nextfloat(y_min), y_max, nextfloat(z_min), z_max, true)
            if __cxxwrap_start!(itor)
                while true
                    pinfo = __cxxwrap_particle_info(itor)
                    particle = Particle(con, pinfo)
                    if __cxxwrap_compute_cell!(cell, con, itor)
                        val_local += fn((particle, cell))
                    end
                    __cxxwrap_inc!(itor) || break
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

function voronoi_volume(vt::ParallelTessellation)
    return sum(vt; init=0.0) do (part, cell)
        volume(cell)
    end
end
