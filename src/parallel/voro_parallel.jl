
# Data structure for computing the parallel version of voro++ in Julia

using Base.Threads

# Container Dimensions
# Needed for Container constructor
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

struct ExtendedContainer{T<:Container}
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
struct VoronoiTessellation{T<:ExtendedContainer}
    # minimum and maximum dimensions for each coordinate
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    z_min::Float64
    z_max::Float64
    domain::Vector{T} #vector of containers
    skin_distance::Float64
    split_dim::Int32 #number of dimension in which box is split
    owner::Vector{Int32} # vector owner of each particle ID
end

#
# GenenerateContainers
#
# Generate Individual Container Dimensions
#
# Parameters: x_min, x_max, y_min, y_max, z_min, z_max
#             d_skin
#             nblocks_xyz
#
# Output: A Container Dimensions vector (cont_dim)
# Each container represents a cube volume for particles, and only the x coordinate will change, y and z remain with the same size, +/- d_skin
#

function GenerateContainers(
    x_min::Float64,  # global bounds
    x_max::Float64,
    y_min::Float64,
    y_max::Float64,
    z_min::Float64,
    z_max::Float64,
    d_skin::Float64, # d_skin size
    nparticles::Integer, # total no. of particles
    ntasks::Int32 = Int32(1), # default number of tasks
)
    lx, ly, lz = (x_max, y_max, z_max) .- (x_min, y_min, z_min)
    split_dim = argmax((lx, ly, lz))
    dist_range = range(
        (x_min, y_min, z_min)[split_dim]
        ;
        stop=(x_max, y_max, z_max)[split_dim],
        length=ntasks + 1
    )
    ilscale = cbrt(nparticles / (OPT_PART_PER_BLOCK * lx * ly * lz))
    containers = map(1:ntasks) do i
        x_lo, y_lo, z_lo = x_in_lo, y_in_lo, z_in_lo = x_min, y_min, z_min
        x_hi, y_hi, z_hi = x_in_hi, y_in_hi, z_in_hi = x_max, y_max, z_max
        if split_dim == 1
            x_lo = i == 1 ? first(dist_range) : dist_range[i] - d_skin
            x_hi = i + 1 > ntasks ? last(dist_range) : dist_range[i+1] + d_skin
            x_in_lo, x_in_hi = dist_range[i], dist_range[i+1]
        elseif split_dim == 2
            y_lo = i == 1 ? first(dist_range) : dist_range[i] - d_skin
            y_hi = i > ntasks ? last(dist_range) : dist_range[i+1] + d_skin
            y_in_lo, y_in_hi = dist_range[i], dist_range[i+1]
        elseif split_dim == 3
            z_lo = i == 1 ? first(dist_range) : dist_range[i] - d_skin
            z_hi = i > ntasks ? last(dist_range) : dist_range[i+1] + d_skin
            z_in_lo, z_in_hi = dist_range[i], dist_range[i+1]
        end
        nx, ny, nz = floor.(
            Int32,
            (x_hi - x_lo, y_hi - y_lo, z_hi - z_lo) .* (ilscale + 1),
        )
        con = Container(
            ;
            bounds = (x_lo, x_hi, y_lo, y_hi, z_lo, z_hi),
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

#
# voronoi_tessellation
#
function parallel_container(con_dims, coords::AbstractArray{<:Real}, ntasks::Integer=2)

    #container dimensions
    x_min, x_max, y_min, y_max, z_min, z_max = con_dims

    # container length in earch coordinate
    lx, ly, lz = Float64.((x_max - x_min, y_max - y_min, z_max - z_min))

    # container volume
    vol = lx * ly * lz

    # number of paticles
    npts = div(length(coords), 3)

    # ?_average
    r_ave = cbrt(vol / npts)
    d_skin = 5 * r_ave

    containers = GenerateContainers(x_min, x_max, y_min, y_max, z_min, z_max, d_skin, npts, Int32(ntasks))
    owner = zeros(Int32, npts) # task that owns each particle
    split_dim = Int32(argmax((lx, ly, lz)))

    tessellation = VoronoiTessellation(
        x_min, x_max, y_min, y_max, z_min, z_max,
        containers,
        d_skin,
        split_dim,
        owner
    )

    # x dimension is there are more than one task
    dr = (lx, ly, lz)[split_dim] / ntasks
    min_coord = (x_min, y_min, z_min)[split_dim]

    # particle
    p = 0

    # sort particles into bins

    #vector of vectors(number of task)
    workload = [Int32[] for _ in 1:ntasks]
    #ghosts = [Int32[] for _ in 1:ntasks]
    for ind in firstindex(coords):3:lastindex(coords)
        p += 1
        x, y, z = coords[ind], coords[ind+1], coords[ind+2]
        coord = (x, y, z)[split_dim]
        i_con::Int32, offset = fldmod1(coord - min_coord, dr)
        if i_con == 0
            @info "" p coord
        end

        owner[p] = i_con
        if i_con > 1 && offset < d_skin
            #push!(ghosts[i_con-1], p)
            push!(workload[i_con-1], p)
        elseif i_con < ntasks && dr - offset < d_skin
            #push!(ghosts[i_con+1], p)
            push!(workload[i_con+1], p)
        end
        push!(workload[i_con], p)

    end


    Threads.@threads for i_con in 1:ntasks
        work = workload[i_con]
        con = tessellation.domain[i_con].con
        for i_part in work
            ind = (i_part - 1) * 3
            x, y, z = coords[ind+1], coords[ind+2], coords[ind+3]
            add_point!(con, i_part, x, y, z)
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
function Base.foreach(fn, vt::VoronoiTessellation)
    @sync for con_id in eachindex(vt.domain)
        @spawn let
            pid = Ref(Int32(0))
            x = Ref(0.0)
            y = Ref(0.0)
            z = Ref(0.0)
            r = Ref(0.0)
            c = VoronoiCell()
            con = vt.domain[con_id].con
            cla = Container_Iterator(con)
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
    mapreduce(fn, op, vt::VoronoiTessellation; init)

Apply function `fn` to each Voronoi cell in `vt`, and then reduce the result using the
    binary function `op`. `op` is assumed to be commutative and transitive, and the order of
    application of `fn` must not have effect on the result. `init` must be neutral element
    for `op` and has to be provided as it defines the type of output.
"""
function Base.mapreduce(fn, op, vt::VoronoiTessellation; init)
    locals = Channel{typeof(init)}(length(vt.domain))
    @sync for con_id in eachindex(vt.domain)
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
    sum(fn, vt::VoronoiTessellation; init)

Sum the results of calling function `fn` on each cell in `vt`. The order of application of
    `fn` must not have effect on the result. `init` must be zero and has to be provided as
    it defines the type of output.
"""
function Base.sum(fn, vt::VoronoiTessellation; init)
    locals = Channel{typeof(init)}(length(vt.domain))
    @sync for con_id in eachindex(vt.domain)
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
