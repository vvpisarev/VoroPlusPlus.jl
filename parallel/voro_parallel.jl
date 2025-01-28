
# Data structure for computing the parallel version of voro++ in Julia

using Base.Threads
using BenchmarkTools

# Container Dimensions
# Needed for Container constructor
struct ContainerDim
    # minimum and maximum dimensions for each coordinate
    x_min::Real
    x_max::Real
    y_min::Real
    y_max::Real
    z_min::Real
    z_max::Real
    # number of blocks for each dimension
    n_x::Integer
    n_y::Integer
    n_z::Integer
end


# Voronoi Tesselation
struct VoronoiTessellation
    # minimum and maximum dimensions for each coordinate
    x_min::Real
    x_max::Real
    y_min::Real
    y_max::Real
    z_min::Real
    z_max::Real
    domain::Vector{Container} #vector of containers
    skin_distance::Float64
    split_dim::Int32 #number of tasks
    split_bounds::Vector{ContainerDim} #vector of container dimensions
    owner::Vector{Int32} # vector owner of each particle ID
end

#
# GenenerateContainerDims
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

function GenerateContainerDims(
    x_min::Float64,  # global bounds
    x_max::Float64,
    y_min::Float64,
    y_max::Float64,
    z_min::Float64,
    z_max::Float64,
    d_skin::Float64, # d_skin size
    nblocks_xyz::Int32 = Int32(6), # n_blocks, same for x, y and z
    ntasks::Int32 = Int32(1), # default number of tasks
)
 
    # Container Dimensions vector (con_dims)
    con_dims = Vector{ContainerDim}(undef, ntasks)

    x_range = range(x_min; stop=x_max, length=ntasks + 1)
    for i in 1:ntasks
        x_lo = i == 1 ? first(x_range) : x_range[i] - d_skin
        x_hi = i == ntasks ? last(x_range) : x_range[i+1] + d_skin
        con_dims[i] = ContainerDim(x_lo, x_hi, y_min, y_max, z_min, z_max, nblocks_xyz, nblocks_xyz, nblocks_xyz)
    end

    return con_dims
end


#
# GenerateContainers!
#
# Parameters:
#           v_coords: vector with the dimensions for each container
#
# Output: A Container vector of ntasks elements, each container is empty
#
function GenerateContainers!(v_coords::Vector{ContainerDim})

    containers = [
        Container(
            ;
            bounds = (con_dim.x_min , con_dim.x_max, con_dim.y_min, con_dim.y_max, con_dim.z_min, con_dim.z_max),
            nblocks = (con_dim.n_x, con_dim.n_y, con_dim.n_z),
            periodic = (false, false, false),
            particles_per_block = 8,
            )
        for con_dim in v_coords
    ]

    return containers
end


#
# GenerateParticles!
#
# Parameters:
#
# Output: A ParticleCoords array of nparticles*3 elements
#
function GenerateParticles!(nparticles::Integer, x_min::Real, x_max::Real, y_min::Real, y_max::Real, z_min::Real, z_max::Real)

    particles = Array{Float64}([])

    for _ in 1:nparticles
        x = x_min + rand() * (x_max - x_min)
        y = y_min + rand() * (y_max - y_min)
        z = z_min + rand() * (z_max - z_min)
        push!(particles, x, y, z)
    end

    return particles

end


#
# voronoi_tessellation
#
function voronoi_tessellation(con_dims, coords::AbstractArray{<:Real}, ids=1:div(length(coords), 3), ntasks::Integer=2)

    #container dimensions
    x_min, x_max, y_min, y_max, z_min, z_max = con_dims
    
    # container length in earch coordinate
    lx, ly, lz = Float64.((x_max - x_min, y_max - y_min, z_max - z_min))
    
    # container volume
    vol = lx * ly * lz
    
    # number of paticles
    npts = length(ids)
    
    # ?_average
    r_ave = cbrt(vol / npts)
    d_skin = 5 * r_ave
    
    con_dims = GenerateContainerDims(x_min, x_max, y_min, y_max, z_min, z_max, d_skin, Int32(6), Int32(ntasks))
    containers = GenerateContainers!(con_dims)
    owner = zeros(Int32, npts) # task that owns each particle
    
    tessellation = VoronoiTessellation(
        x_min, x_max, y_min, y_max, z_min, z_max,
        containers,
        d_skin,
        ntasks,
        con_dims,
        owner
    )
    
    # x dimension is there are more than one task
    dx = lx / ntasks
    
    # particle
    p = 0

    # sort particles into bins

    #vector of vectors(number of task)
    workload = [Int32[] for _ in 1:ntasks]
    #ghosts = [Int32[] for _ in 1:ntasks]
    for ind in firstindex(coords):3:lastindex(coords)
        p += 1
        x, y, z = coords[ind], coords[ind+1], coords[ind+2]
        i_con::Int32, x_rel = fldmod1(x - x_min, dx)
        
        owner[p] = i_con
        if i_con > 1 && x_rel < d_skin
            #push!(ghosts[i_con-1], p)
            push!(workload[i_con-1], p)
        elseif i_con < ntasks && dx - x_rel < d_skin
            #push!(ghosts[i_con+1], p)
            push!(workload[i_con+1], p)
        end
        push!(workload[i_con], p)
        
    end
    
    
    Threads.@threads for i_con in 1:ntasks
        work = workload[i_con]
        con = tessellation.domain[i_con]
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
        draw_particles(tessellation.domain[i_con], "parallel/particles_$(i_con).gnu")
        draw_cells_gnuplot(tessellation.domain[i_con], "parallel/voro_cells_$(i_con).gnu")

    end
    
    return tessellation
end


#
#
function TestVoroTesselation(
    nparticles::Integer,
    x_min::Real, 
    x_max::Real, 
    y_min::Real, 
    y_max::Real, 
    z_min::Real, 
    z_max::Real, 
)

    # number tasks
    ntasks = Threads.nthreads()
    # particle's array
    a_p = GenerateParticles!(nparticles, x_min, x_max, x_min, y_max, z_min, z_max)
    a_p_range = 1:div(length(a_p), 3)
   
    tessellation = voronoi_tessellation((x_min, x_max, y_min, y_max, z_min, z_max), a_p, a_p_range , ntasks)
    #println(typeof(tessellation))

    # Simple Test for Volumes
    cvol = (
        (tessellation.x_max - tessellation.x_min)*
        (tessellation.y_max - tessellation.y_min)*
        (tessellation.z_max - tessellation.z_min)
    )

    vvol = GetVoroVolume(tessellation)

    message = "Container volume: $(cvol)\n"*
                "Voronoi volume: $(vvol)\n"*
                "Difference: $(cvol - vvol)\n"*
                "Is approx?: $(isapprox(vvol, cvol; atol=1e-8))\n" 
    
    println(message)

end


############################################################

function GetVoroVolume(vt::VoronoiTessellation)
    
    partial_vol = [Float64[] for _ in 1:length(vt.domain)]
    
    c = VoronoiCell()

    pid = Ref(Int32(0))
	x = Ref(0.0)
	y = Ref(0.0)
	z = Ref(0.0)
	r = Ref(0.0)

    for (con_id, con) in enumerate(vt.domain)

        cla = Container_Iterator(con)

        if ( start!(cla) )
            if ( compute_cell!(c, con, cla) )
                pos(cla, pid, x, y, z, r);
                if vt.owner[pid[]] == con_id
                    push!(partial_vol[con_id], volume(c))
                end
            end
            while ( next!(cla) )
                if ( compute_cell!(c, con, cla) )
                    pos(cla, pid, x, y, z, r);
                    if vt.owner[pid[]] == con_id
                        push!(partial_vol[con_id], volume(c))
                    end
                end
            end
        end
    end

    t_vol = 0
    for v in partial_vol
        t_vol += sum(v)
    end

    return t_vol

end


############################################################

#TestVoroTesselation(50000, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0)
@elapsed TestVoroTesselation(500000, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0)
#@benchmark TestVoroTesselation(80000, 0.0, 3.0, 0.0, 3.0, 0.0, 3.0)
