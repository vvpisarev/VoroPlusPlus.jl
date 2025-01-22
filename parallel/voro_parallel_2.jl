
########################################################################
# Data structure for computing the parallel version of voro++ in Julia #
########################################################################

using Base.Threads
using BenchmarkTools

# Container Dimensions
# Needed for Container constructor
struct ContainerDim
    x_min::Real # minimum and maximum dimensions for each coordinate for container
    x_max::Real
    y_min::Real
    y_max::Real
    z_min::Real
    z_max::Real
    n_x::Integer # number of blocks for each dimension
    n_y::Integer
    n_z::Integer
end

# Voronoi Tesselation
struct VoronoiTessellation
    x_min::Real # minimum and maximum dimensions for each coordinate
    x_max::Real
    y_min::Real
    y_max::Real
    z_min::Real
    z_max::Real
    domain::Vector{Container} # vector of containers
    skin_distance::Float64 # overlapping distance between neighboring containers
    split_dim::Int32 # number of tasks
    split_bounds::Vector{ContainerDim} # vector of container dimensions
    owner::Vector{Int32} # vector owner of each particle ID
end

#
# SetGeneralParameters!
#
function SetGeneralParameters(
    nparticles::Integer, 
    d_skin::Real, 
    x_min::Real, 
    x_max::Real, 
    y_min::Real, 
    y_max::Real, 
    z_min::Real, 
    z_max::Real, 
    ntasks::Integer
)

    return Dict(
        "particles" => nparticles, 
        "skin" => d_skin, 
        "x_min" => x_min,
        "x_max" => x_max,
        "y_min" => y_min,
        "y_max" => y_max,
        "z_min" => z_min,
        "z_max" => z_max,
        "tasks" => ntasks
    )
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
# Time Complexity: O(1)
# Constant time depends on the number of tasks, a small number
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
# Time Complexity: O(1)
# Constant time depends on the number of tasks, a small number
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
# Time Complexity: O(n)
# Constant time depends on the number of particles n
#
function GenerateParticles!(
    nparticles::Int32, 
    x_min::Float64, 
    x_max::Float64, 
    y_min::Float64, 
    y_max::Float64, 
    z_min::Float64, 
    z_max::Float64)

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
# BaseTesselation
#
# Output: An Empty Tessellation for Computation
#
function BaseTessellation(
    x_min, 
    x_max, 
    y_min, 
    y_max, 
    z_min, 
    z_max,
    npts,
    con_dims,
    containers,
    owner,
    ntasks)
    
    # container length in each coordinate
    lx, ly, lz = Float64.((x_max - x_min, y_max - y_min, z_max - z_min))
    
    # container volume
    vol = lx * ly * lz
    
    # skin distance
    r_ave = cbrt(vol / npts)
    d_skin = 5 * r_ave
    
    tessellation = VoronoiTessellation(
        x_min, x_max, y_min, y_max, z_min, z_max,
        containers,
        d_skin,
        ntasks,
        con_dims,
        owner
    )

    return tessellation
end
    
#
# ComputeTessellation!
#
# Parameters: 
#            tessellation: Empty Tessallation
#            coords: Particle's Vector
#
# Output: Computed Tessellation
#
function ComputeTessellation!(tessellation::VoronoiTessellation, coords::AbstractArray{<:Real})

    # x dimension if there are more than one task
    dx = (tessellation.x_max - tessellation.x_min) / tessellation.split_dim
    
    # particle
    p = 0

    # sort particles into bins

    #vector of vectors(number of task)
    workload = [Int32[] for _ in 1:tessellation.split_dim]
    for ind in firstindex(coords):3:lastindex(coords)
        p += 1
        x, y, z = coords[ind], coords[ind+1], coords[ind+2]
        i_con::Int32, x_rel = fldmod1(x - tessellation.x_min, dx)
        
        owner[p] = i_con
        if i_con > 1 && x_rel < tessellation.skin_distance
            push!(workload[i_con-1], p)
        elseif i_con < tessellation.split_dim && dx - x_rel < tessellation.skin_distance
            push!(workload[i_con+1], p)
        end
        push!(workload[i_con], p)
        
    end
    
    
    Threads.@threads for i_con in 1:tessellation.split_dim
        work = workload[i_con]
        con = tessellation.domain[i_con]
        for i_part in work
            ind = (i_part - 1) * 3
            x, y, z = coords[ind+1], coords[ind+2], coords[ind+3]
            add_point!(con, i_part, x, y, z)
        end

        # After particle assigment, each thread calculates tesselation and generates files for particles and cells
        draw_particles(tessellation.domain[i_con], "parallel/particles_$(i_con).gnu")
        draw_cells_gnuplot(tessellation.domain[i_con], "parallel/voro_cells_$(i_con).gnu")

    end
    
end


#
# Test Volume
#
function TestVolume(tessellation::VoronoiTessellation)

    # Simple Test for Volumes
    cvol = (
        (tessellation.x_max - tessellation.x_min)*
        (tessellation.y_max - tessellation.y_min)*
        (tessellation.z_max - tessellation.z_min)
    )

    partial_vol = [Float64[] for _ in 1:length(tessellation.domain)]
    
    c = VoronoiCell()

    pid = Ref(Int32(0))
	x = Ref(0.0)
	y = Ref(0.0)
	z = Ref(0.0)
	r = Ref(0.0)

    for (con_id, con) in enumerate(tessellation.domain)

        cla = Container_Iterator(con)

        if ( start!(cla) )
            if ( compute_cell!(c, con, cla) )
                pos(cla, pid, x, y, z, r);
                if tessellation.owner[pid[]] == con_id
                    push!(partial_vol[con_id], volume(c))
                end
            end
            while ( next!(cla) )
                if ( compute_cell!(c, con, cla) )
                    pos(cla, pid, x, y, z, r);
                    if tessellation.owner[pid[]] == con_id
                        push!(partial_vol[con_id], volume(c))
                    end
                end
            end
        end
    end

    vvol = 0
    for v in partial_vol
        vvol += sum(v)
    end

    #vvol = sum(map(sum, partial_vol))


    message = "Container volume: $(cvol)\n"*
                "Voronoi volume: $(vvol)\n"*
                "Difference: $(cvol - vvol)\n"*
                "Is approx?: $(isapprox(vvol, cvol; atol=1e-8))\n" 
    
    println(message)

end


############################################################


############################################################


settings = SetGeneralParameters(500000, 0.5, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0, Threads.nthreads())

containerdimensions = GenerateContainerDims(
    settings["x_min"],
    settings["x_max"],
    settings["y_min"],
    settings["y_max"],
    settings["z_min"],
    settings["z_max"],
    settings["skin"],
    Int32(6),
    Int32(settings["tasks"])
)

containers = GenerateContainers!(containerdimensions)

particles = GenerateParticles!(
    Int32(settings["particles"]),
    settings["x_min"],
    settings["x_max"],
    settings["y_min"],
    settings["y_max"],
    settings["z_min"],
    settings["z_max"]
)

len_particles = 1:div(length(particles), 3)

npts = length(len_particles)

owner = zeros(Int32, npts)

tessellation = BaseTessellation(
    settings["x_min"], 
    settings["x_max"], 
    settings["y_min"], 
    settings["y_max"], 
    settings["z_min"], 
    settings["z_max"],
    npts,
    containerdimensions,
    containers,
    owner,
    settings["tasks"]
)

@elapsed ComputeTessellation!(tessellation, particles)

#TestVolume(tessellation)