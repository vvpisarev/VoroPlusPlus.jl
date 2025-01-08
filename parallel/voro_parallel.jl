

# Data structure for computing the parallel version of voro++ in Julia

using Base.Threads
using BenchmarkTools

# Container Dimension
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
    skin_distance::Real
    split_dim::Integer #number of threads
    split_bounds::Vector{ContainerDim} #vector of container dimensions
    owner::Vector{Int32} # vector of "owning" domains for each particle
end

#
# gen_container_dims
#
# Parameters:
#
# Output: A ContainerDim vector of nthr elements
# Each container represents a cube volume for particles, and only the x coordinate will change, y and z remains with the same size
#

function gen_container_dims(
    x_min::Real,  # global bounds
    x_max::Real,
    y_min::Real,
    y_max::Real,
    z_min::Real,
    z_max::Real,
    d_skin::Real,
    nblocks_xyz::Integer = Int32(6),
    ntasks::Integer = 1,
)
    
    coordinates = Vector{ContainerDim}(undef, ntasks)

    x_range = range(Float64(x_min); stop=Float64(x_max), length=ntasks + 1)
    for i in 1:ntasks
        x_lo = i == 1 ? first(x_range) : x_range[i] - d_skin
        x_hi = i == ntasks ? last(x_range) : x_range[i+1] + d_skin
        coordinates[i] = ContainerDim(x_lo, x_hi, y_min, y_max, z_min, z_max, nblocks_xyz, nblocks_xyz, nblocks_xyz)
    end

    return coordinates
end


##################################################################################

#
# GenerateContainers!
#
# Parameters:
# v_coords: vector with the coordinates for each container
# nthr: number of threads
#
# Output: A Container vector of nthr elements, each container is empty
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

    particles = Array{Real}([])

    for _ in 1:nparticles
        x = x_min + rand() * (x_max - x_min)
        y = y_min + rand() * (y_max - y_min)
        z = z_min + rand() * (z_max - z_min)
        push!(particles, x)
        push!(particles, y)
        push!(particles, z)
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
    
    con_dims = gen_container_dims(x_min, x_max, y_min, y_max, z_min, z_max, d_skin, 6, ntasks)
    containers = GenerateContainers!(con_dims)
    owner = zeros(Int32, npts) # task that owns each particle
    
    tessellation = VoronoiTessellation(
        x_min, x_max, y_min, y_max, z_min, z_max,
        containers,
        d_skin,
        1,
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
    for ind in firstindex(coords):3:lastindex(coords)
        p += 1
        x, y, z = coords[ind], coords[ind+1], coords[ind+2]
        i_con::Int32, x_rel = fldmod1(x - x_min, dx)
        
        owner[p] = i_con
        if i_con > 1 && x_rel < d_skin
            push!(workload[i_con-1], p)
        elseif i_con < ntasks && dx - x_rel < d_skin
            push!(workload[i_con+1], p)        
        end
        push!(workload[i_con], p)
        
    end
    
    
    Threads.@threads for _ in 1:ntasks
        work = workload[Threads.threadid()]
        con = tessellation.domain[Threads.threadid()]
        for i_part in work
            ind = (i_part - 1) * 3
            x, y, z = coords[ind+1], coords[ind+2], coords[ind+3]
            add_point!(con, i_part, x, y, z)
        end

        # Simple Test for Volumes
        cvol = (
            (tessellation.split_bounds[Threads.threadid()].x_max - tessellation.split_bounds[Threads.threadid()].x_min)*
            (tessellation.split_bounds[Threads.threadid()].y_max - tessellation.split_bounds[Threads.threadid()].y_min)*
            (tessellation.split_bounds[Threads.threadid()].z_max - tessellation.split_bounds[Threads.threadid()].z_min)
        )
        vvol = sum(volume, tessellation.domain[Threads.threadid()])
        message = "Container volume thread $(Threads.threadid()): $(cvol)\n"*
                    "Voronoi volume thread $(Threads.threadid()): $(vvol)\n"*
                    "Difference thread $(Threads.threadid()): $(cvol - vvol)\n"*
                    "Is approx thread $(Threads.threadid())?: $(isapprox(vvol, cvol; atol=1e-8))\n" 
        
        println(message)

        # After particle assigment, each thread calculates tesselation and generates files for particles and cells

        draw_particles(tessellation.domain[Threads.threadid()], "parallel/particles_$(Threads.threadid()).gnu")
        draw_cells_gnuplot(tessellation.domain[Threads.threadid()], "parallel/voro_cells_$(Threads.threadid()).gnu")

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

    # number of threads / tasks
    ntasks = Threads.nthreads()
    # particle's array
    a_p = GenerateParticles!(nparticles, x_min, x_max, x_min, y_max, z_min, z_max)
    a_p_range = 1:div(length(a_p), 3)
   
    t = voronoi_tessellation((x_min, x_max, y_min, y_max, z_min, z_max), a_p, a_p_range , ntasks)
    println(typeof(t))

end


####################################################################################################

TestVoroTesselation(1000, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0)


