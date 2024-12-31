

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

# Particle Coordinates
struct ParticleCoords
    # x, y, z coordinates for each particle
    xp::Real
    yp::Real
    zp::Real
    # particle id
    idp::Integer
end

# Voronoi Tesselation
struct VoronoiTessellation
    domain::Vector{Container} #vector of containers
    skin_distance::Real
    split_dim::Integer #number of threads
    split_bounds::Vector{ContainerDim} #vector of container dimensions
end

#
# SetContainerCoords!
#
# Parameters:
# offset: min point for each coordinate
# size_xyz: size for each coordinate
# nblocks_xyz: number of block division for each dimension
# nthr: number of threads
#
# Output: A ContainerDim vector of nthr elements
# Each container represents a cube volume for particles, and only the x coordinate will change, y and z remains with the same size
#
#=function SetContainerCoords!(offset::Float64 = Float64(-1), size_xyz::Float64 = Float64(1), nblocks_xyz::Int32 = Int32(6), nthr::Int32 = Int32(1))
    
    coordinates = Vector{ContainerDim}(undef, nthr)
    x_c = offset
    for i in 1:nthr
        coordinates[i] = ContainerDim(x_c, x_c + size_xyz, offset, offset + size_xyz, offset, offset + size_xyz, nblocks_xyz, nblocks_xyz, nblocks_xyz)
        x_c += size_xyz
    end

    return coordinates
end=#


##################################################################################

function GenerateContainerDims(
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
# n_particles: number of random particles
# offset: min point for each coordinate, as containers
# size_xyz: size for each coordinate, as cotainers
# nthr: number of threads
#
# Output: A ParticleCoords vector of nthr elements
#
function GenerateParticles!(nparticles::Integer, x_min::Real, x_max::Real, y_min::Real, y_max::Real, z_min::Real, z_max::Real)

    particles = Vector{ParticleCoords}(undef, nparticles)

    for i in 1:nparticles
        x = x_min + rand() * (x_max - x_min)
        y = y_min + rand() * (y_max - y_min)
        z = z_min + rand() * (z_max - z_min)
        particles[i] = ParticleCoords( x, y, z, i-1 )
    end

    return particles

end

#
# ComputeTessellation!
#
# Parameters:
# n_particles: number of random particles, parameter for function GenerateParticles!
# offset: min point for each coordinate, parameter for SetContainerCoords! and for GenerateParticles!
# size_xyz: size for each coordinate, used for Cotainers and Particles
# nblocks_xyz: number of block division for each dimension, parameter for SetContainerCoords!
# skin_d: skin distance
#
# Output: Files with particles and cells assigned to each container
#

function ComputeTessellation!(
    n_particles::Integer,  
    x_min::Real, 
    x_max::Real, 
    y_min::Real, 
    y_max::Real, 
    z_min::Real, 
    z_max::Real, 
    d_skin::Real, 
    nblocks::Integer
    )

    # number of threads
    nthr = Threads.nthreads()
    # coordinates vector
    con_coords = GenerateContainerDims(x_min, x_max, y_min, y_max, z_min, z_max, d_skin, nblocks, nthr)
    # containers vector
    containers = GenerateContainers!( con_coords )

    # tessellation object
    tessellation = VoronoiTessellation(
        containers,
        d_skin,
        nthr,
        con_coords
    )

    # particles vector
    v_particles = GenerateParticles!( n_particles, x_min, x_max, y_min, y_max, z_min, z_max )

    # parallel execution
    Threads.@threads for _ in 1:Threads.nthreads()

        # Iteration through parcticle vector
        # each thread checks particles from vector
        # if particle's x coordinte is inside container's dimension, assigned to each thread, add this particle to corresponding container
        for i in v_particles
            if ( i.xp > tessellation.split_bounds[Threads.threadid()].x_min ) && ( i.xp < tessellation.split_bounds[Threads.threadid()].x_max )
                add_point!(tessellation.domain[Threads.threadid()], i.idp, i.xp, i.yp, i.zp)
            end
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
    
end

##################################################################################


# Original Design

#= struct VoronoiTessellation
    domain::Vector{Container}
    skin_distance::Float64
    split_dim::Int8
    split_bounds::Vector{Float64}
end =#

#= function ComputeTessellation!(
    vt::VoronoiTessellation,
    #coords::AbstractMatrix{<:Real}
    coords::Vector{ContainerD}
)
end =#

#Base.iterate(vt::VoronoiTessellation)


#################################################################################


# Final Version

ComputeTessellation!(40, 0.0, 8.0, 0.0, 2.0, 0.0, 2.0, 0.01, 6)
