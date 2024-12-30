

# Data structure for computing the parallel version of voro++ in Julia

using Base.Threads
using BenchmarkTools

# Container Dimension
struct ContainerDim
    # minimum and maximum dimensions for each coordinate
    x_min::Float64 
    x_max::Float64
    y_min::Float64
    y_max::Float64
    z_min::Float64
    z_max::Float64
    # number of blocks for each dimension
    n_x::Int32
    n_y::Int32
    n_z::Int32
end

# Particle Coordinates
struct ParticleCoords
    # x, y, z coordinates for each particle
    xp::Float64
    yp::Float64
    zp::Float64
    # particle id
    idp::Int32
end

# Voronoi Tesselation
struct VoronoiTessellation
    domain::Vector{Container} #vector of containers
    skin_distance::Float64
    split_dim::Int32 #number of threads
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
function SetContainerCoords!(offset::Float64 = Float64(-1), size_xyz::Float64 = Float64(1), nblocks_xyz::Int32 = Int32(6), nthr::Int32 = Int32(1))
    
    coordinates = Vector{ContainerDim}(undef, nthr)
    x_c = offset
    for i in 1:nthr
        coordinates[i] = ContainerDim(x_c, x_c + size_xyz, offset, offset + size_xyz, offset, offset + size_xyz, nblocks_xyz, nblocks_xyz, nblocks_xyz)
        x_c += size_xyz
    end

    return coordinates
end

#
# GenerateContainers!
#
# Parameters:
# v_coords: vector with the coordinates for each container
# nthr: number of threads
#
# Output: A Container vector of nthr elements, each container is empty
#
function GenerateContainers!(v_coords::Vector{ContainerDim}, nthr)

    containers = Vector{Container}(undef, nthr)

    for i in 1:nthr
        containers[i] = Container(
        ;
        bounds = (v_coords[i].x_min , v_coords[i].x_max, v_coords[i].y_min, v_coords[i].y_max, v_coords[i].z_min, v_coords[i].z_max),
        nblocks = (v_coords[i].n_x, v_coords[i].n_y, v_coords[i].n_z),
        periodic = (false, false, false),
        particles_per_block = 8,
        )
    end

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
function GenerateParticles!(nparticles::Int32, offset::Float64, size_xyz::Float64, nthr::Int32)

    particles = Vector{ParticleCoords}(undef, nparticles)

    for i in 1:nparticles
        x = offset + rand() * (offset + size_xyz) * nthr
        y = offset + rand() * (offset + size_xyz)
        z = offset + rand() * (offset + size_xyz)
        particles[i] = ParticleCoords(x, y, z, Int32(i-1))
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

function ComputeTessellation!(n_particles::Int32, offset::Float64, size_xyz::Float64, nblocks::Int32, skin_d::Float64)

    # number of threads
    nthr = Threads.nthreads()
    # coordinates vector
    con_coords = SetContainerCoords!(offset, size_xyz, nblocks, Int32(nthr))
    # containers vector
    containers = GenerateContainers!( con_coords, nthr)

    # tessellation object
    tessellation = VoronoiTessellation(
        containers,
        skin_d,
        nthr,
        con_coords
    )

    # particles vector
    v_particles = GenerateParticles!( n_particles, offset, size_xyz, Int32(nthr) )

    # parallel execution
    Threads.@threads for _ in 1:Threads.nthreads()

        # iteration through parcticle vector
        # each thread checks particles from vector
        # if particle's x coordinte is inside container's dimension, assigned to each thread, add this particle to corresponding container
        for i in v_particles
            if ( i.xp >= tessellation.split_bounds[Threads.threadid()].x_min ) && ( i.xp <= tessellation.split_bounds[Threads.threadid()].x_max )
                add_point!(tessellation.domain[Threads.threadid()], i.idp, i.xp, i.yp, i.zp)
            end
        end

        # after particle assigment, each thread calculates tesselation and generates files for particles and cells
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


# Tests

#SetCoords!()
#SetContainerCoords!(Float64(0), Float64(2), Float64(6), Int32(3))
#GenerateContainers!( SetContainerCoords!(Float64(0), Float64(2), Int32(6), Int32(3)) , Int32(3))
#GenerateParticles!( Int32(5), 0.0, 2.0, Int32(3) )


#################################################################################


# Final Version

# Parameters: 
# ComputeTessellation!(number of particles, dimension offset for each coordinate, size for each coordinate, number of block for each coordinate)

ComputeTessellation!(Int32(80), 0.0, 2.0, Int32(6), 1.0)
#@btime ComputeTessellation!(Int32(80), 0.0, 2.0, Int32(6), 1.0)
#@benchmark ComputeTessellation!(Int32(80), 0.0, 2.0, Int32(6), 1.0)

