"""
    AbstractContainer

Supertype for wrappers of C++ container types.
"""
abstract type AbstractContainer end

"""
    AbstractTessellation

Abstract type for Voronoi tessellations
"""
abstract type AbstractTessellation{T} end

"""
    ContainerIterationOrder

Abstract type to specify whether iteration proceeds in a specific order.
"""
abstract type ContainerIterationOrder end

"""
    AbstractVoronoiCell

Abstract type for storing Voronoi cell data.
"""
abstract type AbstractVoronoiCell end

"""
    AbstractWall

Abstract type for wall classes.
"""
abstract type AbstractWall end

"""
    RawContainer

Wrapper for Voro++ `container` type.
"""
RawContainer

"""
    RawContainerPoly

Wrapper for Voro++ `container_poly` type.
"""
RawContainerPoly

"""
    InsertionOrder()

A data structure to store the order of inserted particles into a Voronoi container. Wraps
    `voro::particle_order` type.
"""
InsertionOrder

struct CxxBoxBounds
    ax::Float64
    ay::Float64
    az::Float64
    bx::Float64
    by::Float64
    bz::Float64
end

struct CxxPBC
    px::Bool
    py::Bool
    pz::Bool
end

struct CxxParticleInfo
    x::Float64
    y::Float64
    z::Float64
    r::Float64
    id::Int32
end

struct CxxVec3D
    x::Float64
    y::Float64
    z::Float64
end

struct CxxLoopIndices
    i::Int32
    j::Int32
    k::Int32
    ijk::Int32
    q::Int32
end
