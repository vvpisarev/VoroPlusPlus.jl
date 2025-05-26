"""
    AbstractContainer

Abstract type for Voronoi tessellation containers.
"""
abstract type AbstractContainer end

"""
    AbstractRawContainer

Supertype for direct wrappers of C++ types.
"""
abstract type AbstractRawContainer end

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
    RawContainer

Wrapper for Voro++ `container` type.
"""
RawContainer

"""
    RawContainerPoly

Wrapper for Voro++ `container_poly` type.
"""
RawContainerPoly
