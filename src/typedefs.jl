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
