struct UnspecifiedOrder<:ContainerIterationOrder end

mutable struct Container{C<:AbstractRawContainer, O<:ContainerIterationOrder}<:AbstractContainer
    con::C
    ord::O
end

"""
    container(; bounds, nblocks, periodic=(false, false, false), particles_per_block=8, ordering)

Allocate space for a container of Voronoi cells.
# Keywords
* `bounds`: limits of the bounding box `((xmin, ymin, zmin), (xmax, ymax, zmax))`
* `nblocks`: numbers of computation blocks along each axis
* `periodic::NTuple{3,Bool}`: periodicity in each axis. Default: `(false, false, false)`
* `particles_per_block::Integer`: initially allocate memory for this many particles
    per block. Default: 8
* `ordering`: `UnspecifiedOrder()` or `InsertionOrder()`. Default: `UnspecifiedOrder()`.
"""
function container(
    ;
    bounds::Tuple{NTuple{3,Integer},NTuple{3,Integer}},
    nblocks::NTuple{3,Real},
    periodic::NTuple{3,Bool}=(false, false, false),
    particles_per_block::Integer=8,
    ordering::ContainerIterationOrder=UnspecifiedOrder(),
)
    ((ax, ay, az), (bx, by, bz)) = bounds
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.((ax, bx, ay, by, az, bz))
    n_x, n_y, n_z = Int32.(nblocks)
    px, py, pz = periodic
    ppb = Int32(particles_per_block)

    con = RawContainer(
        x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, px, py, pz, ppb,
    )
    return Container(con, ordering)
end

"""
    polydisperse_container(; bounds, nblocks, periodic=(false, false, false), particles_per_block=8, ordering)

Allocate space for a container of Voronoi cells for polydisperse particles.
# Keywords
* `bounds`: limits of the bounding box `(xmin, xmax, ymin, ymax, zmin, zmax)`
* `nblocks`: numbers of computation blocks along each axis
* `periodic::NTuple{3,Bool}`: periodicity in each axis. Default: `(false, false, false)`
* `particles_per_block::Integer`: initially allocate memory for this many particles
    per block. Default: 8
* `ordering`: `UnspecifiedOrder()` or `InsertionOrder()`. Default: `UnspecifiedOrder()`.
"""
function polydisperse_container(
    ;
    bounds::Tuple{NTuple{3,Integer},NTuple{3,Integer}},
    nblocks::NTuple{3,Real},
    periodic::NTuple{3,Bool}=(false, false, false),
    particles_per_block::Integer=8,
    ordering::ContainerIterationOrder=UnspecifiedOrder(),
)
    ((ax, ay, az), (bx, by, bz)) = bounds
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.((ax, bx, ay, by, az, bz))
    n_x, n_y, n_z = Int32.(nblocks)
    px, py, pz = periodic
    ppb = Int32(particles_per_block)

    con = RawContainerPoly(
        x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, px, py, pz, ppb,
    )
    return Container(con, ordering)
end

"""
    __raw(con::AbstractContainer)

Return the underlying raw container. For raw containers, return the container itself.
"""
__raw(con::AbstractRawContainer) = con
__raw(con::Container) = con.con

@doc """
    __cxxwrap_add_point!(con::RawContainer[, ord::InsertionOrder], id::Int32, x::Float64, y::Float64, z::Float64)

Wrapper for `con.put([ord,] id, x, y, z)`.
""" __cxxwrap_add_point!(::RawContainer)

@doc """
    __cxxwrap_add_point!(con::RawContainerPoly[, ord::InsertionOrder], id::Int32, x::Float64, y::Float64, z::Float64, r::Float64)

Wrapper for `con.put([ord,] id, x, y, z, r)`.
""" __cxxwrap_add_point!(::RawContainerPoly)

"""
    __add_point!(con::AbstractRawContainer, ord::ContainerIterationOrder, id, x, y, z[, r])

Add point respecting `ord` and return the container.
"""
function __add_point!(
    con::RawContainer,
    ::UnspecifiedOrder,
    id::Int32,
    x::Float64,
    y::Float64,
    z::Float64
)
    __cxxwrap_add_point!(con, id, x, y, z)
    return con
end

function __add_point!(
    con::RawContainer,
    ord::InsertionOrder,
    id::Int32,
    x::Float64,
    y::Float64,
    z::Float64
)
    __cxxwrap_add_point!(con, ord, id, x, y, z)
    return con
end

function __add_point!(
    con::RawContainerPoly,
    ::UnspecifiedOrder,
    id::Int32,
    x::Float64,
    y::Float64,
    z::Float64,
    r::Float64,
)
    __cxxwrap_add_point!(con, id, x, y, z, r)
    return con
end

function __add_point!(
    con::RawContainerPoly,
    ord::InsertionOrder,
    id::Int32,
    x::Float64,
    y::Float64,
    z::Float64,
    r::Float64,
)
    __cxxwrap_add_point!(con, ord, id, x, y, z, r)
    return con
end

"""
    add_point!(con::Container, id, (x, y, z[, r]))

Add a point to a container `con`.
# Arguments:
* `con::Container`
* `id::Integer`: must be convertible to `Int32`
* `(x, y, z)::Real`: coordinates of the particle to insert
* `r::Real`: radius (only for polydisperse containers)
"""
@propagate_inbounds function add_point!(con::Container{<:RawContainer}, id::Integer, pt)
    @boundscheck if length(pt) != 3
        throw(ArgumentError("Can only add 3-dimensional points to a VoroPlusPlus Container"))
    end
    x, y, z = pt
    __add_point!(con.con, con.ord, Int32(id), Float64.((x, y, z))...)
    return con
end

@propagate_inbounds function add_point!(con::Container{<:RawContainerPoly}, id::Integer, pt)
    @boundscheck if length(pt) != 4
        throw(ArgumentError("Can only add 3-dimensional points and radius to a polydisperse VoroPlusPlus Container"))
    end
    x, y, z, r = pt
    add_point!(con.con, con.ord, Int32(id), Float64.((x, y, z, r))...)
    return con
end

"""
    bounds(con::AbstractContainer)

Return tuple `((xmin, ymin, zmin), (xmax, ymax, zmax))`.
"""
function bounds(con::AbstractContainer)
    xlo, ylo, zlo, xhi, yhi, zhi = __cxxwrap_bounds(__raw(con))
    return (xlo, ylo, zlo), (xhi, yhi, zhi)
end

"""
    clear!(con::AbstractContainer)

Delete all data in container and return the container.
"""
function Base.empty!(con::AbstractContainer)
    __cxxwrap_clear!(__raw(con))
    return con
end

"""
   periodic(con::AbstractContainer)

Return periodicity flags in X, Y, Z directions.
"""
function periodic(con::AbstractContainer)
    return Bool.(__cxxwrap_periodic(__raw(con)))
end
