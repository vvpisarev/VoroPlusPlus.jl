"""
    UnspecifiedOrder

A singleton type to denote the lack of a specific order of iteration over a Voronoi
    container.
"""
struct UnspecifiedOrder<:ContainerIterationOrder end

mutable struct Container{C<:AbstractRawContainer, O<:ContainerIterationOrder}<:AbstractContainer
    const con::C
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
    bounds::Tuple{NTuple{3,Real},NTuple{3,Real}},
    nblocks::NTuple{3,Integer},
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
    voronoi_tessellation(pos; [id,] bounds, periodic=(false, false, false), ordering)

Allocate space for a container of Voronoi cells and add coordinates stored in `pos`.
# Arguments
* `pos`: vector of positions
# Keywords
* `id`: use those IDs instead of indices of `pos`
* `bounds`: limits of the bounding box `((xmin, ymin, zmin), (xmax, ymax, zmax))`
* `periodic::NTuple{3,Bool}`: periodicity in each axis. Default: `(false, false, false)`
* `ordering`: `UnspecifiedOrder()` or `InsertionOrder()`. Default: `UnspecifiedOrder()`.
"""
function voronoi_tessellation(
    pos::AbstractVector
    ;
    bounds,
    id::AbstractVector{<:Integer}=OneTo(length(pos)),
    periodic=(false, false, false),
    ordering::ContainerIterationOrder=UnspecifiedOrder(),
)
    eachindex(bounds) == OneTo(2) &&
    eachindex(bounds[1]) == eachindex(bounds[2]) == OneTo(3) ||
    throw(DimensionMismatch("Bounds must be a 2-tuple of length-3 arrays or tuples"))

    eachindex(periodic) == OneTo(3) ||
    throw(DimensionMismatch("Periodic flags must be a length-3 boolean array or tuple"))

    eachindex(pos) == eachindex(id) ||
    throw(DimensionMismatch("Position and ID vectors must have the same dimensions"))

    ((ax, ay, az), (bx, by, bz)) = bounds
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.((ax, bx, ay, by, az, bz))
    px, py, pz = (Bool(p) for p in periodic)
    ppb = Int32(8)

    ilscale = cbrt(length(pos) / (OPT_PART_PER_BLOCK * dx * dy * dz))
    nx, ny, nz = floor.(Int32, (dx, dy, dz) .* ilscale .+ 1)

    rcon = RawContainer(
        x_min, x_max, y_min, y_max, z_min, z_max, nx, ny, nz, px, py, pz, ppb,
    )
    con = Container(rcon, ordering)

    for (pid, p) in zip(id, pos)
        add_point!(con, pid, p)
    end
    return con
end

function voronoi_tessellation(
    pos::AbstractVector{T}
    ;
    bounds,
    id::AbstractVector{<:Integer}=OneTo(length(pos)),
    periodic=(false, false, false),
    ordering::ContainerIterationOrder=UnspecifiedOrder(),
) where {T<:Union{SVector{3,<:Real},NTuple{3,Real}}}
    eachindex(bounds) == OneTo(2) &&
    eachindex(bounds[1]) == eachindex(bounds[2]) == OneTo(3) ||
    throw(DimensionMismatch("Bounds must be a 2-tuple of length-3 arrays or tuples"))

    eachindex(periodic) == OneTo(3) ||
    throw(DimensionMismatch("Periodic flags must be a length-3 boolean array or tuple"))

    eachindex(pos) == eachindex(id) ||
    throw(DimensionMismatch("Position and ID vectors must have the same dimensions"))

    ((ax, ay, az), (bx, by, bz)) = bounds
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.((ax, bx, ay, by, az, bz))
    px, py, pz = (Bool(p) for p in periodic)
    ppb = Int32(8)

    ilscale = cbrt(length(pos) / (OPT_PART_PER_BLOCK * dx * dy * dz))
    nx, ny, nz = floor.(Int32, (dx, dy, dz) .* ilscale .+ 1)

    rcon = RawContainer(
        x_min, x_max, y_min, y_max, z_min, z_max, nx, ny, nz, px, py, pz, ppb,
    )
    con = Container(rcon, ordering)

    @inbounds for (pid, p) in zip(id, pos)
        add_point!(con, pid, p)
    end
    return con
end

function voronoi_tessellation(
    pos::AbstractVector{<:Number}
    ;
    bounds,
    id::AbstractVector{<:Integer}=eachindex(1:3:length(pos)),
    periodic=(false, false, false),
    ordering::ContainerIterationOrder=UnspecifiedOrder(),
)
    eachindex(bounds) == OneTo(2) &&
    eachindex(bounds[1]) == eachindex(bounds[2]) == OneTo(3) ||
    throw(DimensionMismatch("Bounds must be a 2-tuple of length-3 arrays or tuples"))

    eachindex(periodic) == OneTo(3) ||
    throw(DimensionMismatch("Periodic flags must be a length-3 boolean array or tuple"))

    length(pos) == length(id) * 3 ||
    throw(DimensionMismatch("Position and ID vectors must have the same dimensions"))

    ((ax, ay, az), (bx, by, bz)) = bounds
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.((ax, bx, ay, by, az, bz))
    px, py, pz = (Bool(p) for p in periodic)
    ppb = Int32(8)

    ilscale = cbrt(length(pos) / (OPT_PART_PER_BLOCK * dx * dy * dz))
    nx, ny, nz = floor.(Int32, (dx, dy, dz) .* ilscale .+ 1)

    rcon = RawContainer(
        x_min, x_max, y_min, y_max, z_min, z_max, nx, ny, nz, px, py, pz, ppb,
    )
    con = Container(rcon, ordering)

    @inbounds for (k, pid) in enumerate(id)
        x, y, z = @views pos[begin+k-1:begin+k+2]
        add_point!(con, pid, (x, y, z))
    end
    return con
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
    bounds::Tuple{NTuple{3,Real},NTuple{3,Real}},
    nblocks::NTuple{3,Integer},
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

__raw_type(::Type{Container{C,O}}) where {C,O} = C
__raw_type(::Type{C}) where {C<:AbstractRawContainer} = C

@doc """
    __cxxwrap_put!(con::RawContainer[, ord::InsertionOrder], id::Int32, x::Float64, y::Float64, z::Float64)

Wrapper for `con.put([ord,] id, x, y, z)`.
""" __cxxwrap_put!(::RawContainer)

@doc """
    __cxxwrap_put!(con::RawContainerPoly[, ord::InsertionOrder], id::Int32, x::Float64, y::Float64, z::Float64, r::Float64)

Wrapper for `con.put([ord,] id, x, y, z, r)`.
""" __cxxwrap_put!(::RawContainerPoly)

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
    __cxxwrap_put!(con, id, x, y, z)
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
    __cxxwrap_put!(con, ord, id, x, y, z)
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
    __cxxwrap_put!(con, id, x, y, z, r)
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
    __cxxwrap_put!(con, ord, id, x, y, z, r)
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
@propagate_inbounds function add_point!(con::Container{<:RawContainer}, id::Integer, pos)
    @boundscheck if length(pos) != 3
        throw(ArgumentError("Can only add 3-dimensional points to a VoroPlusPlus Container"))
    end
    x, y, z = pos
    __add_point!(con.con, con.ord, Int32(id), Float64.((x, y, z))...)
    return con
end

@propagate_inbounds function add_point!(con::Container{<:RawContainerPoly}, id::Integer, posr)
    @boundscheck if length(posr) != 4
        throw(ArgumentError("Can only add 3-dimensional points and radius to a polydisperse VoroPlusPlus Container"))
    end
    x, y, z, r = posr
    __add_point!(con.con, con.ord, Int32(id), Float64.((x, y, z, r))...)
    return con
end

@propagate_inbounds function add_point!(
    con::Container{<:RawContainerPoly}, id::Integer, pos, r::Real
)
    @boundscheck if length(pos) != 3
        throw(ArgumentError("Can only add 3-dimensional points and radius to a polydisperse VoroPlusPlus Container"))
    end
    x, y, z = pos
    __add_point!(con.con, con.ord, Int32(id), Float64.((x, y, z, r))...)
    return con
end

"""
    bounding_box(con::AbstractContainer)

Return tuple `((xmin, ymin, zmin), (xmax, ymax, zmax))`.
"""
function bounding_box(con::AbstractContainer)
    xlo, ylo, zlo, xhi, yhi, zhi = __cxxwrap_bounds(__raw(con))
    return (xlo, ylo, zlo), (xhi, yhi, zhi)
end

"""
    empty!(con::AbstractContainer)

Delete all data in container and return the container.
"""
function Base.empty!(con::AbstractRawContainer)
    __cxxwrap_clear!(con)
    return con
end

function Base.empty!(con::Container{C,O}) where {C,O}
    __cxxwrap_clear!(__raw(con))
    con.ord = O()
    return con
end

"""
    periodicity(con::AbstractContainer)

Return periodicity flags in X, Y, Z directions.
"""
function periodicity(con::AbstractContainer)
    return Bool.(__cxxwrap_periodic(__raw(con)))
end

ordering(con::Container) = con.ord

ordering(con::AbstractRawContainer) = UnspecifiedOrder()
