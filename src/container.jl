struct UnspecifiedOrder<:ContainerIterationOrder end

mutable struct Container{C<:AbstractRawContainer, O<:ContainerIterationOrder}
    con::C
    ord::O
end

"""
    container(; bounds, nblocks, periodic=(false, false, false), particles_per_block=8, ordering)

Allocate space for a container of Voronoi cells.
# Keywords
* `bounds`: limits of the bounding box `(xmin, xmax, ymin, ymax, zmin, zmax)`
* `nblocks`: numbers of computation blocks along each axis
* `periodic::NTuple{3,Bool}`: periodicity in each axis. Default: `(false, false, false)`
* `particles_per_block::Integer`: initially allocate memory for this many particles
    per block. Default: 8
* `ordering`: `UnspecifiedOrder()` or `InsertionOrder()`. Default: `UnspecifiedOrder()`.
"""
function container(
    ;
    bounds::NTuple{6,Real},
    nblocks::NTuple{3,Integer},
    periodic::NTuple{3,Bool}=(false, false, false),
    particles_per_block::Integer=8,
    ordering::ContainerIterationOrder=UnspecifiedOrder(),
)
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.(bounds)
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
    bounds::NTuple{6,Real},
    nblocks::NTuple{3,Integer},
    periodic::NTuple{3,Bool}=(false, false, false),
    particles_per_block::Integer=8,
    ordering::ContainerIterationOrder=UnspecifiedOrder(),
)
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.(bounds)
    n_x, n_y, n_z = Int32.(nblocks)
    px, py, pz = periodic
    ppb = Int32(particles_per_block)

    con = RawContainerPoly(
        x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, px, py, pz, ppb,
    )
    return Container(con, ordering)
end

function __add_point!(
    con::RawContainer,
    ::UnspecifiedOrder,
    id::Int32,
    x::Float64,
    y::Float64,
    z::Float64
)
    return __cxxwrap_add_point!(con, id, x, y, z)
end

function __add_point!(
    con::RawContainer,
    ord::InsertionOrder,
    id::Int32,
    x::Float64,
    y::Float64,
    z::Float64
)
    return __cxxwrap_add_point!(con, ord, id, x, y, z)
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
    return __cxxwrap_add_point!(con, id, x, y, z)
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
    return __cxxwrap_add_point!(con, ord, id, x, y, z)
end

@propagate_inbounds function add_point!(con::Container{RawContainer}, id::Integer, pt)
    @boundscheck if length(pt) != 3
        throw(ArgumentError("Can only add 3-dimensional points to a VoroPlusPlus Container"))
    end
    x, y, z = pt
    __add_point!(con.con, con.ord, Int32(id), Float64.((x, y, z))...)
    return con
end

@propagate_inbounds function add_point!(con::Container{RawContainerPoly}, id::Integer, pt)
    @boundscheck if length(pt) != 4
        throw(ArgumentError("Can only add 3-dimensional points and radius to a polydisperse VoroPlusPlus Container"))
    end
    x, y, z, r = pt
    add_point!(con.con, Int32(id), Float64.((x, y, z, r))...)
    return con
end