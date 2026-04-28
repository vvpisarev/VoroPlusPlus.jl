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
* `periodic`: periodicity in each axis. Default: `(false, false, false)`
* `particles_per_block::Integer`: initially allocate memory for this many particles
    per block. Default: 8
* `ordering`: `UnspecifiedOrder()` or `InsertionOrder()`. Default: `UnspecifiedOrder()`.
"""
function container(
    ;
    bounds,
    nblocks,
    periodic=(false, false, false),
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

    dx, dy, dz = (x_max, y_max, z_max) .- (x_min, y_min, z_min)
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

    dx, dy, dz = (x_max, y_max, z_max) .- (x_min, y_min, z_min)
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
    pos::AbstractVector{<:Real}
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

    dx, dy, dz = (x_max, y_max, z_max) .- (x_min, y_min, z_min)
    ilscale = cbrt(length(id) / (OPT_PART_PER_BLOCK * dx * dy * dz))
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

"""
    ordering(con::Container)

Return the object representing ordering of the container.
"""
ordering(con::Container) = con.ord

ordering(::AbstractRawContainer) = UnspecifiedOrder()

"""
    isinside(pos, con::AbstractContainer)

Test if a given vector lies within the container bounds and any walls.
"""
function isinside(pos, con::AbstractContainer)
    if eachindex(pos) != OneTo(3)
        throw(ArgumentError("Can only test 3D vectors in container"))
    end
    x, y, z = pos
    xx, yy, zz = Float64.((x, y, z))
    return Bool(__cxxwrap_isinside(__raw(con), xx, yy, zz))
end

"""
    total_particles(con::AbstractContainer)

Return the total number of points in the container.
"""
function total_particles(con::AbstractContainer)
    return total_particles(__raw(con))
end

"""
    compute_all_cells(con::AbstractContainer)

Computes all of the Voronoi cells in the container, but does nothing with the output. It is
    useful for measuring the pure computation time of the Voronoi algorithm, without any
    additional calculations such as volume evaluation or cell output.
"""
function compute_all_cells(con::AbstractContainer)
    return compute_all_cells(__raw(con))
end

"""
    nearest_particle(con::AbstractContainer, pos)

Return the parameter of the particle nearest to the position `pos`.
"""
function nearest_particle(con::AbstractContainer, pos)
    if eachindex(pos) != OneTo(3)
        throw(ArgumentError("Can only test 3D vectors in container"))
    end
    x, y, z = pos
    xx, yy, zz = Float64.((x, y, z))
    found, id, px, py, pz = __cxxwrap_find_cell(__raw(con), xx, yy, zz)
    return (Bool(found), Particle(id, (px, py, pz)))
end

"""
    draw_domain_pov(path::AbstractString, con::AbstractContainer)

Export the bounding box of the container in POV-Ray format.
"""
draw_domain_pov

"""
    draw_domain_gnuplot(file, con::AbstractContainer)

Export the bounding box of the container in Gnuplot format. `file` can be a path or an `IO`
    object.
"""
draw_domain_gnuplot

function draw_domain_gnuplot(f, con::AbstractContainer)
    __draw_domain_gnuplot(f, __raw(con))
end

function __draw_domain_gnuplot(path::AbstractString, con::AbstractRawContainer)
    open(path, "w") do io
        __draw_domain_gnuplot(io, con)
    end
end

function __draw_domain_gnuplot(io::IO, con::AbstractRawContainer)
    fmt1 = Format("%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n")
    fmt2 = Format("%g %g %g\n%g %g %g\n\n%g %g %g\n%g %g %g\n\n")
    ax, ay, az, bx, by, bz = __cxxwrap_bounds(con)
    format(io, fmt1, ax, ay, az, bx, ay, az, bx, by, az, ax, by, az)
    format(io, fmt1, ax, ay, az, ax, ay, bz, bx, ay, bz, bx, by, bz)
    format(io, fmt2, ax, by, bz, ax, ay, bz, ax, by, az, ax, by, bz)
    format(io, fmt2, bx, ay, az, bx, ay, bz, bx, by, az, bx, by, bz)
end

"""
    draw_gnuplot(prefix::AbstractString, con::AbstractContainer; domain::Bool=false, cells::Bool=true, particles::Bool=true)

Export container data to text files in gnuplot format. If `prefix` does not have asterisks,
    then files with domain, cells and particles data are named `prefix_domain.gnu`,
    `prefix_cells.gnu` and `prefix_pts.gnu`. If it contains one asterisk, it is replaced by
    `domain`, `cells` and `pts`, respectively. E.g. if `prefix = "/home/user/voronoi.*.dat`,
    then domain would be written to `/home/user/voronoi.domain.dat` etc.
"""
function draw_gnuplot(
    fmask::AbstractString, con::AbstractContainer,
    ;
    domain::Bool=false, cells::Bool=true, particles::Bool=true,
)
    astcnt = count('*', fmask)
    if astcnt > 1
        throw(
            ArgumentError(
                "Maximum one asterisk is allowed in output name template , got \"" * fmask * "\""
            )
        )
    end
    raw_con = __raw(con)

    if domain
        dom_path = astcnt > 0 ? replace(fmask, '*' => "domain") : fmask * "_domain.gnu"
        open(dom_path, "w") do io
            __draw_domain_gnuplot(io, raw_con)
        end
    end

    pt_fmt = Format("%d %g %g %g\n")
    if cells && particles
        cells_path = astcnt > 0 ? replace(fmask, '*' => "cells") : fmask * "_cells.gnu"
        open(cells_path, "w") do cells_io
            #cells_file = Libc.FILE(cells_io)
            pts_path = astcnt > 0 ? replace(fmask, '*' => "pts") : fmask * "_pts.gnu"
            open(pts_path, "w") do pts_io
                for (pt, cell) in con
                    if_valid(cell) do vc
                        (; id, pos) = pt
                        dx, dy, dz = pos
                        draw_gnuplot(cells_io, vc, pos)
                        #__cxxwrap_draw_gnuplot(cells_file, vc, dx, dy, dz)
                        format(pts_io, pt_fmt, id, dx, dy, dz)
                    end
                end
            end
            #close(cells_file)
        end
    elseif cells
        cells_path = astcnt > 0 ? replace(fmask, '*' => "cells") : fmask * "_cells.gnu"
        open(cells_path, "w") do cells_io
            #cells_file = Libc.FILE(cells_io)
            for (pt, cell) in con
                if_valid(cell) do vc
                    (; id, pos) = pt
                    dx, dy, dz = pos
                    draw_gnuplot(cells_io, vc, pos)
                    #__cxxwrap_draw_gnuplot(cells_file, vc, dx, dy, dz)
                end
            end
            #close(cells_file)
        end
    elseif particles
        pts_path = astcnt > 0 ? replace(fmask, '*' => "pts") : fmask * "_pts.gnu"
        open(pts_path, "w") do pts_io
            for pt in eachparticle(con)
                (; id, pos) = pt
                dx, dy, dz = pos
                format(pts_io, pt_fmt, id, dx, dy, dz)
            end
        end
    end
end

"""
    draw_cells_pov(path::AbstractString, con::AbstractContainer)

Export the cells geometry in POV-Ray format.
"""
function draw_cells_pov(path::AbstractString, con::AbstractContainer)
    draw_cells_pov(__raw(con), path)
    return nothing
end

"""
    draw_cells_gnuplot(path::AbstractString, con::AbstractContainer)

Export the cells geometry in Gnuplot format.
"""
function draw_cells_gnuplot(path::AbstractString, con::AbstractContainer)
    draw_cells_gnuplot(__raw(con), path)
    return nothing
end

"""
    draw_particles_pov(path::AbstractString, con::AbstractContainer)

Export the particles in POV-Ray format.
"""
function draw_particles_pov(path::AbstractString, con::AbstractContainer)
    draw_particles_pov(__raw(con), path)
    return nothing
end

"""
    draw_particles(path::AbstractString, con::AbstractContainer)

Export the particles information in text format.
"""
function draw_particles(path::AbstractString, con::AbstractContainer)
    draw_paticles(__raw(con), path)
    return nothing
end

"""
    draw_pov(prefix::AbstractString, con::AbstractContainer; domain::Bool=false, cells::Bool=true, particles::Bool=true)

Export container data to text files in POV-Ray format. If `prefix` does not have asterisks,
    then files with domain, cells and particles data are named `prefix_domain.pov`,
    `prefix_cells.pov` and `prefix_pts.pov`. If it contains one asterisk, it is replaced by
    `domain`, `cells` and `pts`, respectively. E.g. if `prefix = "/home/user/voronoi.*.pov`,
    then domain would be written to `/home/user/voronoi.domain.pov` etc.
"""
function draw_pov(
    fmask::AbstractString, con::AbstractContainer,
    ;
    domain::Bool=false, cells::Bool=true, particles::Bool=true,
)
    astcnt = count('*', fmask)
    if astcnt > 1
        throw(
            ArgumentError(
                "Maximum one asterisk is allowed in output name template , got \"" * fmask * "\""
            )
        )
    end
    raw_con = __raw(con)

    if domain
        dom_path = astcnt > 0 ? replace(fmask, '*' => "domain") : fmask * "_domain.pov"
        draw_domain_pov(raw_con, dom_path)
    end

    if cells
        cells_path = astcnt > 0 ? replace(fmask, '*' => "cells") : fmask * "_cells.gnu"
        draw_cells_pov(raw_con, cells_path)
    end
    if particles
        pts_path = astcnt > 0 ? replace(fmask, '*' => "pts") : fmask * "_pts.gnu"
        draw_particles_pov(raw_con, pts_path)
    end
end
