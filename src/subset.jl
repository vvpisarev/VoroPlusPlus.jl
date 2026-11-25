#
# Region-restricted iteration over container subsets (`voro::c_loop_subset`).
#
# Note on the upstream bug (.known_issues #18): `c_loop_subset::setup_common()` wraps the
#     starting block index (`di = i = step_mod(ci, nx)`) while `next_block()` terminates the
#     x-sweep on the *absolute* index `bi`, so any region whose lower block index is negative
#     (a query crossing a MIN face of a periodic axis) silently returns almost nothing and
#     can even segfault. The `reset_to_*!` methods below apply the verified Julia-side
#     workaround: the region is translated by whole periods along every periodic axis so
#     that its lower edge starts at or above the box minimum, which forces the correct
#     positive-overflow path of `next_block()`. The shift re-expresses the same periodic
#     query, so results are unchanged if upstream ever fixes `c_loops.cc`.

"""
    SubsetSphere

Internal description of a spherical query region after the periodic image shift has been
    applied; stores exactly the parameters passed to the C++ setup call.
"""
struct SubsetSphere
    center::SVector{3,Float64}
    r::Float64
    bounds_test::Bool
end

"""
    SubsetBox

Internal description of a box query region after the periodic image shift has been applied;
    stores exactly the parameters passed to the C++ setup call.
"""
struct SubsetBox
    lower::SVector{3,Float64}
    upper::SVector{3,Float64}
    bounds_test::Bool
end

const SubsetRegion = Union{SubsetSphere, SubsetBox}

"""
    Subset{C<:AbstractContainer}

A stateful, region-restricted iterator over particles of a container. Iterating a `Subset`
    yields `particle => cell` pairs (`particle` is a `Particle`, `cell` a
    `CheckedVoronoiCell`) for the particles located inside the current query region, using
    the same cell-reuse conventions as iteration over a whole container (see
    "Data sharing between iterations").

Create with [`sphere_region`](@ref) or [`box_region`](@ref); change the region with
    `VoroPlusPlus.reset_to_sphere!` or [`reset_to_box!`](@ref); iterate the region with the
    object itself, [`eachparticle`](@ref), or [`eachcell`](@ref).

!!! note "Position convention on periodic containers"
    Particles are always reported at their stored (primary-domain) coordinates, never at
    the periodic image which satisfied the region. Use
    `VoroPlusPlus.minimum_image_distance` / `VoroPlusPlus.minimum_image_displacement` for
    any distance-based post-processing relative to the query.

!!! note "Cells"
    Cells computed through a `Subset` are the genuine tessellation cells of the primary
    domain: a neighbour across a periodic face comes back with its own in-box cell.

!!! note "Single-pass, stateful"
    A `Subset` wraps a single C++ loop object. Starting a `for ... in s` pass re-applies
    the current region, so passes can be repeated and a region change between passes is
    picked up. Two simultaneous iterations over the same `Subset` are not supported; create
    a second `Subset` instead. Region queries inherit the general thread-safety caveats of
    containers. Walls do not restrict subset iteration.

!!! note "Periodic containers"
    `c_loop_subset` is broken upstream for regions starting below the minimum face of a
    periodic axis (silent undercount, possible segfault). `Subset` transparently shifts the
    query by whole periods into a working image and guards against degenerate inputs
    (non-finite queries, empty boxes, and regions spanning one or more whole periods all
    throw `ArgumentError`). Non-periodic axes are clipped to the container bounds by the
    upstream code.
"""
mutable struct Subset{C<:AbstractContainer}
    con::C
    loop::ContainerSubsetIterator
    region::Union{Nothing,SubsetRegion}
    shift::SVector{3,Float64}
    bounds_test::Bool
    cell::Union{Nothing,CheckedVoronoiCell}
end

"""
    __check_subset_container(con)

Throw `ArgumentError` unless `con` is a container class for which a `c_loop_subset`
    constructor exists in the C++ wrapper (`Container` or `ContainerPoly`).
"""
function __check_subset_container(con::AbstractContainer)
    if !(con isa Union{Container,ContainerPoly})
        throw(
            ArgumentError(
                "Region-restricted iteration is only implemented for `Container` and " *
                "`ContainerPoly` (the C++ `c_loop_subset` has no constructor for " *
                "$(nameof(typeof(con)))-backed containers)"
            )
        )
    end
    return con
end

"""
    Subset(con::AbstractContainer)
    Subset(tess::Tessellation)

Create a region-restricted iterator over the container of `tess`. No query region is set
    yet; use `sphere_region`/`box_region` or `VoroPlusPlus.reset_to_sphere!`/`reset_to_box!`
    before iterating.
"""
function Subset(con::C) where {C<:AbstractContainer}
    __check_subset_container(con)
    return Subset{C}(
        con,
        ContainerSubsetIterator(con),
        nothing,
        zero(SVector{3,Float64}),
        true,
        nothing,
    )
end

Subset(tess::Tessellation) = Subset(__raw(tess))

"""
    bounding_box(con::AbstractContainer)

Return tuple `((xmin, ymin, zmin), (xmax, ymax, zmax))` for a raw container.
"""
function bounding_box(con::AbstractContainer)
    (; ax, ay, az, bx, by, bz) = __cxxwrap_bounds(con)
    return (ax, ay, az), (bx, by, bz)
end

"""
    periodicity(con::AbstractContainer)

Return periodicity flags in X, Y, Z directions for a raw container.
"""
function periodicity(con::AbstractContainer)
    (; px, py, pz) = __cxxwrap_periodic(con)
    return px, py, pz
end

"""
    minimum_image_displacement(con, a, b)

Return the displacement `a - b` wrapped into the minimum image along every periodic axis
    of the container (identity on non-periodic axes). `a` and `b` are 3-element positions.
"""
function minimum_image_displacement(con::Union{AbstractContainer,<:Tessellation}, a, b)
    length(a) == length(b) == 3 ||
        throw(ArgumentError("Can only compute displacements of 3D vectors"))
    d = MVector{3,Float64}(Float64.(Tuple(a)) .- Float64.(Tuple(b)))
    per = periodicity(con)
    (xl, yl, zl), (xh, yh, zh) = bounding_box(con)
    lens = (xh - xl, yh - yl, zh - zl)
    @inbounds for k in 1:3
        if per[k]
            d[k] -= lens[k] * round(d[k] / lens[k])
        end
    end
    return SVector(d)
end

"""
    minimum_image_distance(con, a, b)

Return the distance between positions `a` and `b` using the minimum-image convention on
    every periodic axis of the container.
"""
function minimum_image_distance(con::Union{AbstractContainer,<:Tessellation}, a, b)
    d = minimum_image_displacement(con, a, b)
    return sqrt(d[1]^2 + d[2]^2 + d[3]^2)
end

"""
    __subset_point(pos, what)

Validate `pos` as a finite 3-element point and return it as a 3-tuple of `Float64`.
"""
function __subset_point(pos, what)
    if length(pos) != 3
        throw(ArgumentError("$what must be a 3-element point, got $(length(pos)) elements"))
    end
    x, y, z = pos
    xx, yy, zz = Float64.((x, y, z))
    if !all(isfinite, (xx, yy, zz))
        throw(ArgumentError("$what must have finite coordinates, got $pos"))
    end
    return (xx, yy, zz)
end

"""
    __subset_axes(con)

Return `(periodicity, lower_bounds, box_lengths)` as 3-element tuples of
    `Bool`/`Float64`.
"""
function __subset_axes(con::AbstractContainer)
    per = periodicity(con)
    lo, hi = bounding_box(con)
    return per, lo, (hi[1] - lo[1], hi[2] - lo[2], hi[3] - lo[3])
end

"""
    __subset_nearbox(value, lo_k, L_k, what)

Throw `ArgumentError` unless `value` lies within 1e7 periods of the lower face of axis `k`
    of the container. Keeps the query far away from the `double`-to-`int` conversions in
    `voro::c_loop_subset::step_int` (undefined behavior beyond the range of `int`).
"""
function __subset_nearbox(value, lo_k, L_k, what)
    if !(-1e7 * L_k <= value - lo_k <= 1e7 * L_k)
        throw(
            ArgumentError(
                "$what is more than 1e7 box lengths away from the container bounds " *
                "along axis with lower bound $lo_k and length $L_k; the query region " *
                "must stay numerically close to the container"
            )
        )
    end
    return value
end

"""
    VoroPlusPlus.reset_to_sphere!(s::Subset, center, r; bounds_test=true)

Set the query region of `s` to the sphere with center `center` and radius `r`, returning
    `s`. `bounds_test=true` (default) reports only particles inside the sphere; `false`
    reports the superset of particles in all computation blocks overlapping the sphere.

On periodic axes the query is transparently re-expressed as an equivalent periodic image
    (see the `Subset` docstring), so query centers up to 1e7 box lengths away from the
    container behave as their image. The radius must satisfy `2r < L` along every periodic
    axis of length `L`; a region spanning a whole period is rejected with `ArgumentError`
    because the underlying loop cannot represent it without duplicating particles. Throws
    `ArgumentError` on non-finite inputs, `r <= 0`, unsupported container classes, and
    regions farther than 1e7 box lengths from the container (the underlying C++ loop
    converts region coordinates to block indices and has undefined behavior for huge
    magnitudes).
"""
function reset_to_sphere!(s::Subset, center, r; bounds_test::Bool=true)
    con = s.con
    (qx, qy, qz) = __subset_point(center, "Sphere center")
    if !(r isa Real)
        throw(ArgumentError("Sphere radius must be a real number, got $(typeof(r))"))
    end
    rr = Float64(r)
    isfinite(rr) || throw(ArgumentError("Sphere radius must be finite, got $r"))
    rr > 0 || throw(ArgumentError("Sphere radius must be positive, got $r"))

    per, lo, lens = __subset_axes(con)
    c = MVector{3,Float64}((qx, qy, qz))
    shift = zeros(MVector{3,Float64})
    for k in 1:3
        Lk = lens[k]
        if per[k]
            if !(2 * rr < Lk)
                throw(
                    ArgumentError(
                        "Sphere diameter $(2 * rr) exceeds the periodic box length $Lk along " *
                        "axis $k; a region query " *
                        "on a periodic container requires 2r < L."
                    )
                )
            end
            __subset_nearbox(c[k], lo[k], Lk, "Sphere center")
            # Pre-wrap the center into the primary domain, then shift by whole periods so
            # that the region starts at or above the lower face (upstream workaround, see
            # the header comment of this file).
            ck = lo[k] + mod(c[k] - lo[k], Lk)
            n = max(0.0, ceil((lo[k] - (ck - rr)) / Lk))
            shift[k] = n * Lk
            c[k] = ck + shift[k]
        else
            __subset_nearbox(c[k] - rr, lo[k], Lk, "Sphere region (center - radius)")
            __subset_nearbox(c[k] + rr, lo[k], Lk, "Sphere region (center + radius)")
        end
    end

    s.region = SubsetSphere(SVector(c), rr, bounds_test)
    s.shift = SVector(shift)
    s.bounds_test = bounds_test
    __cxxwrap_setup_sphere!(s.loop, c[1], c[2], c[3], rr, bounds_test)
    return s
end

"""
    reset_to_box!(s::Subset, lower, upper; bounds_test=true)

Set the query region of `s` to the axis-aligned box `[lower, upper]`, returning `s`.
    `bounds_test=true` (default) reports only particles inside the box; `false` reports the
    superset of particles in all computation blocks overlapping the box.

Requires `all(lower .< upper)`. On periodic axes the box must be strictly shorter than the
    period along that axis; otherwise an `ArgumentError` is thrown. As for spheres, the
    query is transparently re-expressed as an equivalent periodic image, and non-periodic
    axes are clipped to the container bounds.
"""
function reset_to_box!(s::Subset, lower, upper; bounds_test::Bool=true)
    con = s.con
    (lx, ly, lz) = __subset_point(lower, "Box lower bound")
    (ux, uy, uz) = __subset_point(upper, "Box upper bound")
    l = SVector{3,Float64}((lx, ly, lz))
    u = SVector{3,Float64}((ux, uy, uz))
    if !all(u .> l)
        throw(ArgumentError("Box must be non-empty: need lower < upper on each axis"))
    end

    per, lo, lens = __subset_axes(con)
    offset = zeros(MVector{3,Float64})
    for k in 1:3
        Lk = lens[k]
        if per[k]
            if !(u[k] - l[k] < Lk)
                throw(
                    ArgumentError(
                        "Box width $(u[k] - l[k]) along axis $k equals or exceeds the " *
                        "periodic box length $Lk; a region may not span a whole period."
                    )
                )
            end
            __subset_nearbox(l[k], lo[k], Lk, "Box lower bound")
            # Translate both faces by the same number of periods bringing `lower` into
            # the primary domain (see the header comment of this file).
            offset[k] = lo[k] + mod(l[k] - lo[k], Lk) - l[k]
        else
            __subset_nearbox(l[k], lo[k], Lk, "Box lower bound")
            __subset_nearbox(u[k], lo[k], Lk, "Box upper bound")
        end
    end
    l = l + SVector(offset)
    u = u + SVector(offset)

    s.region = SubsetBox(l, u, bounds_test)
    s.shift = SVector(offset)
    s.bounds_test = bounds_test
    __cxxwrap_setup_box!(s.loop, l[1], u[1], l[2], u[2], l[3], u[3], bounds_test)
    return s
end

"""
    sphere_region(tess, center, r; bounds_test=true) -> Subset

Return a [`Subset`](@ref) of the container of `tess` (a `Tessellation` or a raw container)
    restricted to the sphere with center `center` and radius `r`. See
    `VoroPlusPlus.reset_to_sphere!` for restrictions and conventions.
"""
function sphere_region(tess, center, r; bounds_test::Bool=true)
    return reset_to_sphere!(Subset(tess), center, r; bounds_test=bounds_test)
end

"""
    box_region(tess, lower, upper; bounds_test=true) -> Subset

Return a [`Subset`](@ref) of the container of `tess` (a `Tessellation` or a raw container)
    restricted to the axis-aligned box `[lower, upper]`. See [`reset_to_box!`](@ref) for
    restrictions and conventions.
"""
function box_region(tess, lower, upper; bounds_test::Bool=true)
    return reset_to_box!(Subset(tess), lower, upper; bounds_test=bounds_test)
end

@doc """
    __cxxwrap_setup_sphere!(loop::ContainerSubsetIterator, x, y, z, r, bounds_test)

Raw binding of `voro::c_loop_subset::setup_sphere`. **Bypasses** all validation and the
    periodic minimum-face workaround implemented by [`Subset`](@ref); unsafe on periodic
    containers (silent undercount or crash). Escape hatch for `src/parallel/` code only.
""" __cxxwrap_setup_sphere!(::ContainerSubsetIterator)

@doc """
    __cxxwrap_setup_box!(loop::ContainerSubsetIterator, xmin, xmax, ymin, ymax, zmin, zmax, bounds_test)

Raw binding of `voro::c_loop_subset::setup_box`. **Bypasses** all validation and the
    periodic minimum-face workaround implemented by [`Subset`](@ref); unsafe on periodic
    containers (silent undercount or crash). Escape hatch for `src/parallel/` code only.
""" __cxxwrap_setup_box!(::ContainerSubsetIterator)

"""
    __subset_setup!(s::Subset)

Re-apply the stored (already shifted) region to the C++ loop object.
"""
function __subset_setup!(s::Subset)
    region = s.region
    if region isa SubsetSphere
        c = region.center
        __cxxwrap_setup_sphere!(s.loop, c[1], c[2], c[3], region.r, region.bounds_test)
    elseif region isa SubsetBox
        l, u = region.lower, region.upper
        __cxxwrap_setup_box!(
            s.loop, l[1], u[1], l[2], u[2], l[3], u[3], region.bounds_test,
        )
    else
        throw(
            ArgumentError(
                "No query region set; use sphere_region/box_region or " *
                "VoroPlusPlus.reset_to_sphere!/reset_to_box! before iterating"
            )
        )
    end
    return s
end

"""
    __subset_compute_cell!(cell, con, loop)

Compute the cell of the particle currently under `loop` (a `ContainerSubsetIterator`) by
    decoding the true block from the loop's `ijk` counter, instead of handing the loop
    object to `compute_cell!`. Passing a `c_loop_subset` to `voro`'s `compute_cell(c, vl)`
    is unsafe on periodic containers: `next_block()` leaves `vl.i` at the *absolute*
    (unwrapped) block index while `vl.ijk` wraps, and the Voronoi searcher uses `vl.i` as
    wrapped block coordinates — corrupting the search (silent invalid cells) or reading
    out of bounds (segfault). See docs/src/iteration.md and .known_issues #18.
"""
function __subset_compute_cell!(cell::CheckedVoronoiCell, con::AbstractContainer, loop)
    (; ijk, q) = __cxxwrap_loop_state(loop)
    cell._valid = Bool(__cxxwrap_compute_cell!(cell._cell, con, Int32(ijk), Int32(q)))
    return cell
end

"""
    __subset_current_cell!(s::Subset)

Compute the cell of the particle at the current loop position, lazily allocating the
    reused scratch `CheckedVoronoiCell` if needed.
"""
function __subset_current_cell!(s::Subset)
    cell = s.cell
    if cell === nothing
        cell = CheckedVoronoiCell(VoronoiCell(s.con), false)
        s.cell = cell
    end
    return __subset_compute_cell!(cell, s.con, s.loop)
end

function __subset_particle(s::Subset, loop = s.loop)
    (; x, y, z, r, id) = __cxxwrap_particle_info(loop)
    return Particle(s.con, (id, x, y, z, r))
end

function Base.iterate(s::Subset)
    __subset_setup!(s)
    !__start!(s.loop) && return nothing
    cell = __subset_current_cell!(s)
    return (__subset_particle(s) => cell), (__next!(s.loop), s.loop, cell)
end

function Base.iterate(s::Subset, (hasnext, loop, cell))
    hasnext || return nothing
    cell = __subset_compute_cell!(cell, s.con, loop)
    return (__subset_particle(s, loop) => cell), (__next!(loop), loop, cell)
end

function Base.iterate(c::EachCell{<:Subset})
    s = c.con
    __subset_setup!(s)
    !__start!(s.loop) && return nothing
    cell = __subset_current_cell!(s)
    return cell, (__next!(s.loop), s.loop, cell)
end

function Base.iterate(c::EachCell{<:Subset}, (hasnext, loop, cell))
    hasnext || return nothing
    s = c.con
    cell = __subset_compute_cell!(cell, s.con, loop)
    return cell, (__next!(loop), loop, cell)
end

function Base.iterate(c::EachParticle{<:Subset})
    s = c.con
    __subset_setup!(s)
    !__start!(s.loop) && return nothing
    return __subset_particle(s), (__next!(s.loop), s.loop)
end

function Base.iterate(c::EachParticle{<:Subset}, (hasnext, loop))
    hasnext || return nothing
    __subset_particle(c.con, loop), (__next!(loop), loop)
end

"""
    eachcell(s::Subset)

Return an iterable producing only the (reused) `CheckedVoronoiCell`s of particles in the
    query region of `s`.
"""
function eachcell(s::Subset)
    return EachCell(s)
end

"""
    eachparticle(s::Subset)

Return an iterable producing only the `Particle`s in the query region of `s`.
"""
function eachparticle(s::Subset)
    return EachParticle(s)
end

function Base.isdone(::Subset, (hasnext,))
    return !hasnext
end

particle_type(::Type{<:Subset{C}}) where {C} = particle_type(C)
particle_type(::S) where {S<:Subset} = particle_type(S)

Base.eltype(::Type{<:Subset{C}}) where {C} = Pair{particle_type(C), CheckedVoronoiCell}

Base.IteratorSize(::Type{<:Subset}) = Base.SizeUnknown()