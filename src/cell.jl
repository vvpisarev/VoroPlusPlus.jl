"""
    VoronoiCell

A wrapper for the `voronoicell_neighbor` class representing Voronoi cells with neighbor
    tracking.
"""
VoronoiCell

"""
    VoronoiCell(max_r_sq)

Create a new Voronoi cell with `max_r_sq` as maximum squared radius.
"""
VoronoiCell(::Float64)

"""
    VoronoiCell(con::Container)

Create a new Voronoi cell based on the size of `con`.
"""
VoronoiCell(con::Container) = VoronoiCell(__raw(con))

"""
    CheckedVoronoiCell

Voronoi cell with added validity flag.
"""
struct CheckedVoronoiCell<:AbstractVoronoiCell
    cell::VoronoiCellAllocated
    valid::Bool
end

function CheckedVoronoiCell(con::AbstractRawContainer, itor)
    vc = VoronoiCell(con)
    valid = __cxxwrap_compute_cell!(vc, con, itor)
    return CheckedVoronoiCell(vc, valid)
end

CheckedVoronoiCell(con::AbstractContainer, itor) = CheckedVoronoiCell(__raw(con), itor)

__raw(vc::CheckedVoronoiCell) = vc.cell
__raw(vc::VoronoiCell) = vc

isvalid(vc::CheckedVoronoiCell) = vc.valid
isvalid(::VoronoiCell) = true

"""
    if_valid(fn, vc::CheckedVoronoiCell, default=nothing)

If `vc` valid flag is `true`, apply function `fn` to `vc`, else return `default`.

# Examples

    if_valid(cell, 0.0) do vc
        volume(vc)
    end
"""
function if_valid(fn, vc::CheckedVoronoiCell, default=nothing)
    if isvalid(vc)
        return fn(__raw(vc))
    else
        return default
    end
end

"""
    if_valid(fn, vc::CheckedVoronoiCell, handler::Function)

If `vc` valid flag is `true`, apply function `fn` to `vc`, else call `handler()` which must
    be a zero-argument function.

# Examples

    if_valid(cell, ()->error("Cannot perform operation on a non-valid Voronoi cell")) do vc
        volume(vc)
    end
"""
function if_valid(fn, vc::CheckedVoronoiCell, handler::Function)
    if isvalid(vc)
        return fn(__raw(vc))
    else
        return handler()
    end
end

"""
    if_valid(fn, vc::VoronoiCell, default)

Apply function `fn` to `vc`. `default` is ignored, it is there only to unify the interface
    with `CheckedVoronoiCell`.
"""
function if_valid(fn, vc::VoronoiCell, default=nothing)
    return fn(vc)
end

"""
    volume(vc::AbstractVoronoiCell)

Return the volume of the Voronoi cell. For invalid cells, return 0.0.
"""
function volume(vc::CheckedVoronoiCell)
    isvalid(vc) ? volume(vc.cell) : 0.0
end

"""
    number_of_faces(vc::AbstractVoronoiCell)

Return the number of faces of the Voronoi cell. For invalid cells, return zero.
"""
function number_of_faces(vc::CheckedVoronoiCell)
    isvalid(vc) ? number_of_faces(vc.cell) : zero(Int32)
end

"""
    number_of_edges(vc::AbstractVoronoiCell)

Return the number of edges of the Voronoi cell. For invalid cells, return zero.
"""
function number_of_edges(vc::CheckedVoronoiCell)
    isvalid(vc) ? number_of_edges(vc.cell) : zero(Int32)
end

"""
    centroid(vc::AbstractVoronoiCell)

Return the centroid vector of a valid Voronoi cell or a zero vector if the cell is invalid.
"""
function centroid(vc::AbstractVoronoiCell)
    if_valid(vc, SVector(0.0, 0.0, 0.0)) do vc
        return SVector{3,Float64}(__cxxwrap_centroid(vc))
    end
end

"""
    voronoicell_box((xmin, ymin, zmin), (xmax, ymax, zmax))

Create a Voronoi cell initialized as rectangular cuboid with provided lower and higher
    bounds.
"""
function voronoicell_box((xlo, ylo, zlo), (xhi, yhi, zhi))
    xmin, xmax, ymin, ymax, zmin, zmax = Float64.((xlo, xhi, ylo, yhi, zlo, zhi))
    v = VoronoiCell((xmax - xmin)^2 + (ymax - ymin)^2 + (zmax - zmin)^2)
    __cxxwrap_init!(v, xmin, xmax, ymin, ymax, zmin, zmax)
    return v
end

"""
    voronoicell_tetrahedron(v1, v2, v3, v4)

Create a Voronoi cell initialized as tetrahedron with given vertices.
"""
function voronoicell_tetrahedron((u1, v1, w1), (u2, v2, w2), (u3, v3, w3), (u4, v4, w4))
    x1, y1, z1 = Float64.((u1, v1, w1))
    x2, y2, z2 = Float64.((u2, v2, w2))
    x3, y3, z3 = Float64.((u3, v3, w4))
    x4, y4, z4 = Float64.((u4, v4, w4))
    xmax, xmin = extrema((x1, x2, x3, x4))
    ymax, ymin = extrema((y1, y2, y3, y4))
    zmax, zmin = extrema((z1, z2, z3, z4))
    v = VoronoiCell((xmax - xmin)^2 + (ymax - ymin)^2 + (zmax - zmin)^2)
    __cxxwrap_init_tetrahedron!(v, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
    return v
end

"""
    voronoicell_octahedron(r::Real)

Create a Voronoi cell shaped as octahedron with vertices at `(r, 0, 0)`, `(-r, 0, 0)`,
    `(0, r, 0)`, `(0, -r, 0)`, `(0, 0, r)`, and `(0, 0, -r)`.
"""
function voronoicell_octahedron(r::Real)
    d = Float64(r)
    v = VoronoiCell(4*d*d)
    __cxxwrap_init_octahedron!(v, d)
    return v
end

"""
    reset_to_box!(vc::VoronoiCell, (u1, v1, w1), (u2, v2, w2))

Reshape a Voronoi cell into a rectangular cuboid with provided lower and higher bounds.
"""
function reset_to_box!(vc::VoronoiCell, (u1, v1, w1), (u2, v2, w2))
    xmin, xmax, ymin, ymax, zmin, zmax = Float64.((u1, u2, v1, v2, w1, w2))
    __cxxwrap_init!(vc, xmin, xmax, ymin, ymax, zmin, zmax)
    return vc
end

"""
    reset_to_tetrahedron!(vc::VoronoiCell, v1, v2, v3, v4)

Reshape a Voronoi cell into a tetrahedron with given vertices.
"""
function reset_to_tetrahedron!(
    vc::VoronoiCell,
    (u1, v1, w1),
    (u2, v2, w2),
    (u3, v3, w3),
    (u4, v4, w4)
)
    x1, y1, z1 = Float64.((u1, v1, w1))
    x2, y2, z2 = Float64.((u2, v2, w2))
    x3, y3, z3 = Float64.((u3, v3, w4))
    x4, y4, z4 = Float64.((u4, v4, w4))
    __cxxwrap_init_tetrahedron!(vc, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4)
    return vc
end


"""
    reset_to_octahedron!(vc::VoronoiCell, r::Real)

Reshape a Voronoi cell into an octahedron with vertices at `(r, 0, 0)`, `(-r, 0, 0)`,
    `(0, r, 0)`, `(0, -r, 0)`, `(0, 0, r)`, and `(0, 0, -r)`.
"""
function reset_to_octahedron!(vc::VoronoiCell, r::Real)
    d = Float64(r)
    __cxxwrap_init_octahedron!(vc, d)
    return vc
end

"""
    cut_by_particle_position!(vc::VoronoiCell, pos)

Cut cell `vc` by a particle located at `pos` relative to the cell center.
"""
function cut_by_particle_position!(vc::VoronoiCell, pos)
    if length(pos) == 3
        x, y, z = pos
        u, v, w = Float64.((x, y, z))
        __cxxwrap_nplane!(vc, u, v, w, zero(Int32))
    elseif length(pos) == 4
        x, y, z, rsq = pos
        u, v, w, dsq = Float64.((x, y, z, rsq))
        __cxxwrap_nplane!(vc, u, v, w, dsq, zero(Int32))
    end
    return vc
end

function compute_cell!(
    vc::VoronoiCell, con::AbstractRawContainer, itr
)
    cell_is_valid = __cxxwrap_compute_cell!(vc, con, itr)
    return CheckedVoronoiCell(vc, cell_is_valid)
end

function compute_cell!(vc::CheckedVoronoiCell, con::AbstractRawContainer, itr)
    compute_cell!(vc.cell, con, itr)
end

function compute_cell!(vc::AbstractVoronoiCell, con::AbstractContainer, itr)
    return compute_cell!(vc, __raw(con), itr)
end

#############################

__get_nu(vc::VoronoiCell) = reinterpret(Ptr{Int32}, __cxxwrap_get_nu(vc))
__get_pts(vc::VoronoiCell) = reinterpret(Ptr{Float64}, __cxxwrap_get_pts(vc))
__get_ed(vc::VoronoiCell) = reinterpret(Ptr{Ptr{Int32}}, __cxxwrap_get_ed(vc))
__get_ne(vc::VoronoiCell) = reinterpret(Ptr{Ptr{Int32}}, __cxxwrap_get_ne(vc))

function __vertex_ordering(vc::VoronoiCell)
    nu = __get_nu(vc)
    total_vertices = __get_p(vc)
    return unsafe_wrap(Array, nu, (total_vertices,); own=false)
end

function __vertex_positions(vc::VoronoiCell)
    pts = __get_pts(vc)
    total_vertices = __get_p(vc)
    return unsafe_wrap(Array, pts, (4*total_vertices,); own=false)
end

"""
    __cycle_up(nu::UnsafeIndexable, a::Integer, p::Integer)

Return `vc.cycle_up(a-1, p-1) + 1` if `nu` is `UnsafeIndexable(__get_nu(vc))`.
"""
function __cycle_up(nu::UnsafeIndexable, a::Integer, p::Integer)
    return a == nu[p] ? one(a) : a+one(a)
end

function __reset_edges!(
    vc::VoronoiCell,
    p=__get_p(vc),
    nu=UnsafeIndexable(__get_nu(vc)),
    ed=UnsafeIndexable(__get_ed(vc))
)
    for i in OneTo(p), j in OneTo(nu[i])
        ed_ij = ed[i, j]
        if ed_ij >= zero(ed_ij)
            error("Edge reset routine found a previously untested edge")
        else
            ed[i, j] = -ed_ij - true
        end
    end
    return vc
end

"""
    __test_vector!(v::Vector, T::DataType)

Test if a vector is suitable for storing datatype `T`.
"""
function __test_vector(v::Vector{<:Number})
    try
        push!(v, Float64(pi))
        pop!(v)
    catch
        error("Cannot store Float64 items in the vector")
    end
end

function __test_vector(v::Vector{T}) where {T}
    try
        val = T(SVector{3,Float64}(pi, pi, pi))
        push!(v, val)
        pop!(v)
    catch
        error("Cannot store SVector{3,Float64} items in the vector")
    end
end

"""
    translate!(vc::AbstractVoronoiCell, d)

Translates the vertices of the Voronoi cell by a given vector.
"""
@propagate_inbounds function translate!(vc::AbstractVoronoiCell, d)
    @boundscheck if eachindex(d) != OneTo(3)
        throw(DimensionMismatch("Translation vector must be a length-3 array or tuple."))
    end
    if_valid(vc, vc) do vc
        @inbounds x, y, z = d[1], d[2], d[3]
        __cxxwrap_translate!(vc, x, y, z)
        return vc
    end
end

"""
    number_of_vertices(vc::VoronoiCell)

Return the number of vertices of the cell.
"""
function number_of_vertices(vc::VoronoiCell)
    return __get_p(vc)
end

function vertex_positions!(pos::StdVector{Float64}, vc::VoronoiCell)
    __cxxwrap_vertices!(pos, vc)
    return pos
end

@propagate_inbounds function vertex_positions!(
    pos::StdVector{Float64}, vc::VoronoiCell, offset
)
    @boundscheck if eachindex(offset) != OneTo(3)
        throw(DimensionMismatch("Offset for vertex positions must be a length-3 array or tuple"))
    end
    dr = offset[1], offset[2], offset[3]
    __cxxwrap_vertices!(pos, vc, Float64.(dr)...)
    return pos
end

function vertex_positions!(pos::Vector{<:Number}, vc::VoronoiCell)
    vp_raw = __vertex_positions(vc)
    len = __get_p(vc)
    if length(pos) != 3 * len
        resize!(pos, 3 * len)
    end
    @inbounds for k in 0:len-1
        p = k << 2
        @views pos[begin+3*k:begin+3*k+2] .= 0.5 .* vp_raw[p+1:p+3]
    end
    return pos
end

"""
    vertex_positions!(pos::AbstractVector, vc::AbstractVoronoiCell[, offset])

Fill vector `pos` with the positions of cell vertices, optionally shifted by `offset`.

If `eltype(pos)` is numeric, then 3 vector elements are stored per vertex. Otherwise, 1
    vector is stored per vertex.

`offset` must be a length-3 array or tuple.
"""
@propagate_inbounds function vertex_positions!(
    pos::AbstractVector{<:Number}, vc::VoronoiCell, offset
)
    @boundscheck if eachindex(offset) != OneTo(3)
        throw(DimensionMismatch("Offset for vertex positions must be a length-3 array or tuple"))
    end
    dr = offset[1], offset[2], offset[3]
    vp_raw = __vertex_positions(vc)
    len = __get_p(vc)
    if length(pos) != 3 * len
        resize!(pos, 3 * len)
    end
    @inbounds for k in 0:len-1
        p = k << 2
        @views pos[begin+3*k:begin+3*k+2] .= 0.5 .* vp_raw[p+1:p+3] .+ dr
    end
    return pos
end

function vertex_positions!(pos::AbstractVector, vc::VoronoiCell)
    vp_raw = __vertex_positions(vc)
    len = __get_p(vc)
    if length(pos) != len
        resize!(pos, len)
    end
    @inbounds for k in 0:len-1
        p = k << 2
        pos[begin+k] = 0.5 .* (vp_raw[p+1], vp_raw[p+2], vp_raw[p+3])
    end
    return pos
end

@propagate_inbounds function vertex_positions!(pos::AbstractVector, vc::VoronoiCell, offset)
    @boundscheck if eachindex(offset) != OneTo(3)
        throw(DimensionMismatch("Offset for vertex positions must be a length-3 array or tuple"))
    end
    vp_raw = __vertex_positions(vc)
    len = __get_p(vc)
    if length(pos) != len
        resize!(pos, len)
    end
    @inbounds for k in 0:len-1
        p = k << 2
        pos[begin+k] = 0.5 .* (vp_raw[p+1], vp_raw[p+2], vp_raw[p+3]) .+ offset
    end
    return pos
end

@propagate_inbounds function vertex_positions!(
    pos::AbstractArray{<:Number}, vc::VoronoiCell
)
    vp_raw = __vertex_positions(vc)
    len = __get_p(vc)
    @boundscheck if length(pos) != 3 * len
        throw(DimensionMismatch("The length of the output array does not match the size of vertex array"))
    end
    for k in 0:len-1
        @views pos[begin+3*k:begin+3*k+2] .= 0.5 .* vp_raw[4*k+1:4*k+3]
    end
    return pos
end

@propagate_inbounds function vertex_positions!(
    pos::AbstractArray{<:Number}, vc::VoronoiCell, offset
)
    @boundscheck if eachindex(offset) != OneTo(3)
        throw(DimensionMismatch("Offset for vertex positions must be a length-3 array or tuple"))
    end
    dr = offset[1], offset[2], offset[3]
    vp_raw = __vertex_positions(vc)
    len = __get_p(vc)
    @boundscheck if length(pos) != 3 * len
        throw(DimensionMismatch("The length of the output array does not match the size of vertex array"))
    end
    for k in 0:len-1
        @views pos[begin+3*k:begin+3*k+2] .= 0.5 .* vp_raw[4*k+1:4*k+3] .+ dr
    end
    return pos
end

@propagate_inbounds function vertex_positions!(
    pos::AbstractVector, vc::CheckedVoronoiCell, offset...
)
    if_valid(vc, isempty(pos) ? pos : empty!(pos)) do cell
        vertex_positions!(pos, cell, offset...)
    end
end

@propagate_inbounds function vertex_positions(
    ::Type{Vector{T}}, vc::AbstractVoronoiCell, offset...
) where {T}
    pos = T[]
    if_valid(vc, pos) do vc
        vertex_positions!(pos, vc, offset...)
    end
end

@propagate_inbounds function vertex_positions(
    ::Type{Vector}, vc::AbstractVoronoiCell, offset...
)
    pos = SVector{3,Float64}[]
    if_valid(vc, pos) do vc
        vertex_positions!(pos, vc, offset...)
    end
end

"""
    vertex_positions(vc::AbstractVoronoiCell[, offset])

Return the positions of cell vertices, optionally shifted by `offset`.

The returned object is a `Vector{SVector{3,Float64}}`.

`offset` must be a length-3 array or tuple.
"""
@propagate_inbounds function vertex_positions(vc::AbstractVoronoiCell, offset...)
    pos = SVector{3,Float64}[]
    if_valid(vc, pos) do vc
        vertex_positions!(pos, vc, offset...)
    end
end

"""
    vertex_positions(::Type{<:Matrix}, vc::AbstractVoronoiCell[, offset])

Return the positions of cell vertices, optionally shifted by `offset`, as a 3xN numeric
    matrix.

`offset` must be a length-3 array or tuple.
"""
@propagate_inbounds function vertex_positions(
    ::Type{Matrix{T}}, vc::AbstractVoronoiCell, offset...
) where {T}
    len = isvalid(vc) ? __get_p(vc) : zero(Int32)
    pos = Matrix{T}(undef, 3, len)
    if_valid(vc, pos) do vc
        vertex_positions!(pos, vc, offset...)
        return pos
    end
end

@propagate_inbounds function vertex_positions(
    ::Type{Matrix}, vc::AbstractVoronoiCell, offset...
)
    pos = Float64[]
    if_valid(vc, pos) do vc
        vertex_positions!(pos, vc, offset...)
    end
    return reshape(pos, 3, :)
end

"""
    get_neighbors!(v::AbstractVector, vc::AbstractVoronoiCell)

Fill `v` with neighbor IDs of the cell `vc`.
"""
function get_neighbors!(v::AbstractVector{<:Integer}, vc::AbstractVoronoiCell)
    if_valid(vc, empty!(v)) do vc
        p = __get_p(vc)
        nu = UnsafeIndexable(__get_nu(vc))
        ed = UnsafeIndexable(__get_ed(vc))
        ne = UnsafeIndexable(__get_ne(vc))
        for i in one(p)+true:p
            nu_i = nu[i]
            for j in OneTo(nu_i)
                k = ed[i,j] + true
                if k > zero(k)
                    push!(v, ne[i,j])
                    ed[i,j] = -k
                    l = __cycle_up(nu, ed[i,j + nu_i]+true, k)
                    while true
                        m = ed[k,l] + true
                        ed[k,l] = -m
                        l = __cycle_up(nu, ed[k,l+nu[k]]+true, m)
                        k = m
                        k == i && break
                    end
                end
            end
        end
        __reset_edges!(vc, p, nu, ed)
        return v
    end
end

function append_neighbors!(v::AbstractVector{<:Integer}, vc::AbstractVoronoiCell)
    if_valid(vc, v) do vc
        p = __get_p(vc)
        nu = UnsafeIndexable(__get_nu(vc))
        ed = UnsafeIndexable(__get_ed(vc))
        ne = UnsafeIndexable(__get_ne(vc))
        for i in one(p)+true:p
            nu_i = nu[i]
            for j in OneTo(nu_i)
                k = ed[i,j] + true
                if k > zero(k)
                    push!(v, ne[i,j])
                    ed[i,j] = -k
                    l = __cycle_up(nu, ed[i,j + nu_i]+true, k)
                    while true
                        m = ed[k,l] + true
                        ed[k,l] = -m
                        l = __cycle_up(nu, ed[k,l+nu[k]]+true, m)
                        k = m
                        k == i && break
                    end
                end
            end
        end
        __reset_edges!(vc, p, nu, ed)
        return v
    end
end

function get_neighbors!(v::StdVector{Int32}, vc::VoronoiCell)
    __cxxwrap_get_neighbors!(v, vc)
    return v
end

"""
    get_normals!(v::Vector, vc::AbstractVoronoiCell)

Fill `v` with the normals to the faces of `vc`.

If `v` is a numeric vector, then normals are stored as 3 consecutive items. Otherwise,
    normals are stored as `SVector{3,Float64}`.
"""
get_normals!

function get_normals!(v::Vector, vc::VoronoiCell)
    empty!(v)
    __test_vector(v)
    p = __get_p(vc)
    nu = UnsafeIndexable(__get_nu(vc))
    ed = UnsafeIndexable(__get_ed(vc))
    pts = UnsafeIndexable(__get_pts(vc))
    for i in one(p)+true:p, j in OneTo(nu[i])
        k = ed[i, j] + true
        if k > zero(k)
            __append_normal!(v, vc, ed, nu, pts, i, j, k)
        end
    end
    __reset_edges!(vc, p, nu, ed)
    return v
end

function __append_normal!(v::Vector{<:Number}, vc, ed, nu, pts, i, j, k)
    orig_len = length(v)
    ed[i, j] = -k
    l = __cycle_up(nu, ed[i, nu[i]+j]+true, k)
    n_ed = one(k)
    S = zero(SMatrix{3, 3, Float64, 9})
    xc = pts[4 * k - 3]
    yc = pts[4 * k - 2]
    zc = pts[4 * k - 1]
    push!(v, xc, yc, zc)
    while true
        n_ed += true
        m = ed[k, l] + true
        ed[k, l] = -m
        ux = pts[4 * m - 3]
        uy = pts[4 * m - 2]
        uz = pts[4 * m - 1]
        xc += ux
        yc += uy
        zc += uz
        push!(v, ux, uy, uz)
        l = __cycle_up(nu, ed[k, nu[k]+l]+true, m)
        k = m
        k == i && break
    end
    com = SVector(xc, yc, zc) / n_ed
    # Prepare gyration tensor
    for p in orig_len:3:lastindex(v)-3
        u = SVector(v[p+1], v[p+2], v[p+3]) - com
        S += u .* u'
    end
    resize!(v, orig_len)
    vals, vecs = eigen(Symmetric(S) / n_ed)
    if 1.0 - vals[1] / vals[2] > eps() * 10
        nrm = vecs[:, 1]
        if dot(nrm, com) < 0
            nrm = -nrm
        end
        push!(v, nrm...)
    else
        push!(v, 0.0, 0.0, 0.0)
    end
    return v
end

function __append_normal!(v::Vector, vc, ed, nu, pts, i, j, k)
    ed[i, j] = -k
    l = __cycle_up(nu, ed[i, nu[i]+j]+true, k)
    n_ed = one(k)
    S = zero(SMatrix{3, 3, Float64, 9})
    xc = pts[4 * k - 3]
    yc = pts[4 * k - 2]
    zc = pts[4 * k - 1]
    push!(v, SVector(xc, yc, zc))
    while true
        n_ed += true
        m = ed[k, l] + true
        ed[k, l] = -m
        ux = pts[4 * m - 3]
        uy = pts[4 * m - 2]
        uz = pts[4 * m - 1]
        xc += ux
        yc += uy
        zc += uz
        push!(v, SVector(ux, uy, uz))
        l = __cycle_up(nu, ed[k, nu[k]+l]+true, m)
        k = m
        k == i && break
    end
    com = SVector(xc, yc, zc) / n_ed
    # Prepare gyration tensor
    for p in lastindex(v)-n_ed+1:lastindex(v)
        vpx, vpy, vpz = v[p]
        u = SVector(vpx, vpy, vpz) - com
        S += u .* u'
    end
    resize!(v, length(v) - n_ed)
    vals, vecs = eigen(Symmetric(S) / n_ed)
    if 1.0 - vals[1] / vals[2] > eps() * 10
        nrm = vecs[:, 1]
        if dot(nrm, com) < 0
            nrm = -nrm
        end
        push!(v, nrm)
    else
        push!(v, SVector(0.0, 0.0, 0.0))
    end
    return v
end

function get_normals!(v::AbstractVector, vc::CheckedVoronoiCell)
    if_valid(vc, empty!(v)) do vc
        get_normals!(v, vc)
    end
end

"""
    normals(vc::AbstractVoronoiCell)

Return the normals to the faces of `vc` as `Vector{SVector{3,Float64}}`.
"""
normals(vc::AbstractVoronoiCell) = get_normals!(SVector{3,Float64}[], vc)

function get_face_perimeters!(v::AbstractVector{<:Number}, vc::AbstractVoronoiCell)
    if_valid(vc, empty!(v)) do vc
        p = __get_p(vc)
        nu = UnsafeIndexable(__get_nu(vc))
        ed = UnsafeIndexable(__get_ed(vc))
        pts = UnsafeIndexable(__get_pts(vc))
        for i in one(p)+true:p, j in OneTo(nu[i])
            k = ed[i, j] + true
            if k > zero(k)
                dx = pts[(k<<2)-3] - pts[(i<<2)-3]
                dy = pts[(k<<2)-2] - pts[(i<<2)-2]
                dz = pts[(k<<2)-1] - pts[(i<<2)-1]
                perim = hypot(dx, dy, dz)
                ed[i, j] = -k
                l = __cycle_up(nu, ed[i, nu[i]+j]+true, k)
                while true
                    m = ed[k, l] + true
                    dx = pts[(m<<2)-3] - pts[(k<<2)-3]
                    dy = pts[(m<<2)-2] - pts[(k<<2)-2]
                    dz = pts[(m<<2)-1] - pts[(k<<2)-1]
                    perim += hypot(dx, dy, dz)
                    ed[k, l] = -m
                    l = __cycle_up(nu, ed[k, nu[k]+l]+true, m)
                    k = m
                    k == i && break
                end
                push!(v, 0.5 * perim)
            end
        end
        __reset_edges!(vc, p, nu, ed)
        return v
    end
end

function get_face_perimeters!(v::StdVector{Float64}, vc::AbstractVoronoiCell)
    if_valid(vc, empty!(v)) do cell
        __cxxwrap_face_perimeters!(v, cell)
        return v
    end
end

function face_perimeters(vc::AbstractVoronoiCell)
    v = Float64[]
    if_valid(vc, v) do cell
        get_face_perimeters!(v, vc)
    end
end

function get_face_vertices!(v::AbstractVector{<:Integer}, vc::AbstractVoronoiCell)
    if_valid(vc, empty!(v)) do vc
        vp = 1
        p = __get_p(vc)
        nu = UnsafeIndexable(__get_nu(vc))
        ed = UnsafeIndexable(__get_ed(vc))

        for i in one(p)+true:p, j in OneTo(nu[i])
            k = ed[i, j] + true
            if k > zero(k)
                push!(v, false)
                push!(v, i-true)
                ed[i, j] = -k
                l = __cycle_up(nu, ed[i, nu[i]+j]+true, k)
                while true
                    push!(v, k-true)
                    m = ed[k, l] + true
                    ed[k, l] = -m
                    l = __cycle_up(nu, ed[k, nu[k]+l]+true, m)
                    k = m
                    k == i && break
                end
                vn = length(v)
                v[vp] = vn - vp
                vp = vn + true
            end
        end
        __reset_edges!(vc, p, nu, ed)
        return v
    end
end

function get_face_vertices!(v::StdVector{Int32}, vc::AbstractVoronoiCell)
    if_valid(vc, empty!(v)) do cell
        __cxxwrap_face_vertices!(v, cell)
        return v
    end
end

function get_face_orders!(v::AbstractVector{<:Integer}, vc::AbstractVoronoiCell)
    if_valid(vc, empty!(v)) do vc
        p = __get_p(vc)
        nu = UnsafeIndexable(__get_nu(vc))
        ed = UnsafeIndexable(__get_ed(vc))
        for i in one(p)+true:p, j in OneTo(nu[i])
            k = ed[i, j] + true
            if k > zero(k)
                q = one(k)
                ed[i, j] = -k
                l = __cycle_up(nu, ed[i, nu[i]+j]+true, k)
                while true
                    q += true
                    m = ed[k, l] + true
                    ed[k, l] = -m
                    l = __cycle_up(nu, ed[k, nu[k]+l]+true, m)
                    k = m
                    k == i && break
                end
                push!(v, q)
            end
        end
        __reset_edges!(vc, p, nu, ed)
        return v
    end
end

function get_face_orders!(v::StdVector{Int32}, vc::AbstractVoronoiCell)
    if_valid(vc, empty!(v)) do cell
        __cxxwrap_face_orders!(v, vc)
        return v
    end
end

function face_orders(vc::AbstractVoronoiCell)
    v = Int[]
    if_valid(vc, v) do vc
        get_face_orders!(v, vc)
    end
end

function draw_gnuplot(io::IOStream, vc::VoronoiCell, disp = (0.0, 0.0, 0.0))
    _x, _y, _z = disp
    dx, dy, dz = Float64.((_x, _y, _z))
    file = Libc.FILE(io)
    ccall(
        (:draw_gnuplot_voronoicell, "libvoro++wrap"),
        Cvoid,
        (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Ptr{Libc.FILE}),
        vc.cpp_object, dx, dy, dz, file,
    )
    close(file)
end

function draw_gnuplot(path::AbstractString, vc::VoronoiCell, disp = (0.0, 0.0, 0.0))
    open(path, "w+") do io
        draw_gnuplot(io, vc, disp)
    end
    # _x, _y, _z = disp
    # dx, dy, dz = Float64.((_x, _y, _z))
    # draw_gnuplot(vc, dx, dy, dz, path)
end

"""
    draw_gnuplot(output, vc::VoronoiCell, disp=(0.0, 0.0, 0.0))

Output cell Outputs the edges of the Voronoi cell in gnuplot format to an output stream or
    to a file given by name. `disp` is a displacement to add to cell position.
"""
function draw_gnuplot(io::IO, vc::VoronoiCell, (dx, dy, dz) = (0.0, 0.0, 0.0))
    p = __get_p(vc)
    nu = UnsafeIndexable(__get_nu(vc))
    ed = UnsafeIndexable(__get_ed(vc))
    pts = UnsafeIndexable(__get_pts(vc))

    fmt = Format("%g %g %g\n")
    for i in one(p)+true:p, j in OneTo(nu[i])
        k = ed[i, j] + true
        if k > zero(k)
            format(io, fmt, 0.5 * pts[i<<2-3] + dx, 0.5 * pts[i<<2-2] + dy, 0.5 * pts[i<<2-1] + dz)
            l, m = i, j
            while true
                ed[k, ed[l, nu[l]+m]+true] = -l
                ed[l, m] = -k
                l = k
                format(io, fmt, 0.5 * pts[k<<2-3] + dx, 0.5 * pts[k<<2-2] + dy, 0.5 * pts[k<<2-1] + dz)
                search_edge = false
                m = one(m)
                while m <= nu[l]
                    k = ed[l, m] + true
                    if k > zero(k)
                        search_edge = true
                        break
                    end
                    m += one(m)
                end
                search_edge || break
            end
            print(io, "\n\n")
        end
    end
    __reset_edges!(vc, p, nu, ed)
end


function output_vertex_orders(path::AbstractString, vc::VoronoiCell)

    mode = "w+"
    if isfile(path)
        mode = "a+"
    end

    io = open(path, mode)
    file = Libc.FILE(io)
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "Vertex orders       : ", file)
    ccall(
        (:output_vertex_orders_vorocell, "libvoro++wrap"),
        Cvoid,
        (Ptr{Cvoid}, Ptr{Libc.FILE}),
        vc.cpp_object, file,
    )
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "\n", file)
    close(file)
end


function output_face_perimeters(path::AbstractString, vc::VoronoiCell)

    mode = "w+"
    if isfile(path)
        mode = "a+"
    end

    io = open(path, mode)
    file = Libc.FILE(io)
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "Face perimeters     : ", file)
    ccall(
        (:output_face_perimeters_vorocell, "libvoro++wrap"),
        Cvoid,
        (Ptr{Cvoid}, Ptr{Libc.FILE}),
        vc.cpp_object, file,
    )
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "\n", file)
    close(file)
end

######

function output_face_freq_table(path::AbstractString, vc::VoronoiCell)

    mode = "w+"
    if isfile(path)
        mode = "a+"
    end

    io = open(path, mode)
    file = Libc.FILE(io)
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "Face freq. table    : ", file)
    ccall(
        (:output_face_freq_table_vorocell, "libvoro++wrap"),
        Cvoid,
        (Ptr{Cvoid}, Ptr{Libc.FILE}),
        vc.cpp_object, file,
    )
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "\n", file)
    close(file)
end

function output_face_orders(path::AbstractString, vc::VoronoiCell)

    mode = "w+"
    if isfile(path)
        mode = "a+"
    end

    io = open(path, mode)
    file = Libc.FILE(io)
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "Face orders         : ", file)
    ccall(
        (:output_face_orders_vorocell, "libvoro++wrap"),
        Cvoid,
        (Ptr{Cvoid}, Ptr{Libc.FILE}),
        vc.cpp_object, file,
    )
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "\n", file)
    close(file)
end

function output_face_areas(path::AbstractString, vc::VoronoiCell)

    mode = "w+"
    if isfile(path)
        mode = "a+"
    end

    io = open(path, mode)
    file = Libc.FILE(io)
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "Face areas          : ", file)
    ccall(
        (:output_face_areas_vorocell, "libvoro++wrap"),
        Cvoid,
        (Ptr{Cvoid}, Ptr{Libc.FILE}),
        vc.cpp_object, file,
    )
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "\n", file)
    close(file)
end

function output_normals(path::AbstractString, vc::VoronoiCell)

    mode = "w+"
    if isfile(path)
        mode = "a+"
    end

    io = open(path, mode)
    file = Libc.FILE(io)
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "Face normals        : ", file)
    ccall(
        (:output_normals_vorocell, "libvoro++wrap"),
        Cvoid,
        (Ptr{Cvoid}, Ptr{Libc.FILE}),
        vc.cpp_object, file,
    )
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "\n", file)
    close(file)
end

function output_face_vertices(path::AbstractString, vc::VoronoiCell)

    mode = "w+"
    if isfile(path)
        mode = "a+"
    end

    io = open(path, mode)
    file = Libc.FILE(io)
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "Face vertices       : ", file)
    ccall(
        (:output_face_vertices_vorocell, "libvoro++wrap"),
        Cvoid,
        (Ptr{Cvoid}, Ptr{Libc.FILE}),
        vc.cpp_object, file,
    )
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "\n", file)
    close(file)
end
