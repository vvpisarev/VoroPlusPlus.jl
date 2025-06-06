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
    VoronoiCell(con::AbstractRawContainer)

Create a new Voronoi cell based on the size of `con`.
"""
VoronoiCell(::AbstractRawContainer)

"""
    CheckedVoronoiCell

Voronoi cell with added validity flag.
"""
struct CheckedVoronoiCell<:AbstractVoronoiCell
    cell::VoronoiCellAllocated
    valid::Bool
end

function CheckedVoronoiCell(con::AbstractContainer, itor)
    raw_con = __raw(con)
    vc = VoronoiCell(raw_con)
    return compute_cell!(vc, raw_con, itor)
end

__raw(vc::CheckedVoronoiCell) = vc.cell
__raw(vc::VoronoiCell) = vcat

isvalid(vc::CheckedVoronoiCell) = vc.valid
isvalid(::VoronoiCell) = true

"""
    if_valid(fn, vc::CheckedVoronoiCell, default=missing)

If `vc` valid flag is `true`, apply function `fn` to `vc`, else return `default`.
"""
function if_valid(fn, vc::CheckedVoronoiCell, default=missing)
    if isvalid(vc)
        return fn(__raw(vc))
    else
        return default
    end
end

"""
    if_valid(fn, vc::VoronoiCell, default)

Apply function `fn` to `vc`. `default` is ignored, it is there only to unify the interface
    with `CheckedVoronoiCell`.
"""
function if_valid(fn, vc::VoronoiCell, default=missing)
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

function voronoicell_octahedron(r::Real)
    d = Float64(r)
    v = VoronoiCell(4*d*d)
    __cxxwrap_init_octahedron!(v, d)
    return v
end

function reset_to_box!(vc::VoronoiCell, (u1, v1, w1), (u2, v2, w2))
    xmin, xmax, ymin, ymax, zmin, zmax = Float64.((u1, u2, v1, v2, w1, w2))
    __cxxwrap_init!(vc, xmin, xmax, ymin, ymax, zmin, zmax)
    return vc
end

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
        __cxxwrap_plane!(vc, u, v, w)
    elseif length(pos) == 4
        x, y, z, rsq = pos
        u, v, w, dsq = Float64.((x, y, z, rsq))
        __cxxwrap_plane!(vc, u, v, w, dsq)
    end
    return vc
end

function compute_cell!(
    vc::VoronoiCell, con::AbstractRawContainer, itr
)
    cell_is_valid = convert(Bool, __cxxwrap_compute_cell!(vc, con, itr))
    return CheckedVoronoiCell(vc, cell_is_valid)
end

function compute_cell!(vc::CheckedVoronoiCell, con::AbstractRawContainer, itr)
    compute_cell!(vc.cell, con, itr)
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

function __reset_edges!(
    vc::VoronoiCell,
    p=__get_p(vc),
    nu=UnsafeIndexable(__get_nu(vc)),
    ed=UnsafeIndexable(__get_ed(vc))
)
    for i in OneTo(p)
        for j in OneTo(nu[i])
            ed_ij = ed[i, j]
            if ed_ij >= zero(ed_ij)
                error("Edge reset routine found a previously untested edge")
            else
                ed[i, j] = -ed_ij - true
            end
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

function number_of_vertices(vc::VoronoiCell)
    return __get_p(vc)
end

function vertex_positions!(pos::StdVector{Float64}, vc::VoronoiCell)
    __cxxwrap_vertices!(pos, vc)
    return pos
end

function vertex_positions!(pos::AbstractVector{<:Number}, vc::VoronoiCell)
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

@propagate_inbounds function vertex_positions!(pos::AbstractArray{<:Number}, vc::VoronoiCell)
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

function vertex_positions!(pos::AbstractVector, vc::CheckedVoronoiCell)
    if_valid(vc, isempty(pos) ? pos : empty!(pos)) do cell
        vertex_positions!(pos, cell)
    end
end

function vertex_positions(::Type{Vector{T}}, vc::AbstractVoronoiCell) where {T}
    pos = T[]
    if_valid(vc, pos) do vc
        vertex_positions!(pos, vc)
    end
end

function vertex_positions(::Type{Vector}, vc::AbstractVoronoiCell)
    pos = SVector{3,Float64}[]
    if_valid(vc, pos) do vc
        vertex_positions!(pos, vc)
    end
end

function vertex_positions(vc::AbstractVoronoiCell)
    pos = SVector{3,Float64}[]
    if_valid(vc, pos) do vc
        vertex_positions!(pos, vc)
    end
end

function vertex_positions(::Type{Matrix{T}}, vc::AbstractVoronoiCell) where {T}
    len = isvalid(vc) ? __get_p(vc) : zero(Int32)
    pos = Matrix{T}(undef, 3, len)
    if_valid(vc, pos) do vc
        vertex_positions!(pos, vc)
        return pos
    end
end

function vertex_positions(::Type{Matrix}, vc::AbstractVoronoiCell)
    pos = Float64[]
    if_valid(vc, pos) do vc
        vertex_positions!(pos, vc)
    end
    return reshape(pos, 3, :)
end

function get_neighbors!(v::Vector{<:Integer}, vc::VoronoiCell)
    empty!(v)
    p = __get_p(vc)
    nu = UnsafeIndexable(__get_nu(vc))
    ed = UnsafeIndexable(__get_ed(vc))
    ne = UnsafeIndexable(__get_ne(vc))
    for i in one(p)+true:p
        nu_i = nu[i]
        for j in OneTo(nu_i)
            k = ed[i,j] + true
            if k >= one(k)
                push!(v, ne[i,j])
                ed[i,j] = -k
                l = __cycle_up(vc, ed[i,j + nu_i], k-true) + true
                while true
                    m = ed[k,l] + true
                    ed[k,l] = -m
                    l = __cycle_up(vc, ed[k,l+nu[k]], m-true) + true
                    k = m
                    k == i && break
                end
            end
        end
    end
    __reset_edges!(vc, p, nu, ed)
    return v
end

function get_neighbors!(v::StdVector{Int32}, vc::VoronoiCell)
    __cxxwrap_get_neighbors!(v, vc)
    return v
end

function get_neighbors!(v::AbstractVector{<:Integer}, vc::CheckedVoronoiCell)
    if_valid(vc, isempty(v) ? v : empty!(v)) do vc
        get_neighbors!(v, vc)
    end
end

function get_normals!(v::Vector, vc::VoronoiCell)
    empty!(v)
    __test_vector(v)
    p = __get_p(vc)
    nu = UnsafeIndexable(__get_nu(vc))
    ed = UnsafeIndexable(__get_ed(vc))
    pts = UnsafeIndexable(__get_pts(vc))
    for i in one(p)+true:p
        for j in Base.OneTo(nu[i])
            k = ed[i, j] + true
            if k >= one(k)
                __append_normal!(v, vc, ed, nu, pts, i, j, k)
            end
        end
    end
    __reset_edges!(vc, p, nu, ed)
    return v
end

function __append_normal!(v::Vector{<:Number}, vc, ed, nu, pts, i, j, k)
    orig_len = length(v)
    ed[i, j] = -k
    l = __cycle_up(vc, ed[i, nu[i]+j], k-true) + true
    n_ed = one(k)
    Sxx = Syy = Szz = Sxy = Sxz = Syz = 0.0
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
        l = __cycle_up(vc, ed[k, nu[k]+l], m-true) + true
        k = m
        k == i && break
    end
    com = SVector(xc, yc, zc) / n_ed
    # Prepare gyration tensor
    for p in orig_len:3:lastindex(v)-3
        ux, uy, uz = (v[p+1], v[p+2], v[p+3]) .- com
        Sxx += ux * ux
        Syy += uy * uy
        Szz += uz * uz
        Sxy += ux * uy
        Sxz += ux * uz
        Syz += uy * uz
    end
    resize!(v, orig_len)
    S = Symmetric(SMatrix{3,3,Float64,9}(Sxx, Sxy, Sxz, Sxy, Syy, Syz, Sxz, Syz, Szz))
    vals, vecs = eigen(S / n_ed)
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
    l = __cycle_up(vc, ed[i, nu[i]+j], k-true) + true
    n_ed = one(k)
    Sxx = Syy = Szz = Sxy = Sxz = Syz = 0.0
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
        l = __cycle_up(vc, ed[k, nu[k]+l], m-true) + true
        k = m
        k == i && break
    end
    com = SVector(xc, yc, zc) / n_ed
    # Prepare gyration tensor
    for p in lastindex(v)-n_ed+1:lastindex(v)
        ux, uy, uz = v[p] .- com
        Sxx += ux * ux
        Syy += uy * uy
        Szz += uz * uz
        Sxy += ux * uy
        Sxz += ux * uz
        Syz += uy * uz
    end
    resize!(v, length(v) - n_ed)
    S = Symmetric(SMatrix{3,3,Float64,9}(Sxx, Sxy, Sxz, Sxy, Syy, Syz, Sxz, Syz, Szz))
    vals, vecs = eigen(S / n_ed)
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
    if_valid(vc, isempty(v) ? v : empty!(v)) do vc
        get_normals!(v, vc)
    end
end

normals(vc::AbstractVoronoiCell) = get_normals!(SVector{3,Float64}[], vc)

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


function draw_gnuplot(io::IO, vc::VoronoiCell, (dx, dy, dz) = (0.0, 0.0, 0.0))
    p = __get_p(vc)
    nu = UnsafeIndexable(__get_nu(vc))
    ed = UnsafeIndexable(__get_ed(vc))
    pts = UnsafeIndexable(__get_pts(vc))
	
    fmt = Format("%g %g %g\n")
    for i in one(p)+true:p
        for j in Base.OneTo(nu[i])
            k = ed[i, j] + true
            if k >= one(k)
                format(io, fmt, 0.5 * pts[i<<2-3] + dx, 0.5 * pts[i<<2-2] + dy, 0.5 * pts[i<<2-1] + dz)
                l, m = i, j
                while true
                    ed[k, ed[l, nu[l]+m]] = -l
                    ed[l, m] = -k
                    l = k
                    format(io, fmt, 0.5 * pts[k<<2-3] + dx, 0.5 * pts[k<<2-2] + dy, 0.5 * pts[k<<2-1] + dz)
                    search_edge = false
                    m = one(m)
                    while m <= nu[l]
                        k = ed[l, m] + true
                        if k >= one(k)
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
