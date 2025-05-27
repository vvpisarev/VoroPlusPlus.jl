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

function VoronoiCell(con::AbstractContainer, itor)
    raw_con = __raw(con)
    vc = VoronoiCell(raw_con)
    compute_cell!(vc, raw_con, itor)
    return vc
end

"""
    CheckedVoronoiCell

Voronoi cell with added validity flag.
"""
struct CheckedVoronoiCell<:AbstractVoronoiCell
    cell::VoronoiCellAllocated
    valid::Bool
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
    return vc, cell_is_valid
end

#############################

__get_nu(vc::VoronoiCell) = reinterpret(Ptr{Int32}, __cxxwrap_get_nu(vc))
__get_pts(vc::VoronoiCell) = reinterpret(Ptr{Float64}, __cxxwrap_get_pts(vc))

function __vertex_ordering(vc::VoronoiCell)
    nu = __get_nu(vc)
    len = __get_current_vertices(vc)
    return unsafe_wrap(Array, nu, (len,); own=false)
end

function __vertex_positions(vc::VoronoiCell)
    pts = __get_pts(vc)
    cur_vertices = __get_current_vertices(vc)
    return unsafe_wrap(Array, pts, (4*cur_vertices,); own=false)
end

function vertex_positions!(pos::AbstractVector{<:Real}, vc::VoronoiCell)
    vp_raw = __vertex_positions(vc)
    len = __get_p(vc)
    if length(pos) != 3 * len
        resize!(pos, 3 * len)
    end
    @inbounds for k in 0:len-1
        @views pos[begin+3*k:begin+3*k+2] .= 0.5 .* vp_raw[4*k+1:4*k+3]
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

@propagate_inbounds function vertex_positions!(pos::AbstractArray{<:Real}, vc::VoronoiCell)
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

function vertex_positions(::Type{Vector{T}}, vc::VoronoiCell) where {T<:Real}
    len = __get_p(vc)
    pos = Vector{T}(undef, 3*len)
    return vertex_positions!(pos, vc)
end

function vertex_positions(::Type{Vector{T}}, vc::VoronoiCell) where {T}
    len = __get_p(vc)
    pos = Vector{T}(undef, len)
    return vertex_positions!(pos, vc)
end

function vertex_positions(::Type{Vector}, vc::VoronoiCell)
    len = __get_p(vc)
    pos = Vector{SVector{3,Float64}}(undef, len)
    return vertex_positions!(pos, vc)
end

function vertex_positions(vc::VoronoiCell)
    len = __get_p(vc)
    pos = Vector{SVector{3,Float64}}(undef, len)
    return vertex_positions!(pos, vc)
end

function vertex_positions(::Type{Matrix{T}}, vc::VoronoiCell) where {T}
    len = __get_p(vc)
    pos = Matrix{T}(undef, 3, len)
    vertex_positions!(pos, vc)
    return pos
end

function vertex_positions(::Type{Matrix}, vc::VoronoiCell)
    len = __get_p(vc)
    pos = Matrix{Float64}(undef, 3, len)
    vertex_positions!(pos, vc)
    return pos
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


function output_vertices(path::AbstractString, vc::VoronoiCell, dx::Float64, dy::Float64, dz::Float64 )

    mode = "w+"
    if isfile(path)
        mode = "a+"
    end

    io = open(path, mode)
    file = Libc.FILE(io)
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "Vertex positions    : ", file)
    ccall(
        (:output_vertices_positions, "libvoro++wrap"), Cvoid,
        (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Ptr{Libc.FILE}),
        vc.cpp_object, dx, dy, dz, file,)
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "\n", file)
    close(file)
end

function output_vertices(path::AbstractString, vc::VoronoiCell )

    mode = "w+"
    if isfile(path)
        mode = "a+"
    end

    io = open(path, mode)
    file = Libc.FILE(io)
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "Vertex positions    : ", file)
    ccall(
        (:output_vertices_nopos, "libvoro++wrap"), Cvoid,
        (Ptr{Cvoid}, Ptr{Libc.FILE}),
        vc.cpp_object, file,
    )
    ccall(:fputs, Cint, (Cstring, Ptr{Libc.FILE}), "\n", file)
    close(file)
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
