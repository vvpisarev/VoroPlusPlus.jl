function VoronoiCell(xlo, xhi, ylo, yhi, zlo, zhi)
    xmin, xmax, ymin, ymax, zmin, zmax = Float64.((xlo, xhi, ylo, yhi, zlo, zhi))
    v = VoronoiCell()
    init!(v, xmin, xmax, ymin, ymax, zmin, zmax)
    return v
end


#############################

function __vertex_ordering(vc::VoronoiCell)
    nu = reinterpret(Ptr{Int32}, __get_nu(vc))
    len = __get_current_vertices(vc)
    return unsafe_wrap(Array, nu, (len,); own=false)
end

function __vertex_positions(vc::VoronoiCell)
    nu = reinterpret(Ptr{Float64}, __get_pts(vc))
    len = __get_current_vertices(vc)
    return unsafe_wrap(Array, nu, (4*len,); own=false)
end

function vertex_positions!(pos::AbstractArray, vc::VoronoiCell)
    vp_raw = __vertex_positions(vc)
    len = num_vertices(vc)
    if length(pos) != 3 * len
        resize!(pos, 3 * len)
    end
    @inbounds for k in 0:len-1
        @views pos[3*k+1:3*k+3] .= vp_raw[4*k+1:4*k+3]
    end
    return pos
end

function vertex_positions(::Type{Vector{T}}, vc::VoronoiCell) where {T}
    len = num_vertices(vc)
    pos = Vector{T}(undef, 3*len)
    return vertex_positions!(pos, vc)
end

function vertex_positions(::Type{Vector}, vc::VoronoiCell)
    len = num_vertices(vc)
    pos = Vector{Float64}(undef, 3*len)
    return vertex_positions!(pos, vc)
end

function vertex_positions(::Type{Matrix{T}}, vc::VoronoiCell) where {T}
    len = num_vertices(vc)
    pos = Matrix{T}(undef, 3, len)
    vertex_positions!(pos, vc)
    return pos
end

function vertex_positions(::Type{Matrix}, vc::VoronoiCell)
    len = num_vertices(vc)
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