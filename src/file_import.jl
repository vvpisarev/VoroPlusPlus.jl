"""
    read_particles(input; bounds, periodic=(false, false, false), ordering)

read particle data from `input` which can be a path or I/O object.
# Keywords
* `bounds`: limits of the bounding box `((xmin, ymin, zmin), (xmax, ymax, zmax))`
* `periodic`: periodicity in each axis. Default: `(false, false, false)`
* `ordering`: `UnspecifiedOrder()` or `InsertionOrder()`. Default: `UnspecifiedOrder()`.
"""
function read_particles(
    input::Union{IO,AbstractString},
    ;
    bounds,
    periodic=(false, false, false),
    ordering::ContainerIterationOrder=UnspecifiedOrder(),
)
    eachindex(bounds) == OneTo(2) &&
    eachindex(bounds[1]) == eachindex(bounds[2]) == OneTo(3) ||
    throw(DimensionMismatch("Bounds must be a 2-tuple of length-3 arrays or tuples"))

    eachindex(periodic) == OneTo(3) ||
    throw(DimensionMismatch("Periodic flags must be a length-3 boolean array or tuple"))
    ((ax, ay, az), (bx, by, bz)) = bounds
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.((ax, bx, ay, by, az, bz))
    px, py, pz = (Bool(p) for p in periodic)
    ppb = Int32(8)

    dx = x_max - x_min
    dy = y_max - y_min
    dz = z_max - z_min

    particles = Particle{Nothing}[]
    for ln in eachline(input)
        spl = eachsplit(ln)
        id, x, y, z = parse.((Int32, Float64, Float64, Float64), spl)
        if (px || (ax <= x <= bx)) && (py || (ay <= y <= by)) && (pz || (az <= z <= bz))
            push!(particles, particle(id, (x, y, z)))
        else
            @warn "Particle out of bounds" x y z
        end
    end

    ilscale = cbrt(length(particles) / (OPT_PART_PER_BLOCK * dx * dy * dz))
    nx, ny, nz = floor.(Int32, (dx, dy, dz) .* ilscale .+ 1)
    raw_con = RawContainer(
        x_min, x_max, y_min, y_max, z_min, z_max, nx, ny, nz, px, py, pz, ppb,
    )
    con = Container(raw_con, ordering)

    @inbounds for p in particles
        add_point!(con, p)
    end
    return con
end

"""
    read_polydisperse_particles(input; bounds, periodic=(false, false, false), ordering)

read particle data from `input` which can be a path or I/O object.
# Keywords
* `bounds`: limits of the bounding box `((xmin, ymin, zmin), (xmax, ymax, zmax))`
* `periodic`: periodicity in each axis. Default: `(false, false, false)`
* `ordering`: `UnspecifiedOrder()` or `InsertionOrder()`. Default: `UnspecifiedOrder()`.
"""
function read_polydisperse_particles(
    input::Union{IO,AbstractString},
    ;
    bounds,
    periodic=(false, false, false),
    ordering::ContainerIterationOrder=UnspecifiedOrder(),
)
    eachindex(bounds) == OneTo(2) &&
    eachindex(bounds[1]) == eachindex(bounds[2]) == OneTo(3) ||
    throw(DimensionMismatch("Bounds must be a 2-tuple of length-3 arrays or tuples"))

    eachindex(periodic) == OneTo(3) ||
    throw(DimensionMismatch("Periodic flags must be a length-3 boolean array or tuple"))
    ((ax, ay, az), (bx, by, bz)) = bounds
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.((ax, bx, ay, by, az, bz))
    px, py, pz = (Bool(p) for p in periodic)
    ppb = Int32(8)

    dx = x_max - x_min
    dy = y_max - y_min
    dz = z_max - z_min

    particles = Particle{Nothing}[]
    for ln in eachline(input)
        spl = eachsplit(ln)
        id, x, y, z, r = parse.((Int32, Float64, Float64, Float64, Float64), spl)
        if (px || (ax <= x <= bx)) && (py || (ay <= y <= by)) && (pz || (az <= z <= bz))
            push!(particles, particle(id, (x, y, z), r))
        else
            @warn "Particle out of bounds" x y z
        end
    end

    ilscale = cbrt(length(particles) / (OPT_PART_PER_BLOCK * dx * dy * dz))
    nx, ny, nz = floor.(Int32, (dx, dy, dz) .* ilscale .+ 1)
    raw_con = RawContainerPoly(
        x_min, x_max, y_min, y_max, z_min, z_max, nx, ny, nz, px, py, pz, ppb,
    )
    con = Container(raw_con, ordering)

    @inbounds for p in particles
        add_point!(con, p)
    end
    return con
end

function read_particles!(
    con::Container{<:RawContainer},
    input::Union{IO,AbstractString},
)
    ((ax, ay, az), (bx, by, bz)) = bounding_box(con)
    px, py, pz = periodicity(con)

    particles = Particle{Nothing}[]
    for ln in eachline(input)
        spl = eachsplit(ln)
        id, x, y, z = parse.((Int32, Float64, Float64, Float64), spl)
        if (px || (ax <= x <= bx)) && (py || (ay <= y <= by)) && (pz || (az <= z <= bz))
            push!(particles, particle(id, (x, y, z)))
        else
            @warn "Particle out of bounds" x y z
        end
    end

    @inbounds for p in particles
        add_point!(con, p)
    end
    return con
end

function read_particles!(
    con::Container{<:RawContainerPoly},
    input::Union{IO,AbstractString},
)
    ((ax, ay, az), (bx, by, bz)) = bounding_box(con)
    px, py, pz = periodicity(con)

    particles = Particle{Float64}[]
    for ln in eachline(input)
        spl = eachsplit(ln)
        id, x, y, z, r = parse.((Int32, Float64, Float64, Float64, Float64), spl)
        if (px || (ax <= x <= bx)) && (py || (ay <= y <= by)) && (pz || (az <= z <= bz))
            push!(particles, particle(id, (x, y, z), r))
        else
            @warn "Particle out of bounds" x y z
        end
    end

    @inbounds for p in particles
        add_point!(con, p)
    end
    return con
end
