"""
    add_wall!(con::AbstractContainer, w::AbstractWall)

Add wall `w` to container `con`.
"""
function add_wall!(con::AbstractContainer, w::AbstractWall)
    __cxxwrap_add_wall!(__raw(con), w)
    return con
end

function point_inside(pt, w::AbstractWall)
    x, y, z = pt
    xx, yy, zz = Float64.((x, y, z))
    return convert(Bool, __cxxwrap_point_inside(w, xx, yy, zz))
end

"""
    point_inside_walls(pt, con::AbstractContainer)

Check if point `pt` is inside all walls added to `con`.
"""
function point_inside_walls(pt, con::AbstractContainer)
    x, y, z = pt
    xx, yy, zz = Float64.((x, y, z))
    return convert(Bool, __cxxwrap_point_inside_walls(__raw(con), xx, yy, zz))
end

"""
    wall_plane((a, b, c), d[, wall_id=-99])

Define a planar wall described by the equation `ax + by + cz = d`. The inside of the wall is
    defined as `ax + by + cz < d`.
"""
function wall_plane(nrm, distance::Real)
    nx, ny, nz = nrm
    a, b, c = Float64.((nx, ny, nz))
    d = Float64(distance)
    return WallPlane(a, b, c, d)
end

function wall_plane(nrm, distance::Real, wid::Integer)
    nx, ny, nz = nrm
    a, b, c = Float64.((nx, ny, nz))
    d = Float64(distance)
    return WallPlane(a, b, c, d, convert(Int32, wid))
end

function wall_sphere(center, radius::Real)
    x, y, z = center
    xc, yc, zc = Float64.((x, y, z))
    r = Float64(radius)
    return WallSphere(xc, yc, zc, r)
end

function wall_sphere(center, radius::Real, wid::Integer)
    x, y, z = center
    xc, yc, zc = Float64.((x, y, z))
    r = Float64(radius)
    return WallSphere(xc, yc, zc, r, wid)
end

function wall_cylinder(origin, axis, radius::Real)
    x, y, z = origin
    xc, yc, zc = Float64.((x, y, z))
    ax, ay, az = axis
    xa, ya, za = Float64.((ax, ay, az))
    r = Float64(radius)
    return WallSphere(xc, yc, zc, xa, ya, za, r)
end

function wall_cylinder(origin, axis, radius::Real, wid::Integer)
    x, y, z = origin
    xc, yc, zc = Float64.((x, y, z))
    ax, ay, az = axis
    xa, ya, za = Float64.((ax, ay, az))
    r = Float64(radius)
    return WallSphere(xc, yc, zc, xa, ya, za, r, wid)
end

function wall_cone(apex, axis, angle::Real)
    x, y, z = apex
    xc, yc, zc = Float64.((x, y, z))
    ax, ay, az = axis
    xa, ya, za = Float64.((ax, ay, az))
    ang = Float64(angle)
    return WallCone(xc, yc, zc, xa, ya, za, ang)
end

function wall_cone(apex, axis, angle::Real, wid::Integer)
    x, y, z = apex
    xc, yc, zc = Float64.((x, y, z))
    ax, ay, az = axis
    xa, ya, za = Float64.((ax, ay, az))
    ang = Float64(angle)
    return WallCone(xc, yc, zc, xa, ya, za, ang, wid)
end
