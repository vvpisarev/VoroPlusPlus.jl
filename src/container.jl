"""
    Container(; bounds, nblocks, periodic=(false, false, false), particles_per_block=8)

Allocate space for a container of Voronoi cells.
# Keywords
* `bounds`: limits of the bounding box `(xmin, xmax, ymin, ymax, zmin, zmax)`
* `nblocks`: numbers of computation blocks along each axis
* `periodic::NTuple{3,Bool}`: periodicity in each axis. Default: `(false, false, false)`
* `particles_per_block::Integer`: initially allocate memory for this many particles
    per block. Default: 8
"""
function Container(
    ;
    bounds::NTuple{6,Real},
    nblocks::NTuple{3,Integer},
    periodic::NTuple{3,Bool}=(false, false, false),
    particles_per_block::Integer=8,
)
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.(bounds)
    n_x, n_y, n_z = Int32.(nblocks)
    px, py, pz = periodic
    ppb = Int32(particles_per_block)

    return Container(
        x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, px, py, pz, ppb,
    )
end

Base.@propagate_inbounds function add_point!(con::Container, id::Integer, pt)
    @boundscheck if length(pt) != 3
        throw(ArgumentError("Can only add 3-dimensional points to a VoroPlusPlus Container"))
    end
    x, y, z = pt
    return add_point!(con, Int32(id), Float64.((x, y, z))...)
end

#############################################################################

function __get_particle_pos(cl::Container_Iterator)
   
    ccall(
        (:get_pos, "libvoro++wrap"),
        Particle_Info,
        (Ptr{Cvoid},),
        cl.cpp_object,
    )

end

function __find_voro_cell(cn::Container, coord = (0.0, 0.0, 0.0))

    _x, _y, _z = coord
    dx, dy, dz = Float64.((_x, _y, _z))
   
    ccall(
        (:find_voro_cell, "libvoro++wrap"),
        Fnd_Voro_Cell,
        (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble),
        cn.cpp_object, dx, dy, dz,
    )

end

###########################################################################


###########   Container_Poly

function Container_Poly(
    ;
    bounds::NTuple{6,Real},
    nblocks::NTuple{3,Integer},
    periodic::NTuple{3,Bool}=(false, false, false),
    particles_per_block::Integer=8,
)
    x_min, x_max, y_min, y_max, z_min, z_max = Float64.(bounds)
    n_x, n_y, n_z = Int32.(nblocks)
    px, py, pz = periodic
    ppb = Int32(particles_per_block)

    return Container_Poly(
        x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z, px, py, pz, ppb,
    )
end

Base.@propagate_inbounds function conp_add_point!(con::Container_Poly, id::Integer, pt)
    @boundscheck if length(pt) != 4
        throw(ArgumentError("Can only add 3-dimensional points and radius to a VoroPlusPlus Container_Poly"))
    end
    x, y, z, r = pt
    return conp_add_point!(con, Int32(id), Float64.((x, y, z, r))...)
end