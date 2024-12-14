
###########   Containter Periodic Poly (conprdply)   ##############

#=* The class constructor sets up the geometry of container.
 * (bx_) The x coordinate of the first unit vector.
 * (bxy_,by_) The x and y coordinates of the second unit vector.
 * (bxz_,byz_,bz_) The x, y, and z coordinates of the third unit vector.
 * (nx_,ny_,nz_) the number of grid blocks in each of the three coordinate directions.
 * init_mem_ the initial memory allocation for each block. 

    double bx_,double bxy_,double by_,double bxz_,double byz_,double bz_,
    int nx_,int ny_,int nz_, 
    int init_mem_
 =#

#function Container_Periodic_Poly(bx_::Float64, bxy_::Float64, by_::Float64, bxz_::Float64, byz_::Float64, bz_::Float64, nx_::Int32, ny_::Int32, nz_::Int32, init_mem_::Int32)
#
#    return Container_Periodic_Poly(bx_, bxy_, by_, bxz_, byz_, bz_, nx_, ny_, nz_, init_mem_)
#    
#end

function Container_Periodic_Poly(
    ;
    v_data::NTuple{6,Real},
    nblocks::NTuple{3,Integer},
    initial_memory::Integer=1,
)
    bx_, bxy_, by_, bxz_, byz_, bz_ = Float64.(v_data)
    nx_, ny_, nz_ = Int32.(nblocks)
    init_mem_ = Int32(initial_memory)

    return Container_Periodic_Poly(
        bx_, bxy_, by_, bxz_, byz_, bz_, nx_, ny_, nz_, init_mem_,
    )
end

#=
Base.@propagate_inbounds function conprdply_add_point!(conprdply::Container_Periodic_Poly, id::Integer, pt)
    @boundscheck if length(pt) != 4
        throw(ArgumentError("Can only add 3-dimensional points and radius to a VoroPlusPlus Container_Periodic_Poly"))
    end
    x, y, z, r = pt
    return conprdply_add_point!(conprdply, Int32(id), Float64.((x, y, z, r))...)
end
=#


function conprdply_compute_ghost_cell(con::Container_Periodic_Poly, c::VoronoiCell, dx::Float64, dy::Float64, dz::Float64, dr::Float64 )

    
    return ccall(
        (:compute_ghost_cell_conprdply, "libvoro++wrap"), Bool,
        (Ptr{Cvoid}, Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble),
        con.cpp_object, c.cpp_object, dx, dy, dz, dr,)
end