
function GenerateParticles!(nparticles::Integer, x_min::Real, x_max::Real, y_min::Real, y_max::Real, z_min::Real, z_max::Real)

    particles = Array{Real}([])

    for _ in 1:nparticles
        x = x_min + rand() * (x_max - x_min)
        y = y_min + rand() * (y_max - y_min)
        z = z_min + rand() * (z_max - z_min)
        push!(particles, x)
        push!(particles, y)
        push!(particles, z)
    end

    return particles

end

function voronoi_tessellation(con_dims, coords::AbstractArray{<:Real}, ids=1:div(length(coords), 3), ntasks::Integer=2)

    #container dimensions
    x_min, x_max, y_min, y_max, z_min, z_max = con_dims
    
    # container length in earch coordinate
    lx, ly, lz = Float64.((x_max - x_min, y_max - y_min, z_max - z_min))
    
    # container volume
    vol = lx * ly * lz
    
    # number of paticles
    npts = length(ids)
    
    # ?_average
    r_ave = cbrt(vol / npts)
    d_skin = 5 * r_ave
    
    con_dims = gen_container_dims(x_min, x_max, y_min, y_max, z_min, z_max, d_skin, ntasks)
    containers = GenerateContainers!(con_dims)
    owner = zeros(Int32, npts) # task that owns each particle
    
    tessellation = VoronoiTessellation(
        x_min, x_max, y_min, y_max, z_min, z_max,
        containers,
        skin_d,
        1,
        con_dims,
        owner
    )
    
    # x dimension is there are more than one task
    dx = lx / ntasks
    
    # particle
    p = 0

    # sort particles into bins

    #vector of vectors(number of task)
    workload = [Int32[] for _ in 1:ntasks]
    for ind in firstindex(coords):3:lastindex(coords)
        p += 1
        x, y, z = coords[ind], coords[ind+1], coords[ind+2]
        i_con::Int32, x_rel = fldmod1(x - x_min, dx)
        
        owner[p] = i_con
        if p > 1 && x_rel < d_skin
            push!(workload[i_con-1], p)
        elseif p < ntasks && dx - x_rel < d_skin # length coords instead ntasks
            push!(workload[i_con+1], p)
        end
        push!(workload[i_con], p)
    end
    
    Threads.@threads for i_con in 1:ntasks
        work = workload[i_con]
        con = tessellation.domain[i_con]
        for i_part in work
            ind = (i_part - 1) * 3
            x, y, z = coords[ind+1], coords[ind+2], coords[ind+3]
            add_point!(con, i_part, x, y, z)
        end
    end
    
    return tessellation
end


#p = GenerateParticles!(10, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0)
#println(p)