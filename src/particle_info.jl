struct Particle{T<:Union{Nothing,Float64}}
    id::Int32
    pos::SVector{3,Float64}
    radius::T
end

function particle(id::Integer, pos; radius::Union{Nothing,Float64}=nothing)
    return Particle(id, pos, radius)
end

function Particle(con::Container{<:RawContainer}, (id, x, y, z))
    return Particle{Nothing}(id, (x, y, z), nothing)
end

function Particle(con::Container{<:RawContainerPoly}, (id, x, y, z, r))
    return Particle{Float64}(id, (x, y, z), r)
end
