struct Particle{T<:Union{Nothing,Float64}}
    id::Int32
    pos::SVector{3,Float64}
    radius::T
end

function particle(id::Integer, pos; radius::Union{Nothing,Float64}=nothing)
    return Particle(id, pos, radius)
end

