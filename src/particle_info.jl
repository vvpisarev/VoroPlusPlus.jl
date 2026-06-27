struct Particle{T<:Union{Nothing,Float64}}
    id::Int32
    pos::SVector{3,Float64}
    radius::T
end

function particle(id::Integer, pos; radius::Union{Nothing,Real}=nothing)
    return particle(id, pos, radius)
end

function particle(id::Integer, pos, ::Nothing)
    return Particle{Nothing}(id, pos, nothing)
end

function particle(id::Integer, pos, radius::Real)
    return Particle{Float64}(id, pos, radius)
end

function Particle(::Container, (id, x, y, z))
    return Particle{Nothing}(id, (x, y, z), nothing)
end

function Particle(::ContainerPoly, (id, x, y, z, r))
    return Particle{Float64}(id, (x, y, z), r)
end

function Particle(::Container, pinfo::CxxParticleInfo)
    (; id, x, y, z) = pinfo
    return Particle{Nothing}(id, (x, y, z), nothing)
end

function Particle(::ContainerPoly, pinfo::CxxParticleInfo)
    (; id, x, y, z, r) = pinfo
    return Particle{Float64}(id, (x, y, z), r)
end

Particle(con::Tessellation, itr...) = Particle(__raw(con), itr...)

@propagate_inbounds function add_point!(con::Tessellation{<:Container}, p::Particle)
    add_point!(con, p.id, p.pos)
end

@propagate_inbounds function add_point!(
    con::Tessellation{<:ContainerPoly}, p::Particle{Float64},
)
    add_point!(con, p.id, p.pos, p.radius)
end

particle_type(::Type{C}) where {C<:Tessellation} = particle_type(__raw_type(C))

particle_type(::Type{<:Container}) = Particle{Nothing}

particle_type(::Type{<:ContainerPoly}) = Particle{Float64}

particle_type(::C) where {C<:Tessellation} = particle_type(C)
