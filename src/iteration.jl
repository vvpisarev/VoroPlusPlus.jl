function __container_iterator(con::AbstractRawContainer, ::UnspecifiedOrder)
    return RawContainerIterator(con)
end

function __container_iterator(con::AbstractRawContainer, ord::InsertionOrder)
    return InsertionOrderIterator(con, ord)
end

function __start!(itr_state)
    return Bool(__cxxwrap_start!(itr_state))
end

function __next!(itr_state)
    return Bool(__cxxwrap_inc!(itr_state))
end

function Base.iterate(con::Container)
    itor = __container_iterator(con.con, con.ord)
    !__start!(itor) && return nothing
    cell = CheckedVoronoiCell(con, itor)
    particle_tup = particle_info(itor)
    particle = Particle(con, particle_tup)
    return (particle => cell), (__next!(itor), itor)
end

function Base.iterate(con::Container, (hasnext, itor))
    hasnext || return nothing
    cell = CheckedVoronoiCell(con, itor)
    particle_tup = particle_info(itor)
    particle = Particle(con, particle_tup)
    return (particle => cell), (__next!(itor), itor)
end

struct Unsafe{C<:AbstractContainer}<:AbstractContainer
    container::C

    Unsafe{C}(con::C) where {C<:AbstractContainer} = new{C}(con)
end

function Unsafe(con::Unsafe)
    return Unsafe(con.container)
end

function Unsafe(con::C) where {C<:AbstractContainer}
    return Unsafe{C}(con)
end

__raw(con::Unsafe) = __raw(con.container)
__raw_type(::Type{Unsafe{C}}) where {C} = __raw_type(C)

ordering(con::Unsafe) = ordering(con.container)

function Base.iterate(ucon::Unsafe)
    con = ucon.container
    itor = __container_iterator(__raw(con), ordering(con))
    !__start!(itor) && return nothing
    cell = CheckedVoronoiCell(con, itor)
    particle_tup = particle_info(itor)
    particle = Particle(con, particle_tup)
    return (particle => cell), (__next!(itor), itor, cell)
end

function Base.iterate(ucon::Unsafe, (hasnext, itor, cell))
    con = ucon.container
    hasnext || return nothing
    cell = compute_cell!(cell, con, itor)
    particle_tup = particle_info(itor)
    particle = Particle(con, particle_tup)
    return (particle => cell), (__next!(itor), itor, cell)
end

struct EachCell{C<:AbstractContainer}
    con::C
end

struct EachParticle{C<:AbstractContainer}
    con::C
end

function eachcell(con::AbstractContainer)
    return EachCell(con)
end

function eachparticle(con::AbstractContainer)
    return EachParticle(con)
end

function Base.iterate(con::EachCell)
    raw_con = __raw(con.con)
    ord = ordering(con.con)
    itor = __container_iterator(raw_con, ord)
    !__start!(itor) && return nothing
    cell = CheckedVoronoiCell(raw_con, itor)
    return cell, (__next!(itor), itor)
end

function Base.iterate(con::EachCell, (hasnext, itor))
    hasnext || return nothing
    raw_con = __raw(con.con)
    cell = CheckedVoronoiCell(raw_con, itor)
    return cell, (__next!(itor), itor)
end

function Base.iterate(con::EachCell{<:Unsafe})
    raw_con = __raw(con.con)
    ord = ordering(con.con)
    itor = __container_iterator(raw_con, ord)
    !__start!(itor) && return nothing
    cell = CheckedVoronoiCell(raw_con, itor)
    return cell, (__next!(itor), itor, cell)
end

function Base.iterate(con::EachCell{<:Unsafe}, (hasnext, itor, cell))
    hasnext || return nothing
    raw_con = __raw(con.con)
    cell = compute_cell!(cell, raw_con, itor)
    return cell, (__next!(itor), itor, cell)
end

function Base.iterate(con::EachParticle)
    raw_con = __raw(con.con)
    ord = ordering(con.con)
    itor = __container_iterator(raw_con, ord)
    !__start!(itor) && return nothing
    particle_tup = particle_info(itor)
    particle = Particle(raw_con, particle_tup)
    return particle, (__next!(itor), itor)
end

function Base.iterate(con::EachParticle, (hasnext, itor))
    hasnext || return nothing
    raw_con = __raw(con.con)
    particle_tup = particle_info(itor)
    particle = Particle(raw_con, particle_tup)
    return particle, (__next!(itor), itor)
end

function Base.isdone(::Union{AbstractContainer,EachCell,EachParticle}, (hasnext,))
    return !hasnext
end

Base.eltype(::Type{<:EachCell}) = CheckedVoronoiCell

Base.eltype(::Type{<:EachParticle{C}}) where {C} = particle_type(C)

function Base.eltype(::Type{C}) where {C<:AbstractContainer}
    PT = particle_type(C)
    return Pair{PT, CheckedVoronoiCell}
end

Base.IteratorSize(::Type{<:AbstractContainer}) = Base.SizeUnknown()

Base.IteratorSize(::Type{<:EachParticle}) = Base.SizeUnknown()

Base.IteratorSize(::Type{<:EachCell}) = Base.SizeUnknown()
