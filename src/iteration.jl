"""
    __container_iterator(con::AbstractRawContainer, ord::ContainerIterationOrder)

Return an iterator over the container according to the ordering represented by `ord`.
"""
function __container_iterator(con::AbstractContainer, ::UnspecifiedOrder)
    return ContainerIterator(con)
end

function __container_iterator(con::AbstractContainer, ord::InsertionOrder)
    return InsertionOrderIterator(con, ord)
end

"""
    __start!(itr_state)

Start the iteration. Return `true` if there is the next iteration element, `false` if
    container is empty.
"""
function __start!(itr_state)
    return Bool(__cxxwrap_start!(itr_state))
end

"""
    __next!(itr_state)

Propagate the iteration state one step forward. Return `true` if there is the next iteration
    element, `false` if iterator is at end.
"""
function __next!(itr_state)
    return Bool(__cxxwrap_inc!(itr_state))
end

function Base.iterate(con::Tessellation)
    itor = __container_iterator(con.con, con.ord)
    !__start!(itor) && return nothing
    cell = CheckedVoronoiCell(con.con, itor)
    particle_tup = __cxxwrap_particle_info(itor)
    particle = Particle(con.con, particle_tup)
    return (particle => cell), (__next!(itor), itor, cell)
end

function Base.iterate(con::Tessellation, (hasnext, itor, cell))
    hasnext || return nothing
    cell = compute_cell!(cell, con, itor)
    particle_tup = __cxxwrap_particle_info(itor)
    particle = Particle(con.con, particle_tup)
    return (particle => cell), (__next!(itor), itor, cell)
end

struct EachCell{C<:Tessellation}
    con::C
end

struct EachParticle{C<:Tessellation}
    con::C
end

"""
    eachcell(con::AbstractContainer)

Return an iterable which produces only the Voronoi cells upon iteration.
"""
function eachcell(con::Tessellation)
    return EachCell(con)
end

"""
    eachcell(con::AbstractContainer)

Return an iterable which produces only the particles data upon iteration.
"""
function eachparticle(con::Tessellation)
    return EachParticle(con)
end

function Base.iterate(con::EachCell)
    raw_con = __raw(con.con)
    ord = ordering(con.con)
    itor = __container_iterator(raw_con, ord)
    !__start!(itor) && return nothing
    cell = CheckedVoronoiCell(raw_con, itor)
    return cell, (__next!(itor), itor, cell)
end

function Base.iterate(con::EachCell, (hasnext, itor, cell))
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
    particle_tup = __cxxwrap_particle_info(itor)
    particle = Particle(raw_con, particle_tup)
    return particle, (__next!(itor), itor)
end

function Base.iterate(con::EachParticle, (hasnext, itor))
    hasnext || return nothing
    raw_con = __raw(con.con)
    particle_tup =__cxxwrap_particle_info(itor)
    particle = Particle(raw_con, particle_tup)
    return particle, (__next!(itor), itor)
end

function Base.isdone(::Union{Tessellation,EachCell,EachParticle}, (hasnext,))
    return !hasnext
end

Base.eltype(::Type{<:EachCell}) = CheckedVoronoiCell

Base.eltype(::Type{<:EachParticle{C}}) where {C} = particle_type(C)

function Base.eltype(::Type{C}) where {C<:Tessellation}
    PT = particle_type(C)
    return Pair{PT, CheckedVoronoiCell}
end

Base.IteratorSize(::Type{<:Tessellation}) = Base.SizeUnknown()

Base.IteratorSize(::Type{<:EachParticle}) = Base.SizeUnknown()

Base.IteratorSize(::Type{<:EachCell}) = Base.SizeUnknown()
