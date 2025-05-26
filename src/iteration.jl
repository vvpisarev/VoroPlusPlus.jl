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
    cell = VoronoiCell(con, itor)
    particle_tup = particle_info(itor)
    particle = Particle(con, particle_tup)
    return (particle => cell), itor
end

function Base.iterate(con::Container, itor)
    !__next!(itor) && return nothing
    cell = VoronoiCell(con, itor)
    particle_tup = particle_info(itor)
    particle = Particle(con, particle_tup)
    return (particle => cell), itor
end

struct Unsafe{T}
    container::T
end

function Base.iterate(ucon::Unsafe{<:AbstractContainer})
    con = ucon.container
    itor = __container_iterator(con.con, con.ord)
    !__start!(itor) && return nothing
    cell = VoronoiCell(con, itor)
    particle_tup = particle_info(itor)
    particle = Particle(con, particle_tup)
    return (particle => cell), (itor, cell)
end

function Base.iterate(ucon::Unsafe{<:AbstractContainer}, (itor, cell))
    con = ucon.container
    !__next!(itor) && return nothing
    compute_cell!(cell, con.con, itor)
    particle_tup = particle_info(itor)
    particle = Particle(con, particle_tup)
    return (particle => cell), (itor, cell)
end
