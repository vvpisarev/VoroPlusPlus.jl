function Base.iterate(con::Container)
    itr_state = Container_Iterator(con)
    !start!(itr_state) && return nothing
    cell = VoronoiCell(con)
    while !compute_cell!(cell, con, itr_state)
        next!(itr_state) || return nothing
    end
    return cell, itr_state
end

function Base.iterate(con::Container, itr_state::Container_Iterator)
    !next!(itr_state) && return nothing
    cell = VoronoiCell(con)
    while !compute_cell!(cell, con, itr_state)
        next!(itr_state) || return nothing
    end
    return cell, itr_state
end

struct Unsafe{T}
    container::T
end

function Base.iterate(ucon::Unsafe{<:Container})
    con = ucon.container
    itr_state = Container_Iterator(con)
    !start!(itr_state) && return nothing
    cell = VoronoiCell(con)
    while !compute_cell!(cell, con, itr_state)
        next!(itr_state) || return nothing
    end
    return cell, (itr_state, cell)
end

function Base.iterate(con::Unsafe{<:Container}, (itr_state, cell))
    !next!(itr_state) && return nothing
    while !compute_cell!(cell, con.container, itr_state)
        next!(itr_state) || return nothing
    end
    return cell, (itr_state, cell)
end
