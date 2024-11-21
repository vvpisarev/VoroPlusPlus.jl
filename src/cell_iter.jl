function Base.iterate(con::Container)
    itr_state = ContainerIterator(con)
    !start!(itr_state) && return nothing
    cell = VoronoiCell(con)
    while !compute_cell!(cell, con, itr_state)
        next!(itr_state) || return nothing
    end
    return cell, itr_state
end

function Base.iterate(con::Container, itr_state::ContainerIterator)
    !next!(itr_state) && return nothing
    cell = VoronoiCell(con)
    while !compute_cell!(cell, con, itr_state)
        next!(itr_state) || return nothing
    end
    return cell, itr_state
end
