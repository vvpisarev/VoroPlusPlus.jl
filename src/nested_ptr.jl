"""
    UnsafeNested{Ptr}

Wrapper for nested pointers to allow index access. The indices are 1-based, i.e. `p[k]`
    in Julia corresponds to `p[k-1]` in C/C++
"""
struct UnsafeIndexable{P<:Ptr}
    ptr::P
end

function UnsafeIndexable(ptr::CxxPtr{T}) where {T}
    return UnsafeIndexable(reinterpret(Ptr{T}, ptr))
end

function Base.getindex(p::UnsafeIndexable, ind::Integer)
    return unsafe_load(p.ptr, ind)
end

function Base.getindex(
    p::UnsafeIndexable{Ptr{T}}, i::Integer, inds...
) where {T<:Union{Ptr,CxxPtr}}
    next_ptr = unsafe_load(p.ptr, i)
    next = UnsafeIndexable(next_ptr)
    return next[inds...]
end

function Base.setindex!(p::UnsafeIndexable, x, ind::Integer)
    return unsafe_store!(p.ptr, x, ind)
end

function Base.setindex!(
    p::UnsafeIndexable{Ptr{T}}, x, i::Integer, inds...
) where {T<:Union{Ptr,CxxPtr}}
    next_ptr = unsafe_load(p.ptr, i)
    next = UnsafeIndexable(next_ptr)
    return next[inds...] = x
end
