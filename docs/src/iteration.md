# Iteration

Iteration over a `Container` gives pairs `particle => cell` where `particle` is a struct
    holding integer ID, position and (optionally) radius, `cell` is a `CheckedVoronoiCell`. 

## `Unsafe` wrapper

By default, a fresh cell is allocated on each iteration, which may be suboptimal in terms 
    of efficiency. If you want to reuse the same `VoronoiCell` object through iteration, 
    there is wrapper `Unsafe{<:AbstractContainer}`.

### Example
```julia
vol = 0.0

for (p, vc) in Unsafe(con)
    vol += volume(vc)
end
```
In this example, `vc` will hold a reference to the same `VoronoiCell` object which will be 
    modified during iteration.

## `eachparticle` and `eachcell`

In some cases, you don't need either cell information or particle information. Then, use 
    `eachparticle` and `eachcell` functions.

### Example
```julia
vol = 0.0

for vc in eachcell(Unsafe(con))
    vol += volume(vc)
end
```