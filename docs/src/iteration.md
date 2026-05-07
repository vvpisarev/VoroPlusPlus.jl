# Iteration

Iteration over a `Container` gives pairs `particle => cell` where `particle` is a struct
    holding integer ID, position and (optionally) radius, `cell` is a `CheckedVoronoiCell`. 

## Data sharing between iterations

`CheckedVoronoiCell` is a mutable container. Upon iteration, the same object is reused 
    throughout a container iteration. So, the actual Voronoi cell object returned by the 
    iteration does not change identity. The delayed computation of properties will result 
    in the property computed for an object state not corresponding to the iteration it 
    has been saved.

### Example
```julia
vol = 0.0

for (p, vc) in container
    vol += volume(vc)
end
```
In this example, `vc` will hold a reference to the same `VoronoiCell` object which will be 
    modified during iteration.

### Getting a copy of a Voronoi cell

If an independent copy of a Voronoi cell is required, there are several options:
1. `copy(cell)` will create an independent copy of `VoronoiCell` or `CheckedVoronoiCell`.
2. `take!(cell)` gets a `VoronoiCell` object from `CheckedVoronoiCell` writing an 
    uninitialized cell into the latter
3. `convert(VoronoiCell, ::CheckedVoronoiCell)` returns a fresh copy, so that `push!`ing
    or assigning a `CheckedVoronoiCell` to an `Array{VoronoiCell}` would automatically 
    create an independent copy of a cell object

That is rarely needed in a typical workflow but may be used if the tessellation analysis
    required comparison between Voronoi cells or other operations involving multiple 
    cells at once.

### Example

1. `[copy(vc) for (p, vc) in container]` is a `Vector{CheckedVoronoiCell}` with independent 
    items
2. `[take!(vc) for (p, vc) in container]` is a `Vector{VoronoiCellAllocated}` with independent 
    items
3. `VoronoiCell[vc for (p, vc) in container]` should contain independent items

## `eachparticle` and `eachcell`

In some cases, you don't need either cell information or particle information. Then, use 
    `eachparticle` and `eachcell` functions.

### Example
```julia
vol = 0.0

for vc in eachcell(container)
    vol += volume(vc)
end
```