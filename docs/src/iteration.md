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
## Region-restricted iteration

To iterate only over particles located inside a spherical or box-shaped region, use
    `sphere_region` and `box_region`. They return a `VoroPlusPlus.Subset`: a stateful
    iterator over a *region* of the container which yields the same
    `particle => cell` pairs as iteration over the whole container and supports
    `eachparticle` and `eachcell` in the same way.

```julia
for (particle, cell) in sphere_region(container, (0.5, 0.5, 0.5), 0.1)
    # cells of particles within 0.1 of (0.5, 0.5, 0.5)
end

s = box_region(container, (0.0, 0.0, 0.0), (0.2, 0.2, 1.0))
ids = [particle.id for particle in eachparticle(s)]
VoroPlusPlus.reset_to_sphere!(s, (0.9, 0.1, 0.1), 0.05)   # reuse for another region
```

The region of an existing `Subset` is changed with `VoroPlusPlus.reset_to_sphere!` or
    `reset_to_box!` (the latter extends the cell reshaping function of the same name).
    Starting a new iteration pass re-applies the current region, so a `Subset` can be
    iterated repeatedly. `bounds_test=false` relaxes the region to the superset of
    particles in all computation blocks overlapping it.

### Conventions on periodic containers

* **Positions**: particles are reported at their stored (primary-domain) coordinates,
    never at the periodic image which satisfied the region. Use
    `VoroPlusPlus.minimum_image_distance` or `VoroPlusPlus.minimum_image_displacement`
    when computing distances from the query.
* **Cells**: cells are the genuine tessellation cells in the primary domain — a neighbor
    across a periodic face comes back with its own in-box cell.
* **Restrictions**: on periodic axes the region must be shorter than one period
    (`2r < L` for spheres, `upper - lower < L` for boxes), otherwise an `ArgumentError`
    is thrown. Non-periodic axes are clipped to the container bounds. Non-finite
    coordinates, non-positive radii, empty boxes and queries very far from the container
    (more than 1e7 box lengths) are rejected with `ArgumentError` before any C++ code is
    called.

!!! note "Upstream workaround"
    The underlying Voro++ `c_loop_subset` loop is broken for regions which start below
    the minimum face of a periodic axis (silent undercount, sometimes a crash), and its
    block indices are corrupted for cells computed through such a loop. `Subset`
    transparently re-expresses every query as a periodic image in a working position and
    computes cells from the true block indices, so results are exact. Raw use of
    `ContainerSubsetIterator` together with the raw
    `VoroPlusPlus.__cxxwrap_setup_sphere!` / `__cxxwrap_setup_box!` bindings bypasses
    these guards and remains unsafe on periodic containers.

```@docs
sphere_region

box_region

VoroPlusPlus.Subset

VoroPlusPlus.reset_to_sphere!

reset_to_box!(::VoroPlusPlus.Subset, ::Any, ::Any)

VoroPlusPlus.minimum_image_displacement

VoroPlusPlus.minimum_image_distance

eachparticle(::VoroPlusPlus.Subset)

eachcell(::VoroPlusPlus.Subset)
```
