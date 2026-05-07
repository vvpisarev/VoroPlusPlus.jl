# Voronoi cells

Typically, the users should not need to instantiate cell objects directly and instead use 
    ones produced by container iteration.

For completeness, the constructor `VoronoiCell` is exposed. It is a wrapper for `voro::voronoicell_neighbor`
    class

```@docs
VoroPlusPlus.VoronoiCell
```

## `CheckedVoronoiCell`

Iteration over a container produces `CheckedVoronoiCell` objects which contain a reference to 
    Voronoi cell object and a boolean flag showing if the cell computation was successful.
    The public API for `CheckedVoronoiCell` is the same as for `VoronoiCell`, but the latter 
    throws an error or returns a default value (if it makes sense) if the flag is `false`.

The object without the boolean flag is obtained as `checked_cell[]` without transfering the 
    ownership or as `take!(checked_cell)` with transfering the ownership.

## Conditional computation of properties

To avoid explicit validity flag check, there is an utility function `if_valid` which allows 
    computations of properties for underlying Voronoi cell if the validity flag is `true` 
    and returning a default value or triggering a fallback function otherwise.

For most functions from this package, the default is already implemented, so that the 
    interface is shared between `VoronoiCell` and `CheckedVoronoiCell`.

```@docs
VoroPlusPlus.if_valid
```

## Initialize a cell with a given shape

```@docs
voronoicell_box

voronoicell_tetrahedron

voronoicell_octahedron
```

## Reshape a Voronoi cell

```@docs
translate!

reset_to_box!

reset_to_tetrahedron!

reset_to_octahedron!

cut_by_particle_position!
```

## Properties of Voronoi cells
```@docs
volume

number_of_faces

number_of_edges

number_of_vertices

centroid

vertex_positions

vertex_positions!

get_normals!

normals
```

## Cell output

```@docs
draw_gnuplot(::IO, ::VoronoiCell)
```

