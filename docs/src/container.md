# Container

Container is a box in which the tessellation is performed.

Voro++ supports periodic and non-periodic mono- and polydisperse containers.

VoroPlusPlus.jl type `Container` encapsulates `container` or `container_poly` object 
from Voro++ with an ordering object.

## Creating a container

### Allocate an empty container
```@docs
VoroPlusPlus.container

VoroPlusPlus.polydisperse_container
```

### Read particle data from file
```@docs
read_particles

read_polydisperse_particles
```

## Particle ordering

Ordering defines in which order the particles and their Voronoi cells are given upon
    iteration over a container.

`UnspecifiedOrder()` means that the order of insertion is not preserved upon iteration.

`InsertionOrder()` means that the iteration produces particles in the same order they've
    been inserted. This type wraps `voro::particle_order` class.

With Voro++, it's possible to insert some particles with order tracking and some without
    it. In Julia, this is not possible through public interface.

```@docs
VoroPlusPlus.InsertionOrder

VoroPlusPlus.UnspecifiedOrder
```

## Geometric properties of a container

```@docs
bounding_box

periodicity
```

## Adding particles to the tessellation

```@docs
add_point!
```

## Output data to a file

```@docs
draw_gnuplot(::Any, ::Container)
```