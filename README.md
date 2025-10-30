# VoroPlusPlus.jl
Julia bindings to [Voro++ library](https://math.lbl.gov/voro++/).

From Voro++ reference manual:
> Voro++ is a software library for carrying out three-dimensional computations of the Voronoi tessellation. 
>    A distinguishing feature of the Voro++ library is that it carries out cell-based calculations, 
>    computing the Voronoi cell for each particle individually, rather than computing the Voronoi 
>    tessellation as a global network of vertices and edges. It is particularly well-suited for applications 
>    that rely on cell-based statistics, where features of Voronoi cells (eg. volume, centroid, number 
>    of faces) can be used to analyze a system of particles.

## Installation

This package depends on the binary of Voro++ library and the `CxxWrap` wrapper. The source code for the
    wrapper is available at https://github.com/vvpisarev/voropp-cxxwrap-julia. The pre-built binaries 
    for Linux-glibc platforms including both libraries are available as the package `voropp_wrapper_jll`.
    For other platforms, manual build of Voro++ and the wrapper library is required (see build instructions
    in the wrapper repo).

This package requires Julia 1.10 or later.

### Installation from a personal package registry

```julia
julia> using Pkg

julia> pkg"registry add https://gitlab.com/pisarevvv/samma-registry.git"

julia> add VoroPlusPlus
```

### Add package by Github link

```julia
julia> using Pkg

julia> pkg"add https://github.com/vvpisarev/voropp_wrapper_jll.jl#voropp_wrapper-v0.1.1+0"

julia> pkg"add https://github.com/vvpisarev/VoroPlusPlus.jl"
```

## Configuration (only for manual wrapper builds)

Before the first use, configure the location of C++ wrapper:
```julia
julia> using VoroPlusPlus

julia> VoroPlusPlus.set_wrapper_path("/path/to/libvoro++wrap/lib")
```
Then, restart Julia and the bindings should work in the new session.