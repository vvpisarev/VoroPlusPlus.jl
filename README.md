# VoroPlusPlus.jl
Julia bindings to Voro++ library

## Installing C++ wrapper library

Download and build the C++ side of this package (https://github.com/vvpisarev/voropp-cxxwrap-julia)
    according to the instructions in README.

## Configuration

Before the first use, configure the location of C++ wrapper:
```julia
julia> using VoroPlusPlus

julia> VoroPlusPlus.set_wrapper_path("/path/to/libvoro++wrap/lib")
```
Then, restart Julia and the bindings should work in the new session.