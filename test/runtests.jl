using VoroPlusPlus

using CxxWrap
using Printf: Format, format
using Random
using StaticArrays
using Test

foreach(
    rm,
    path for path in readdir()
    if endswith(path, r".gnu|.txt|.pov|.custom1|.custom2|.custom3|.points|.vec|.vol")
)

#Implemented
include("basic.jl")
include("cell_statistics.jl")
#include("custom.jl")
#include("degenerate.jl")
#include("extra.jl")

# Need to check
#include("interface.jl")
#include("no_release.jl")

# Under Testing refactoring
#include("particle_info.jl")
