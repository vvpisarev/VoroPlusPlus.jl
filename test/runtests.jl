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

foreach(
    rm,
    joinpath("./image_data", path) for path in readdir("./image_data")
    if endswith(path, r".gnu|.txt|pts.pov|cells.pov|domain.pov")
)

#Implemented
include("basic.jl")
include("cell_statistics.jl")
include("iteration.jl")
include("drawing.jl")
#include("custom.jl")
#include("degenerate.jl")
#include("extra.jl")

# Need to check
#include("interface.jl")
#include("no_release.jl")

# Under Testing refactoring
#include("particle_info.jl")
