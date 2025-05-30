using VoroPlusPlus

using CxxWrap
using Printf: Format, format
using Random
using Test

foreach(rm, filter(endswith(".gnu"), readdir()))
foreach(rm, filter(endswith(".txt"), readdir()))
foreach(rm, filter(endswith(".pov"), readdir()))
foreach(rm, filter(endswith(".custom1"), readdir()))
foreach(rm, filter(endswith(".custom2"), readdir()))
foreach(rm, filter(endswith(".custom3"), readdir()))
foreach(rm, filter(endswith(".points"), readdir()))
foreach(rm, filter(endswith(".vec"), readdir()))
foreach(rm, filter(endswith(".vol"), readdir()))
#Implemented
include("basic.jl")
#include("custom.jl")
#include("degenerate.jl")
#include("extra.jl")

# Need to check
#include("interface.jl")
#include("no_release.jl")

# Under Testing refactoring
#include("centroid.jl")
#include("particle_info.jl")
