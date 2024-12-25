using VoroPlusPlus

using Printf: Format, format
using Test

foreach(rm, filter(endswith(".gnu"), readdir()))
foreach(rm, filter(endswith(".txt"), readdir()))
foreach(rm, filter(endswith(".pov"), readdir()))
foreach(rm, filter(endswith(".custom*"), readdir()))
foreach(rm, filter(endswith(".points"), readdir()))
foreach(rm, filter(endswith(".vec"), readdir()))
foreach(rm, filter(endswith(".vol"), readdir()))
include("basic.jl")
include("custom.jl")
include("degenerate.jl")
include("extra.jl")
include("interface.jl")
include("no_release.jl")
