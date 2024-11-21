using VoroPlusPlus

using Printf: Format, format
using Test

foreach(rm, filter(endswith(".gnu"), readdir()))
include("basic.jl")
