
#!/bin/bash

#=
    exec julia --project --threads=$1 -e 'using VoroPlusPlus; include("parallel/voro_parallel_2.jl")'
=#