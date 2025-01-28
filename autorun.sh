
#!/bin/bash

#=
    exec julia --project -t $1 -e 'using VoroPlusPlus; include("./src/parallel/test/voro_par_test.jl"); @btime sum(volume, tessellation; init=0.0)'
=#