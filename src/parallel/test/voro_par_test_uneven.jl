
# using Base.Threads

using Base.Threads
using BenchmarkTools
using Random

using VoroPlusPlus
#
# GenerateParticles!
#
# Parameters:
#
# Output: A ParticleCoords array of nparticles*3 elements
#
function GenerateParticles(
    nparticles::Integer, x_min::Real, x_max::Real, y_min::Real, y_max::Real, z_min::Real, z_max::Real
    ;
    rng = copy(Random.GLOBAL_RNG)
)

    particles = Float64[]

    for _ in 1:nparticles
        x = x_min + rand(rng) * (x_max - x_min)
        y = y_min + rand(rng) * (y_max - y_min)
        z = z_min + rand(rng) * (z_max - z_min)
        push!(particles, x, y, z)
    end

    return particles

end
#
#
function TestVoroTesselation(
    nparticles::Integer,
    x_min::Real,
    x_max::Real,
    y_min::Real,
    y_max::Real,
    z_min::Real,
    z_max::Real,
    imbalance::Real=2.0,
)

    # number tasks
    ntasks = Threads.nthreads()
    # particle's array
    nleft = floor(Int, nparticles / (1 + imbalance))
    nright = nparticles - nleft
    a_left = GenerateParticles(nleft, x_min, 0.5 * (x_min + x_max), y_min, y_max, z_min, z_max)
    a_right = GenerateParticles(nright, 0.5 * (x_min + x_max), x_max, y_min, y_max, z_min, z_max)
    a_p = append!(a_left, a_right)
    a_p_range = 1:div(length(a_p), 3)

    tessellation = voronoi_tessellation((x_min, x_max, y_min, y_max, z_min, z_max), a_p, a_p_range , ntasks)
    #println(typeof(tessellation))

    # Simple Test for Volumes
    cvol = (
        (tessellation.x_max - tessellation.x_min)*
        (tessellation.y_max - tessellation.y_min)*
        (tessellation.z_max - tessellation.z_min)
    )

    vvol = GetVoroVolume(tessellation)

    message = "Container volume: $(cvol)\n"*
                "Voronoi volume: $(vvol)\n"*
                "Difference: $(cvol - vvol)\n"*
                "Is approx?: $(isapprox(vvol, cvol; atol=1e-8))\n"

    println(message)

end

############################################################


############################################################
#
# SetGeneralParameters!
#
function SetGeneralParameters(
    nparticles::Integer,
    x_min::Real,
    x_max::Real,
    y_min::Real,
    y_max::Real,
    z_min::Real,
    z_max::Real,
    ntasks::Integer
)

    return Dict(
        "particles" => nparticles,
        "x_min" => x_min,
        "x_max" => x_max,
        "y_min" => y_min,
        "y_max" => y_max,
        "z_min" => z_min,
        "z_max" => z_max,
        "tasks" => ntasks
    )
end

function generate_uneven(
    nparticles::Integer,
    x_min::Real,
    x_max::Real,
    y_min::Real,
    y_max::Real,
    z_min::Real,
    z_max::Real,
    imbalance::Real=2.0,
    ;
    rng=copy(Random.GLOBAL_RNG)
)

    # particle's array
    nleft = floor(Int, nparticles / (1 + imbalance))
    nright = nparticles - nleft
    a_left = GenerateParticles(nleft, x_min, 0.5 * (x_min + x_max), y_min, y_max, z_min, z_max; rng)
    a_right = GenerateParticles(nright, 0.5 * (x_min + x_max), x_max, y_min, y_max, z_min, z_max; rng)
    a_p = append!(a_left, a_right)
    return a_p
end


npts = isempty(ARGS) ? 500_000 : parse(Int, ARGS[1])
ntasks = checkbounds(Bool, ARGS, 2) ? parse(Int, ARGS[2]) : Threads.nthreads()

settings = SetGeneralParameters(npts, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0, ntasks)

particles = generate_uneven(
    Int32(settings["particles"]),
    settings["x_min"],
    settings["x_max"],
    settings["y_min"],
    settings["y_max"],
    settings["z_min"],
    settings["z_max"]
)

con_dims = (
    settings["x_min"],
    settings["x_max"],
    settings["y_min"],
    settings["y_max"],
    settings["z_min"],
    settings["z_max"],
)

tessellation = VoroPlusPlus.parallel_container(con_dims, particles, ntasks)
# tessellation = BaseTessellation(
#     settings["x_min"],
#     settings["x_max"],
#     settings["y_min"],
#     settings["y_max"],
#     settings["z_min"],
#     settings["z_max"],
#     npts,
#     owner,
#     settings["tasks"]
# )

@info "nthreads = $(Threads.nthreads())" npts ntasks size(tessellation.domain) @btime sum(volume, $tessellation; init=0.0)
