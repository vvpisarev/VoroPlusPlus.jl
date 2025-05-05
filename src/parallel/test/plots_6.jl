using CairoMakie

fig = Figure()

ax = Axis(
    fig[1, 1]
)

let
    points = NTuple{2, Float64}[]

    for iy in 0:4
        for ix in 0:10
            x = ix + 0.5 * isodd(iy) + 0.15 * randn()
            y = iy * 0.5 * sqrt(3) + 0.15 * randn()

            push!(points, (x, y))
        end
    end

    voronoiplot!(ax, points; color=:white)
    save("figure6.pdf", fig, pdf_version="1.4")
end
