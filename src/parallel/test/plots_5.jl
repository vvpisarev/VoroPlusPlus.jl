using CairoMakie

npart_3M = 3e6

x = [1, 8, 16, 24, 32, 64, 128, 256]

t_a_3M = [51.0, 8.63, 8.9, 8.38, 8.32, 7.15, 7.18, 7.04]

y_a_3M = first(t_a_3M) ./ t_a_3M

fig = Figure()

ax = Axis(
    fig[1, 1]
    ;
    limits=(1, nothing, 1.0, 8.5),
    xlabel="No. of Tasks", ylabel="Speedup with 8 threads",
    yticks=1:8
)

@views scatterlines!(ax, x, y_a_3M;
color=:red, strokecolor=:red, markercolor=:transparent, strokewidth = 1.5, marker=:circle, label="Node A - 3M particles")

@views hlines!(ax, [8.0];
color=:black, linestyle=:dash, label="Ideal speedup")

@views hlines!(ax, [8.0 * 0.9];
color=:blue, linestyle=:dash, label="90% efficiency")

@views hlines!(ax, [8.0 / (4/3)];
color=:blue, linestyle=:dot, label="Expected without load balancing")

axislegend(; position=:rb)

save("figure5.pdf", fig, pdf_version="1.4")
