using CairoMakie

npart = 2e6

x = [1, 2, 4, 6, 8, 12, 16, 24, 32, 40, 48, 64]

y1 = [34.434, 16.680, 8.584, 5.934, 4.479, 2.948, 2.385, 1.717, 1.335, 1.141, NaN, NaN]
y2 = [28.914, 14.995, 7.515, 5.305, 5.018, 2.802, 2.665, 1.599, 1.289, 0.920, 0.752, NaN]
y3 = [35.897, 17.558, 8.937, 6.052, 4.477, 3.062, NaN, 1.516, 1.173, NaN, 0.816, 0.667]

h1 = npart ./ y1 ./ x
h2 = npart ./ y2 ./ x
h3 = npart ./ y3 ./ x

fig = Figure()

ax1 = Axis(
    fig[1, 1]
    ;
    xlabel="No. of Threads", ylabel="Particles/(thread × sec)",
)

ax2 = Axis(
    fig[1, 1]
    ;
    yaxisposition = :right, ylabel="Parallel efficiency",
)

hidespines!(ax2)
hidexdecorations!(ax2)

y1t = (!isnan).(y1)
y2t = (!isnan).(y2)
y3t = (!isnan).(y3)
@views scatterlines!(ax1, x[y1t], h1[y1t], label="cHARISMa Node A")
@views scatterlines!(ax1, x[y2t], h2[y2t], label="cHARISMa Node D")
@views scatterlines!(ax1, x[y3t], h3[y3t], label="cHARISMa Node E")
axislegend(; position=:rb)

save("figure3.pdf", fig, pdf_version="1.4")
