using CairoMakie

x = [1, 2, 4, 6, 8, 12, 16, 24, 32, 40, 48, 64]

y1 = [34.434, 16.680, 8.584, 5.934, 4.479, 2.948, 2.385, 1.717, 1.335, 1.141, NaN, NaN]
y2 = [28.914, 14.995, 7.515, 5.305, 5.018, 2.802, 2.665, 1.599, 1.289, 0.920, 0.752, NaN]
y3 = [35.897, 17.558, 8.937, 6.052, 4.477, 3.062, NaN, 1.516, 1.173, NaN, 0.816, 0.667]

fig = Figure()

ax = Axis(
    fig[1, 1]
    ;
    xscale=log2, yscale=log2,
    xticks=2 .^(0:6),
    yticks=2.0.^(-1:5),
    xlabel="No. of Threads", ylabel="Wall Time (s)",
)
y1t = (!isnan).(y1)
y2t = (!isnan).(y2)
y3t = (!isnan).(y3)
@views scatterlines!(ax, x[y1t], y1[y1t], label="cHARISMa Node A")
@views scatterlines!(ax, x[y2t], y2[y2t], label="cHARISMa Node D")
@views scatterlines!(ax, x[y3t], y3[y3t], label="cHARISMa Node E")
lines!(ax, [1, 64], [25.0, 25/64], label="Ideal scaling")
axislegend()

save("figure1.pdf", fig, pdf_version="1.4")
