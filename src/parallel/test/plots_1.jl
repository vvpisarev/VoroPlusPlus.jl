using CairoMakie

x = [1, 2, 4, 6, 8, 12, 16, 24, 32, 40, 48, 64]

y_a_2M = [34.434, 16.680, 8.584, 5.934, 4.479, 2.948, 2.385, 1.717, 1.335, 1.141, NaN, NaN]
y_d_2M = [28.914, 14.995, 7.515, 5.305, 5.018, 2.802, 2.665, 1.599, 1.289, 0.920, 0.752, NaN]
y_e_2M = [35.897, 17.558, 8.937, 6.052, 4.477, 3.062, NaN, 1.516, 1.173, NaN, 0.816, 0.667]

y_a_20M = [334.322, 167.974, 87.271, 59.152, 43.485, 29.304, 22.023, 16.809, 14.248, 10.674, NaN, NaN]
y_d_20M = [289.480, 147.410, 74.099, 65.041, 39.759, 32.246, 20.988, 14.275, 11.236, 8.172, 7.652, NaN]
y_e_20M = [358.326, 179.526, 89.889, 60.250, 44.907, 30.461, NaN, 15.770, 12.697, NaN, 8.480, 7.460]

fig = Figure()

ax = Axis(
    fig[1, 1]
    ;
    xscale=log2, yscale=log2,
    xticks=2 .^(0:6),
    yticks=2.0.^(-1:2:9),
    xlabel="No. of Threads", ylabel="Wall Time (s)",
)
y1t = (!isnan).(y_a_2M)
y2t = (!isnan).(y_d_2M)
y3t = (!isnan).(y_e_2M)
@views scatterlines!(ax, x[y1t], y_a_2M[y1t];
color=:red, strokecolor=:red, markercolor=:transparent, strokewidth = 1.5, marker=:circle, label="Node A - 2M")
@views scatterlines!(ax, x[y2t], y_d_2M[y2t];
color=:green, strokecolor=:green, markercolor=:transparent, strokewidth = 1.5, marker=:rect, label="Node D - 2M")
@views scatterlines!(ax, x[y3t], y_e_2M[y3t];
color=:blue, strokecolor=:blue, markercolor=:transparent, strokewidth = 1.5, marker=:utriangle, label="Node E - 2M")

@views scatterlines!(ax, x[y1t], y_a_20M[y1t], color=:red, marker=:circle, label="Node A - 20M")
@views scatterlines!(ax, x[y2t], y_d_20M[y2t], color=:green, marker=:rect, label="Node D - 20M")
@views scatterlines!(ax, x[y3t], y_e_20M[y3t], color=:blue, marker=:utriangle, label="Node E - 20M")
lines!(ax, [1, 64], [25.0, 25/64], color=:black, linestyle=:dash, label="Ideal scaling")
lines!(ax, [1, 64], [250.0, 250/64], color=:black, linestyle=:dash)
axislegend(; position=:lb)

save("figure1.pdf", fig, pdf_version="1.4")
