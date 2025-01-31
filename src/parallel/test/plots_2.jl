using Plots

function parpersec(par, v)
    r = []
    for i in v
        if i != 0.0 && par != 0
            push!(r, par / i)
        else
            push!(r, 0.0)
        end
    end
    return r
end

x = [1, 2, 4, 6, 8, 12, 16, 24, 32, 40, 48, 64]

y1 = [34.434, 16.680, 8.584, 5.934, 4.479, 2.948, 2.385, 1.717, 1.335, 1.141, 0.0, 0.0]
y2 = [28.914, 14.995, 7.515, 5.305, 5.018, 2.802, 2.665, 1.599, 1.289, 0.920, 0.752, 0.0]
y3 = [35.897, 17.558, 8.937, 6.052, 4.477, 3.062, 0.0, 1.516, 1.173, 0.0, 0.816, 0.667]

h1 = parpersec(200000, y1)
h2 = parpersec(200000, y2)
h3 = parpersec(200000, y3)
println(h1)
println(h2)
println(h3)

plot(x, [h1 h2 h3], label=["Cluster Type A/B" "Cluster Type C/D" "Cluster Type E"], lw=[2 2 2])
plot!(legend=:outerright, legendcolumns=3)
title!("particles/sec vs number of threads")
xlabel!("No. of Threads")
ylabel!("particles/sec")

