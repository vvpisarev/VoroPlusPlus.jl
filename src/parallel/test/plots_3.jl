using Plots

function parpersecthread(par, th, v)
    r = []
    
    for (g, i) in zip(th, v)
        if i != 0.0 && par != 0 && g != 0
            push!(r, (par / i) / g)
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

k1 = parpersecthread(200000, x, y1)
k2 = parpersecthread(200000, x, y2)
k3 = parpersecthread(200000, x, y3)
println(k1)
println(k2)
println(k3)

plot(x, [k1 k2 k3], label=["Cluster Type A/B" "Cluster Type C/D" "Cluster Type E"], lw=[2 2 2])
plot!(legend=:outerright, legendcolumns=3)
title!("particles/sec/thread vs number of threads")
xlabel!("No. of Threads")
ylabel!("particles/sec/thread")
