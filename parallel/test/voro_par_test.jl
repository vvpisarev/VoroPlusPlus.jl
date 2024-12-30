
# using Base.Threads

using Base.Threads

vorocell_vector = Array{VoronoiCell}(undef, Threads.nthreads())

Threads.@threads for th = 1:Threads.nthreads()

    v = VoronoiCell()
    init!(v,-1, 1, -1, 1, -1, 1)
    vorocell_vector[th] = v

    for _ in 1:250
        x = 2 * rand() - 1
        y = 2 * rand() - 1
        z = 2 * rand() - 1
        rsq = x * x + y * y + z * z
        if rsq > 0.01 && rsq < 1
            r = 1 / sqrt(rsq)
            x *= r
            y *= r
            z *= r
            add_plane!(vorocell_vector[th], x, y, z, 1)
        end
    end

    draw_gnuplot("parallel/single_thread_cell_$(Threads.threadid()).gnu", v)

end

