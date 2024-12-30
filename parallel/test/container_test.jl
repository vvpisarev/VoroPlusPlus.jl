
(x_min, x_max) = (y_min, y_max) = (z_min, z_max) = -1, 1
cvol = (x_max - x_min) * (y_max - y_min) * (z_max - z_min)

n_x = n_y = n_z = 6

nparticles = 50

# Create a container with the geometry given above, and make it
# non-periodic in each of the three coordinates. Allocate space for
# eight particles within each computational block
con = Container(
    ;
    bounds = (x_min, x_max, y_min, y_max, z_min, z_max),
    nblocks = (n_x, n_y, n_z),
    periodic = (false, false, false),
    particles_per_block = 8,
)

# Randomly add particles into the container
for i in 1:nparticles
    x = x_min + rand() * (x_max - x_min)
    y = y_min + rand() * (y_max - y_min)
    z = z_min + rand() * (z_max - z_min)
    add_point!(con, i - 1, x, y, z)
end

# Sum up the volumes, and check that this matches the container volume
#=vvol = sum(volume, con)
@info "Container volume" cvol
@info "Voronoi volume" vvol
@info "Difference" cvol - vvol
@test isapprox(vvol, cvol; atol=1e-8)
@test isapprox(sum(volume, VoroPlusPlus.Unsafe(con)), cvol; atol=1e-8)=#

# Output the particle positions in gnuplot format
draw_particles(con, "parallel/particles.gnu")

# Output the Voronoi cells in gnuplot format
draw_cells_gnuplot(con, "parallel/voro_cells.gnu")
