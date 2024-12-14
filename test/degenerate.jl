
@testset "Degenerate" begin
    
    pi_ = 3.1415926535897932384626433832795
    n = 32
    step = 2*pi_/n
    theta = pi_/4-0.25

	v = VoronoiCell()
	
	init!(v, -1, 1, -1, 1, -1, 1)

	#for phi=0;phi<2*pi_-0.5*step;phi+=step
    for phi in 0:step:(2*pi_-0.5*step)

		x=cos(theta)
        y=cos(phi)*sin(theta)
        z=sin(phi)*sin(theta)
		add_plane!(v, x, y, z, 1)
		add_plane!(v, -x, y, z, 1)
		add_plane!(v, y, x, z, 1)
		add_plane!(v, y, -x, z, 1)
		add_plane!(v, y, z, x, 1)
		add_plane!(v, y, z, -x, 1)
    end

	@test check_relations(v) === nothing
	@test check_duplicates(v) === nothing
	@test draw_gnuplot(v, 0, 0, 0, "degenerate.gnu") === nothing
	@test draw_pov(v, 0, 0, 0, "degenerate_v.pov") === nothing

end



@testset "Degenerate_02" begin
	
	pi_ = 3.1415926535897932384626433832795
	points = 100
	n = 64
	step = 2*pi/n
	theta = 0.04

	v = VoronoiCell()
	n_ = 0

	init!(v, -1, 1, -1, 1, -1, 1)

	while n_<points 
		
		x = 2*rand()-1
		y = 2*rand()-1
		z = 2*rand()-1

		rsq = x*x+y*y+z*z

		if rsq<0.01 || rsq>1
			continue
		end

		r = 1/sqrt(rsq)
		x *= r
		y *= r
		z *= r

		rsq = sqrt(x*x+y*y)
		r = z/rsq

		for phi in rand()*step:step:(2*pi_)
		# for(phi=rnd()*step;phi<2*pi;phi+=step)
			add_plane!(v, x*cos(theta)+sin(theta)*(-y*cos(phi)/rsq-x*r*sin(phi)),
				y*cos(theta)+sin(theta)*(x*cos(phi)/rsq-y*r*sin(phi)),
				z*cos(theta)+sin(theta)*rsq*sin(phi),1);
		end

		n_ += 1

	end

	@test draw_gnuplot(v, 0, 0, 0, "degenerate2.gnu") === nothing
	@test draw_pov(v, 0,0,0,"degenerate2_v.pov") === nothing
	@test draw_pov_mesh(v, 0,0,0,"degenerate2_m.pov") === nothing

end