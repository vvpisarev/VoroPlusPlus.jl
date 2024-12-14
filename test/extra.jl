
@testset "Box Cut" begin
    
    cx=1.5
    cy=1.5
    cz=1.5

	v = VoronoiCell()
	init!(v, -8, 8, -8, 8, -8, 8) === nothing

    for x in (cx-0.5):0.1:(cx+0.55)
        for y in (cy-0.5):0.1:(cy+0.55)
            for z in (cz-0.5):0.1:(cz+0.55)
                add_plane!(v, x, y, z)
            end
        end
    end
        
	@test draw_gnuplot(v, 0, 0, 0, "box_cut.gnu") === nothing

	init!(v, cx-0.5, cx+0.5, cy-0.5, cy+0.5, cz-0.5, cz+0.5)
	@test draw_gnuplot(v, 0, 0, 0, "box_cut.points") === nothing

end


@testset "Cut Region" begin
    
    pi_ = 3.1415926535897932384626433832795
    tolwidth = 1e-7
    theta_step = pi_/200

	v = VoronoiCell()

	#FILE *fp=safe_fopen("cell_cut_region.gnu","w");
    fp = open("cell_cut_region.gnu","w")

	@test init_octahedron(v, 1) === nothing

	add_plane!(v, 0.4, 0.3, 1, 0.1)
	draw_gnuplot(v, 0, 0, 0, "cell.gnu")

    for theta in theta_step*0.5:pi_:theta_step

        phi_step = 2*pi/ ( trunc( Int, 2*pi*sin(theta)/theta_step ) + 1 )

        for phi in phi_step*0.5:2*pi_:phi_step

			x=sin(theta)*cos(phi);
			y=sin(theta)*sin(phi);
			z=cos(theta);

			rmin = 0
            rmax = 1
            
            while ( plane_intersects(v, x, y, z, rmax) )
                rmax *= 2
            end

            while ( rmax-rmin > tolwidth )
                r = ( rmax + rmin )*0.5
				if ( plane_intersects(v, x, y, z, r) )
                    rmin=r
				else 
                    rmax=r
                end
            end

			r = ( rmax + rmin )*0.5
			x *= r
            y *= r
            z *= r

			#fprintf(fp,"%g %g %g\n",x,y,z);
            write(fp, "$(x) $(y) $(z)\n")
        end
    end

    close(fp)

end


@testset "Superellipsoid" begin
    
    #double x,y,z,rsq,r;
	v = VoronoiCell()

	init!(v, -1, 1, -1, 1, -1, 1)

	
	for i in 0:5000:1
		x = 2*rand()-1
		y = 2*rand()-1
		z = 2*rand()-1
		rsq = x*x*x*x+y*y*y*y+z*z*z*z

		if ( rsq>0.01&&rsq<1 )
            r = 1/sqrt( sqrt(rsq) )
			x *= r
            y *= r
            z *= r
			add_plane!( v, x*x*x, y*y*y, z*z*z, x*x*x*x+y*y*y*y+z*z*z*z )
        end
	end

	@test draw_gnuplot(v, 0, 0, 0, "superellipsoid.gnu") === nothing
	@test draw_pov(v, 0, 0, 0, "superellipsoid_v.pov") === nothing
	@test draw_pov_mesh(v, 0, 0, 0, "superellipsoid_m.pov") === nothing

end