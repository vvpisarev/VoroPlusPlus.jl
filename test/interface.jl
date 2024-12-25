@testset "Find Voro Cell" begin

    h = 0.05
    hcube = h*h*h
    particles = 20

	rx = Ref(0.0)
	ry = Ref(0.0)
	rz = Ref(0.0)
	pid = Ref(Int32(0))

	con =  Container(0,1,0,1,0,1,5,5,5,false,false,false,8)

	#for i in 0:1:particles
	#	x = rand()
	#	y = rand()
	#	z = rand()
	#	add_point!(con, i, x, y, z)
    #end

	add_point!(con, 18, 0.890233, 0.348893, 0.0641713)
	add_point!(con, 8, 0.156679, 0.400944, 0.12979)
	add_point!(con, 19, 0.020023, 0.457702, 0.0630958)
	add_point!(con, 6, 0.141603, 0.606969, 0.0163006)
	add_point!(con, 1, 0.79844, 0.911647, 0.197551)
	add_point!(con, 2, 0.335223, 0.76823, 0.277775)
	add_point!(con, 9, 0.108809, 0.998925, 0.218257)
	add_point!(con, 12, 0.493583, 0.972775, 0.292517)
	add_point!(con, 14, 0.400229, 0.891529, 0.283315)
	add_point!(con, 11, 0.296032, 0.637552, 0.524287)
	add_point!(con, 16, 0.0697553, 0.949327, 0.525995)
	add_point!(con, 17, 0.0860558, 0.192214, 0.663227)
	add_point!(con, 0, 0.840188, 0.394383, 0.783099)
	add_point!(con, 3, 0.55397, 0.477397, 0.628871)
	add_point!(con, 13, 0.771358, 0.526745, 0.769914)
	add_point!(con, 5, 0.916195, 0.635712, 0.717297)
	add_point!(con, 10, 0.512932, 0.839112, 0.61264)
	add_point!(con, 7, 0.242887, 0.137232, 0.804177)
	add_point!(con, 4, 0.364784, 0.513401, 0.95223)
	add_point!(con, 15, 0.352458, 0.807725, 0.919026)


	@test draw_particles(con, "find_voro_cell_p.gnu") === nothing

    f1 = open("find_voro_cell.vec","w")

    for x in 0.5*h:h:1 
        for y in 0.5*h:h:1
            if (find_voronoi_cell(con, x, y, 0.5, rx, ry, rz, pid))
                write(f1, "$(x) $(y) $(0.5) $(rx[]-x) $(ry[]-y) $(rz[]-0.5) $(sqrt( (rx[]-x)*(rx[]-x)+(ry[]-y)*(ry[]-y)+(rz[]-0.5)*(rz[]-0.5) ))\n")
            #else
                #fprintf(stderr,"# find_voronoi_cell error for %g %g 0.5\n",x,y)
            end
        end
    end

	close(f1)

	samp_v = Array{Int32}(undef, particles)

	for z in 0.5*h:h:1
		for y in 0.5*h:h:1 
			for x in 0.5*h:h:1
				if ( find_voronoi_cell(con, x, y, z, rx, ry, rz, pid) )
					samp_v[pid[]+1] += 1 
					#if pid[] < 20
						#tmp = samp_v[pid[]+1]
						#tmp =+ 1
						#samp_v[pid[]+1] = tmp
					#end
				#else 
					#fprintf(stderr,"# find_voronoi_cell error for %g %g %g\n",x,y,z);
				end
			end
		end
	end

	f1 = open("find_voro_cell.vol","w")

	cla = Container_Iterator(con)
	c = VoronoiCell()

	x = Ref(0.0)
	y = Ref(0.0)
	z = Ref(0.0)
	r = Ref(0.0)

	if ( start!(cla) )
		if ( compute_cell!(c, con, cla) )
			pos(cla, pid, x, y, z, r);
			write(f1,"$(pid[]) $(x[]) $(y[]) $(z[]) $(volume(c)) $(samp_v[pid[]+1]*hcube)\n")
			draw_gnuplot("find_voro_cell_v.gnu", c, (x[], y[], z[]))
		end
		while ( next!(cla) )
			if ( compute_cell!(c, con, cla) )
				pos(cla, pid, x, y, z, r);
				write(f1,"$(pid[]) $(x[]) $(y[]) $(z[]) $(volume(c)) $(samp_v[pid[]+1]*hcube)\n")
				draw_gnuplot("find_voro_cell_v.gnu", c, (x[], y[], z[]))
			end
		end
	end

	close(f1)

end