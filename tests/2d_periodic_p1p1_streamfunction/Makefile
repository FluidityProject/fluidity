PROJECT = 2d_periodic_p1p1_streamfunction

input: clean
	gmsh -2 -o 2d_square.msh src/2d_square.geo
	../../scripts/gmsh2triangle --2d 2d_square.msh
	
	
	#set up parallel run
	../../bin/periodise 2d_periodic_p1p1_streamfunction.flml
	mpiexec -np 4 ../../bin/flredecomp -i 1 -o 4 2d_periodic_p1p1_streamfunction_periodised 2d_periodic_p1p1_streamfunction_periodised_flredecomp
	spud-set 2d_periodic_p1p1_streamfunction_periodised_flredecomp.flml /simulation_name 2d_periodic_p1p1_streamfunction_parallel

	
clean:
	rm -f *.stat
	rm -f *vtu
	rm -f fluidity* matrixdump*
	rm -f *.ele *.face *.node *.edge
	rm -f *.msh
	rm -f *halo
	rm -f *periodise*