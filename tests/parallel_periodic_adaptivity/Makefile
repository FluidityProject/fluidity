default: input

input: clean
	../../bin/gmsh2triangle --2d src/MMS_C.msh
	../../bin/periodise periodic_adaptivity.flml
	mpiexec -n 3 ../../bin/flredecomp -i 1 -o 3  periodic_adaptivity_periodised periodic_adaptivity_periodised_parallel

clean: clean-mesh clean-run
clean-mesh:
	rm -f *.ele *.edge *.face *.node *.halo 
	rm -f gmsh.log gmsh_err.log
	rm -f periodic_adaptivity_periodised*.flml
clean-run:
	rm -f *.pvtu *.vtu *.detectors *.detectors.dat *.stat
	rm -f matrixdump matrixdump.info
	rm -f fluidity.err-* fluidity.log-*
	rm -rf periodic_adaptivity_?
	rm -rf surface_*
	rm -rf mesh_*
	rm -rf *log-?
	rm -f edge_weights*

