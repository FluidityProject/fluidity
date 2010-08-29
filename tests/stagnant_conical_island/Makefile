input:	clean
	gmsh -2 src/2d_island.geo 
	../../scripts/gmsh2triangle --2d src/2d_island.msh
	../../bin/fldecomp -n 3 2d_island

clean:
	rm -f src/*.msh *.vtu *.pvtu *.stat *.log-? *.err-? *.node *.ele *.halo *.edge *_checkpoint* *.convergence *.detectors* matrixdump*
