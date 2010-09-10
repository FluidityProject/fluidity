input:	clean
	gmsh -2 src/2d_island.geo 
	../../scripts/gmsh2triangle --2d src/2d_island.msh
	../../bin/fldecomp -n 3 2d_island

clean:
	rm -rf *_checkpoint*
	rm -f src/*.msh *.pvtu *.stat *.log-? *.err-? *.node *.ele *.halo *.edge *.convergence *.detectors* matrixdump*
	rm -f */*.vtu
	rm -rf conical_island_*
