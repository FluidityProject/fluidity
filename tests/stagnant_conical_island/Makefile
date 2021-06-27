input:	clean
	gmsh -2 src/2d_island.geo 
	cp src/2d_island.msh .

clean:
	rm -rf *_checkpoint*
	rm -f src/*.msh *.pvtu *.stat *.log-? *.err-? *.msh *.halo *.convergence *.detectors* matrixdump*
	rm -f */*.vtu
	rm -rf conical_island_*
	rm -rf *flredecomp*
