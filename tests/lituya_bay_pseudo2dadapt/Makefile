input:	clean
	gmsh -3 src/inletslice.geo
	cp src/inletslice.msh .


clean:
	rm -f *.vtu *.stat *.log* *.err* src/*.msh *.msh *_checkpoint* *.convergence *.pvtu matrixdump*
