input:	clean
	gmsh -2 src/inletslice.geo
	../../scripts/gmsh2triangle --2d src/inletslice.msh
	fldecomp -n 8 -i triangle inletslice


clean:
	rm -f *.vtu *.stat *.log* *.err* *.msh *.halo *.node *.ele *.edge *_checkpoint* *.convergence *.pvtu matrixdump*
