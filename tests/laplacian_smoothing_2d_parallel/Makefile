input:	clean
	gmsh -2 -clmax 0.05 src/square.geo
clean:
	rm -f *.halo *.msh *.pvtu *.vtu *.stat *.log* matrixdump* *~ fluidity.* *.err* src/*.msh
	rm -rf laplacian_smoothing_2d_?
