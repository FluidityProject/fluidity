input:	clean
	gmsh -1 -clmax 0.01 src/line.geo
clean:
	rm -f *.halo *.pvtu *.vtu *.stat *.log* matrixdump* *~ fluidity.* *.msh *.err* src/*.msh
	rm -rf laplacian_smoothing_1d_?
