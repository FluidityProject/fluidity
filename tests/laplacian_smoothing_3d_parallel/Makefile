input:	clean
	gmsh -3 -clmax 0.03 src/unit_cube.geo
clean:
	rm -f *.pvtu *.msh *.halo *.vtu *.stat *.log* matrixdump* *~ fluidity.* *.err* src/*.msh
	rm -rf laplacian_smoothing_3d_?
