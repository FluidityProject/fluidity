input:	clean
	gmsh -2 -clmax 0.05 src/square.geo

clean:
	rm -f *.pvtu *.msh *.halo *.vtu *.stat *.log* matrixdump* *~ fluidity.* *.err* src/*.msh
	rm -rf lineal_smoothing_2d_?
