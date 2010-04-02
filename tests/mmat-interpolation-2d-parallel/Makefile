input: clean
	gmsh -2 src/2d_square.geo
	../../scripts/gmsh2triangle --2d src/2d_square.msh
	../../bin/fldecomp -n 2  2d_square

clean:
	rm -rf *.stat *.vtu *.pvtu *.node *.ele *.halo *.edge *checkpoint* *convergence* matrixdump* *.log-? *.err-? src/2d_square.msh
