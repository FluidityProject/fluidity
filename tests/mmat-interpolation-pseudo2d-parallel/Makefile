input: clean
	gmsh -3 src/2d_square.geo
	../../bin/gmsh2triangle src/2d_square.msh
	../../bin/fldecomp -n 2  2d_square

clean:
	rm -rf *.stat *.vtu *.pvtu *.node *.ele *.halo *.face *checkpoint* *convergence* *.log-? *.err-? src/2d_square.msh
