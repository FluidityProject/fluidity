input: clean
	gmsh -format msh2 -3 src/pseudo2d_square.geo
	mv src/pseudo2d_square.msh .

clean:
	rm -rf *.halo *.stat *.detectors *.pvtu *.vtu *checkpoint* *convergence* *.log-? *.err-? *.geo *.msh
