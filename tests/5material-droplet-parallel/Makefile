input: clean
	../../bin/gmsh2triangle --2d src/square.msh
	../../bin/fldecomp -n 2  square

clean:
	rm -rf *.halo *.stat *.vtu *.pvtu *.node *.ele *.edge *.log-0 *.err-0 *checkpoint *.convergence
	rm -f *.err-? *.log-?
	rm -rf 5material-droplet_? 5material-droplet_??
