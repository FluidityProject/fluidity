input: clean
	../../bin/gmsh2triangle --2d src/square.msh
	mpiexec -n 2 ../../bin/flredecomp -i 1 -o 2 -v -l 5material-droplet 5material-droplet_flredecomp

clean:
	rm -rf *.halo *.stat *.vtu *.pvtu *.node *.ele *.edge *.log-0 *.err-0 *checkpoint *.convergence
	rm -f *.err-? *.log-?
	rm -rf 5material-droplet_? 5material-droplet_?? \
	matrixdump matrixdump.info
	rm -rf *flredecomp*
