input: clean
	../../scripts/gmsh2triangle --2d src/square.msh
	../../bin/fldecomp -n 2 -i triangle square

clean:
	rm -rf *.stat *.vtu *.pvtu *.node *.ele *.edge *.log-0 *.err-0 *checkpoint *.convergence
