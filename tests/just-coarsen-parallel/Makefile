input: clean
	gmsh -3 -o just-coarsen-parallel.msh just-coarsen-parallel.geo
	../../bin/gmsh2triangle just-coarsen-parallel.msh
	../../bin/fldecomp -n 2  -f just-coarsen-parallel


clean:
	rm -f  *.msh *.ele *.face *.node *.halo *.vtu *.stat matrixdump* *log* *err*
