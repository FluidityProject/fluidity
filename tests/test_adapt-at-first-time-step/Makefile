input: clean
	gmsh -3 ./src/cube-1.geo
	../../scripts/gmsh2triangle ./src/cube-1.msh
	../../bin/fldecomp -i triangle -n 2 cube-1

clean:
	rm -f *.ele *.face *.node *.halo *.vtu *.s *.d.1 *.stat *.log-* *.err-* *.pvtu *~ ./src/*.msh ./src/*~ 
