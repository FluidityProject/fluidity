default: input

input: src/square.tar.bz2
	tar -xvjf src/square.tar.bz2
	../../bin/triangle2gmsh square

clean: clean-mesh clean-run
clean-mesh:
	rm -f square.edge square.ele square.node
clean-run:
	rm -f square_?*.edge square_?*.ele square_?*.halo square_?*.node
	rm -f *.log* *.err* test_flredecomp.flml
	rm -f test_flredecomp_CoordinateMesh_*
