default: input

input: src/square.tar.bz2
	tar -xvjf src/square.tar.bz2
	../../bin/triangle2gmsh square

clean: clean-mesh clean-run
clean-mesh:
	rm -f square.edge square.ele square.node square.msh
clean-run:
	rm -f square_?*.halo square_?*.msh
	rm -f fldecomp.log
