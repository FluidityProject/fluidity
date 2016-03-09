default: input

input: src/cube.tar.bz2
	tar -xvjf src/cube.tar.bz2
	../../bin/triangle2gmsh cube

clean: clean-mesh clean-run
clean-mesh:
	rm -f cube.ele cube.face cube.node cube.msh
clean-run:
	rm -f cube_?*.halo cube_?*.msh
	rm -f fldecomp.log
