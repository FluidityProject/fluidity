default: input

input: src/cube.tar.bz2
	tar -xvjf src/cube.tar.bz2

clean: clean-mesh clean-run
clean-mesh:
	rm -f cube.ele cube.face cube.node
clean-run:
	rm -f cube_?*.ele cube_?*.face cube_?*.halo cube_?*.node
	rm -f fldecomp.log
