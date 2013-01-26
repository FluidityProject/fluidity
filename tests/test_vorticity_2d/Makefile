default: input

input: clean
	gmsh -2 src/square.geo -o square.msh

clean: clean-mesh clean-run
clean-mesh:
	rm -f *.ele *.edge *.face *.node *.msh
clean-run:
	rm -f matrixdump matrixdump.info
	rm -f fluidity.err-* fluidity.log-* *.stat *.vtu
