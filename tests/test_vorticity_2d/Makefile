default: input

input: clean
	tar -xvjf src/square.tar.bz2

clean: clean-mesh clean-run
clean-mesh:
	rm -f *.ele *.edge *.face *.node *.vtu *.s *.stat src/*.msh
	rm -f gmsh.log gmsh_err.log
clean-run:
	rm -f matrixdump matrixdump.info
	rm -f fluidity.err-* fluidity.log-*
