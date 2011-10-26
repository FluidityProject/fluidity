default: input

input: clean
	gmsh -2 square.geo > gmsh.log 2> gmsh_err.log
	../../bin/gmsh2triangle -2 square.msh

clean:
	rm -f *.ele *.edge *.face *.node *.vtu *.s *.stat *.log *.msh
