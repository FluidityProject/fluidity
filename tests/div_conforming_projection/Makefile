default: input

input: clean
	gmsh -3 src/basin1.geo > gmsh.log 2> gmsh_err.log -o basin1.msh
	../../bin/gmsh2triangle basin1.msh

clean:
	rm -f *.vtu *.s *.stat *.log *.msh *.node *.ele *.edge
