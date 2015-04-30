default: input

input: clean
	gmsh -3 src/basin1.geo > gmsh.log 2> gmsh_err.log -o basin1.msh
	gmsh -3 src/basin2.geo > gmsh.log 2> gmsh_err.log -o basin2.msh
	gmsh -3 src/basin3.geo > gmsh.log 2> gmsh_err.log -o basin3.msh

clean:
	rm -f *.vtu *.s *.stat *.log *.msh *.node *.ele *.edge
