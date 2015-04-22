default: input

input: clean
	gmsh -3 src/basin.geo > gmsh.log 2> gmsh_err.log
	cp src/basin.msh .

clean:
	rm -f *.ele *.edge *.face *.node *.vtu *.s *.stat *.log src/*.msh
