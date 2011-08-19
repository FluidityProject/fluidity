default: input

input: clean
	gmsh -2 src/mid-res.geo > gmsh.log 2> gmsh_err.log
	../../scripts/gmsh2triangle src/mid-res.msh
	mv mid-res.ele sphere.ele
	mv mid-res.edge sphere.edge
	mv mid-res.node sphere.node

clean:
	rm -f *.ele *.edge *.face *.node *.vtu *.s *.stat *.log src/*.msh
