default: input

input: clean
	cd src; gmsh -2 sphere.geo
	../../scripts/gmsh2triangle src/sphere.msh

clean:
	rm -f *.ele *.edge *.face *.node *.vtu *.stat *.log
