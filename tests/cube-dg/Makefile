input: clean
	gmsh -3 src/cube-1.geo  -o cube-1.msh
	../../scripts/gmsh2triangle cube-1.msh

clean:
	rm -f cube-1 cube-10 *.ele *.face *.node *.vtu *.s *.d.1 *.stat
