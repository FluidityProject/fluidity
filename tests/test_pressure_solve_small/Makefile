input: clean
	gmsh -3 cube_unstructured_small.geo
	../../bin/gmsh2triangle cube_unstructured_small.msh

clean:
	rm -f *.msh *.ele *.face *.node *.dat
