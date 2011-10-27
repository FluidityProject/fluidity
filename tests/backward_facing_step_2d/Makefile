input:
	gmsh -2 -o unstructured.msh src/unstructured.geo && ../../bin/gmsh2triangle --2d unstructured.msh

clean:
	rm -f *.vtu *.log-0 *.err-0 *.ele *.msh *.edge *.node *.stat
