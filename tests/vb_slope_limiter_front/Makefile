input: default
default: clean
	gmsh -2 slice_unstructured.geo
	../../scripts/gmsh2triangle --2d slice_unstructured.msh
clean: 
	rm -f *.msh *.node *.ele *.face *.edge *.vtu *.log *.stat