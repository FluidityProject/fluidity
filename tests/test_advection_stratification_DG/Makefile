input: default
default: clean
	gmsh -2 slice_unstructured.geo
	../../bin/gmsh2triangle --2d slice_unstructured.msh
	gmsh -2 slice_unstructured_2.geo
	../../bin/gmsh2triangle --2d slice_unstructured_2.msh

clean: 
	rm -f *.msh *.node *.ele *.face *.edge *.vtu *.log *.stat