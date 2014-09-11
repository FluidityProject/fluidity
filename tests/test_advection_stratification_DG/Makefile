input: default
default: clean
	gmsh -2 slice_unstructured.geo
	gmsh2triangle --2d slice_unstructured.msh
	gmsh -2 slice_unstructured_2.geo
	gmsh2triangle --2d slice_unstructured_2.msh

clean: 
	rm -f *.msh *.node *.ele *.face *.edge *.vtu *.log *.stat \
	matrixdump matrixdump.info
