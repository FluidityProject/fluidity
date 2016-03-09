input: default
default: clean
	gmsh -2 slice_unstructured.geo
	gmsh -2 slice_unstructured_2.geo

clean: 
	rm -f *.msh *.node *.ele *.face *.edge *.vtu *.log *.stat \
	matrixdump matrixdump.info
