input: default
default: clean
	gmsh -format msh2 -2 -bin slice_unstructured.geo  -o slice_unstructured.msh

clean: 
	rm -f *.msh *.vtu *.log *.stat \
	matrixdump matrixdump.info
