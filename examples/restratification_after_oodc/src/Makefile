preprocess:
	gmsh -2 circle.geo 
	../../../scripts/gmsh2triangle circle.msh; 
	../../../bin/fldecomp -n 32  circle 

input: preprocess

clean:
	rm -f *.ele *.edge *.face *.node *.halo *.poly *.vtu *.stat *.log-0 *.err-0 */matrixdump*
