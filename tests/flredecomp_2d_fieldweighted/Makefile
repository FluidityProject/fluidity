MODEL=flredecomp-2d-fieldweighted

default: input

input:  clean
	gmsh -2 -optimize src/Subduction_Mesh.geo
	cp src/Subduction_Mesh.msh .
	../../bin/gmsh2triangle --2d -i Subduction_Mesh.msh

clean:
	rm -f  *.ele *.edge *.node *.vtu *.stat *.msh *.detectors fluidity.* Parallel*
	rm -R -f $(MODEL)_*

