MODEL=flredecomp-2d-fieldweighted

default: input

input:  clean
	gmsh -2 -bin -optimize src/Subduction_Mesh.geo -o Subduction_Mesh.msh

clean:
	rm -f *.pvtu *.stat *.msh *.halo *.detectors fluidity.* Parallel* *.pyc

