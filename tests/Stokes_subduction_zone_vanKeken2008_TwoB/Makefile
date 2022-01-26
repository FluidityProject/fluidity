MODEL=vanKekenTwoB

default: input

input:  clean
	gmsh -2 -optimize src/Subduction_Mesh.geo
	cp src/Subduction_Mesh.msh .

clean:
	rm -f  *.vtu *.stat src/*.msh *.msh *halo *.detectors fluidity.* *.dat *Parallel* *.pyc
	rm -R -f $(MODEL)_*

