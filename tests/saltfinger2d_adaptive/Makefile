input: clean
	cp src/unstructured-2.geo .
	gmsh -2 -algo del2d unstructured-2.geo

clean:
	rm -fr saltfinger-2d_* unstructured-2.* *.stat *.msh *.halo
	rm -rf *flredecomp*
