input: clean
	gmsh -3 -o oodc-unstructured.msh oodc-unstructured.geo 

clean:
	rm -f *.msh *.halo *.vtu *.stat matrixdump* *log* *err*
	rm -rf oodc-unstructured-parallel_?*
	rm -rf *flredecomp*
