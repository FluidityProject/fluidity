input:
	##gmsh -3 -o unstructured.msh src/unstructured.geo 
	../../scripts/gmsh2triangle unstructured.msh
	../../bin/fldecomp -n 64 -i triangle -f unstructured

clean:
	rm -f *.log *.vtu *.log-0 *.err-0 *.msh *.stat *.node *.ele *.face *.halo
