input:
	##gmsh -3 -o unstructured.msh src/unstructured.geo 
	../../bin/gmsh2triangle unstructured.msh
	../../bin/fldecomp -n 64  -f unstructured

clean:
	rm -f *.log *.vtu *.log-0 *.err-0 *.msh *.stat *.node *.ele *.face *.halo
