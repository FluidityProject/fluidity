input:
	gmsh -2 -o unstructured.msh src/unstructured.geo && ../../bin/gmsh2triangle --2d unstructured.msh
	../../bin/fldecomp -n 2 unstructured
	cp sam_output/* .

clean:
	rm -f *.stat *vtu *.log-0 *.err-0 *.ele *.node *.edge *.halo *.msh
	rm -rf backward_facing_step_2d_?
