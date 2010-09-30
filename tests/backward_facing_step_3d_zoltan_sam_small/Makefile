input: clean
	gmsh -3 -o unstructured.msh src/unstructured.geo 
	../../scripts/gmsh2triangle unstructured.msh
	../../bin/fldecomp -n 3  -f unstructured

clean:
	rm -rf backward_facing_step_3d_? backward_facing_step_3d_?.pvtu *.log *vtu *.log-0 *.err-0 *.msh *.stat *.node *.ele *.face *.halo
