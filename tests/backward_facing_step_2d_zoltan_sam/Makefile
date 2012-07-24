input:
	gmsh -2 -o unstructured.msh src/unstructured.geo && ../../bin/gmsh2triangle --2d unstructured.msh
	mpiexec -n 2 ../../bin/flredecomp -i 1 -o 2 -v -l backward_facing_step_2d backward_facing_step_2d_flredecomp
	cp -r sam_output/* .

clean:
	rm -f *.stat *vtu *.log-0 *.err-0 *.ele *.node *.edge *.halo *.msh
	rm -rf backward_facing_step_2d_?
	rm -rf backward_facing_step_2d_sam_9*
	rm -rf *flredecomp*
