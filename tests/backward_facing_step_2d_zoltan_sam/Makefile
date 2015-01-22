input:
	gmsh -2 -o unstructured.msh src/unstructured.geo
	cp -r sam_output/* .

clean:
	rm -f *.stat *vtu *.log-0 *.err-0 *.ele *.node *.edge *.halo *.msh
	rm -rf backward_facing_step_2d_?
	rm -rf backward_facing_step_2d_sam_9* \
	matrixdump matrixdump.info
	rm -rf *flredecomp*
