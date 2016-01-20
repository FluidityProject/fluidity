default: input

input: src/square-unstructured.tar.bz2
	tar -xvjf src/square-unstructured.tar.bz2
	../../bin/triangle2gmsh square-unstructured

clean: clean-mesh clean-run
	rm -rf parallel_p0_consistent_interpolation_*
clean-mesh:
	rm -f square-unstructured.edge square-unstructured.ele square-unstructured.node
	rm -f square-unstructured_?*.edge square-unstructured_?*.ele square-unstructured_?*.node square-unstructured_?*.halo
clean-run:
	rm -f parallel_p0_consistent_interpolation_?*.pvtu parallel_p0_consistent_interpolation_?*_?*.vtu parallel_p0_consistent_interpolation.stat fluidity.err-* fluidity.log-* \
	matrixdump matrixdump.info
	rm -rf *flredecomp*
