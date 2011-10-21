FLMLMODEL = parallel_p2_consistent_interpolation
MESH = cube-unstructured

default: input

input: src/$(MESH).geo
	gmsh -3 src/$(MESH).geo
	../../bin/gmsh2triangle src/$(MESH).msh
	../../bin/fldecomp -n 3 $(MESH)

clean: clean-mesh clean-run
	rm -rf parallel_p2_consistent_interpolation_*
clean-mesh:
	rm -f $(MESH).face $(MESH).ele $(MESH).node
	rm -f $(MESH)_?*.face $(MESH)_?*.ele $(MESH)_?*.node $(MESH)_?*.halo
	rm -f src/$(MESH).msh
clean-run:
	rm -f $(FLMLMODEL)_?*.pvtu $(FLMLMODEL)_?*_?*.vtu $(FLMLMODEL).stat fluidity.err-* fluidity.log-*
