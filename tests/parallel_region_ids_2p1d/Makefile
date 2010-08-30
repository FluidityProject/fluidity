FLMLMODEL = parallel_region_ids
MESH = square-unstructured
NP = 3

default: input

input: src/$(MESH).geo
	gmsh -2 src/$(MESH).geo
	../../scripts/gmsh2triangle --2d src/$(MESH).msh
	../../bin/fldecomp -n $(NP) $(MESH)

clean: clean-mesh clean-run
clean-mesh:
	rm -f $(MESH).edge $(MESH).ele $(MESH).node
	rm -f $(MESH)_?*.edge $(MESH)_?*.ele $(MESH)_?*.node $(MESH)_?*.halo
	rm -f src/$(MESH).msh 
clean-run:
	rm -f $(FLMLMODEL)_?*.pvtu $(FLMLMODEL)_?*_?*.vtu $(FLMLMODEL).stat fluidity.err-* fluidity.log-*
