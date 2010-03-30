FLMLMODEL = parallel_p2_consistent_interpolation
MESH = square-unstructured

default: input

input: src/$(MESH).tar.bz2
	tar -xvjf src/$(MESH).tar.bz2
	fldecomp -n 2 $(MESH)

clean: clean-mesh clean-run
clean-mesh:
	rm -f $(MESH).edge $(MESH).ele $(MESH).node
	rm -f $(MESH)_?*.edge $(MESH)_?*.ele $(MESH)_?*.node $(MESH)_?*.halo
clean-run:
	rm -f $(FLMLMODEL)_?*.pvtu $(FLMLMODEL)_?*_?*.vtu $(FLMLMODEL).stat fluidity.err-* fluidity.log-*
