FLMLMODEL = tri_mixing_bin_test_serial
MESH = square

default: input

input: src/$(MESH).msh
	gmsh2triangle -2 src/$(MESH).msh 

clean: clean-mesh clean-run-debug
clean-mesh:
	rm -f $(MESH).edge $(MESH).ele $(MESH).node
clean-run:
	rm -f $(FLMLMODEL)_?*.vtu $(FLMLMODEL).stat fluidity.err-0 fluidity.log-0
clean-run-debug: clean-run
