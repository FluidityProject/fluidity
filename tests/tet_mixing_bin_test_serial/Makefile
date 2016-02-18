FLMLMODEL = tet_mixing_bin_test_serial
MESH = cube

default: input

input: src/$(MESH).msh
	cp src/$(MESH).msh . 

clean: clean-mesh clean-run-debug
clean-mesh:
	rm -f $(MESH).face $(MESH).ele $(MESH).node
clean-run:
	rm -f $(FLMLMODEL)_?*.vtu $(FLMLMODEL).stat fluidity.err-0 fluidity.log-0 \
	matrixdump matrixdump.info
clean-run-debug: clean-run
