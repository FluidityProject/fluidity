FLMLMODEL = line_mixing_bin_test_serial
MESH = line

default: input

input:

clean: clean-mesh clean-run-debug
clean-mesh:
clean-run:
	rm -f $(FLMLMODEL)_?*.vtu $(FLMLMODEL).stat fluidity.err-0 fluidity.log-0
clean-run-debug: clean-run
