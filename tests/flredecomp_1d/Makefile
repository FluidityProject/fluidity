default: input

input:
	../../bin/interval --dx=0.0008 0 1 line

clean: clean-mesh clean-run
clean-mesh:
	rm -f line.bound line.ele line.node line.msh
clean-run:
	rm -f line_?*.bound line_?*.ele line_?*.halo line_?*.node
	rm -f *.log* *.err* test_flredecomp.flml
	rm -f test_flredecomp_CoordinateMesh_*
