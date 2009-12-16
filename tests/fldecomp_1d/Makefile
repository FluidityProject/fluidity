default: input

input: src/line.tar.bz2
	tar -xvjf src/line.tar.bz2

clean: clean-mesh clean-run
clean-mesh:
	rm -f line.bound line.ele line.node
clean-run:
	rm -f line_?*.bound line_?*.ele line_?*.halo line_?*.node
	rm -f fldecomp.log
