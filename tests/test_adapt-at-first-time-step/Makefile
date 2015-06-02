input: clean
	gmsh -3 ./src/cube-1.geo
	cp src/cube-1.msh .

clean:
	rm -f *.ele *.face *.node *.halo *.vtu *.s *.d.1 *.stat *.log-* *.err-* *.pvtu *~ ./src/*.msh ./src/*~  \
	matrixdump matrixdump.info
	rm -rf *flredecomp*
	rm -rf test_adapt-at-first-time-step_[0-9]*
