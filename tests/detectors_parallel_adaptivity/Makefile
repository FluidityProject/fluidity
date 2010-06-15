default: input

input: clean
	../../scripts/gmsh2triangle --2d src/square.msh
	fldecomp -n 2 -v square

clean: clean-mesh clean-run
clean-mesh:
	rm -f *.ele *.edge *.face *.node *.halo *.pvtu *.vtu *.detectors *.detectors.dat *.stat
	rm -f gmsh.log gmsh_err.log
clean-run:
	rm -f matrixdump matrixdump.info
	rm -f fluidity.err-* fluidity.log-*
