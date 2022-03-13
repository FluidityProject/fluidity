input: clean
	gmsh -3 -o cube.msh src/cube.geo

clean:
	rm -f *.log* *.err* smoothed_bcs_parallel.flml
	rm -rf out* *.msh *.halo
	rm -rf matrixdump matrixdump.info
