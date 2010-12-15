input: clean
	gmsh -3 -bin src/standing_wave.geo -o standing_wave.msh
	../../bin/fldecomp -m gmsh -n 3 standing_wave

clean:
	rm -f standing_wave.stat
	rm -f fluidity.err* fluidity.log*
	rm -f *.msh *.halo
	rm -rf velocity_logs/
	rm -f matrixdump.info
	rm -f petsc_readnsolve.log* petsc_readnsolve.err*
	rm -f velocity_matrixdump pressure_matrixdump
	rm -rf standing_wave_[0-9]*
