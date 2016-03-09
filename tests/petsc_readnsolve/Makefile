include ../tools.mk

input: clean
	gmsh -3 -bin src/standing_wave.geo -o standing_wave.msh

clean:
	rm -f standing_wave.stat
	rm -f fluidity.err* fluidity.log*
	rm -f *.msh *.halo
	rm -rf velocity_logs/
	rm -f matrixdump.info
	rm -f petsc_readnsolve.log* petsc_readnsolve.err*
	rm -f velocity_matrixdump pressure_matrixdump
	rm -rf standing_wave_[0-9]*
	rm -rf *flredecomp*

mesh_rename:
	$(SED) --in-place s/failing_velocity_flredecomp_CoordinateMesh/flredecomp_CoordinateMesh/g standing_wave_failing_velocity_flredecomp.flml;
	$(SED) --in-place s/failing_pressure_flredecomp_CoordinateMesh/flredecomp_CoordinateMesh/g standing_wave_failing_pressure_flredecomp.flml;
