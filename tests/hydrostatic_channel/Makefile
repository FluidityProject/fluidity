input: clean
	gmsh -2 -bin src/standing_wave.geo -o standing_wave.msh

clean:
	rm -f *.stat *.detectors *.detectors.dat
	rm -f *.vtu *.pvtu standing_wave_flredecomp.flml
	rm -f fluidity.* flredecomp.*
	rm -f *.msh *.halo
	rm -f matrixdump*
	rm -rf standing_wave_[0-9]*
