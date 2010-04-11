PROJECT = standing_wave

input: clean
	gmsh -2 src/standing_wave.geo -o standing_wave.msh
	../../scripts/gmsh2triangle --2d standing_wave.msh
	../../bin/fldecomp -n 2 standing_wave

clean:
	rm -f *.vtu *.pvtu *.stat
	rm -f fluidity.err-? fluidity.log-?
	rm -f *.ele *.edge *.node *.msh *.halo
	rm -f matrixdump*