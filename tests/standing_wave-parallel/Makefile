PROJECT = standing_wave

input: clean
	gmsh -2 src/standing_wave.geo -o standing_wave.msh
	../../scripts/gmsh2triangle --2d standing_wave.msh
	../../bin/fldecomp -n 2 standing_wave

clean:
	rm -f *.vtu *.pvtu *.stat *.detectors
	rm -f fluidity.err-? fluidity.log-?
	rm -f *.ele *.edge *.node *.msh *.halo *.face
	rm -rf matrixdump* *_checkpoint*
	rm -rf run1_logs/
