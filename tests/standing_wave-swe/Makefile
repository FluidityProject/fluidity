PROJECT = standing_wave

input: clean
	gmsh -2 -bin src/standing_wave.geo -o standing_wave.msh

clean:
	rm -f *.dat *.vtu *.pvtu *.stat *.detectors.h5part
	rm -f fluidity.err-? fluidity.log-?
	rm -f *.ele *.edge *.node *.msh *.halo *.face
	rm -rf matrixdump* *_checkpoint*
