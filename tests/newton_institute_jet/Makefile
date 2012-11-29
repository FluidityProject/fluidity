PROJECT = standing_wave

input: clean
	gmsh -2 -bin src/square.geo -o square.msh

clean:
	rm -f *.dat *.vtu *.pvtu *.stat *.detectors
	rm -f fluidity.err-? fluidity.log-?
	rm -f *.ele *.edge *.node *.msh *.halo *.face
	rm -rf matrixdump* *_checkpoint*
