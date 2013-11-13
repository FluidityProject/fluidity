input: clean
	gmsh -2 -bin src/channel.geo -o channel.msh

clean:
	rm -f *.dat *.vtu *.pvtu *.stat *.detectors
	rm -f fluidity.err-? fluidity.log-?
	rm -f *.ele *.edge *.node *.msh *.halo *.face
	rm -rf matrixdump* *_checkpoint*
