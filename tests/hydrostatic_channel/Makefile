PROJECT = standing_wave

input: clean
	gmsh -2 -bin src/standing_wave.geo -o standing_wave.msh

clean:
	rm -f $(PROJECT).stat
	rm -f *.vtu *.detectors
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.msh
