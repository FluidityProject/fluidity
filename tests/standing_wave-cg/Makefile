PROJECT = standing_wave

input: clean
	gmsh -3 src/standing_wave.geo
	mv src/standing_wave.msh .
	../../scripts/gmsh2triangle standing_wave.msh

clean:
	rm -f $(PROJECT) $(PROJECT)0 $(PROJECT).fld
	rm -f $(PROJECT).s $(PROJECT).stat
	rm -f *.vtu *.detectors
	rm -f fluidity.err-0 fluidity.log-0
	rm -f *.ele *.face *.node *.msh
